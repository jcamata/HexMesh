close all
clear all
clc

% set bounds (in degrees and minutes)
% [ minlat° minlat' maxlat° maxlat' minlon° minlon' maxlon° maxlon']
% NOTE THAT YOU SHOULD CONSIDER A LARGER AREA THAN WHAT YOU INTEND TO MESH
% l = [40 00 52 00   -5 00    8 00]; % France
% l = [43 40 44 00    5 30    6 00]; % Cadarache, France
% l = [47 00 48 00   -4 00   -3 00]; % Belle-Ile, France
% l = [38 00 39 00   20 00   21 00]; % Kefalonia, Greece
%l = [37 10 37 40  138 15  138 55]; % Kashiwazaki, Japan
% l = [18 30 21 00 -157 00 -154 00]; % Mauna Loa, Hawai
l = [38 00 38 32 20 10 20 56]; % test small Kefalonia, Greece
% l = [38 03 38 31 20 18 20 50]; % test small Kefalonia, Greece
% l = [20 00 21 00 -156 00 -155 00]; % test small Mauna Loa, Hawai

% choose output directory
outdir = '.';

%choose output name
output_name = 'Kefalonia';

% characteristic length over which details of the coastline are removed
% put H<=0 if you want no smoothing (this is very expensive because many
% small elements will be created)
H = .05; % in units of lon/lat

% minimum water depth
minwater = -50;

% set the path for wget function
wget_path = '/opt/local/bin';
% set the path for dos2unix function
dos2unix_path = '/opt/local/bin';
% set the path for stl2gts function
stl2gts_path = '/opt/local/bin';
setenv('PATH', [getenv('PATH') ':',wget_path,':',dos2unix_path,':',stl2gts_path]);

% after the STL files are generated, you should run in Terminal
% >> dos2unix bathy.stl
% >> stl2gts < bathy.stl > bathy.gts
% >> dos2unix topo.stl
% >> stl2gts < topo.stl > topo.gts
% and move the GTS files in the ./input directory for use by HexMesh

%--------------------------------------------------------------------------
% DO NOT MODIFY BELOW THIS LINE
%--------------------------------------------------------------------------
% todo
% 1) vérifier que les interpolations des NaN sont bien faites
% 2) changer la lecture des données pour que le rétrécissement soit bien
%    pris en compte dans HexMesh
% 3) refaire le changement de z dans la bathymétrie

% transform coordinates
latbnds = l(1):(l(3)+1);
latcrop = [l(1)+l(2)/60 l(3)+l(4)/60];
lonbnds = l(5):(l(7)+1);
loncrop = [l(5)+l(6)/60 l(7)+l(8)/60];

% get topography and plot
topo = topographySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
    'crop', [latcrop loncrop]);
nmax = 1e4;
if length(topo.lon)>nmax || length(topo.lat)>nmax
    error('the STL file will be too large')
end

% get bathymetry and plot
bathy = bathymetrySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
    'crop', [latcrop loncrop]);

% get coastlines (for now only treating Ocean and Land)
[x,y] = ndgrid(latbnds(1:end-1), lonbnds(1:end-1));
x = x(:); y = y(:);
water = struct('Ocean',{[]},'River',{[]},'Lake',{[]}, 'Land',{[]},'Isle',{[]});
for i1 = 1:length(x)
    wat = swbd_shore( [ y(i1) y(i1)+1 x(i1) x(i1)+1 ], outdir );
    if isempty(wat)
        ocean = [ y(i1) y(i1) y(i1)+1 y(i1)+1 y(i1);
            x(i1) x(i1)+1 x(i1)+1 x(i1) x(i1) ];
        water.Ocean = [water.Ocean; {ocean}];
    else
        water.Ocean = [water.Ocean; wat.Ocean];
        %        water.River = [water.River; wat.River];
        %        water.Lake = [water.Lake; wat.Lake];
        water.Land = [water.Land; wat.Land];
        %        water.Isle = [water.Isle; wat.Isle];
    end
end

% smooth coastlines
h = figure;
plotCoastline( water.Ocean, 'k-', h );
plotCoastline( water.Land, 'b-', h );
water = smoothCurve(water,H);
plotCoastline( water.Ocean, 'r--o', h );
plotCoastline( water.Land, 'r--o', h );

%% bathy 1 -> only the bathymetry and the 
% vertical elements in the coastlines

% integrate bathymetry and coastlines
zz = 0;
bathy1 = altimetryCoastline( bathy, water, 3, zz);
[dist,distc] = distanceCoastline( bathy1, water );
z = bathy1.Points(:,3);

% replace positive values in water by default negative value
ind = dist<0 & z>=0;
z(ind) = minwater;

 % replace the topography for a constant z
ind = dist>0;
z(ind) = 1*max(max(topo.z));
% remove elements on land (depending on value at center of element)
ind = distc>0;
X = [ bathy1.Points(:,1:2) z ];
%T = cleanFlatT( bathy.ConnectivityList(:,:), X, 1e-10 );
T = cleanFlatT( bathy1.ConnectivityList(~ind,:), X, 1e-10 );
bathy1 = triangulation( T, X(:,1), X(:,2), X(:,3) );
figure; trisurf( bathy1 );

% add vertical elements to make sure the STL crosses the z=0 surface
altz = double(1*max(max(topo.z)));
ind = abs(bathy1.Points(:,3))<1e-8;
bnd = freeBoundary(bathy1);
bnd = bnd(all(ind(bnd),2),:);
[bndnodes,~,indnodes] = unique(bnd);
indnodes = reshape(indnodes, size(bnd)) + size(bathy1.Points,1);
newnodes = [bathy1.Points(bndnodes,1:2) altz*ones(length(bndnodes),1)];
Nodes = [bathy1.Points; newnodes];
newelts = [ [bnd indnodes(:,1)];
    [bnd(:,2) indnodes(:,2) indnodes(:,1)] ];
newelts = cleanFlatT( newelts, Nodes, 1e-10 );
Elts = [bathy1.ConnectivityList; newelts];
bathy1 = triangulation( Elts, Nodes );
figure(30);hold all; trisurf( bathy1 );
%% bathy2 -> only the topography with a larger z
% it allow the material identification in HexMesh

% integrate bathymetry and coastlines
zz = double(1.5*max(max(topo.z)));
bathy2 = altimetryCoastline( bathy, water, 3, zz);
[dist,distc] = distanceCoastline( bathy2, water );
z = bathy2.Points(:,3);

% replace positive values in water by default negative value
ind = dist<0 & z>=0;
z(ind) = minwater;

 % replace the topography for a constant z
ind = dist>0;
z(ind) = 1.5*max(max(topo.z));
% remove elements on land (depending on value at center of element)
ind = distc>0;
X = [ bathy2.Points(:,1:2) z ];
%T = cleanFlatT( bathy.ConnectivityList(:,:), X, 1e-10 );
T = cleanFlatT( bathy2.ConnectivityList(ind,:), X, 1e-10 );
bathy2 = triangulation( T, X(:,1), X(:,2), X(:,3) );
figure; trisurf( bathy2 );

% add vertical elements to make sure the STL crosses the z=0 surface
altz = double(1*max(max(topo.z)));
ind = abs(bathy2.Points(:,3))<1e-8;
bnd = freeBoundary(bathy2);
bnd = bnd(all(ind(bnd),2),:);
[bndnodes,~,indnodes] = unique(bnd);
indnodes = reshape(indnodes, size(bnd)) + size(bathy2.Points,1);
newnodes = [bathy2.Points(bndnodes,1:2) altz*ones(length(bndnodes),1)];
Nodes = [bathy2.Points; newnodes];
newelts = [ [bnd indnodes(:,1)];
    [bnd(:,2) indnodes(:,2) indnodes(:,1)] ];
newelts = cleanFlatT( newelts, Nodes, 1e-10 );
Elts = [bathy2.ConnectivityList; newelts];
bathy2 = triangulation( Elts, Nodes );
figure(30);hold all; trisurf( bathy2 );

%bulid the bathy
bathy.Points = [bathy1.Points; bathy2.Points];
bathy.ConnectivityList = [bathy1.ConnectivityList; bathy2.ConnectivityList + size(bathy1.Points,1)];

%% integrate topography with coastline
topo = altimetryCoastline( topo, water ,4, 0 );
[ dist, distc ] = distanceCoastline( topo, water);
z = topo.Points(:,3);

% replace all values in water by zero
ind = dist<=0;
z(ind) = 0;

% at distance H around the coastline, put all negative z on land to zero
ind = dist>0 & dist<H & z<0;
z(ind) = 0;
X = [ topo.Points(:,1:2) z ];
T = cleanFlatT( topo.ConnectivityList, X, 1e-10 );
topo = triangulation( T, X(:,1), X(:,2), X(:,3) );
figure; trisurf( topo ); shading flat;

%% write stl files
% try
%     system('wget http://mooring.ucsd.edu/software/matlab/mfiles/toolbox/data/stl/write_stl.m');
% catch
%     error('please try to dowload write_stl function')
% end
% write topography STL file
if ~isempty(topo)
    [xtopo,ytopo] = lonlat2m(topo.Points(:,1),topo.Points(:,2));
    xtopo = xtopo - min(xtopo); ytopo = ytopo - min(ytopo);
    write_stl( fullfile(outdir,[output_name,'_topo.stl']), ...
        [xtopo ytopo topo.Points(:,3)], topo.ConnectivityList');
end

% write bathymetry STL file
if ~isempty(bathy)
    [xbathy,ybathy] = lonlat2m(bathy.Points(:,1),bathy.Points(:,2));
    xbathy = xbathy - min(xbathy); ybathy = ybathy - min(ybathy);
    write_stl( fullfile(outdir,[output_name,'_bathy.stl']), ...
        [xbathy ybathy bathy.Points(:,3)], bathy.ConnectivityList');
end
% convert and move to input folder
system(['dos2unix ',output_name,'_topo.stl']);
system(['dos2unix ',output_name,'_bathy.stl']);
system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
system(['mv ',output_name,'_topo.stl ../input/.']);
system(['mv ',output_name,'_bathy.stl ../input/.']);