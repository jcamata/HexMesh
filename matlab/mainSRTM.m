close all
clear all

% set bounds (in degrees and minutes)
% [ minlat° minlat' maxlat° maxlat' minlon° minlon' maxlon° maxlon']
% NOTE THAT YOU SHOULD CONSIDER A LARGER AREA THAN WHAT YOU INTEND TO MESH
% l = [40 00 52 00   -5 00    8 00]; % France
% l = [43 48 44 05    5 30    6 00]; % Cadarache, France
% l = [47 00 48 00   -4 00   -3 00]; % Belle-Ile, France
l = [38 00 39 00   20 00   21 00]; % Kefalonia, Greece
% l = [37 10 37 40  138 15  138 55]; % Kashiwazaki, Japan
% l = [18 30 21 00 -157 00 -154 00]; % Mauna Loa, Hawai
% l = [38 40 38 60   20 40   20 60]; % test small Kefalonia, Greece
% l = [20 00 21 00 -156 00 -155 00]; % test small Mauna Loa, Hawai

% choose output directory
outdir = '.';

% characteristic length over which details of the coastline are removed
% put H<=0 if you want no smoothing (this is very expensive because many
% small elements will be created)
H = 0.1; % in units of lon/lat

% minimum water depth
minwater = -10;

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

% transform coordinates
latbnds = l(1):(l(3)+1);
latcrop = [l(1)+l(2)/60 l(3)+l(4)/60];
lonbnds = l(5):(l(7)+1);
loncrop = [l(5)+l(6)/60 l(7)+l(8)/60];

% get topography and plot
topographySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
                     'crop', [latcrop loncrop]);
outtopo = topographySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
                     'crop', [latcrop loncrop]);

% get bathymetry and plot
bathymetrySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
                     'crop', [latcrop loncrop]);
outbathy = bathymetrySRTM(latbnds, lonbnds, outdir, 'interp', 'merge', ...
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

% integrate bathymetry and coastlines
[ trib, z, dist, distc, xc ] = altimetryCoastline( outbathy, water, 3 );

% replace positive values in water by default negative value
ind = dist<0 & z>=0;
z(ind) = minwater;

% remove nodes on land (depending on value at center of element) and nodes
% that are outside the bounding box
ind = (xc(:,1)<loncrop(1)) | (xc(:,1)>loncrop(2)) | ...
      (xc(:,2)<latcrop(1)) | (xc(:,2)>latcrop(2));
ind = ind | (distc>0);
if any(~ind)
    trib = triangulation( trib.ConnectivityList(~ind,:), ...
                                   trib.Points(:,1), trib.Points(:,2), z );

% add vertical elements to make sure the STL crosses the z=0 surface
    altz = 1000;
    ind = abs(trib.Points(:,3))<1e-8;
    bnd = freeBoundary(trib);
    bnd = bnd(all(ind(bnd),2),:);
    [bndnodes,~,indnodes] = unique(bnd);
    indnodes = reshape(indnodes, size(bnd)) + size(trib.Points,1);
    newnodes = [trib.Points(bndnodes,1:2) altz*ones(length(bndnodes),1)];
    newelts1 = [bnd indnodes(:,1)];
    newelts2 = [bnd(:,2) indnodes(:,2) indnodes(:,1)];
    Elts = [trib.ConnectivityList; newelts1; newelts2];
    Nodes = [trib.Points; newnodes];
    clear trib
    trib = triangulation( Elts, Nodes );
    figure; trisurf( trib );
else
    trib = [];
end

% integrate topography with coastline
[ trit, z, dist, distc, xc ] = altimetryCoastline( outtopo, water );

% replace all values in water by zero
ind = dist<=0;
z(ind) = 0;

% at distance H around the coastline, put all negative z on land to zero
ind = dist>0 & dist<H & z<0;
z(ind) = 0;

% remove nodes that are outside the bounding box
ind = (xc(:,1)<loncrop(1)) | (xc(:,1)>loncrop(2)) | ...
      (xc(:,2)<latcrop(1)) | (xc(:,2)>latcrop(2));
if any(~ind)
    trit = triangulation( trit.ConnectivityList(~ind,:), ...
                                   trit.Points(:,1), trit.Points(:,2), z );
%    figure; trisurf( trit ); shading flat;
else
    trit = [];
end

return

% write topography STL file
nmax = 1e4;
if length(outtopo.lon)<nmax && length(outtopo.lat)<nmax
    if ~isempty(trit)
        [xtopo,ytopo] = lonlat2m(trit.Points(:,1),trit.Points(:,2));
        write_stl( fullfile(outdir,'topo.stl'), ...
                  [xtopo ytopo trit.Points(:,3)], trit.ConnectivityList');
    end
else
    error('the STL file is too large')
end

% write bathymetry STL file
if ~isempty(trib)
    [xbathy,ybathy] = lonlat2m(trib.Points(:,1),trib.Points(:,2));
    write_stl( fullfile(outdir,'bathy.stl'), ...
                 [xbathy ybathy trib.Points(:,3)], trib.ConnectivityList');
end            
             