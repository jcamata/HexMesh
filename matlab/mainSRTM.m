close all
clear all

% set bounds (in degrees and minutes)
% [ minlat° minlat' maxlat° maxlat' minlon° minlon' maxlon° maxlon']
% l = [40 00 52 00  -5 00   8 00]; % France
l = [43 48 44 05   4 30   5 00]; % Cadarache, France
% l = [38 00 38 35  20 00  21 00]; % Kefalonia, Greece
% l = [37 10 37 40 138 15 138 55]; % Kashiwazaki, Japan

% choose output directory
outdir = '.';

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

% get coastlines and rivers
water = swbd_shore( [ loncrop latcrop ], outdir );

% interpolate coastlines
[lon,lat] = ndgrid(outbathy.lon(:),outbathy.lat(:));
z = double(outbathy.z');
if ~isempty(water)
    for i1 = 1:length(water.Ocean)
        lon1 = water.Ocean{i1}(1,1:end-1)';
        lat1 = water.Ocean{i1}(2,1:end-1)';
        ind = lon1>=loncrop(1) & lon1<=loncrop(2) ...
            & lat1>=latcrop(1) & lat1<=latcrop(2);
        lon = [lon(:); water.Ocean{i1}(1,ind)'];
        lat = [lat(:); water.Ocean{i1}(2,ind)'];
        z = [z(:); zeros(sum(ind),1)];
    end
else
    disp('no water body information found')
end
tri = delaunay(lon,lat);

% transform heights in bathymetry to +9999
in = false(size(lon));
for i1 = 1:length(water.Ocean)
    in = in | inpolygon(lon,lat,water.Ocean{i1}(1,:),water.Ocean{i1}(2,:));
end
z(~in) = 9999;
% figure; trimesh(tri,lon,lat,z)


% transform angles to meters
nmax = 1e4;
if length(outtopo.lon)<nmax && length(outtopo.lat)<nmax
    [xtopo,ytopo] = ndgrid(outtopo.lon,outtopo.lat);
    [xtopo,ytopo] = lonlat2m(xtopo,ytopo);
    [xbathy,ybathy] = lonlat2m(lon,lat);
else
    error('the STL file is too large')
end

% write STL files
write_stl( fullfile(outdir,'topo.stl'), xtopo, ytopo, double(outtopo.z') );
write_stl( fullfile(outdir,'bathy.stl'), [xbathy ybathy z], tri' );








