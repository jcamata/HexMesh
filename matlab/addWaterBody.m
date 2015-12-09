function [lon,lat,z] = addWaterBody(water,lon,lat,z,fz,depth,loncrop,latcrop)
for i1 = 1:length(water)
    % create nodes for the water body borders
    lont = water{i1}(1,:)';
    latt = water{i1}(2,:)';
    % create nodes at the center of the water triangles
    tri = delaunay(lont,latt);
    lonc = mean(lont(tri),2);
    latc = mean(latt(tri),2);
    % create final structure
    lon = [lon(:); lont; lonc];
    lat = [lat(:); latt; latc];
    z = [z(:); fz(lont,latt); fz(lonc,latc)-depth];
    % cleanup
    ind = lon<=loncrop(2) & lon>=loncrop(1) & lat<=latcrop(2) & lat>=latcrop(1);
    lon = lon(ind);
    lat = lat(ind);
    z = z(ind);
end

