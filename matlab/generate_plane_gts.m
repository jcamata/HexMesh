function generate_plane_gts(filename, xmin, xmax, ymin, ymax, z_level, nx, ny)
% GENERATE_PLANE_GTS Generates a planar GTS surface at a given z-elevation.
%
% Usage:
%   generate_plane_gts(filename, xmin, xmax, ymin, ymax, z_level, nx, ny)
%
% Example:
%   generate_plane_gts('../input/flat_layer_z50000.gts', -10000, 125000, -10000, 125000, 50000, 20, 20);

if nargin < 7, nx = 20; end
if nargin < 8, ny = 20; end

fprintf('Generating planar GTS surface: %s at z = %.1f m...\n', filename, z_level);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[X, Y] = meshgrid(x, y);

n_vertices = nx * ny;
vertices = zeros(n_vertices, 3);
for j = 1:ny
    for i = 1:nx
        vid = (j-1)*nx + i;
        vertices(vid, :) = [X(j, i), Y(j, i), z_level];
    end
end

% Build triangles (2 triangles per grid quad)
triangles = [];
for j = 1:(ny-1)
    for i = 1:(nx-1)
        v1 = (j-1)*nx + i;
        v2 = (j-1)*nx + (i+1);
        v3 = j*nx + (i+1);
        v4 = j*nx + i;
        
        triangles = [triangles; v1, v2, v3; v1, v3, v4];
    end
end

% Build unique edges
edge_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
edges = [];
tri_edges = zeros(size(triangles, 1), 3);

for t = 1:size(triangles, 1)
    tv = triangles(t, :);
    pairs = [tv(1), tv(2); tv(2), tv(3); tv(3), tv(1)];
    for p = 1:3
        va = min(pairs(p, 1), pairs(p, 2));
        vb = max(pairs(p, 1), pairs(p, 2));
        key = sprintf('%d_%d', va, vb);
        if isKey(edge_map, key)
            eid = edge_map(key);
        else
            eid = size(edges, 1) + 1;
            edges = [edges; va, vb];
            edge_map(key) = eid;
        end
        tri_edges(t, p) = eid;
    end
end

% Write GTS file
fid = fopen(filename, 'w');
if fid == -1
    error('Could not open file %s for writing.', filename);
end

fprintf(fid, '%d %d %d GtsSurface GtsFace GtsEdge GtsVertex\n', ...
    size(vertices, 1), size(edges, 1), size(triangles, 1));

for v = 1:size(vertices, 1)
    fprintf(fid, '%.6f %.6f %.6f\n', vertices(v, 1), vertices(v, 2), vertices(v, 3));
end

for e = 1:size(edges, 1)
    fprintf(fid, '%d %d\n', edges(e, 1), edges(e, 2));
end

for t = 1:size(triangles, 1)
    fprintf(fid, '%d %d %d\n', tri_edges(t, 1), tri_edges(t, 2), tri_edges(t, 3));
end

fclose(fid);
fprintf('Successfully generated %s (%d vertices, %d edges, %d triangles).\n', ...
    filename, size(vertices, 1), size(edges, 1), size(triangles, 1));

end
