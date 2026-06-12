% Define a flat plane (1x1 in x-y) at a given z-height
function generatePlaneSTL(filename, z_height)
    % Vertices of a 1x1 square in x-y plane, at height z
    vertices = [
        0 0 z_height;  % Vertex 1
        1 0 z_height;  % Vertex 2
        1 1 z_height;  % Vertex 3
        0 1 z_height;  % Vertex 4
    ];
    
    % Faces (two triangles forming a square)
    faces = [
        1 2 3;  % First triangle
        1 3 4;  % Second triangle
    ];
    
    % Write to STL file
    stlwrite(filename, faces, vertices);
    fprintf('Generated: %s (z = %.1f)\n', filename, z_height);
end

% Generate the two STL files
generatePlaneSTL('plane_z10.stl', 10);  % Plane at z = 10
generatePlaneSTL('plane_z5.stl', 5);    % Plane at z = 5

function stlwrite(filename, faces, vertices)
    fid = fopen(filename, 'w');
    fprintf(fid, 'solid plane\n');
    for i = 1:size(faces, 1)
        tri = vertices(faces(i,:), :);
        fprintf(fid, 'facet normal 0 0 1\n');
        fprintf(fid, '  outer loop\n');
        fprintf(fid, '    vertex %.6f %.6f %.6f\n', tri(1,:));
        fprintf(fid, '    vertex %.6f %.6f %.6f\n', tri(2,:));
        fprintf(fid, '    vertex %.6f %.6f %.6f\n', tri(3,:));
        fprintf(fid, '  endloop\n');
        fprintf(fid, 'endfacet\n');
    end
    fprintf(fid, 'endsolid plane\n');
    fclose(fid);
end