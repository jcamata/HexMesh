function data = read_pillow_file(filename)
% Lê o arquivo pillow com a nova linha de edges (com \n no final)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Não consegui abrir: %s', filename);
    end

    data = struct('octree_id', {}, 'edges', {}, 'pillows', {});
    idx = 1;

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(strtrim(line)), continue; end

        % Cabeçalho
        tokens = sscanf(line, 'Tenho %d verticies no pillow do octree %d');
        if length(tokens) ~= 2, continue; end

        num_pillows = tokens(1);
        octree_id   = tokens(2);

        % ==================== LISTA DE EDGES ====================
        fgetl(fid);                          % pula "Lista de edges"
        edges_line = fgetl(fid);             % ex: "0 3 5 7 11"  ou "" se nenhum
        fgetl(fid);                          % pula "Fim Lista de edges"

        edges = sscanf(edges_line, '%d')'; 

        % ==================== PILLOWS ====================
        %pillows = struct('id',cell(num_pillows,1),'pa',{},'pb',{},'x',{},'y',{},'z',{},'elements',{});
        pillows = struct('id', cell(num_pillows,1), ...
                         'pa', cell(num_pillows,1), ...
                         'pb', cell(num_pillows,1), ...
                         'x',  cell(num_pillows,1), ...
                         'y',  cell(num_pillows,1), ...
                         'z',  cell(num_pillows,1), ...
                         'elements', cell(num_pillows,1));
        num_pillows = tokens(1);
        octree_id   = tokens(2);

        pillows = struct('id', cell(num_pillows,1), ...
                         'pa', cell(num_pillows,1), ...
                         'pb', cell(num_pillows,1), ...
                         'x',  cell(num_pillows,1), ...
                         'y',  cell(num_pillows,1), ...
                         'z',  cell(num_pillows,1), ...
                         'elements', cell(num_pillows,1));

        for i = 1:num_pillows
            % Pillow id
            line = fgetl(fid);
            pillows(i).id = sscanf(line, 'Pillow id:%d');

            % pa pb
            line = fgetl(fid);
            t = sscanf(line, 'pa:%d pb:%d');
            pillows(i).pa = t(1);
            pillows(i).pb = t(2);

            % x y z
            line = fgetl(fid);
            t = sscanf(line, 'x:%d y:%d z:%d');
            pillows(i).x = t(1);
            pillows(i).y = t(2);
            pillows(i).z = t(3);

            % Número de elementos
            line = fgetl(fid);
            num_elem = sscanf(line, 'Tenho %d elementos');

            pillows(i).elements = struct('id', cell(num_elem,1), 'faces', cell(num_elem,1));

            for j = 1:num_elem
                % Element id
                line = fgetl(fid);
                pillows(i).elements(j).id = sscanf(line, 'Element:%d');

                % Número de faces
                line = fgetl(fid);
                num_faces = sscanf(line, 'Tenho %d faces');

                faces = zeros(num_faces, 1);
                for k = 1:num_faces
                    line = fgetl(fid);
                    faces(k) = sscanf(line, ' face %d');
                end
                pillows(i).elements(j).faces = faces;
            end
        end


        % Salva este octree
        data(idx).octree_id = octree_id;
        data(idx).edges     = edges;
        data(idx).pillows   = pillows;
        idx = idx + 1;
    end

    fclose(fid);
    fprintf('Pronto! %d octree(s) lidos com sucesso.\n', length(data));
end