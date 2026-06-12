% check_conformity.m
% Lê os arquivos de debug do pillow e reporta:
%   - pillow_postpass    : quantos nós foram corrigidos pelo post-pass
%   - pillow_residual    : quantos nós ficaram sem pillow do lado oposto
%   - pillow_conformity  : quantas faces não-conformes restaram

mpi_size = 1;
mpi_rank = 0;
tag = sprintf('%04d_%04d', mpi_size, mpi_rank);

%% --- Post-pass ---
pp_file = sprintf('pillow_postpass_%s.txt', tag);
fid = fopen(pp_file, 'r');
if fid == -1
    fprintf('[postpass] Arquivo não encontrado: %s\n', pp_file);
else
    % Última linha contém o total
    lines = {};
    while ~feof(fid)
        l = fgetl(fid);
        if ischar(l), lines{end+1} = l; end %#ok<AGROW>
    end
    fclose(fid);
    if ~isempty(lines)
        total_line = lines{end};
        n = sscanf(total_line, 'PostPass total: %d');
        fprintf('[postpass] Nós corrigidos pelo post-pass: %d\n', n);
        if n > 0
            fprintf('  (ver %s para detalhes por elemento)\n', pp_file);
        end
    end
end

%% --- Residual ---
res_file = sprintf('pillow_residual_%s.txt', tag);
fid = fopen(res_file, 'r');
if fid == -1
    fprintf('[residual] Arquivo não encontrado: %s\n', res_file);
else
    lines = {};
    while ~feof(fid)
        l = fgetl(fid);
        if ischar(l), lines{end+1} = l; end %#ok<AGROW>
    end
    fclose(fid);
    if ~isempty(lines)
        total_line = lines{end};
        n = sscanf(total_line, 'Residual: %d');
        if n == 0
            fprintf('[residual] OK — sem nós com lado de pillow faltando\n');
        else
            fprintf('[residual] AVISO: %d nó(s) sem pillow do lado oposto\n', n);
            fprintf('  (ver %s para detalhes)\n', res_file);
        end
    end
end

%% --- Conformidade geométrica ---
conf_file = sprintf('pillow_conformity_%s.txt', tag);
fid = fopen(conf_file, 'r');
if fid == -1
    fprintf('[conformity] Arquivo não encontrado: %s\n', conf_file);
else
    % Última linha: "Conformity check: N non-conforming face(s)"
    lines = {};
    face_a = [];
    face_b = [];
    while ~feof(fid)
        l = fgetl(fid);
        if ~ischar(l), continue; end
        lines{end+1} = l; %#ok<AGROW>

        % Parsear pares de faces não-conformes para plotagem
        t = sscanf(l, 'face iel_a=%d face_a=%d nodes_a=%d %d %d %d');
        if numel(t) == 6
            face_a(end+1,:) = t'; %#ok<AGROW>
        end
        t = sscanf(l, '     iel_b=%d face_b=%d nodes_b=%d %d %d %d');
        if numel(t) == 6
            face_b(end+1,:) = t'; %#ok<AGROW>
        end
    end
    fclose(fid);

    if ~isempty(lines)
        total_line = lines{end};
        n = sscanf(total_line, 'Conformity check: %d');
        if n == 0
            fprintf('[conformity] OK — malha conforme (0 faces não-conformes)\n');
        else
            fprintf('[conformity] ERRO: %d face(s) não-conforme(s)\n', n);
            fprintf('  Elementos envolvidos (iel_a):\n');
            for k = 1:size(face_a,1)
                fprintf('    iel_a=%d face=%d  <->  iel_b=%d face=%d\n', ...
                    face_a(k,1), face_a(k,2), face_b(k,1), face_b(k,2));
            end
            fprintf('  (ver %s para nós completos)\n', conf_file);
        end
    end
end

fprintf('\nResumido:\n');
fprintf('  Use mpi_size=%d mpi_rank=%d — altere no topo do script se necessário.\n', mpi_size, mpi_rank);
