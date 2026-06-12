ccc

% ============================================================
% DIAGNÓSTICO AUTORITATIVO
% ============================================================
% premesh_conformity: conformidade ANTES do pillow (após MovingNodes)
% pillow_conformity : conformidade APÓS o pillow + post-pass
%
% Se premesh > 0 → bug em MovingNodes, não no pillow.
% Se premesh == 0 e pillow > 0 → bug no próprio Pillowing.
% Se ambos == 0 → pillow OK; não-conformidade visível vem de
%                outra etapa (PML, otimizador, elementos de transição).

tag = '0001_0000';

checks = { ...
    sprintf('premesh_conformity_%s.txt', tag), '[PRE-PILLOW ] conformidade', 'Conformity check: %d'; ...
    sprintf('pillow_conformity_%s.txt',  tag), '[POS-PILLOW ] conformidade', 'Conformity check: %d'; ...
    sprintf('node_consistency_%s.txt',   tag), '[CONSISTÊNCIA] nós',         'Node consistency: %d'; ...
};

for ii = 1:size(checks,1)
    fname = checks{ii,1};
    label = checks{ii,2};
    fmt   = checks{ii,3};
    fid = fopen(fname, 'r');
    if fid == -1
        fprintf('%s  arquivo não encontrado: %s\n', label, fname);
    else
        txt = {};
        while ~feof(fid), l = fgetl(fid); if ischar(l), txt{end+1} = l; end; end
        fclose(fid);
        n = sscanf(txt{end}, fmt);
        if isempty(n), n = -1; end
        if n == 0
            fprintf('%s  OK — 0\n', label);
        else
            fprintf('%s  PROBLEMA: %d\n', label, n);
        end
    end
end
fprintf('\n');
fprintf('Diagnóstico: se PRE-PILLOW ou CONSISTÊNCIA > 0 → bug em MovingNodes.\n');
fprintf('             se apenas POS-PILLOW > 0           → bug no Pillowing.\n');
fprintf('             se todos 0 → pillow OK; gap visível vem de outra etapa.\n\n');
fprintf('(os avisos abaixo são estado por-octree ANTES do post-pass — apenas informativos)\n\n');

% ============================================================
% Exemplo de uso
data = read_pillow_file('pillow_0001_0000.txt');

%%
count = 1;
for id = 1:numel(data)
    flaaag = 0;
    for ip = 1:numel(data(id).pillows)
        if (data(id).pillows(ip).pa == -1 || data(id).pillows(ip).pb == -1)
            fprintf("Octree %d, pillow %d\n",id-1,ip-1);
            fprintf("\n");
            display(data(id).edges);
            fprintf("\n");
            defe(count,:) = [id;numel(data(id).edges)];
            flaaag = 1;
            %keyboard
        end
    end
    if flaaag
        count = count+1;
    end
    ned(id) = numel(data(id).edges);
end
plot(ned)
hold on
plot(defe(:,1),defe(:,2),'or')
return
% Ver o que tem no octree 0
oct = data(1);                     % primeiro octree
p   = oct.pillows(1);              % primeiro pillow

fprintf('Octree %d tem %d pillows\n', oct.octree_id, length(oct.pillows));
fprintf('Pillow id=%d  pa=%d pb=%d  pos=(%d,%d,%d)\n', ...
    p.id, p.pa, p.pb, p.x, p.y, p.z);

% Ver os elementos e faces do primeiro pillow
for el = 1:length(p.elements)
    fprintf('  Elemento %d tem %d faces: ', p.elements(el).id, length(p.elements(el).faces));
    disp(p.elements(el).faces');
end