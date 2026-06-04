ccc

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