ccc

fid = fopen("Octree_0001_0000.txt");
Oct.mat = [];
Oct.ele = [];
Oct.id = [];
while ~feof(fid)
    line = fgetl(fid);
    aux = sscanf(line,'Octree %d element %d mat %d');
    Oct(aux(1)+1).id = aux(1);
    Oct(aux(1)+1).mat = [Oct(aux(1)+1).mat aux(3)];
    Oct(aux(1)+1).ele = [Oct(aux(1)+1).ele aux(2)];
end
fclose(fid);
%%
id2f = 3428

for ioc = 1:numel(Oct)
if(any(Oct(ioc).ele == id2f))
    Oct(ioc).ele
end
end

