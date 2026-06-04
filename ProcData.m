ccc

%d = h5read("/Users/lac/Documents/Git/HexMesh/SEM/traces/capteurs.0001.h5",'/UU_0000');
d = h5read("/Users/lac/Documents/Git/HexMesh/ARGOSTOLI/traces/capteurs.0000.h5",'/UU_0000');

pos = h5read("/Users/lac/Documents/Git/HexMesh/SEM/traces/capteurs.0001.h5",'/UU_0000_pos');

var = h5read("/Users/lac/Documents/Git/HexMesh/SEM/traces/capteurs.0001.h5",'/Variables');

%%
 
figure
plot(d(1,:),d(2,:))
hold on
plot(d(1,:),d(3,:))
plot(d(1,:),d(4,:))

 
figure
plot(d(1,:),d(5,:))
hold on
plot(d(1,:),d(6,:))
plot(d(1,:),d(7,:))

figure
plot(d(1,:),d(8,:))
