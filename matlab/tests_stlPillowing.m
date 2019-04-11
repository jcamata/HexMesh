close all
clear variables
clc

teste1 = 0;
teste2 = 1;
teste3 = 0;
teste4 = 0;
teste5 = 0;

% set the path for wget function
wget_path = '/opt/local/bin';
% set the path for dos2unix function
dos2unix_path = '/opt/local/bin';
% set the path for stl2gts function
stl2gts_path = '/opt/local/bin';
setenv('PATH', [getenv('PATH') ':',wget_path,':',dos2unix_path,':',stl2gts_path]);

%% test 1
if(teste1)
    output_name = 'teste1';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 2 2];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 1]*5;
    y = [0 1]*5;
    z = [1 1; 1 1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
end


%% test 2
if(teste2)
    output_name = 'teste2';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 2 2];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [1 1]*1;
    y = [0 1]*5;
    z = [-1 1; -1 1]*10;
    write_stl([output_name,'_bathy.stl'],x,y,z)
end
%% test 3
if(teste3)
    output_name = 'teste3';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 2 2];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 1]*5;
    y = [1 1]*1;
    z = [-1 -1; 1 1]*10;
    write_stl([output_name,'_bathy.stl'],x,y,z)
end
%% test 4
if(teste4)
    output_name = 'teste4';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 0.25 0.5 1]*5;
    y = [0.2 0.21 0.22 1]*5;
    z = [-6 -6 -6 -6; 1 1 1 1; 1 1 1 1; 1.1 1.1 1.1 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
end
%% test 5
if(teste5)
    output_name = 'teste5';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 1]*5;
    y = [0 1]*5;
    z = [-5 -5; 1.1 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
end


%%
system(['dos2unix ',output_name,'_topo.stl']);
system(['dos2unix ',output_name,'_bathy.stl']);
system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
system(['mv ',output_name,'_topo.stl ../input/.']);
system(['mv ',output_name,'_bathy.stl ../input/.']);