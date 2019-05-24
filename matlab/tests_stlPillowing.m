close all
clear variables
clc

teste1 = 0; %z
teste2 = 0; %x
teste3 = 0; %y

% cortes tipo "2"
teste4 = 0;
teste5 = 0;
teste6 = 0;
teste7 = 0;

%cortes tipo "3"
teste8  = 1;
teste9  = 0;
teste10 = 0;
teste11 = 0;
teste12  = 0;
teste13 = 0;
teste14 = 0;
teste15 = 0;

%cortes tipo "4"
teste16  = 0;
teste17 = 0;
teste18 = 0;
teste19 = 0;
teste20  = 0;
teste21  = 0;
teste22 = 0;
teste23 = 0;

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
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
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
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
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
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
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
    x = [0 1]*5;
    y = [0 1]*5;
    z = [-5 -5; 1.1 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
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
    z = [-5 1.1; -5 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
end
%% test 6
if(teste6)
    output_name = 'teste6';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 1]*5;
    y = [0.9 0]*5;
    z = [-5 -5; 1.1 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
end
%% test 7
if(teste7)
    output_name = 'teste7';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [1 0]*5;
    y = [0 1]*5;
    z = [-5 1.1; -5 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
end
%% test 8
if(teste8)
    output_name = 'teste8';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [0 1]*5;
    y = [-0.5 0.5]*5;
    z = [-7 -5; -3 1.1]*1.5+8;
    write_stl([output_name,'_bathy.stl'],x,y,z)

    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
end
%% test 9
if(teste9)
    output_name = 'teste9';
    
    %topo
    x = [0 1]*5;
    y = [0 1]*5;
    z = [2 2; 3 3];
    write_stl([output_name,'_topo.stl'],x,y,z)
    
    %bathy
    x = [1 0]*5;
    y = [0 1]*5;
    z = [-5 1.1; -5 1.1]*1.5;
    write_stl([output_name,'_bathy.stl'],x,y,z)
    
    system(['dos2unix ',output_name,'_topo.stl']);
    system(['dos2unix ',output_name,'_bathy.stl']);
    system(['stl2gts <',output_name,'_topo.stl > ../input/',output_name,'_topo.gts']);
    system(['stl2gts <',output_name,'_bathy.stl > ../input/',output_name,'_bathy.gts']);
    system(['mv ',output_name,'_topo.stl ../input/.']);
    system(['mv ',output_name,'_bathy.stl ../input/.']);
end
%%





return
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