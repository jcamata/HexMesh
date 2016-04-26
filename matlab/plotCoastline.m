function h = plotCoastline( c, type, h )
% PLOTCOASTLINE to plot the coastline
%
% syntax h = plotCoastline( h, c )
%
% c: coastline
% type: plotting options
% h: handle to the figure where the coastline should be plotted

% initialization
if nargin<3
    h = figure;
else
    h = figure(h);
end

% loop on coastlines
for i1 = 1:length(c)    
    hold on; plot(c{i1}(1,:),c{i1}(2,:),type);
end