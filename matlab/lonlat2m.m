function [x,y] = lonlat2m(x,y,ymoyen)
% LONLAT2M transforms longitude/latitude (in degrees) coordinates 
% to meters
%
% syntax: [x,y] = lonlat2m(lon,lat)

% Earth radius
R = 6371;
% distance corresponding to 1 degree of latitude
nm = pi*R/180*1000;
% consideration of average y rather local y (to obtain a square mesh)
if nargin>2
    ym = ymoyen;
else
    ym = y;
end

% change of coordinates
x = x.*cos(ym/(180*pi))*nm;
y = y*nm;
