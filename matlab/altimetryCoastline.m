function dt = altimetryCoastline( Xb, Cl, lp ,zz)
% BATHYMETRYCOASTLINE to integrate the information from bathymetry and 
% coastlines
% 
% syntax: [X,z0] = bathyTopoCoastline( Xb, Cl )
%
% Xb structured array of information on bathymetry
% Cl structured array of information on Coastlines
% lp handle to a figure for plotting the resulting constraints
%
% dt: delaunay triangulation of the mesh 

% initialization
if nargin<3
    lp = 0;
end
if lp>0
    figure(lp);
end
    
% mesh based on bathymetry only
[lon,lat] = ndgrid(Xb.lon(:),Xb.lat(:));
lon = lon(:);
lat = lat(:);
z = double(Xb.z'); 
z = z(:);
cx = zeros(0,1);
cy = zeros(0,1);

% bounding box
box = [min(lon) max(lon) min(lat) max(lat)];

% add nodes based on Coastlines
[lon,lat,z,cx,cy] = integrateCoastlines(Cl,lon,lat,z,cx,cy,lp,zz);

% construct delaunay triangulation
dt = delaunayTriangulation(lon,lat,[cx cy]);

% remove elements that are completely outside the bounding box
indin = (dt.Points(:,1)>=box(1)) & (dt.Points(:,1)<=box(2)) & ...
        (dt.Points(:,2)>=box(3)) & (dt.Points(:,2)<=box(4));
indin = indin(dt.ConnectivityList);
indon = (dt.Points(:,1)==box(1)) | (dt.Points(:,1)==box(2)) | ...
        (dt.Points(:,2)==box(3)) | (dt.Points(:,2)==box(4));
indon = indon(dt.ConnectivityList);
ind = any( indin & ~indon, 2) | all( indin, 2 );

% construct initial map of altitudes
z = scatteredInterpolant( lon, lat, z );
z = z(dt.Points(:,1),dt.Points(:,2));

% create final triangulation
X = [ dt.Points(:,1:2) z ];
T = cleanFlatT( dt.ConnectivityList(ind,:), X, 1e-10 );
dt = triangulation( T, X(:,1), X(:,2), X(:,3) );

%==================================
% FUNCTION INTEGRATECOASTLINES
function [lon,lat,z,cx,cy] = integrateCoastlines(Cl,lon,lat,z,cx,cy,lp,zz)
if ~isempty(Cl)
    for i1 = 1:length(Cl.Ocean)
        [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl.Ocean{i1},lon,lat,z,cx,cy,lp,zz);
    end    
    for i1 = 1:length(Cl.Land)
        [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl.Land{i1},lon,lat,z,cx,cy,lp,zz);
    end    
else
    disp('no water body information found')
end

%==================================
% FUNCTION INTEGRATESINGLECOASTLINE
function [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl,lon,lat,z,cx,cy,lp,zz)

% initialization
lon1 = Cl(1,1:end-1)';
lat1 = Cl(2,1:end-1)';

% pre-existing bounding box (erase anything outside of that)
mx = min(lon);
Mx = max(lon);
my = min(lat);
My = max(lat);
ind = (lon1>=mx) & (lon1<=Mx) & (lat1>=my) & (lat1<=My);

% add one node on each cut end to keep the segments complete
ind = ind | (ind([end 1:end-1]) & ~ind([2:end 1])) | ...
           (~ind([end 1:end-1]) &  ind([2:end 1]));

% bounding box for the current coastline (remove nodes that are between two
% elements on the perimeter)
mx1 = min(lon1);
Mx1 = max(lon1);
my1 = min(lat1);
My1 = max(lat1);
gerr=1e-5;
Ndp = abs(lon1-mx1)<gerr | abs(lon1-Mx1)<gerr | abs(lat1-my1)<gerr | abs(lat1-My1)<gerr;
Eltp = Ndp & Ndp([2:end 1]); % elements along the perimeter
rmNd = Eltp([end 1:end-1]) & Eltp;
ind = ind & ~rmNd;

% use diff to separate between different segments of the coastline
indd = find(diff(ind)==1)+1;
if ind(1); indd = [1;indd]; end
indf = find(diff(ind)==-1);
if ind(end); indf = [indf;length(ind)]; end

% create the segment to constrain the triangulation
cx1 = zeros(0,1);
cy1 = zeros(0,1);
for i1 = 1:length(indd)
    cx1 = [cx1(:); (indd(i1):(indf(i1)-1))'];
    cy1 = [cy1(:); ((indd(i1)+1):indf(i1))'];
end

% renumbering to account for nodes that will not be kept
for i1=length(lon1):-1:1
    if ~ind(i1)
        cx1(cx1>i1) = cx1(cx1>i1)-1;
        cy1(cy1>i1) = cy1(cy1>i1)-1;
    end
end
lon1 = lon1(ind);
lat1 = lat1(ind);

% plot constraint network
if lp>0
    for i1 = 1:length(cx1)
        hold on; plot( lon1([cx1(i1) cy1(i1)]), lat1([cx1(i1) cy1(i1)]), 'g--x' );
    end
end

% integrate with previous values
Nl = length(lon(:));
lon = [lon(:); lon1];
lat = [lat(:); lat1];
z = [z(:); double(zz)*ones(size(lon1))];
cx = [cx(:); Nl+cx1 ];
cy = [cy(:); Nl+cy1 ];
