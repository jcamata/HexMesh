function [dt,z,dist,distc,xc] = altimetryCoastline( Xb, Cl, lp )
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
% z0: altitude at the points of the delaunaytriangulation
% dist: distance of each node of dt to the closest coastline (positive on 
%       land
% distc: same as dist, for the center of the elements in dt
% xc: position of the center of the elements in tri

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

% add nodes based on Coastlines
[lon,lat,z,cx,cy] = integrateCoastlines(Cl,lon,lat,z,cx,cy,lp);

% construct delaunay triangulation
dt = delaunayTriangulation(lon,lat,[cx cy]);

% construct initial map of altitudes
z = scatteredInterpolant( lon, lat, z );
z = z(dt.Points(:,1),dt.Points(:,2));

% compute position of element centers
xc = incenter(dt);
Nc = size(xc,1);

% construct function distance to coastline (positive on land) for all the
% nodes, and also the center of the elements
distOcean = makeLS( [[dt.Points(:,1) dt.Points(:,2)]; xc], Cl.Ocean );
distLand = -makeLS( [[dt.Points(:,1) dt.Points(:,2)]; xc], Cl.Land );
dist = min(abs([distOcean distLand]),[],2);
ind = dist==abs(distOcean);
dist(ind) = distOcean(ind);
dist(~ind) = distLand(~ind);
distc = dist(end-Nc+1:end);
dist = dist(1:end-Nc);

%==================================
% FUNCTION INTEGRATECOASTLINES
function [lon,lat,z,cx,cy] = integrateCoastlines(Cl,lon,lat,z,cx,cy,lp)
if ~isempty(Cl)
    for i1 = 1:length(Cl.Ocean)
        [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl.Ocean{i1},lon,lat,z,cx,cy,lp);
    end    
    for i1 = 1:length(Cl.Land)
        [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl.Land{i1},lon,lat,z,cx,cy,lp);
    end    
else
    disp('no water body information found')
end

%==================================
% FUNCTION INTEGRATESINGLECOASTLINE
function [lon,lat,z,cx,cy] = integrateSingleCoastline(Cl,lon,lat,z,cx,cy,lp)

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
z = [z(:); zeros(size(lon1))];
cx = [cx(:); Nl+cx1 ];
cy = [cy(:); Nl+cy1 ];
%==================================
% FUNCTION MAKELS to construct the level-set distance function
function LS = makeLS( X, polis )
% LS = makeLS( X, polis )
%
% Esta funcion arma un campo level set poniendo en cada nodo 
% la distancia con signo a las rectas descritan en polis.
%
% INPUT
%	X           Posiciones nodales donde se tiene que calcular distancias.
%	polis{i}    cell de poligonos
%
% OUTPUT
%  LS          level set nodal

% R. Cottereau 04/2008

N = size(X,1);
LS = Inf*ones(N,1);

for i1 = 1:length(polis)
    LS = min( [LS DistanceToPolygon( X, polis{i1} )], [], 2 );
end

for i1 = 1:length(polis)
    ind = inpolygon( X(:,1), X(:,2), polis{i1}(1,:), polis{i1}(2,:) );
    LS( ind ) = -LS( ind );
end

LS = real( LS );

%==========================================================================
% Distance from all the points in X to polygon poli1
function dist1 = DistanceToPolygon( X, poli1 )

dist1 = Inf * ones(size(X,1),1);
for i1=1:size(poli1,2)-1
    if (poli1(1,i1)==poli1(1,i1+1))&&(poli1(2,i1)==poli1(2,i1+1))
        dist2 = DistanceToPoint( X, poli1(:,i1) );
    else
        dist2 = DistanceToSegment( X, poli1(:,i1), poli1(:,i1+1) );
    end
    dist1 = min( [dist1 dist2], [], 2 );
end

%==========================================================================
% Distance from all the points in X to segment defined by X1 and X2
% depending on the scalar product s=(X-X1,X2-X1), the point X is either
% closest to X2 (s>1), to X1 (s<0), or to the segment itself (0<s<1)
function dist1 = DistanceToSegment( X, X1, X2)

v1 = [X(:,1)-X1(1) X(:,2)-X1(2)];
%s1 = v1.^2;
v2 = X2'-X1';
s2 = v2*v2';
s = (v1*v2') / s2;
ind0 = find( (s>0) & (s<1) );
ind1 = find( s >= 1 );
ind2 = find( s <= 0 );

dist1 = zeros(size(X,1),1);
dist1( ind0 ) = sqrt( v1(ind0,1).^2 + v1(ind0,2).^2 - s(ind0).^2*s2 );
dist1( ind1 ) = sqrt( (X(ind1,1)-X2(1)).^2 + (X(ind1,2)-X2(2)).^2 );
dist1( ind2 ) = sqrt( v1(ind2,1).^2 + v1(ind2,2).^2 );

%==========================================================================
% Distance from all the points in X to point defined by X1 
function dist1 = DistanceToPoint( X, X1)
x = X(:,1)-X1(1) ;
y = X(:,2)-X1(2);
dist1 = sqrt( x.^2 + y.^2 );
