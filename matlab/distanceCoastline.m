function [dist,distc] = distanceCoastline(dt,Cl)
% DISTANCECOASTLINE to compute the distance from the elements in a
% triangulation to a coastline
%
% syntax: [dist,distc] = distanceCoastline(dt)
%
% dt: delaunay triangulation of the mesh 
% Cl structured array of information on Coastlines
%
% dist: distance of each node of dt to the closest coastline (positive on 
%       land
% distc: same as dist, for the center of the elements in dt

% compute position of element centers
xc = incenter(dt);
Nc = size(xc,1);
  
% construct function distance to coastline (positive on land) for all the
% nodes, and also the center of the elements
distOcean = makeLS( [[dt.Points(:,1) dt.Points(:,2)]; xc(:,1:2)], Cl.Ocean );
distLand = -makeLS( [[dt.Points(:,1) dt.Points(:,2)]; xc(:,1:2)], Cl.Land );
dist = min(abs([distOcean distLand]),[],2);
ind = dist==abs(distOcean);
dist(ind) = distOcean(ind);
dist(~ind) = distLand(~ind);
distc = dist(end-Nc+1:end);
dist = dist(1:end-Nc);

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
