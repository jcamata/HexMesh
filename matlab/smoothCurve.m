function c = smoothCurve(c,h)
% SMOOTHCURVE to smoothen a curve on a caracteristic length h
%
% syntax: c = smoothCurve(c,h)
%
% c: 2*N closed curve of (x,y) pairs of points
% h: caracteristic size of details that should be removed (larger h means
%    more smoothing)
% lout: logical indicating that the interior

% loop on all curves
if ~isempty(c)
    for i1 = 1:length(c.Ocean)
        c.Ocean{i1} = smoothSingleCurve(c.Ocean{i1},h);
    end
    c.Ocean = c.Ocean(~cellfun('isempty',c.Ocean));
    for i1 = 1:length(c.Land)
        c.Land{i1} = smoothSingleCurve(c.Land{i1},h);
    end
    c.Land = c.Land(~cellfun('isempty',c.Land));
end

% FUNCTION SMOOTHSINGLECURVE
function c = smoothSingleCurve(c,h)

% constants
N = size(c,2);
cx = c(1,:);
cy = c(2,:);
gerr = 1e-8;

% initialization
ind = false(1,N);

% identify and keep sides of the box
mincx = min(cx);
mincy = min(cy);
maxcx = max(cx);
maxcy = max(cy);
indb = (abs(cx-mincx)<gerr) | (abs(cy-mincy)<gerr) | ...
       (abs(cx-maxcx)<gerr) | (abs(cy-maxcy)<gerr );
ind( indb ) = true;

% identify and keep points separated by distance h
cumd = cumsum(sqrt(diff(cx).^2+diff(cy).^2));
Nh = ceil(cumd(end)/h);
for i1 = 1:Nh
    ind( find(abs(cumd-(i1-1)*h)==min(abs(cumd-(i1-1)*h)),1) ) = true;
end

% identify corners that are too sharp and remove
% TODO

% output
c = c(:,ind);

% remove completely if there are less than two nodes or if area is too
% small
if (size(c,2)<3) || (polyarea(c(1,:),c(2,:))<h^2)
    c = [];
    return
end
