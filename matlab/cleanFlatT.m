function T = cleanFlatT( T, X, gerr )
% get rid of elements that have area zero
d = size(X,2);
if d==2
    xi = X(:,1); yi = X(:,2);
    S = polyarea( xi(T)', yi(T)' )';
elseif d==3
    S = (cross(X(T(:,2),:)-X(T(:,1),:),X(T(:,3),:)-X(T(:,1),:)));
    S = sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2);
else
    error('unknown dimension type')
end
ind = S > gerr;
T = T(ind,:);
