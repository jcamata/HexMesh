function [T,X] = cleanT(T,X,gerr)
% get rid of elements that are repeated
[~,ind] = unique( sort(T,2) , 'rows', 'first', 'legacy' );
T = T(ind,:);
% get rid of elements that have area zero
xi = X(:,1); yi = X(:,2);
ind = polyarea( xi(T)', yi(T)' )' > gerr;
T = T(ind,:);
% get rid of nodes that are not used in T
indX = unique( T, 'legacy' );
for i1 = 1:length(indX)
    T( T==indX(i1) ) = i1;
end
X = X(indX,:);
