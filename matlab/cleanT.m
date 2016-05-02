function [T,X] = cleanT(T,X,gerr)
% get rid of elements that are repeated
[~,ind] = unique( sort(T,2) , 'rows', 'first', 'legacy' );
T = T(ind,:);
% get rid of elements that have area zero
T = cleanFlatT(T,X,gerr);
% get rid of nodes that are not used in T
indX = unique( T, 'legacy' );
for i1 = 1:length(indX)
    T( T==indX(i1) ) = i1;
end
X = X(indX,:);
