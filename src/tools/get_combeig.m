function [A allind] = get_combeig(V, p, u)
% -------------------------------------------------
% function A = get_combeig(V, p, u)
% 
% This function gives an array (A) of nchoosek(p,u) 
% matrices each one being a permutation of the columns
% of matrix V.
% -------------------------------------------------
comb = nchoosek(p, u);
ind = 1 : p;
allind = nchoosek(ind, u);
for i = 1 : comb,
    A(i, :, :) = V(:, allind(i, : ));
end
