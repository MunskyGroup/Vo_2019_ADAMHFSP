function y = marginal(w, nfsp, dim)
% find the 1-d marginal distribution of the FSP solution. This requires 
% the FSP to be hyper-rectangular.
%
% Arguments:
% ---------
% 
% w: the column vector storing the solution. The indexing is big-endian:
% (0,0,..,0) -> (0,0,..,1) -> (0,0,..,2) ... 
%
% nfsp: maximum number of molecules in each dimension. 
%
% dim: the index of the species/dimension of which we compute the marginal
%      distribution.

N = length(nfsp);
arg = num2cell(nfsp+1);
tensor = reshape(w, arg{1:N});
pd = unique([dim (1:N)], 'stable');
tensor = permute(tensor, pd);
tensor = reshape(tensor, nfsp(dim)+1, prod( nfsp( pd(2:end) ) +1 ));
y = sum(tensor, 2);
end