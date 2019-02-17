function X = fsp_get_states( nmax )
% Given the max molecular counts in the FSP, generate all states in the
% hyper-cubic FSP.
%
% Input:
% ====
% nmax: vector of length N for max molecule populations.
%
% Output:
% =====
%
% X( 1:nst, 1:N ) : array of states, nst = number of states, N = number of
% species.

N = length( nmax );

X = cell(1,N);
Xvec = cell(1,N);
for i = 1:N
    Xvec{i} = (0:nmax(i));
end
[X{1:N}] = ndgrid(Xvec{1:N}); % Find all states in the N-d hyper-rectangle ...
for i = 1:N
    X{i} = reshape(X{i}, numel(X{i}), 1);
end
X = cell2mat(X);
end
