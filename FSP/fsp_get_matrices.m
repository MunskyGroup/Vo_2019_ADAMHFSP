function Aterms = fsp_get_matrices(S, ind_prop, nmax)
% Generate the parameter-independent terms to form the FSP matrix on a
% hyper-rectangle.
%
% Arguments:
% ---------
% S         : stoichiometry matrix of size N x M, where N is the number of species,
%             M the number of reactions.
%
% ind_prop  : function handle to compute the parameter-independent part of
%             the propensities.
%
% nmax      : vector of maximum numbers of molecules of all species.
%
%
%% Initialize variables
N = size(S,1);   % Number of species.
M = size(S,2);   % Number of reactions.

%% Compute Inf. Gen. Matrix. on the box defined by the first N constraints
Nst = prod(nmax+1);                     % Total number of states.

% Find the propensities of all reactions at all states in the
% hyper-rectangle ...
X = cell(1,N);
Xvec = cell(1,N);
for i = 1:N
    Xvec{i} = (0:nmax(i));
end
[X{1:N}] = ndgrid(Xvec{1:N}); % Find all states in the N-d hyper-rectangle ...
for i = 1:N
    X{i} = reshape(X{i}, numel(X{i}), 1);
end
props = ind_prop(X{1:N});
%keyboard
Aterms = cell(1,M);
for mu = 1:M
    % transform mu^th stoichiometry vector into 1D coordinates ...
    stoich_1D = S(1,mu);
    for i = 2:N
        stoich_1D = stoich_1D + S(i,mu)*prod(nmax(1:i-1)+1);
    end

    X_new = cell(1,N);
    for i = 1:N
        X_new{i} = X{i} + S(i,mu);
    end
    vaild_constraints = (cell2mat(X_new) >= zeros(Nst,N)) & (cell2mat(X_new) <= repmat(nmax, Nst, 1));
    props_keep = props(:,mu); props_keep(min(vaild_constraints,[],2)~=1) = 0;
    %keyboard
    % adding the matrix terms corresponding to the mu-th reaction
    Aterms{mu} = spdiags(props_keep,-stoich_1D,Nst,Nst)-spdiags(props(:,mu),0,Nst,Nst);
end
end
