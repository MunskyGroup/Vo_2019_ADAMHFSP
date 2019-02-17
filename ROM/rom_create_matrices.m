function reduced_matrices = rom_create_matrices( t_array, Acell, P0, PHIcell)
%reduced_matrices = rom_create_matrices( t_array, Acell, P0, PHIcell)
%generates matrices that represent the reduced system by projecting the
%full dynamics p' = Ap to a sequence of reduced order models q' = B{k}*q over
% sub-intervals [t{k}, t{k+1}], K = 0, 1, .., nt.
%
% Input:
% ======
% t_array: array of time nodes t{1}, .., t{nt} that delimit the subdivision
% of the global time interval [0,T].
%
% Acell: cell containing the parameter-independent matrices A{1}, .., A{M},
% assuming that the system operator is parameter-separable.
%
% P0: initial system state.
%
% PHIcell: cell containing the bases of the reduced subspaces, PHIcell{k}
% is the basis used for the interval [t_{k-1}, t_{k}].
%
% Output:
% ======
% reduced_matrices: struct with the following fields
%
%   B: Mx(nt) cell array, where B(k, j) = PHIcell{j}'*Acell{k}*PHIcell{j}
%
%   q0: initial state of the reduced model on the interval [0, t{1}],
%   computed from P0 and the first reduced basis.
%
%   T: cell array of basis transition matrices to transform the result at
%   the end of [t{k}, t{k+1}] to the start of [t{k+1},t{k+2}].
%
%   err0: initial error made by projecting P0 to the subspace spanned by
%   PHIcell{1}

nt = length(t_array); M = length(Acell);

Bcell = cell( M, nt );

for it = 1:nt
    for k = 1:M
        Bcell{k, it} = PHIcell{it}'*(Acell{k}*PHIcell{it});
    end
end

q0 = PHIcell{1}'*P0;
Tcell = cell( nt-1, 1 );
for it = 1:nt-1
    Tcell{it} = PHIcell{it+1}'*PHIcell{it};
end

err0 = norm( PHIcell{1}*q0 - P0, 2 );

reduced_matrices = struct('B', {Bcell}, 'T', {Tcell}, 'q0', {q0}, 'err0', err0);
end
