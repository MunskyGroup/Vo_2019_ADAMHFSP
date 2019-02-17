function   [PHIcell, reduced_matrices_new] = rom_update_matrices( t_array, Acell, P0, ...
    reduced_matrices_old, ...
    PHIold, PHIadd )
%rom_update_matrices Update the reduced model with information from the new basis
%vectors. This is a technical routine to be used with rom_expand.
%
% [PHIcell, reduced_matrices_new] = rom_update_matrices( t_array, Acell, P0, ...
%     reduced_matrices_old, ...
%     PHIold, PHIadd )
%
% Arguments:
% ===========
%
% t_array (in): one-dimensional array of endpoints of the reduced basis set's
% time subintervals.
%
% Acell (in): one-dimensional cell that contains the parameter-independent term
% in the affine decomposition of the FSP matrix. These matrices must have
% the same size.
%
% P0 (in): starting initial probability distribution. P0 must have the same
% length as the number of columns of each matrix in Acell.
%
% reduced_matrices_old (in) and reduced_matrices_new (out) are structs with the following
% fields:
%
%   B: Mx(nt) cell array, where B(k, j) = PHIcell{j}'*Acell{k}*PHIcell{j}
%
%   q0: initial state of the reduced model on the interval [0, t{1}],
%   computed from P0 and the first reduced basis.
%
%   T: cell array of basis transition matrices to transform the result at
%   the end of [t{k}, t{k+1}] to the start of [t{k+1},t{k+2}].
%
% PHIold (in): one-dimensional cell that contains the set of reduced bases to be
% updated. PHIold{j} is the orthogonal basis matrix used to project the
% full model at the j-th time subinterval.
%
% PHIadd (in): one-dimensional cell that contains the orthogonal bases to
% be concatenated to the bases in PHIold. PHIadd{j} will be concatenated to
% PHIold{j}.
%
% PHIcell (out): updated basis set obtained by concatenating bases in
% PHIold and PHIadd.

nt = length( t_array );
M = length( Acell );

PHIcell = cell(nt, 1);
Bnew = cell(M, nt);
Tnew = cell(nt-1, 1);

for it = 1:nt        
    if (isempty(PHIadd{it})==0)
        % Generate new PHI matrices
        PHIcell{it} = [PHIold{it} PHIadd{it}];
    
        % Generate the new B matrices
        for i = 1:M
            B1 = reduced_matrices_old.B{i, it}; % upper left
            B2 = PHIadd{it}'*(Acell{i}*PHIadd{it}); % lower right
            B3 = PHIold{it}'*(Acell{i}*PHIadd{it}); % upper right
            B4 = (PHIadd{it}'*Acell{i})*PHIold{it}; % lower left
            Bnew{i, it} = [B1 B3; B4 B2];
        end              
        
        if (it < nt)
            % Generate new T matrices
            T1 = reduced_matrices_old.T{it};
            T2 = PHIadd{it+1}'*PHIadd{it};
            T3 = PHIold{it+1}'*PHIadd{it};
            T4 = PHIadd{it+1}'*PHIold{it};
            Tnew{it} = [T1 T3; T4 T2];
        end
    else
        PHIcell{it} = PHIold{it};
        for i = 1:M
            Bnew{i, it} = reduced_matrices_old.B{i, it};
        end
        if (it < nt)
            % Generate new T matrices
            T1 = reduced_matrices_old.T{it};
            T2 = PHIadd{it+1}'*PHIadd{it};
            T3 = PHIold{it+1}'*PHIadd{it};
            T4 = PHIadd{it+1}'*PHIold{it};
            Tnew{it} = [T1 T3; T4 T2];
        end
    end
end
%Update q0
q0 = PHIcell{1}'*P0;

err0 = norm( PHIcell{1}*q0 - P0, 2 );

reduced_matrices_new = struct('B', {Bnew}, 'q0', {q0}, 'T', {Tnew}, 'err0', err0);

end
