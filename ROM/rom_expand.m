function [new_model, flag] = rom_expand(params, old_model, full_model, rom_get_basis)
%rom_expand Extend the parametric reduced model for the FSP using the
%projection-based reduced model at the specifiec parameters. The
%implementation of this procedure is somewhat generic to allow for
%different ways of generating the reduced models.
%
% The full FSP model is assumed to depend "linearly" on the parameters
%
% A(theta) = theta_1*A_1 + theta_2*A_2 + .. + theta_M*A_M
%
% new_model = rom_expand(params, old_model, full_model, rom_get_basis)
% expands old_model by concatenating the basis matrices with a reduced
% basis set generated at the parameter params.
%
% Both new_model and old_model are structs that store information about the
% reduced dynamics of the full FSP model, with the following fields:
%
%   nmax: row vector of maximum number of molecules of the chemical
%   species. For example, nmax = [100 200] means the full FSP contains
%   states in {0,...,100}x{0,...,200}
%
%   tbasis_array: row vectors containing the time mesh for the reduced
%   model.
%
%   PHIcell: cell containing the set of orthogonal bases that generate the reduced model.
%
%   matrices: struct that store the matrices necessary for the fast
%   assembly of the reduced ODE system. See rom_create_matrices for
%   details.
%
% Input:
% ======
% params: column vector of parameter combination at which to collect the
% additional reduced basis set.
%
% old_model: struct for the old model that we wish to extend (see details
% above).
%
% full_model: struct that stores information about the full model, with the
% following fields
%     Acell: cell to store the set of parameter-dependent FSP matrices that
%     constitute the full model matrix.
%
%     P0: initial full FSP probability vector.
%
%     nmax: maximum number of molecules in the full FSP model.
%
% rom_get_basis: function to generate the reduced basis set at the
% parameter params. This function need to have the following interface
%
%   PHI = rom_get_basis(param, T, Afull_cell, P0, max_basis_dim, tol, varargin)
%
% Output:
% =======
% new_model: struct to store the extended model (see the detailed description above).
%
% flag: 0 if the model is updated normally, 1 if the model is not updated
% due to basis sizes exceeding maximum capacity

new_model = old_model;

max_basis_size = 500;

Local_bases = rom_get_basis(params, old_model.tbasis_array, full_model.Acell, full_model.P0);

flag = 0;

if (isempty(old_model.PHIcell))
    new_model.PHIcell = Local_bases;
    new_model.matrices = rom_create_matrices(old_model.tbasis_array, full_model.Acell, full_model.P0, new_model.PHIcell);
else
    % Compute the part of the new POD basis that is not already included in
    % the current basis
    maxAdd = 0;
    
    for it = 1:length(old_model.PHIcell)
        if (size(old_model.PHIcell{it}, 2) < max_basis_size)
            Local_bases{it} = orth_reduce( old_model.PHIcell{it}, Local_bases{it} );
            maxAdd = max(maxAdd, size(Local_bases{it}, 2));
        else
            Local_bases{it} = zeros(size(old_model.PHIcell{it},1), 0);
        end
    end
    
    if (maxAdd == 0) 
        flag = 1;
    end
    
    fprintf('Adding max %d vectors into the basis. Updating the reduced model ... \n', maxAdd);
    [ PHIcell, new_rom_matrices ] = rom_update_matrices( old_model.tbasis_array, full_model.Acell, ...
        full_model.P0, ...
        old_model.matrices,...
        old_model.PHIcell, ...
        Local_bases );
    
    new_model.PHIcell = PHIcell;
    new_model.matrices = new_rom_matrices;
end
end
