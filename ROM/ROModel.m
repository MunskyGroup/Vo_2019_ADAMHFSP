classdef ROModel
    %ROMMODEL Class for the reduced-order models of the full FSP.    
    % We use a time subdivision scheme, which splits the full time inerval
    % [0, tfinal] to subintervals [t(j), t(j+1)]. Each of these
    % subintervals corresponds to a different orthogonal basis that the
    % full FSP dynamic is projected onto (see our paper for more detail).
    % Contact:
    % ========
    % Huy Vo. huydvo@colostate.edu
    properties        
        nmax % size of the FSP, nmax(i) is the maximum number of molecules of species i.
        t_out % one-dimensional scalar array of output times of the reduced model solutions.
        tbasis_array % one-dimensional scalar array of subinterval endpoints for the reduced basis set.
        PHIcell % cell array of reduced bases
        matrices % cell array to store the reduced model matrices        
        max_basis_size = 500; % maximum basis size allowed
        basis_tol = 1.0e-8; % tolerance for the basis construction routine
    end
    
    methods
        function obj = ROModel(nmax, t_out, t_basis)
            %ROMMODEL Construct an instance of ROModel
            % nmax(in): size of the FSP, nmax(i) is the maximum number of molecules of species i.
            % t_out(in): one-dimensional scalar array of time points the model is supposed to output solutions at.
            % t_basis(in): endpoints of the subintervals for the basis set.
            obj.nmax = nmax;
            obj.tbasis_array = t_basis;
            obj.t_out = t_out;
            obj.PHIcell = [];
            obj.matrices = [];
        end
        
        function P_out = solve(obj, params, in_log10_space)
            %SOLVE solves the reduced FSP model at the user-input parameters and
            %output the solution vectors at the model's output times.
            %
            % Syntax:
            % ======
            % P_out = obj.solve(params)
            %
            % P_out = solve(obj, params)
            %
            % Arguments:
            % ==========
            % obj (in): the FullFSPModel object (if not using the 'dot' syntax)
            %
            % params (in): parameter values to solve the reduced model at.
            %
            % in_log10_space (in, optional): true if params is in
            % log10-transformed space (default is false).
            %
            % P_out (out): solutions at the time points stored in
            % obj.t_array. (note: these vectors have the same size as the
            % full model)
            if (nargin < 3)
                in_log10_space = false;
            end
            if (in_log10_space)
                params = 10.^params;
            end
            P_out = rom_solve(params, obj.tbasis_array, obj.PHIcell, obj.t_out, obj.matrices);
        end
    end 
    methods(Static)
        function [rom_new, flag] = update(log_params, rom_old, FSP_model, custom_local_basis_fun)
            % UPDATE reduced model update.
            % 
            % Syntax:
            % =======
            % [rom_new, flag] = update(log_params, rom_old, FSP_model, custom_local_basis_fun)
            %
            % Arguments:
            % =========
            % log_params (in): parameters at which to collect the local reduced
            % moodel (in log10 space).
            % 
            % rom_old (in): the ROModel object to update.
            %
            % FSP_model (in): the FullFSPModel object, storing information
            % about the full parameter-dependent FSP model.
            %
            % custom_local_basis_fun (in, optional): function that builds a
            % reduced basis set for the FSP dynamics at a specific
            % point in parameter space.
            %
            % rom_new (out): the updated ROModel object.
            % 
            % flag (out): output flag
            %       flag = 0 if the basis update is successful.
            %       flag = 1 if all bases exceed the maximum allowable size
            %       (stored in rom_old.max_basis_size).
            if (nargin < 4)
                custom_local_basis_fun = @(par, tbasis_array, Acell, P0) krylov_get_basis(par, tbasis_array, Acell, P0, rom_old.max_basis_size, rom_old.basis_tol);
            end            
            [rom_new, flag] = rom_expand(10.^log_params, rom_old, FSP_model, custom_local_basis_fun);
        end
    end
end

