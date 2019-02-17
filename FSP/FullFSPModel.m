classdef FullFSPModel
    %FULLFSPMODEL Class for storing the full FSP model. 
    % We assume that the model has an affine dependence on the parameters:
    % A(theta) = c_1(theta)*A_1 + .. + c_M(theta)*A_M
    % where theta is the model parameter, c_1, .., c_M are scalar functions and A_1, .., A_M are constant
    % matrices. Each A_i corresponds to a reaction in the reaction network.
    % This is equivalent to requiring that each propensity function of the
    % CME could be factored into 
    %        a_i(x; theta) = c_i (theta) * b_i(x)
    %
    % Contact:
    % ========
    % Huy Vo. huydvo@colostate.edu
    properties
        Acell % one-dimensional cell array to store the constant matrices in the affine decomposition (see the class description above).
        t_array % one-dimensional scalar array of output times of the FSP solutions.
        P0 % initial probability vector
        nmax % size of the FSP, nmax(i) is the maximum number of molecules of species i.
    end
    
    methods
        function obj = FullFSPModel(t_array, S, ind_prop, nmax, x0)
            %FULLFSPMODEL Construct an instance of FullFSPModel.
            %   Arguments:
            %   =========
            %   t_array : one-dimensional scalar array of output times of the FSP solutions.
            %   S : stoichiometry matrix of size (number of
            %   species)x(number of reactions)
            %   ind_prop : function handle to compute the
            %   parameter-independent part of the propensities (i.e., the
            %   b_i(x) factors in the class description). For example (for a model with two species x,y):
            %       ind_prop = @(x,y) [ones(length(x),1), 1./(1+y), x, ones(length(y),1), 1./(1+x.^2), y];
            %   nmax : size of the FSP, nmax(i) is the maximum number of molecules of species i.
            %   x0 : initial state.
            Acell = fsp_get_matrices(S, ind_prop, nmax); % parameter-independent terms of the FSP matrix
            P0 = zeros(prod(nmax+1),1);
            
            if (length(x0) > 1)
                for i = 1:length(x0)
                    arg{i} = x0(:,i) + 1;
                end
                n0 = sub2ind(nmax+1, arg{1:end});
            else
                n0 = x0+1;
            end
            
            P0(n0) = 1;
            
            obj.t_array = t_array;
            obj.Acell  = Acell;
            obj.P0 = P0;
            obj.nmax = nmax;
        end
        
        function P_out = solve(obj, params, in_log10_space)
            %SOLVE solves the FSP model at the user-input parameters and
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
            % params (in): parameter values to solve the FSP at.
            %
            % in_log10_space (in, optional): true if params is in
            % log10-transformed space (default is false).
            %
            % P_out (out): solutions at the time points stored in
            % obj.t_array.
            if (nargin < 3)
                in_log10_space = false;
            end
            if (in_log10_space)
                params = 10.^params;
            end
            P_out = fsp_solve(params, obj.t_array, obj.Acell, obj.P0);
        end
    end
end

