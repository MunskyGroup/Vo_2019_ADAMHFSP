function P_out = fsp_solve(parameters, t_array, Acell, P0)
%
% Arguments:
% ---------
%
% parameters: vector of parameters.
%
% t_array: vector of time points.
%
% Acell: cell array of size 1 x M, where M is the number of reactions containing
%        the parameter-independent terms of the FSP matrix.
%
% P0: initial probability vector;
%
% Output:
% ------
% P_out: cell array of size 1 x nt. P_out{k} is the approximate solution at
%        time t_array(k). P_out{k} has the same length as the full FSP
%        solution.

M = size(Acell,2);
nt = length(t_array);

A_now = parameters(1)*Acell{1};
for i = 2:M
    A_now = A_now + parameters(i)*Acell{i};
end

[P_out,~,~] = expvt(t_array, A_now, P0, 1.0e-14, 30);

P_out = mat2cell(P_out, length(P0), ones(1, nt));

end
