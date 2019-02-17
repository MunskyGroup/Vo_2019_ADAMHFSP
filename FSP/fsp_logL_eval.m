function LL = fsp_logL_eval(params0, full_model, data_cell, unobserved_species)
%fsp_logL_eval Summary of this function goes here
%   Detailed explanation goes here
if (~isempty(find(params0<0, 1)))
    LL = -Inf;
    return
end
P_full_cell = full_model.solve(params0);
LL = loglikelihood(data_cell, P_full_cell, full_model.nmax, unobserved_species);
end
