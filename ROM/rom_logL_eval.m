function LL = rom_logL_eval(params, rom, data_cell, unobserved_species)
%rom_logL_eval Evaluate an approximation to the log-likelihood of the
%CME model using Krylov basis.
if (~isempty(find(params<0,1)))
    LL = -Inf;    
    return
end
P_cell = rom.solve(params);
LL = loglikelihood_rom(data_cell, P_cell, rom.nmax, unobserved_species);
end
