function datcell = data_gen( stoich_mat, propensity, x0, t_array, n_sample, unobservable)
% datcell = data_gen( stoich_mat, propensity, x0, t_array, n_sample,
% unobservable) synthesizes smFISH data at specific measurement time points
% (for numerical testing purpose).
datcell = cell( length(t_array ), 1 );
n_species = length(x0);

for it = 1:length( t_array )

  t_ssa = t_array(it);


  ssa_samples = zeros(n_sample, n_species);
  for i = 1:n_sample
    ssa_samples(i, :) = ssa(t_ssa, x0, stoich_mat, propensity);
  end

  ssa_samples = ssa_samples(:, setdiff((1:n_species), unobservable));

  [observed_states, ~, i_observed] = unique(ssa_samples, 'rows');
  n_observed = size(observed_states, 1);

  [wbin, jbin] = hist(i_observed, unique(i_observed));

  observed_freq = zeros(n_observed, 1);
  observed_freq(jbin) = wbin;

  datcell{it} = [observed_states observed_freq];
end

end
