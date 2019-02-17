function end_state = ssa(t_f, x0, stoich_mat, propensity)
% SSA stochastic simulation algorithm.
%
% t_f(in): final time.
%
% x0 (in): starting intial vector.
%
% stoich_mat (in): stoichiometry matrix of size (n.o.reactions) x
% (n.o.species).
%
% propensity (in): function handle, propensity(x) = [a_1(x) a_2(x) .. a_M(x)]
% where a_k(x) is the k-th propensity function evaluated at x.
%
% Reference:
% =========
% ?Gillespie, D. T. (1977). Exact Stochastic Simulation of Coupled Chemical Reactions. J. Phys. Chem., 81(25), 2340?2361. 
%
  x = x0;
  t = 0;
  while (t < t_f)
    a = propensity(x);
    a0 = sum(a);

    u1 = rand();
    u2 = rand();

    % Update time
    tau = min(log(1/u1)/a0, t_f - t);
    t = t + tau;

    % Update state
    r = find(cumsum(a) > u2*a0, 1, 'first');
    x = x + stoich_mat(r, :);
  end
  end_state = x;
end
