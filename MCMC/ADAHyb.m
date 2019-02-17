function [sampling_result, perf_data, approx_model] = ADAHyb( params0, ...
    lprior_eval, ... 
    full_ll_eval, full_model, ... 
    approx_ll_eval, approx_model, ... 
    approx_update, tol, ... 
    phase1_length, phase2_length, var, ichange,...
    n_non_adaptive_cov, epsi)
%ADAHYB a hybrid chain that spends its first portion like ADAMH and its
%remaining portion using exclusively the reduced model-based likelihood
%(that is built by
%ADAMH in the first phase).
%
% Syntax:
% =======
%   [sampling_result, perf_data, approx_model] = ADAHyb( params0, ...
%     lprior_eval, ... 
%     full_ll_eval, full_model, ... 
%     approx_ll_eval, approx_model, ... 
%     approx_update, tol, ... 
%     phase1_length, phase2_length, var, ichange,...
%     n_non_adaptive_cov, epsi)
%
% Arguments:
% =========
%
%   params0 (in): chain starting point
%
%   lprior_eval (in): function handle to evaluate the log-prior, for example
%                 lprior_eval = @(x, args) -x(1)^2 - x(2)^2;
%
%   full_ll_eval (in): function handle to compute the exact log-likelihood, for
%   example
%                  full_ll_eval = @(x, args) f(x);
%
%   full_model (in): data structure needed for full log-likelihood evaluations.
%
%   approx_ll_eval (in): function to compute the approximate
%   log-likelihood, having the same syntax as full_ll_eval.
%
%   approx_model (in): data structure of the surrogate model.
%
%   approx_update (in): function to perform updates on approx_model, need
%   to be in the form
%           new_approx_model = approx_update(x_update, approx_model_old, full_model)
%   here x_update is the parameter sample used to update the approximate
%   model.
%
%   tol (in): tolerance threshold for model update sensitivity. Update is
%   considered if |approx_ll_eval(x) - full_ll_eval(x)| >
%   tol*|full_ll_eval(x)|
%
%   phase1_length (in): number of MCMC iterations using ADAMH transition kernel.
%
%   phase2_length (in): number of MCMC interations using AMH and
%   approximate posterior from the reduced model output from phase 1.
%
%   var (in): initial variance of the Gaussian proposal density.
%
%   ichange (in, optional): indices of the components to be perturbed in the MCMC.
%
%   n_non_adaptive_cov (in, optional): number of iterations where the proposal covariance
%   is not updated.
%
%   epsi (in, optional): parameter in the covariance update formula.
%
%   sampling_result (out): struct with phase1 and phase2 storing the
%   outputs of the ADAMH phase and the ROM-based AMH phase. (see the
%   documentation of ADAMH and AMH for the structures of the outputs).
%
%   perf_data (out): struct with phase1 storing performance information of
%   the ADAMH phase and phase2 storing the information of the ROM-based AMH
%   phase.
%
%   approx_model (out): the approximate model used in the second phase.
%
%
%Contact:
%=======
%   Huy Vo. huydvo@colostate.edu
%%
[damh_result, damh_perf, approx_model] = ...
    ADAMH( params0, ...
    lprior_eval,... % function to evaluate log(prior)
    full_ll_eval, full_model, ... % function and data for full model likelihood evaluation
    approx_ll_eval, approx_model, ... % function and data for approximate likelihood evaluation
    approx_update, tol, ... % function to update the approximation
    phase1_length, var, ichange,...
    n_non_adaptive_cov, epsi);
%% Run MCMC with pure ROM evaluations
phase2_start = damh_result.samples(:,end);
[rom_amh_result, rom_amh_perf] = AMH(phase2_start, ...
    lprior_eval, ...
    approx_ll_eval, approx_model,...
    phase2_length, var, ichange, n_non_adaptive_cov, epsi);
%%
sampling_result = struct('phase1', damh_result, 'phase2', rom_amh_result);
perf_data = struct('phase1', damh_perf, 'phase2', rom_amh_perf);
end

