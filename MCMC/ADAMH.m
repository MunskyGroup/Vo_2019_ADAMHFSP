function [sampling_result, perf_data, approx_model ] = ...
    ADAMH( sample0, ...
    lprior_eval, ...
    full_ll_eval, full_model, ...
    approx_ll_eval, approx_model, ...
    approx_update, tol, ...
    nchain, var, ichange,...
    n_non_adaptive_cov, epsi, adapt_prob)

% Adaptive Delayed Acceptance Metropolis with diminishing model
% adaptations.
%
%   Syntax:
%   =======
%
%   [sampling_result, perf_data, approx_model ] = ...
%     ADAMH( sample0, ...
%     lprior_eval, ...
%     full_ll_eval, full_model, ...
%     approx_ll_eval, approx_model, ...
%     approx_update, tol, ...
%     nchain, var, ichange,...
%     n_non_adaptive_cov, epsi)
%
%   Arguments:
%   ==========
%
%   sample0 (in): chain starting point
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
%   nchain (in): number of MCMC iterations.
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
%   adapt_prob (in, optional): vector of at least nchain elements,
%   the ADAMH algorithm will ignore model updates at step i with
%   probability 1 - adapt_prob(i).
%
%   sampling_result (out): struct to stores the posterior chain with the
%   following fields
%           samples: array of size d x nchain storing the accepted samples.
%
%           L: unormalized log-posterior value of the accepted samples.
%
%           L_approx: approximate value of the unormalized log-posterior
%
%           value of the accepted samples, computed via approximate_model.
%
%           rejection: number of MCMC rejections.
%
%           cov_matrix: final covariance matrix used by the adaptive
%           Gaussian proposal.
%
%           proposals: array of size d x nchain storing values of the
%           value proposed at each iteration.
%
%   perf_data (out): struct to store other information about the MCMC run
%   with the fields
%           samples_update: array to store the sample values at which the
%           local reduced models are collected to build the global model.
%
%           cput: 1xnchain vector to store the CPU time taken.
%
%           n_full_eval: number of full model evaluations.
%
%           n_approx_eval: number of approximate model evaluations.
%
%           n_first_reject: number of proposals rejected by the first
%           (cheap) stage.
%
%           n_second_reject: number of proposals rejected by the second
%           (expensive) stage. For efficiency, this value should be as
%           small as possible.
%
%           n_promote: number of proposals promoted to the second-stage
%           evaluation.
%
%           full_eval_time: 1 x nchain array to store the time taken to
%           evaluate the full models.
%
%           approx_eval_time: 1 x nchain array to store the time taken to
%           evaluate the reduced models.
%
%           approx_update_time: 1 x nchain array to store the time taken to
%           update the model.
%
%   approx_model (out): final approximate model output by the ADAMH.
%
%   Contact:
%   =======
%   Huy Vo. huydvo@colostate.edu
if (nargin < 11 || isempty(ichange))
    ichange = 1:length(sample0);
end
if (nargin < 12 || isempty(n_non_adaptive_cov))
    n_non_adaptive_cov = 1;
end
if (nargin < 13 || isempty(epsi))
    epsi = 1.0e-6;
end
if (nargin < 14 || isempty(adapt_prob))
    adapt_prob = 2.^(-(1:nchain)/10000);
end

np = length(sample0); % number of parameter components
samples = zeros(np, nchain);
proposals = zeros(np, nchain);
rejection = 0;

% Initial proposal covariance matrix
cov_matrix = var*eye(length(ichange));
% coefficient for the adaptive covariance matrix
sd = 2.4^2/length(ichange);

% arrays to store the exact and approximate log-posteriors of accepted parameters
Lposterior = zeros(nchain, 1);
Lposterior_approx = zeros(nchain, 1);

% record performance info
n_full_eval = 1;
n_approx_eval = 0;
n_first_reject = 0;
n_second_reject = 0;
n_promote = 0;
cputime = zeros(nchain, 1);
samples_update = sample0;
full_eval_time = zeros(nchain, 1); % record the time of the full eval part in each iteration
approx_eval_time = zeros(nchain, 1); % record the time of the approx eval part in each iteration
approx_update_time = zeros(nchain, 1); % record the time to update the approx model in each iteration

Lposterior(1) = full_ll_eval(sample0, full_model) + lprior_eval(sample0);

t_update = tic;
approx_model_zero = approx_model;
approx_model = approx_update(sample0, approx_model, full_model);
approx_update_time(1) = toc(t_update);

Lposterior_approx(1) = approx_ll_eval(sample0, approx_model) + lprior_eval(sample0);

solver_time_ratio = 1000;
mcmc_time = tic;
stop_update_flag = 0;
L_min = Lposterior(1);
samples(:,1) = sample0;
for iter = 1:nchain-1
    disp('-------');
    disp(['i=' num2str(iter)]);
    disp(['L=', num2str(Lposterior(iter))]);
    disp(['L_approx=', num2str(Lposterior_approx(iter))]);
    disp(['params=', num2str(samples(:,iter)')]);  % display in linear space
    disp(['rejection rate=', num2str(rejection*100/iter) '%']);
    disp(['full/approx time ratio=', num2str(solver_time_ratio)]);
    disp(['false promotion rate=', num2str(100*n_second_reject/n_promote)]);
    disp(['reject with fsp=', num2str(n_second_reject)]);
    % Compute the covariance matrix
    if (iter <= n_non_adaptive_cov)
        cov_matrix = var*eye(length(ichange));
    else
        cov_matrix = sd*cov(samples(ichange, 1:iter)') + sd*epsi*eye(length(ichange));
    end
    
    %% First stage acceptance
    promote = false;
    Y = samples(:, iter);
    
    % Generate a new parameter guess here ...
    pcandidate = Y;
    pcandidate(ichange) = mvnrnd( Y(ichange)', cov_matrix ) ;
    proposals(:, iter) = pcandidate;
    
    % Evaluate the log-likelihood of the candidate here ...
    t_approx = tic;
    L_reduced_candidate = approx_ll_eval(pcandidate, approx_model);
    lprior_candidate = lprior_eval(pcandidate);
    approx_eval_time(iter) = toc(t_approx);
    
    lposterior_reduced_candidate = L_reduced_candidate + lprior_candidate;
    n_approx_eval = n_approx_eval + 1;
    
    r1 = rand();
    
    if log(r1) < (lposterior_reduced_candidate - Lposterior_approx(iter))
        promote = true;
        Y = pcandidate;
    end
    
    %% Second stage acceptance
    if ( promote )
        n_promote = n_promote+1;
        % Evaluate the full likelihood
        t_full = tic;
        lposterior_candidate = full_ll_eval(Y, full_model) + ...
            + lprior_eval(Y);
        n_full_eval = n_full_eval + 1;
        full_eval_time(iter) = full_eval_time(iter) + toc(t_full);
        solver_time_ratio = min([full_eval_time(iter)/approx_eval_time(iter)]);
        % log of acceptance probability for the second stage
        beta = lposterior_candidate - Lposterior(iter) + (Lposterior_approx(iter) - lposterior_reduced_candidate);
        
        r1 = rand();
        if (log(r1) < beta)
            samples(:,iter +1) = Y;
            Lposterior(iter+1) = lposterior_candidate;
            Lposterior_approx(iter+1) = lposterior_reduced_candidate;
            
            % Reset the basis to discard those built in low posterior density
            % region
            if (Lposterior(iter+1) - L_min >= 20)
                fprintf('Reset the reduced model\n');
                approx_model = approx_model_zero;
                [approx_model, stop_update_flag] = approx_update(Y, approx_model, full_model);
                samples_update = [samples_update Y];
                approx_update_time(iter) = approx_update_time(iter) + toc(t_update);
                Lposterior_approx(iter+1) = approx_ll_eval(samples(:,iter+1), approx_model)+lprior_eval(samples(:,iter+1));
                L_min = Lposterior(iter+1);
            else
                % Expand the basis if the relative error in log-likelihood is high
                rel_error = abs(lposterior_reduced_candidate - lposterior_candidate)/abs(lposterior_candidate);
                if (rel_error > tol && stop_update_flag == 0)
                    %         if ((error_indicator > tol) && (i_adaptive == true))
                    r1 = rand();
                    if (r1 < adapt_prob(iter))
                        % decide to update the reduced model using the information from
                        % the accepted parameter
                        t_update = tic;
                        [approx_model, stop_update_flag] = approx_update(Y, approx_model, full_model);
                        samples_update = [samples_update Y];
                        approx_update_time(iter) = approx_update_time(iter) + toc(t_update);
                        Lposterior_approx(iter+1) = approx_ll_eval(samples(:,iter+1), approx_model)+lprior_eval(samples(:,iter+1));
                        
                        X_lr = [ones(1,1) Y'];
                        y_lr = Lposterior(iter+1) - Lposterior_approx(iter+1);
                    end
                end
            end
        else
            n_second_reject = n_second_reject + 1;
            samples(:,iter+1) = samples(:,iter);
            Lposterior(iter+1) = Lposterior(iter);
            Lposterior_approx(iter+1) = Lposterior_approx(iter);
            rejection = rejection+1;
        end
    else
        n_first_reject = n_first_reject + 1;
        samples(:,iter+1) = samples(:,iter);
        Lposterior(iter+1) = Lposterior(iter);
        Lposterior_approx(iter+1) = Lposterior_approx(iter);
        rejection = rejection+1;
    end
    
    cputime(iter) = toc(mcmc_time);
end
sampling_result = struct('samples', samples, 'L', Lposterior, 'L_approx', Lposterior_approx, 'rejection', rejection, 'cov_matrix', cov_matrix,...
    'proposals', proposals);
perf_data = struct( 'samples_update', samples_update, 'cput', cputime, ...
    'n_full_eval', n_full_eval,'n_approx_eval', n_approx_eval, 'n_first_reject', n_first_reject,...
    'n_second_reject', n_second_reject, 'n_promote', n_promote, ...
    'full_eval_time', full_eval_time, 'approx_eval_time', approx_eval_time, 'approx_update_time', approx_update_time);
end
