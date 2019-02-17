function [sampling_result, perf_data] = AMH(sample0, ...
    lprior_eval, ...
    full_ll_eval, full_model,...
    nchain, var, ichange, n_non_adaptive, epsi,...
    printing, fid)
% An implementation of the adaptive Metropolis algorithm.
%
% Syntax:
% =======
%
% function [sampling_result, perf_data] = AMH(sample0, ...
%     lprior_eval, ...
%     full_ll_eval, full_model,...
%     nchain, var, ichange, n_non_adaptive, epsi,...
%     printing, fid)
%
% Arguments:
% =========
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
%   full_model (in): data structure needed for the log-likelihood evaluations.
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
%   printing (in, optional): if set to 1 the chain information will be
%   printed, if set to 0 no printing. (If left blank, we assume the default
%   value is 1 ).
%
%   fid (in, optional): if specified, id of the file to print the chain
%   information. Otherwise, all printings are on screen.
%
%   sampling_result (out): struct to stores the posterior chain with the
%   following fields
%           samples: array of size d x nchain storing the accepted samples.
%           L: unormalized log-posterior value of the accepted samples.
%           n_reject: number of MCMC rejections.
%           cov_matrix: final covariance matrix used by the adaptive
%           Gaussian proposal.
%           proposals: array of size d x nchain storing values of the
%           value proposed at each iteration.
%
%   perf_data (out): struct to store other information about the MCMC run
%   with the fields
%           cput: 1xnchain vector to store the CPU time taken.
%
% Reference:
% =========
% Haario, H., Saksman, Ee., & Tamminen, J. (2001). An Adaptive Metropolis Algorithm. Bernoulli, 7(2), 223.
%  https://doi.org/10.2307/3318737
%
% Contact (about this implementation):
% ===================================
% Huy Vo. huydvo@colostate.edu

if (nargin < 7 || isempty(ichange))
    ichange = 1:length(sample0);
end
if (nargin < 8 || isempty(n_non_adaptive))
    n_non_adaptive = 10;
end
if (nargin < 9 || isempty(epsi))
    epsi = 1.0e-6;
end
if (nargin < 10 || isempty(printing))
    printing = 1;
    fid = [];
end
if (nargin < 11)
    fid = [];
end


np = length(sample0);
params = zeros(np, nchain);
proposals = zeros(np, nchain);
sd = 2.4^2/length(ichange);
cov_matrix = var*eye(length(ichange));

params(:,1) = sample0;
Lposterior = zeros(nchain,1);

cputime = zeros(nchain,1);

% Evaluate the log-likelihood of the first parameter here ...
Lposterior(1) = full_ll_eval(sample0, full_model) + lprior_eval(sample0);

rejection = 0;

tic
for i = 1:nchain-1
    if (printing)
        if (isempty(fid))
            disp('-------');
            disp(['i=' num2str(i)]);
            disp(['L=', num2str(Lposterior(i))]);
            disp(['params=', num2str(params(:,i)')]);
            disp(['rejection rate=', num2str(rejection*100/i) '%']);
        else
            fprintf(fid,'-------\n');
            fprintf(fid,'i= %d \n',i);
            fprintf(fid, 'L= %.2e \n', Lposterior(i));
            fprintf(fid, 'params= %s \n', sprintf('%.2e , ',params(:,i)'));
            fprintf(fid, 'rejection rate= %.2f per cent \n', rejection*100/i);
        end
    end
    
    % Compute the covariance matrix
    if (i <= n_non_adaptive)
        cov_matrix = var*eye(length(ichange));
    else
        cov_matrix = sd*cov(params(ichange, 1:i)') + sd*epsi*eye(length(ichange));
    end
    
    % Generate a new parameter guess here ...
    pcandidate = params(:,i);
    pcandidate(ichange) = mvnrnd( params(ichange,i)', cov_matrix );
    proposals(:,i) = pcandidate;
    
    % Evaluate the log-likelihood of the candidate here ...
    lposterior_candidate = full_ll_eval(pcandidate, full_model) + lprior_eval(pcandidate);
    
    % Decide whether to accept or reject the candidate with the Metropolis criteria
    if log(rand()) > (lposterior_candidate-Lposterior(i))
        Lposterior(i+1) = Lposterior(i);
        params(:,i+1) = params(:,i);
        rejection = rejection+1;
    else
        Lposterior(i+1) = lposterior_candidate;
        params(:,i+1) = pcandidate;
    end
    cputime(i) = toc;
end

sampling_result = struct('samples', params, 'L', Lposterior, 'n_reject', rejection, 'cov_mat', cov_matrix, ...
    'proposals', proposals);
perf_data = struct('cput', cputime);
end
