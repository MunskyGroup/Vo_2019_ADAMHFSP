clear;
%%
rng(0, 'twister');
%% Parameters of the numerical algorithmss

% Parameters for the initial likelihood optimization
Npop = 200; % population size for each generation of genetic algorithm
NGen = 5; % number of generations in the genetic algorithm
NEvalLocalOpt = 500; % maximum number of evaluations for fmincon local search

% Parameters for the MCMC runs
nchain = 100000; % number of MCMC iterations
nchain1 = 10000; % number of MCMC iterations in the first/ADAMH phase of the Hybrid algorithm
nchain2 = 90000; % number of MCMC iterations in the second/approximate phase of the Hybrid algorithm

n_non_adaptive_cov = 10; % number of steps where the Gaussian proposal covariance is not updated
epsi = 1.0e-6; % parameter in the covariance update formula
var = 0.00001; % initial proposal is Gaussian with mean zero and covariance var*eye(d)
tol_damh = 1.0e-4; % log-likelihood approximation error relative tolerance to trigger model update in DAMH

max_basis_size = 100;
basis_tol = 1.0e-6;
%% Parameters for the full model
model_name = 'toggle';

S = [1 0;1 0;-1 0; 0 1; 0 1; 0 -1]'; % stoichiometry matrix

nmax = [ 100 100 ];
% time points of interest
T_data = [2 6 8]*3600;
nb = 10;
T_basis = linspace(max(T_data)/nb, max(T_data), nb);
T_basis = sort( unique( [ T_basis T_data ] ), 'ascend' );

ayx = 2.6e-3; axy = 6.1e-3; nyx = 3; nxy = 2.1; % we assume these parameters are known
%the parameter combination to generate the 'data'
true_params = [2.2e-3, 1.7e-2, 3.8e-4, 6.8e-5, 1.6e-2, 3.8e-4]';
params_dim = length(true_params);
% indices of parameters that will vary
ichange = [ 1 2 3 4 5 6];
unobserved_species = [];
ind_prop = @(x,y) [ones(length(x),1), 1./(1+ayx*(y.^nyx)), x, ones(length(y),1), 1./(1+axy*(x.^nxy)), y]; % the parameter-independent part of the propensities
x0 = [0 0];

%Define the prior
params_min = [1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6]';
params_max = [1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]';
lprior_eval = @(log_params) log(prod(10.^log_params >= params_min)*prod(10.^log_params <= params_max));
save([model_name '_basic_info.mat']);
%% Generate data
nsample = 500;
prop = @(x) true_params'.*ind_prop(x(1), x(2));
data = data_gen( S', prop, x0, T_data, nsample, [] );

rng_state = rng;
save([model_name '_synthetic_data.mat'], 'data', 'rng_state');
%% Define the full model-based loglikelihood evaluation
load([model_name '_synthetic_data.mat'], 'data', 'rng_state');
rng(rng_state);
full_model = FullFSPModel(T_data, S, ind_prop, nmax, x0);
full_eval = @(log_params, full_model) fsp_logL_eval(10.^log_params, full_model, data, unobserved_species);
%% Define the reduced moodel-based loglikelihood evaluation
reduced_model0 = ROModel(nmax, T_data, T_basis);
reduced_model0.max_basis_size = max_basis_size;
reduced_model0.basis_tol = basis_tol;
approx_eval = @(log_params, reduced_model) rom_logL_eval(10.^log_params, reduced_model, data, unobserved_species);
%% Run opt algorithm with a small number of generations to find a good starting point for the chain
opt_fun = @(log_param) -full_eval(log_param, full_model);

%Step 1: GA with custom mutation function
Constrs = struct('LB',log10(params_min'),'UB',log10(params_max'));
Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,Constrs);

ga_opts = gaoptimset('Display','iter','useparallel',0,...
    'MutationFcn',Mut_Fun,...
    'PopulationSize',Npop);
ga_opts = gaoptimset(ga_opts,'EliteCount',1,'Generations',NGen);
tic
[params_ga, L_ga, ga_flag, ga_output, ga_population] = ga(opt_fun, params_dim, [], [], [], [], log10(params_min'), log10(params_max'), [], ga_opts);
ga_time = toc;

% Step 2: Local refinement with fmincon
fmincon_opt = optimoptions('fmincon','Display', 'iter', 'MaxFunctionEvaluations', NEvalLocalOpt);
[params_opt, L_opt, flag] = fmincon(opt_fun, params_ga, [], [], [], [], log10(params_min'), log10(params_max'), [], fmincon_opt);

rng_state = rng;
save([model_name '_opt_results.mat'],  'rng_state', 'params_opt', 'L_opt', 'params_ga', 'L_ga', 'ga_population');
%% Run the MCMC algorithms
load([model_name '_opt_results.mat']);
load([model_name '_synthetic_data.mat'], 'data');
params0 = params_opt';
%% Run adaptive MH with full evaluations
rng(rng_state);
[full_amh_result, full_amh_perf] = AMH(params0, ...
    lprior_eval, ...
    full_eval, full_model,...
    nchain, var, ichange, n_non_adaptive_cov, epsi);
rng_state = rng;
save([model_name '_amh_results.mat'], 'full_amh_result', 'full_amh_perf', 'rng_state', '-v7.3');
%% Run ADAMH with reduced model
load([model_name '_amh_results.mat'], 'rng_state');
rng(rng_state);
[damh_result, damh_perf, approx_model] = ...
    ADAMH( params0, ...
    lprior_eval,... % function to evaluate log(prior)
    full_eval, full_model, ... % function and data for full model likelihood evaluation
    approx_eval, reduced_model0, ... % function and data for approximate likelihood evaluation
    @ROModel.update, tol_damh, ... % function to update the approximation
    nchain, var, ichange,...
    n_non_adaptive_cov, epsi);
rng_state = rng;
save([model_name '_damh_results.mat'], 'damh_result', 'damh_perf', 'rng_state', '-v7.3');
save([model_name '_approx_model.mat'], 'approx_model');
%% Run hybrid ADAMH-ABC chain
load([model_name '_damh_results.mat'], 'approx_model', 'rng_state');
rng(rng_state);
[dahyb_result, dahyb_perf, hyb_approx_model] = ...
    ADAHyb( params0, ...
    lprior_eval,... % function to evaluate log(prior)
    full_eval, full_model, ... % function and data for full model likelihood evaluation
    approx_eval, reduced_model0, ... % function and data for approximate likelihood evaluation
    @ROModel.update, tol_damh, ... % function to update the approximation
    nchain1, nchain2, var, ichange,...
    n_non_adaptive_cov, epsi);
rng_state = rng;
save([model_name '_dahyb_results.mat'], 'dahyb_result', 'dahyb_perf', 'rng_state', '-v7.3');
save([model_name '_dahyb_approx_model.mat'], 'hyb_approx_model');
