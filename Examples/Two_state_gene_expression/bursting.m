clear;
%%
rng(0, 'twister');
%% Parameters of the numerical algorithmss

% Parameters for the initial likelihood optimization
Npop = 100; % population size for each generation of genetic algorithm
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

% Parameters for reduced basis routines
max_basis_size = 500; % maximum number of basis vectors per subinterval for the reduced model
basis_tol = 1.0e-8; % Krylov basis construction stopping tolerance
%% Parameters for the full model
model_name = 'bursting';
S = [1 0;-1 0;0 1;0 -1]'; % stoichiometry matrix
nmax = [1 1100];
nb = 100;
T_basis = linspace(1/nb, 1, nb);
T_data = linspace(0.1, 1, 10);
T_basis = sort( unique( [ T_basis T_data ] ), 'ascend' );
true_params = [0.5, 0.8, 1000, 1]';
params_dim = length(true_params);
ind_prop = @(gene, rna) [gene==0, gene==1, gene==1, rna]; % the parameter-independent part of the propensities
% indices of parameters that will vary
ichange = [ 1 2 3 4];
% gene state is not observable
unobserved_species = [1];
x0 = [1 0];
% Define the prior
params_min = [1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6]';
params_max = [10, 10, 10000, 10]';
lprior_eval = @(log_params) log(prod(10.^log_params >= params_min)*prod(10.^log_params <= params_max));
save([model_name '_basic_info.mat']);
%% Generate data
nsample = 200;
true_prop = @(X) true_params'.*ind_prop(X(1), X(2));
data = data_gen( S', true_prop, x0, T_data, nsample, unobserved_species );
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

%Step 1: GA with Brian's mutation function
Constrs = struct('LB',log10(params_min'),'UB',log10(params_max'));
Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,Constrs);

ga_opts = gaoptimset('Display','iter','useparallel',0,...
    'MutationFcn',Mut_Fun,...
    'PopulationSize',Npop);
ga_opts = gaoptimset(ga_opts,'EliteCount',1,'Generations',5);
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
load([model_name '_damh_results.mat'], 'rng_state');
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