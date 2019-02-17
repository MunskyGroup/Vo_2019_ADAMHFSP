clear

addpath mcmcstat
include_all

model = 'hog3d';
methods = {'ADAMH-FSP-Krylov', 'Hybrid'};
param_names = {'k_{\\text{gene}}^{+}', 'k_{\\text{gene}}^{-}', 'k_r', '\\gamma_{\\text{nuc}}', 'k_{\\text{trans}}', '\\gamma_{\\text{cyt}}'};
param_names_plot = {'k_{gene}^{+}', 'k_{gene}^{-}', 'k_r', '\gamma_{nuc}', 'k_{trans}', '\gamma_{cyt}'};

flag = plot_figures(model, methods, param_names, param_names_plot);
flag = analyze_perf_hog3d(model, methods, param_names, param_names_plot);


