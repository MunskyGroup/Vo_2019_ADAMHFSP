clear;

model = 'bursting';
methods = {'AMH-FSP', 'ADAMH-FSP-Krylov', 'Hybrid'};
param_names = {'k_{\\text{on}}', 'k_{\\text{off}}', 'k_r', '\\gamma'};
param_names_plot = {'k_{on}', 'k_{off}', 'k_r', '\gamma'};

flag = plot_figures(model, methods, param_names, param_names_plot);
flag = analyze_perf(model, methods, param_names, param_names_plot);

