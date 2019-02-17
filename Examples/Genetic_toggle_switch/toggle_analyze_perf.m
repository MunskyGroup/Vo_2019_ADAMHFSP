clear;

model = 'toggle';
methods = {'AMH-FSP', 'ADAMH-FSP-Krylov', 'Hybrid'};
param_names = {'k_{0X}', 'k_{1X}', '\\gamma_X', 'k_{0Y}', 'k_{1Y}', '\\gamma_Y'};
param_names_plot = {'k_{0X}', 'k_{1X}', '\gamma_X', 'k_{0Y}', 'k_{1Y}', '\gamma_Y'};

flag = plot_figures(model, methods, param_names, param_names_plot);
flag = analyze_perf(model, methods, param_names, param_names_plot);

