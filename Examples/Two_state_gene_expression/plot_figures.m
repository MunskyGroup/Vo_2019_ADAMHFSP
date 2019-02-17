function flag = plot_figures(model, methods, param_names, param_names_plot)
%% Setting sizes
width = 2.5;     % Width in inches
height = 2;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
%%
damh_color = [34, 139, 34]/256; %Forest Green
amh_color = [220, 20, 60]/256; %Crimson
amh_rom_color = [128, 0, 128]/256; %Purple
lr_damh_color = [0, 0, 34]/256;
lr_amh_rom_color = [50, 0, 0]/256;
gp_damh_color = [220, 220, 60]/256;

ls = {'-', '--', '-.', '-.', ':', '-.'};
markers = {'none', 'none', 'none', 'none', 'none', 'o'};
true_color = [0.1, 0.6, 0.8];
%%
stats = cell(size(methods));

load([model '_basic_info.mat']);


if (~exist('unobserved_species', 'var'))
    unobserved_species = [];
end

for i = 1:length(methods)
    method = methods{i};
    if (contains(string(method), 'ADAMH-FSP-Krylov'))
        load([model '_damh_results.mat'], 'damh_result', 'damh_perf');
        samples = damh_result.samples(:,:)';
        L = damh_result.L(:);
        cput = max(damh_perf.cput);
        cumcput = damh_perf.cput(1:end-1);
        n_full_eval = damh_perf.n_full_eval;
        n_reject = damh_result.rejection;
        n_full_reject = damh_perf.n_second_reject;
        color = damh_color;
    elseif (contains(string(method), 'Hybrid'))
        load([model '_dahyb_results.mat']);
        samples = [dahyb_result.phase1.samples(:,:)'; dahyb_result.phase2.samples(:,:)'];
        L = dahyb_result.phase2.L(:);
        cput = max(dahyb_perf.phase2.cput) + max(dahyb_perf.phase1.cput);
        cumcput = [dahyb_perf.phase1.cput(1:end-1); max(dahyb_perf.phase1.cput)+ dahyb_perf.phase2.cput(1:end-1)];
        n_full_eval = 0;
        n_reject = dahyb_result.phase2.n_reject;
        n_full_reject = 0;
        color = lr_amh_rom_color;
    elseif (contains(string(method), 'AMH-FSP'))
        load([model '_amh_results.mat']);
        samples = full_amh_result.samples(:,:)';
        L = full_amh_result.L(:);
        cumcput = full_amh_perf.cput(1:end-1);
        cput = max(full_amh_perf.cput);
        n_full_eval = length(full_amh_result.samples);
        n_reject = full_amh_result.n_reject;
        n_full_reject = n_reject;
        color = amh_color;
    end
    
    % use the parameter average
    mESS = multiESS(samples(:, :));
    
    perf = struct('samples', samples, 'mESS', mESS, 'cput', cput, 'cumcput', cumcput,...
        'cpu_per_eff_sample', cput/mESS, 'n_full_evals', n_full_eval, ...
        'n_reject', n_reject, 'n_full_reject', n_full_reject, 'color',color);
    stats{i} = perf;
    
end
%% Plot all posterior marginal
figure(1);
% Set size for display
defpos = get(0,'defaultFigurePosition');
set(1,'Position', [defpos(1) defpos(2) width*100*ceil(params_dim/2), height*100*3]);

% Set the size for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width*ceil(params_dim/2), height*3];
set(1, 'PaperPosition', defsize);

for k = 1:params_dim
    for i = 1:length(methods)
        hold on;
        subplot(2, ceil(params_dim/2), k);
        [f, x] = ksdensity(stats{i}.samples(:,k)');
        mplot(i) = plot(x, f, 'Color', stats{i}.color, 'LineStyle', ls{i}, 'Marker', markers{i});
        xlabel(['$\log_{10}(' param_names_plot{k} ')$'], 'Interpreter', 'latex');
    end    
    h = gca;
    line([log10(true_params(k)) log10(true_params(k))], h.YLim, 'Color', true_color, 'LineStyle', '-');
    axis tight;
    legend(mplot, methods, 'Location', 'northoutside');
end
print([model '_marginals.eps'], '-depsc');
%% Plot the ACFs
figure(2);
% Set size for display
defpos = get(0,'defaultFigurePosition');
set(2,'Position', [defpos(1) defpos(2) width*100*ceil(params_dim/2), height*100*3]);

% Set the size for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width*ceil(params_dim/2), height*3];
set(2, 'PaperPosition', defsize);

for k = 1:params_dim
    for i = 1:length(methods)
        hold on;
        subplot(2, ceil(params_dim/2), k);
        y = acf(stats{i}.samples(:,k), 500);
        mplot(i) = plot(y, 'Color', stats{i}.color, 'LineStyle', ls{i}, 'Marker', markers{i});
        xlabel('Lag');
        ylabel('Autocorrelation');
        title(['$\log_{10}(' param_names_plot{k} ')$'], 'Interpreter', 'latex');
    end        
    axis tight;
    legend(mplot, methods, 'Location', 'northoutside');
end
print([model '_acfs.eps'], '-depsc');
%% Plot the relative error in the DAMH log-posterior
figure(4);
rel_error = abs(damh_result.L_approx - damh_result.L)./abs(damh_result.L);
histogram(rel_error, 'Normalization', 'probability');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('Relative error in log-likelihood approximation');
ylabel('Fraction');
% axis tight;
print([model '_lpos_error.eps'], '-depsc');
%% Plot where the reduced basis was taken
figure(5);
ip = [1 2];
scatter(damh_result.samples(ip(1),:), damh_result.samples(ip(2),:), 'b'); hold on;
scatter(damh_perf.samples_update(ip(1),:), damh_perf.samples_update(ip(2),:), 'r', 'filled');
xlabel(['$\log_{10}(' param_names_plot{ip(1)} ' )$'], 'Interpreter', 'latex');
ylabel(['$\log_{10}(' param_names_plot{ip(2)} ' )$'], 'Interpreter', 'latex');
leg = legend('Accepted', 'Reduced basis build');
set(leg, 'Location', 'northoutside');
print([model '_rb_points.eps'], '-depsc');
%% Plot the cpu time
figure(6);
for i = 1:length(methods)
    hold on;
    plot(stats{i}.cumcput, 'Color', stats{i}.color, 'LineStyle', ls{i}, 'Marker', markers{i});
end
leg = legend(methods);
set(leg, 'Location', 'northwest');
xlabel('Iteration');
ylabel('CPU time (sec)');
axis tight;
print([model '_cpu_time.eps'], '-depsc');
%% Scatter plot for the log-posterior true values and approximations
figure(7);
scatter(damh_result.L, damh_result.L_approx, 2);
xlabel('Full FSP');
ylabel('Reduced model');
axis tight;
R = corrcoef(damh_result.L, damh_result.L_approx);
leg = legend(['corr = ' num2str(R(1,2))]);
set(leg, 'Location', 'northwest');
print([model '_lpos.eps'], '-depsc');
%%
flag = 0;
end
