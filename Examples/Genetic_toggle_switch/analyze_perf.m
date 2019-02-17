function flag = analyze_perf(model, methods, param_names, param_names_plot)
%%
stats = cell(size(methods));
i=1;

load([model '_basic_info.mat']);

if (~exist('unobserved_species', 'var'))
    unobserved_species = [];
end

i_burn_in = 1;
for method = methods
    if (contains(string(method), 'DAMH'))
        load([model '_damh_results.mat'], 'damh_result', 'damh_perf');
        load([model '_approx_model.mat'], 'approx_model');
        samples = damh_result.samples(:,i_burn_in:end)';
        L = damh_result.L(i_burn_in:end);
        cput = max(damh_perf.cput);
        n_full_eval = damh_perf.n_full_eval;
        n_reject = damh_result.rejection;
        n_full_reject = damh_perf.n_second_reject;        
    elseif (contains(string(method), 'AMH-FSP'))
        load([model '_amh_results.mat']);
        samples = full_amh_result.samples(:,i_burn_in:end)';
        L = full_amh_result.L(i_burn_in:end);
        cput = max(full_amh_perf.cput);
        n_full_eval = length(full_amh_result.samples);
        n_reject = full_amh_result.n_reject;
        n_full_reject = n_reject;        
    elseif (contains(string(method), 'Hybrid'))
        load([model '_dahyb_results.mat']);
        samples = [dahyb_result.phase1.samples(:,:)'; dahyb_result.phase2.samples(:,:)'];
        L = dahyb_result.phase2.L(:);
        cput = max(dahyb_perf.phase2.cput) + max(dahyb_perf.phase1.cput);
        n_full_eval = dahyb_perf.phase1.n_full_eval;
        n_reject = dahyb_result.phase1.rejection + dahyb_result.phase2.n_reject;
        n_full_reject = dahyb_perf.phase1.n_second_reject;        
    end
    mESS = multiESS(samples);
    
    nchain = size(samples, 1);
    % use the parameter average
    map_est = mean(samples);
    std_est = std(samples);
    map_L = mean(L);
    [~, geweke_p] = geweke(samples);
    iact_val = iact(samples);
    
    perf = struct('mESS', mESS, 'geweke', geweke_p, 'iact', iact_val, 'cput', cput,...
        'cpu_per_eff_sample', cput/mESS, 'n_full_evals', n_full_eval, ...
        'map_est', map_est, 'std', std_est, 'map_L', map_L, 'n_reject', n_reject, 'n_full_reject', n_full_reject);
    stats{i} = perf;
    i = i + 1;
end
%% Number of basis per interval in reduced model
max_num_basis = 1;
for i = 1:length(approx_model.PHIcell)
    max_num_basis = max([max_num_basis size(approx_model.PHIcell{i}, 2)]);
end
write_file([model '_max_num_basis.tex'], '%d', max_num_basis);
%% Number of samples for basis builds
n_basis_builds = size(damh_perf.samples_update, 2);
write_file([model '_num_model_updates.tex'], '%d', n_basis_builds);
%% Compute savings
second_stage_accept_percent = 100*(damh_perf.n_full_eval - damh_perf.n_second_reject)/damh_perf.n_full_eval;
write_file([model '_second_stage_accept_percent.tex'], '%.2f', second_stage_accept_percent);

dahyb_cput_saving_percent = 100*(max(full_amh_perf.cput) - max(dahyb_perf.phase1.cput) - max(dahyb_perf.phase2.cput))/max(full_amh_perf.cput);
write_file([model '_dahyb_cput_saving_percent.tex'], '%.2f', dahyb_cput_saving_percent);

damh_cput_saving_percent = 100*(max(full_amh_perf.cput) - max(damh_perf.cput))/max(full_amh_perf.cput);
write_file([model '_damh_cput_saving_percent.tex'], '%.2f', damh_cput_saving_percent);

damh_full_eval_saving_percent = 100*(length(damh_result.L) - damh_perf.n_full_eval)/length(damh_result.L);
write_file([model '_damh_full_eval_saving_percent.tex'], '%.2f', damh_full_eval_saving_percent);
%% Compute median and average of the relative errors
rel_error = abs(damh_result.L_approx - damh_result.L)./abs(damh_result.L);
median_rel_error = median(rel_error);
avg_rel_error = mean(rel_error);
write_file([model '_lpos_median_rel_err.tex'], 'publication', median_rel_error);
write_file([model '_lpos_avg_rel_err.tex'], 'publication', avg_rel_error);
%% Figure out the when the last basis update occurs
[lia, locb] = ismember(damh_perf.samples_update', damh_result.proposals', 'rows');
last_update_iter = max(locb);
last_update_portion = where_is_it(last_update_iter, nchain);
write_file([model '_last_update_iter.tex'], '%d', last_update_iter);
write_file([model '_last_update_portion.tex'], '%s', last_update_portion);
%% Make a table for CPU time
fid = fopen([model '_table.tex'], 'w+');
fprintf(fid, '\\begin{tabular}{l  c c c c c c} \n');
fprintf(fid, '\\toprule \n');
fprintf(fid, [' & mESS & \\thead{CPU \\\\ time \\\\ ($\\sec$)} & \\thead{$\\frac{\\text{CPU time}}{\\text{mESS}}$ ($\\sec$)} & \\thead{Number of \\\\ full evaluations} &'...
               ' \\thead{Number of \\\\ rejections} & \\thead{Number of \\\\ rejections by \\\\ full FSP} \\\\ \n']);
fprintf(fid, '\\midrule \n');
for i = 1:length(stats)
    fprintf(fid, ' %s & %.2f & %.2f & %.2f & %d & %d & %d \\\\ \n',...
        string(methods{i}), ...
        stats{i}.mESS,...
        stats{i}.cput,...
        stats{i}.cpu_per_eff_sample,...
        stats{i}.n_full_evals,...
        stats{i}.n_reject,...
        stats{i}. n_full_reject);
end
fprintf(fid, '\\bottomrule \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
%% The table for the prior bounds
n_params = length(true_params);
fid = fopen([model '_prior.tex'], 'w+');
fprintf(fid, ['\\begin{tabular}{l' repmat('c', 1, n_params) '}\n']);
fprintf(fid, '\\toprule \n');
fprintf(fid, 'Parameter ');
for i = 1:n_params
    fprintf(fid, ['& $' param_names{i} '$ ']);
end
fprintf(fid, '\\\\ \n');
fprintf(fid, '\\midrule \n');
fprintf(fid, ['Lower bound' num2str(params_min', '& %1.2e') '\\\\ \n']);
fprintf(fid, ['Upper bound' num2str(params_max', '& %1.2e') '\\\\ \n']);
fprintf(fid, '\\bottomrule \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
%% The table to compare true and identified parameters
fid = fopen([model '_param_id.tex'], 'w+');
fprintf(fid, '\\begin{tabular}{l');
for i = 1:length(methods)
    fprintf(fid, ' c c');
end
fprintf(fid, ' c } \n');

fprintf(fid, '\\toprule \n');

fprintf(fid, 'Parameter ');
for i = 1:length(methods)
    fprintf(fid, [' & \\multicolumn{2}{c}{' methods{i} '} ']);
end
fprintf(fid, ' & True \\\\ \n');

for i = 1:length(methods)
    fprintf(fid, '& mean & std ');
end
fprintf(fid, ' & $ $ \\\\ \n');

fprintf(fid, '\\midrule \n');

for j = 1:n_params
    fprintf(fid, ['$\\log_{10}(' param_names{j} ')$']);
    for i = 1:length(methods)
        fprintf(fid, ['&' num2str(stats{i}.map_est(j), '%1.2e') '&' ...
        num2str(stats{i}.std(j), '%1.2e') ]);
    end

    fprintf(fid, ['&' num2str(log10(true_params(j)), '%1.2e') '\\\\ \n']);
end
fprintf(fid,'\\bottomrule  \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
%% The table of geweke p-value and IACT
fid = fopen([model '_diagnostics.tex'], 'w+');
fprintf(fid, '\\begin{tabular}{l');
for i = 1:length(methods)
    fprintf(fid, ' c c ');
end
fprintf(fid, '} \n');

fprintf(fid, '\\toprule \n');

fprintf(fid, 'Parameter ');
for i = 1:length(methods)
    fprintf(fid, [' & \\multicolumn{2}{c}{' methods{i} '} ']);
end
fprintf(fid, '\\\\ \n');

for i = 1:length(methods)
    fprintf(fid, '& Geweke & IACT ');
end
fprintf(fid, ' \\\\ \n');

fprintf(fid, '\\midrule \n');

for j = 1:n_params
    fprintf(fid, ['$\\log_{10}(' param_names{j} ')$']);
    for i = 1:length(methods)
        fprintf(fid, ['&' num2str(stats{i}.geweke(j), '%1.2e') '&' num2str(stats{i}.iact(j), '%1.2e') ]);
    end

    fprintf(fid, '\\\\ \n');
end
fprintf(fid,'\\bottomrule  \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
%% CPU time break down for ADAMH
total_time_fsp = sum(damh_perf.full_eval_time);
total_time_rom = sum(damh_perf.approx_eval_time);
total_time_update = sum(damh_perf.approx_update_time);
total_time = max(damh_perf.cput);

fid = fopen([model '_time_breakdown.tex'], 'w+');
fprintf(fid, '\\begin{tabular}{l c c} \n');
fprintf(fid, '\\toprule \n');
fprintf(fid, 'Component & Time occupied ($\\sec$) & Fraction of total time (per cent) \\\\ \n');
fprintf(fid, '\\midrule \n');
fprintf(fid, '%s & %.2f & %.2f \\\\', 'Full FSP Evaluation', total_time_fsp, 100*total_time_fsp/total_time);
fprintf(fid, '%s & %.2f & %.2f \\\\', 'Reduced Model Evaluation', total_time_rom, 100*total_time_rom/total_time);
fprintf(fid, '%s & %.2f & %.2f \\\\', 'Reduced Model Update', total_time_update, 100*total_time_update/total_time);
fprintf(fid, '%s & %.2f & %.2f \\\\', 'Total', total_time, 100);
fprintf(fid, '\\bottomrule \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

%% Average solve time of full and reduced model

% process the chain results
accept = damh_result.samples(:, 2:end)==damh_result.proposals(:,1:end-1);
accept = prod(accept, 1);
i_accept = find(accept);

% solve time at the accepted parameters and their frequencies
full_solve_time = damh_perf.full_eval_time(i_accept);
approx_solve_time = damh_perf.approx_eval_time(i_accept);
freq = [diff(i_accept) ( length(accept) - i_accept(end)+1)];

avg_fsp_time = dot(full_solve_time, (freq/sum(freq)));
avg_rom_time = dot(approx_solve_time, (freq/sum(freq)));
%%
fid = fopen([model '_avg_solver_time.tex'], 'w+');
fprintf(fid, '\\begin{tabular}{|l|c|} \n');
fprintf(fid, '\\hline \n');
fprintf(fid, 'Model & Time ($\\sec$) \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid, '%s & %.2f \\\\', 'Full FSP', avg_fsp_time);
fprintf(fid, '%s & %.2f \\\\', 'Reduced Model', avg_rom_time);
fprintf(fid, '\\hline \n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
%%
write_file([model '_avg_full_solve_cput.tex'], '%.2f', avg_fsp_time);
write_file([model '_avg_rom_solve_cput.tex'], '%.2f', avg_rom_time);
%% Average speedup (very rough estimate)
avg_speedup_factor = avg_fsp_time/avg_rom_time;
avg_speedup_percentage = 100*(avg_fsp_time - avg_rom_time)/avg_fsp_time;
fid = fopen([model '_avg_speedup_factor.tex'], 'w+');
fprintf(fid, '%.2f', avg_speedup_factor);
fclose(fid);
fid = fopen([model '_avg_speedup_percentage.tex'], 'w+');
fprintf(fid, '%.2f', avg_speedup_percentage);
fclose(fid);
%%
close all
flag = 0;
end

function write_file(filename, format, value)
    fid = fopen(filename, 'w+');
    if (strcmp('publication', format))
        exponent = floor(log10(value));
        digits = value/(10^exponent);
        fprintf(fid, "%.2f \\times 10^{ %d }", digits, exponent);
    elseif (strcmp('%d',format) && (0<= value && value <= 10))
        format = '%s';
        switch value
            case 0
                value = "zero";
            case 1
                value = "one";
            case 2
                value = "two";
            case 3
                value = "three";
            case 4
                value = "four";
            case 5
                value = "five";
            case 6
                value = "six";
            case 7 
                value = "seven";
            case 8
                value = "eight";
            case 9
                value = "nine";
            case 10
                value = "ten";
        end
        fprintf(fid, format, value);
    else
        fprintf(fid, format, value);  
    end
    fclose(fid);    
end

function s = where_is_it(m, n)
% check whether m belong to the first tenth, fifth, third or half of the
% sequence 1, 2, .., n
r = m/n;
if (r <= 0.1)
    s = "first tenth";
    return;
elseif (r <= 0.2)
    s = "first fifth";
    return;
elseif (r <= 1/3)
    s = "first third";
    return;
elseif (r <= 0.5)
    s = "first half";
    return;
else
    s = "second half";
    return;
end
end