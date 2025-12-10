%% Load and Visualize Optimization Results
% This script loads saved optimization results and provides various
% post-processing capabilities including:
% - Summary statistics
% - Regenerate all plots
% - Animate any optimization run
% - Compare different runs
%
% Usage:
%   1. Run this script directly (uses default 'optimization_results.mat')
%   2. Or call functions directly:
%      - display_summary('optimization_results.mat')
%      - animate_run('optimization_results.mat', run_idx)
%      - regenerate_plots('optimization_results.mat')
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% Configuration
results_file = 'optimization_results.mat';

%% Load Results
fprintf('=============================================================\n');
fprintf('  POST-PROCESSING: Optimization Results Viewer\n');
fprintf('=============================================================\n\n');

if ~exist(results_file, 'file')
    error('Results file not found: %s\nRun stress_optimization_auxiliary_holes.m first.', results_file);
end

fprintf('Loading results from: %s\n', results_file);
data = load(results_file);
all_results = data.all_results;
fprintf('Results loaded successfully.\n\n');

%% Display Summary
display_summary(all_results);

%% Interactive Menu
while true
    fprintf('\n--- POST-PROCESSING MENU ---\n');
    fprintf('1. Display summary statistics\n');
    fprintf('2. Regenerate all plots\n');
    fprintf('3. Animate best run\n');
    fprintf('4. Animate specific run\n');
    fprintf('5. Compare all runs\n');
    fprintf('6. Export optimal design parameters\n');
    fprintf('7. Re-run stress analysis on optimal design\n');
    fprintf('8. Save animation as GIF\n');
    fprintf('0. Exit\n');
    fprintf('\n');

    choice = input('Enter choice (0-8): ');

    switch choice
        case 0
            fprintf('Exiting post-processor.\n');
            break;
        case 1
            display_summary(all_results);
        case 2
            regenerate_plots(all_results);
        case 3
            animate_run(all_results, all_results.global_best_idx);
        case 4
            fprintf('Available runs: 1 to %d\n', length(all_results.runs));
            run_idx = input('Enter run index to animate: ');
            if run_idx >= 1 && run_idx <= length(all_results.runs)
                animate_run(all_results, run_idx);
            else
                fprintf('Invalid run index.\n');
            end
        case 5
            compare_all_runs(all_results);
        case 6
            export_optimal_design(all_results);
        case 7
            rerun_stress_analysis(all_results);
        case 8
            fprintf('Available runs: 1 to %d (best run: %d)\n', length(all_results.runs), all_results.global_best_idx);
            run_idx = input('Enter run index to save as GIF (0 for best): ');
            if run_idx == 0
                run_idx = all_results.global_best_idx;
            end
            if run_idx >= 1 && run_idx <= length(all_results.runs)
                save_animation_gif(all_results, run_idx);
            else
                fprintf('Invalid run index.\n');
            end
        otherwise
            fprintf('Invalid choice. Please enter 0-8.\n');
    end
end

%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function display_summary(all_results)
%DISPLAY_SUMMARY Display optimization summary statistics

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf('                    OPTIMIZATION SUMMARY\n');
    fprintf('=============================================================\n');
    fprintf('Timestamp: %s\n', datestr(all_results.timestamp));
    fprintf('\n');

    % Optimal design (4 variables: x_aux, a, b, theta_R)
    x_opt = all_results.global_best_x;
    fprintf('OPTIMAL DESIGN (4 variables: y_aux=0, mirror symmetric):\n');
    fprintf('  x_aux:   %.3f mm\n', x_opt(1));
    fprintf('  y_aux:   0 mm (fixed at centerline)\n');
    fprintf('  a:       %.3f mm\n', x_opt(2));
    fprintf('  b:       %.3f mm\n', x_opt(3));
    fprintf('  theta:   %.3f rad (%.1f deg) [both holes, visually mirrored]\n', x_opt(4), rad2deg(x_opt(4)));
    fprintf('  Aspect ratio: %.2f\n', x_opt(2)/x_opt(3));
    fprintf('\n');

    % Performance metrics
    summary = all_results.summary;
    fprintf('PERFORMANCE:\n');
    fprintf('  Baseline max stress:     %.2f MPa (Kt = %.3f)\n', summary.baseline_stress, summary.baseline_Kt);
    fprintf('  Optimal max stress:      %.2f MPa (Kt = %.3f)\n', summary.optimal_stress, summary.optimal_Kt);
    fprintf('  Stress reduction:        %.1f%%\n', summary.stress_reduction_percent);
    fprintf('\n');

    % Computation statistics
    fprintf('COMPUTATION:\n');
    fprintf('  Number of multi-starts:  %d\n', summary.num_starts);
    fprintf('  Best run (start point):  %d\n', all_results.global_best_idx);
    fprintf('  Total function evals:    %d\n', summary.total_function_evals);
    fprintf('  Total optimization time: %.1f seconds\n', summary.total_time_seconds);
    fprintf('\n');

    % Run statistics
    fvals = zeros(length(all_results.runs), 1);
    for i = 1:length(all_results.runs)
        fvals(i) = all_results.runs{i}.fval;
    end
    fprintf('RUN STATISTICS:\n');
    fprintf('  Best objective:   %.2f MPa\n', min(fvals));
    fprintf('  Worst objective:  %.2f MPa\n', max(fvals));
    fprintf('  Mean objective:   %.2f MPa\n', mean(fvals));
    fprintf('  Std deviation:    %.2f MPa\n', std(fvals));
    fprintf('=============================================================\n');
end

function regenerate_plots(all_results)
%REGENERATE_PLOTS Regenerate all visualization plots from saved results

    fprintf('Regenerating plots...\n');

    params = all_results.params;
    baseline = all_results.baseline;
    x_opt = all_results.global_best_x;

    % Expand optimal design to 6-var format
    x_opt_full = [x_opt(1); 0; x_opt(2); x_opt(3); x_opt(4); x_opt(4)];

    % --- Figure 1: Geometry Comparison ---
    figure('Name', 'Geometry Comparison', 'Position', [50, 50, 1400, 500]);

    subplot(1, 3, 1);
    [g_base, ~] = create_baseline_geometry(params);
    pdegplot(g_base, 'FaceLabels', 'off');
    axis equal;
    axis([-1.1*params.plate.length, 1.1*params.plate.length, ...
          -1.1*params.plate.width, 1.1*params.plate.width]);
    xlabel('x (mm)'); ylabel('y (mm)');
    title('Baseline (Central Hole Only)');
    grid on;

    % Initial guess from best run
    x0_best = all_results.runs{all_results.global_best_idx}.x0;
    x0_full = [x0_best(1); 0; x0_best(2); x0_best(3); x0_best(4); x0_best(4)];

    subplot(1, 3, 2);
    try
        [g_init, ~] = create_geometry(x0_full, params);
        pdegplot(g_init, 'FaceLabels', 'off');
        hold on;
        plot(x0_best(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(-x0_best(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        hold off;
    catch
        text(0.5, 0.5, 'Geometry creation failed', 'HorizontalAlignment', 'center');
    end
    axis equal;
    axis([-1.1*params.plate.length, 1.1*params.plate.length, ...
          -1.1*params.plate.width, 1.1*params.plate.width]);
    xlabel('x (mm)'); ylabel('y (mm)');
    title('Initial Guess');
    grid on;

    subplot(1, 3, 3);
    try
        [g_opt, ~] = create_geometry(x_opt_full, params);
        pdegplot(g_opt, 'FaceLabels', 'off');
        hold on;
        plot(x_opt(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(-x_opt(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        % Draw orientation arrows
        theta_R = x_opt(4);
        a = x_opt(2);
        quiver(x_opt(1), 0, a*cos(theta_R)*0.8, a*sin(theta_R)*0.8, 'r', 'LineWidth', 2, 'AutoScale', 'off');
        quiver(-x_opt(1), 0, -a*cos(-theta_R)*0.8, a*sin(-theta_R)*0.8, 'b', 'LineWidth', 2, 'AutoScale', 'off');
        hold off;
    catch
        text(0.5, 0.5, 'Geometry creation failed', 'HorizontalAlignment', 'center');
    end
    axis equal;
    axis([-1.1*params.plate.length, 1.1*params.plate.length, ...
          -1.1*params.plate.width, 1.1*params.plate.width]);
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('Optimal (%.1f%% reduction)', all_results.summary.stress_reduction_percent));
    grid on;

    % --- Figure 2: Convergence History ---
    figure('Name', 'Convergence History', 'Position', [50, 50, 1200, 500]);

    subplot(1, 2, 1);
    hold on;
    for i = 1:length(all_results.runs)
        if ~isempty(all_results.runs{i}.history.fval)
            hist_fval = all_results.runs{i}.history.fval;
            if i == all_results.global_best_idx
                plot(1:length(hist_fval), hist_fval, 'b-', 'LineWidth', 2);
            else
                plot(1:length(hist_fval), hist_fval, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            end
        end
    end
    yline(baseline.max_stress, 'r--', 'LineWidth', 1.5, 'Label', 'Baseline');
    hold off;
    xlabel('Function Evaluation');
    ylabel('Max von Mises Stress (MPa)');
    title('Convergence History (All Runs)');
    legend('Best Run', 'Other Runs', 'Baseline', 'Location', 'northeast');
    grid on;

    subplot(1, 2, 2);
    best_hist = all_results.runs{all_results.global_best_idx}.history;
    if ~isempty(best_hist.fval)
        plot(1:length(best_hist.fval), best_hist.fval, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        yline(baseline.max_stress, 'r--', 'LineWidth', 1.5, 'Label', 'Baseline');
        [min_fval, min_idx] = min(best_hist.fval);
        plot(min_idx, min_fval, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
        hold off;
    end
    xlabel('Function Evaluation');
    ylabel('Max von Mises Stress (MPa)');
    title(sprintf('Best Run Convergence (Start %d)', all_results.global_best_idx));
    legend('Objective', 'Baseline', 'Optimal', 'Location', 'northeast');
    grid on;

    % --- Figure 3: Final Objective Distribution ---
    figure('Name', 'Final Objective Distribution', 'Position', [50, 50, 800, 400]);

    fvals = zeros(length(all_results.runs), 1);
    for i = 1:length(all_results.runs)
        fvals(i) = all_results.runs{i}.fval;
    end

    subplot(1, 2, 1);
    bar(fvals);
    hold on;
    bar(all_results.global_best_idx, fvals(all_results.global_best_idx), 'g');
    yline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Run Index');
    ylabel('Final Objective (MPa)');
    title('Final Objective by Run');
    legend('Other Runs', 'Best Run', 'Baseline', 'Location', 'northeast');
    grid on;

    subplot(1, 2, 2);
    histogram(fvals, 10);
    hold on;
    xline(all_results.global_best_fval, 'g-', 'LineWidth', 2);
    xline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Final Objective (MPa)');
    ylabel('Count');
    title('Objective Distribution');
    legend('Distribution', 'Best', 'Baseline', 'Location', 'northeast');
    grid on;

    fprintf('Plots regenerated.\n');
end

function animate_run(all_results, run_idx)
%ANIMATE_RUN Animate a specific optimization run

    if run_idx < 1 || run_idx > length(all_results.runs)
        fprintf('Invalid run index: %d\n', run_idx);
        return;
    end

    fprintf('Animating run %d...\n', run_idx);

    run_result = all_results.runs{run_idx};
    params = all_results.params;
    baseline = all_results.baseline;

    history = run_result.history;
    if isempty(history.x) || size(history.x, 1) < 2
        fprintf('Insufficient history for animation in run %d.\n', run_idx);
        return;
    end

    n_frames = min(30, size(history.x, 1));
    frame_indices = round(linspace(1, size(history.x, 1), n_frames));

    fig = figure('Name', sprintf('Design Evolution - Run %d', run_idx), ...
                 'Position', [50, 50, 1200, 500]);

    for frame = 1:length(frame_indices)
        idx = frame_indices(frame);
        x_frame = history.x(idx, :)';
        fval_frame = history.fval(idx);

        % Expand 4-var to 6-var: [x_aux, y_aux=0, a, b, theta_R, theta_L=theta_R]
        x_full = [x_frame(1); 0; x_frame(2); x_frame(3); x_frame(4); x_frame(4)];

        clf(fig);

        subplot(1, 2, 1);
        try
            [g_frame, ~] = create_geometry(x_full, params);
            pdegplot(g_frame, 'FaceLabels', 'off');
            hold on;
            plot(x_frame(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            plot(-x_frame(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
            hold off;
        catch
            text(0.5, 0.5, 'Geometry failed', 'HorizontalAlignment', 'center');
        end
        axis equal;
        axis([-1.1*params.plate.length, 1.1*params.plate.length, ...
              -1.1*params.plate.width, 1.1*params.plate.width]);
        xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('Design at Iteration %d', idx));
        grid on;

        % Display current design variables
        text(-params.plate.length*0.9, -params.plate.width*0.8, ...
             sprintf('x_{aux}=%.1f, a=%.1f, b=%.1f, \\theta=%.1f°', ...
                     x_frame(1), x_frame(2), x_frame(3), rad2deg(x_frame(4))), ...
             'FontSize', 9);

        subplot(1, 2, 2);
        plot(1:idx, history.fval(1:idx), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        yline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
        plot(idx, fval_frame, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        xlabel('Function Evaluation');
        ylabel('Max von Mises Stress (MPa)');
        title(sprintf('\\sigma_{max} = %.1f MPa', fval_frame));
        xlim([0, size(history.x, 1) + 1]);
        ylim([min(history.fval)*0.95, max([max(history.fval), baseline.max_stress])*1.05]);
        legend('Objective', 'Baseline', 'Current', 'Location', 'northeast');
        grid on;

        sgtitle(sprintf('Run %d - Frame %d/%d', run_idx, frame, n_frames));
        drawnow;
        pause(0.15);
    end

    fprintf('Animation complete.\n');
end

function save_animation_gif(all_results, run_idx)
%SAVE_ANIMATION_GIF Save optimization run animation as GIF file

    if run_idx < 1 || run_idx > length(all_results.runs)
        fprintf('Invalid run index: %d\n', run_idx);
        return;
    end

    run_result = all_results.runs{run_idx};
    params = all_results.params;
    baseline = all_results.baseline;

    history = run_result.history;
    if isempty(history.x) || size(history.x, 1) < 2
        fprintf('Insufficient history for animation in run %d.\n', run_idx);
        return;
    end

    % GIF settings
    gif_filename = sprintf('optimization_animation_run%d.gif', run_idx);
    delay_time = 0.2;  % seconds between frames
    n_frames = min(40, size(history.x, 1));
    frame_indices = round(linspace(1, size(history.x, 1), n_frames));

    fprintf('Saving animation to %s (%d frames)...\n', gif_filename, n_frames);

    fig = figure('Name', sprintf('Design Evolution - Run %d', run_idx), ...
                 'Position', [50, 50, 1200, 500], 'Color', 'white');

    for frame = 1:length(frame_indices)
        idx = frame_indices(frame);
        x_frame = history.x(idx, :)';
        fval_frame = history.fval(idx);

        % Expand 4-var to 6-var
        x_full = [x_frame(1); 0; x_frame(2); x_frame(3); x_frame(4); x_frame(4)];

        clf(fig);

        % Left subplot: Geometry
        subplot(1, 2, 1);
        try
            [g_frame, ~] = create_geometry(x_full, params);
            pdegplot(g_frame, 'FaceLabels', 'off');
            hold on;
            plot(x_frame(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            plot(-x_frame(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
            hold off;
        catch
            text(0.5, 0.5, 'Geometry failed', 'HorizontalAlignment', 'center');
        end
        axis equal;
        axis([-1.1*params.plate.length, 1.1*params.plate.length, ...
              -1.1*params.plate.width, 1.1*params.plate.width]);
        xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('Design at Iteration %d', idx));
        grid on;

        % Display current design variables
        text(-params.plate.length*0.9, -params.plate.width*0.8, ...
             sprintf('x_{aux}=%.1f, a=%.1f, b=%.1f, \\theta=%.1f°', ...
                     x_frame(1), x_frame(2), x_frame(3), rad2deg(x_frame(4))), ...
             'FontSize', 9);

        % Right subplot: Convergence
        subplot(1, 2, 2);
        plot(1:idx, history.fval(1:idx), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        yline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
        plot(idx, fval_frame, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        xlabel('Function Evaluation');
        ylabel('Max von Mises Stress (MPa)');
        title(sprintf('\\sigma_{max} = %.1f MPa', fval_frame));
        xlim([0, size(history.x, 1) + 1]);
        ylim([min(history.fval)*0.95, max([max(history.fval), baseline.max_stress])*1.05]);
        legend('Objective', 'Baseline', 'Current', 'Location', 'northeast');
        grid on;

        sgtitle(sprintf('Run %d - Frame %d/%d', run_idx, frame, n_frames));
        drawnow;

        % Capture frame for GIF
        frame_data = getframe(fig);
        im = frame2im(frame_data);
        [imind, cm] = rgb2ind(im, 256);

        % Write to GIF
        if frame == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end

        fprintf('  Frame %d/%d saved\n', frame, n_frames);
    end

    close(fig);
    fprintf('Animation saved to: %s\n', gif_filename);
end

function compare_all_runs(all_results)
%COMPARE_ALL_RUNS Display detailed comparison of all optimization runs

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf('                    ALL RUNS COMPARISON\n');
    fprintf('=============================================================\n');
    fprintf('%-5s | %-10s | %-10s | %-10s | %-8s | %-10s | %s\n', ...
            'Run', 'Objective', 'x_aux', 'a', 'b', 'theta_R', 'ExitFlag');
    fprintf('%s\n', repmat('-', 1, 75));

    for i = 1:length(all_results.runs)
        run = all_results.runs{i};
        x = run.x_opt;
        marker = '';
        if i == all_results.global_best_idx
            marker = ' *BEST*';
        end
        fprintf('%-5d | %-10.2f | %-10.2f | %-10.2f | %-8.2f | %-10.1f | %d%s\n', ...
                i, run.fval, x(1), x(2), x(3), rad2deg(x(4)), run.exitflag, marker);
    end

    fprintf('%s\n', repmat('-', 1, 75));
    fprintf('Baseline stress: %.2f MPa\n', all_results.baseline.max_stress);
    fprintf('=============================================================\n');
end

function export_optimal_design(all_results)
%EXPORT_OPTIMAL_DESIGN Export optimal design parameters to a text file

    x_opt = all_results.global_best_x;
    summary = all_results.summary;

    filename = 'optimal_design_parameters.txt';
    fid = fopen(filename, 'w');

    fprintf(fid, 'OPTIMAL AUXILIARY HOLE DESIGN PARAMETERS\n');
    fprintf(fid, '=========================================\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));

    fprintf(fid, 'DESIGN VARIABLES (4 total):\n');
    fprintf(fid, '  x_aux   = %.4f mm (distance from center)\n', x_opt(1));
    fprintf(fid, '  y_aux   = 0 mm (fixed at centerline)\n');
    fprintf(fid, '  a       = %.4f mm (semi-major axis)\n', x_opt(2));
    fprintf(fid, '  b       = %.4f mm (semi-minor axis)\n', x_opt(3));
    fprintf(fid, '  theta_R = %.4f rad = %.2f deg\n', x_opt(4), rad2deg(x_opt(4)));
    fprintf(fid, '  theta_L = %.4f rad = %.2f deg (= theta_R, visually mirrored)\n\n', x_opt(4), rad2deg(x_opt(4)));

    fprintf(fid, 'DERIVED QUANTITIES:\n');
    fprintf(fid, '  Aspect ratio (a/b): %.3f\n', x_opt(2)/x_opt(3));
    fprintf(fid, '  Ellipse area: %.3f mm^2\n\n', pi * x_opt(2) * x_opt(3));

    fprintf(fid, 'PERFORMANCE:\n');
    fprintf(fid, '  Baseline stress:  %.4f MPa\n', summary.baseline_stress);
    fprintf(fid, '  Optimal stress:   %.4f MPa\n', summary.optimal_stress);
    fprintf(fid, '  Stress reduction: %.2f%%\n\n', summary.stress_reduction_percent);

    fprintf(fid, 'PLATE GEOMETRY:\n');
    params = all_results.params;
    fprintf(fid, '  Plate half-width:  %.1f mm\n', params.plate.width);
    fprintf(fid, '  Plate half-length: %.1f mm\n', params.plate.length);
    fprintf(fid, '  Central hole radius: %.1f mm\n', params.centralHole.radius);

    fclose(fid);

    fprintf('Optimal design exported to: %s\n', filename);
end

function rerun_stress_analysis(all_results)
%RERUN_STRESS_ANALYSIS Re-run FEA on optimal design with detailed output

    fprintf('Re-running stress analysis on optimal design...\n');

    params = all_results.params;
    x_opt = all_results.global_best_x;

    % Expand to 6-var format
    x_opt_full = [x_opt(1); 0; x_opt(2); x_opt(3); x_opt(4); x_opt(4)];

    % Run analysis
    result = run_detailed_analysis_internal(x_opt_full, params);

    fprintf('\n');
    fprintf('DETAILED STRESS ANALYSIS RESULTS:\n');
    fprintf('  Max von Mises stress: %.4f MPa\n', result.max_stress);
    fprintf('  Location of max stress: (%.2f, %.2f) mm\n', result.loc_max(1), result.loc_max(2));
    fprintf('  Stress concentration factor Kt: %.4f\n', result.Kt);
    fprintf('  Mesh nodes: %d\n', size(result.mesh.Nodes, 2));
    fprintf('  Mesh elements: %d\n', size(result.mesh.Elements, 2));

    % Plot stress distribution
    figure('Name', 'Stress Analysis - Optimal Design', 'Position', [50, 50, 1200, 500]);

    subplot(1, 2, 1);
    pdeplot(result.mesh, 'XYData', result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(result.loc_max(1), result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal; colorbar;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('von Mises Stress (Max = %.1f MPa)', result.max_stress));

    subplot(1, 2, 2);
    pdeplot(result.mesh, 'XYData', result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(result.loc_max(1), result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    xlim([-100, 100]); ylim([-60, 60]);
    colorbar;
    xlabel('x (mm)'); ylabel('y (mm)');
    title('Zoomed View Near Holes');

    fprintf('Analysis complete.\n');
end

%% ========================================================================
%  GEOMETRY FUNCTIONS
%  ========================================================================

function [g, status] = create_baseline_geometry(params)
%CREATE_BASELINE_GEOMETRY Create plate with central hole only

    L = params.plate.length;
    W = params.plate.width;
    R_c = params.centralHole.radius;

    R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];
    C1 = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

    gdm = [R1, C1];
    ns = char('R1', 'C1')';
    [g, ~] = decsg(gdm, 'R1 - C1', ns);
    status = 0;
end

function [g, status] = create_geometry(x, params)
%CREATE_GEOMETRY Create plate with central and auxiliary holes
%   x = [x_aux, y_aux, a, b, theta_R, theta_L] (6 variables)

    L = params.plate.length;
    W = params.plate.width;
    R_c = params.centralHole.radius;
    n_pts = params.mesh.n_ellipse_pts;

    x_aux = x(1); y_aux = x(2);
    a = x(3); b = x(4);
    theta_R = x(5); theta_L = x(6);

    R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];
    C_central = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

    theta_pts = linspace(0, 2*pi, n_pts + 1);
    theta_pts = theta_pts(1:end-1);

    x_local = a * cos(theta_pts);
    y_local = b * sin(theta_pts);

    cos_R = cos(theta_R); sin_R = sin(theta_R);
    x_right = cos_R * x_local - sin_R * y_local + x_aux;
    y_right = sin_R * x_local + cos_R * y_local + y_aux;

    cos_L = cos(theta_L); sin_L = sin(theta_L);
    x_left = -(cos_L * x_local - sin_L * y_local) - x_aux;
    y_left = sin_L * x_local + cos_L * y_local + y_aux;

    P_right = [2; n_pts; x_right(:); y_right(:)];
    P_left = [2; n_pts; x_left(:); y_left(:)];

    max_len = max([length(R1), length(C_central), length(P_right), length(P_left)]);
    R1 = [R1; zeros(max_len - length(R1), 1)];
    C_central = [C_central; zeros(max_len - length(C_central), 1)];
    P_right = [P_right; zeros(max_len - length(P_right), 1)];
    P_left = [P_left; zeros(max_len - length(P_left), 1)];

    gdm = [R1, C_central, P_right, P_left];
    ns = char('R1', 'C1', 'P_right', 'P_left')';
    [g, ~] = decsg(gdm, 'R1 - C1 - P_right - P_left', ns);
    status = 0;
end

function result = run_detailed_analysis_internal(x, params)
%RUN_DETAILED_ANALYSIS_INTERNAL Run FEA with full output

    [g, ~] = create_geometry(x, params);

    model = femodel(AnalysisType="structuralStatic", Geometry=g);
    model.MaterialProperties = materialProperties(...
        YoungsModulus=params.material.E, PoissonsRatio=params.material.nu);
    model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
    model.EdgeBC(4) = edgeBC(XDisplacement=0);
    model.VertexBC(1) = vertexBC(YDisplacement=0);

    min_feature = min([x(3), x(4), params.constraints.min_gap_central]);
    Hmax = min_feature / params.mesh.Hmax_factor;
    model = generateMesh(model, Hmax=Hmax);

    R = solve(model);
    stress = R.Stress;
    vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

    [max_stress, idx] = max(vonMises);
    nodes = R.Mesh.Nodes;

    result.max_stress = max_stress;
    result.Kt = max_stress / params.load.tension;
    result.vonMises = vonMises;
    result.mesh = R.Mesh;
    result.loc_max = nodes(:, idx);
    result.stress = stress;
    result.Displacement = R.Displacement;
    result.R = R;
end
