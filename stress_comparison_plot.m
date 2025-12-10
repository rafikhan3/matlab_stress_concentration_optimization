%% Stress Comparison Plot: Initial vs Optimal Design
% This standalone script generates a side-by-side von Mises stress comparison
% between the initial design parameters and the optimal design parameters.
%
% The plot includes:
%   - Baseline (central hole only)
%   - Initial design (starting guess from best optimization run)
%   - Optimal design (converged solution)
%
% Usage:
%   1. Run this script after completing optimization
%   2. Requires: optimization_results.mat (from stress_optimization_auxiliary_holes.m)
%
% Output:
%   - Figure with 3 subplots showing von Mises stress distribution
%   - Optionally exports to PNG for presentation use
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% Configuration
results_file = 'optimization_results.mat';
export_figure = true;  % Set to true to auto-export PNG
export_filename = 'stress_comparison.png';
export_resolution = 300;  % DPI

%% Load Results
fprintf('=============================================================\n');
fprintf('  STRESS COMPARISON: Initial vs Optimal Design\n');
fprintf('=============================================================\n\n');

if ~exist(results_file, 'file')
    error('Results file not found: %s\nRun stress_optimization_auxiliary_holes.m first.', results_file);
end

fprintf('Loading results from: %s\n', results_file);
data = load(results_file);
all_results = data.all_results;
fprintf('Results loaded successfully.\n\n');

%% Extract Data
params = all_results.params;
baseline = all_results.baseline;
best_idx = all_results.global_best_idx;

% Check if stored results contain mesh/stress data
has_stored_results = isfield(all_results, 'optimal_result') && ...
                     isfield(all_results.optimal_result, 'mesh') && ...
                     isfield(all_results.optimal_result, 'vonMises');

if has_stored_results
    fprintf('Using stored stress results from optimization run.\n');
    optimal_result = all_results.optimal_result;

    % Check for initial result
    if isfield(all_results, 'initial_result') && ...
       isfield(all_results.initial_result, 'mesh')
        initial_result = all_results.initial_result;
        has_initial = true;
    else
        has_initial = false;
    end
else
    fprintf('Stored results incomplete. Re-running stress analysis...\n');
    has_initial = false;

    % Re-run analysis for optimal design
    x_opt = all_results.global_best_x;
    x_opt_full = [x_opt(1); 0; x_opt(2); x_opt(3); x_opt(4); x_opt(4)];
    optimal_result = run_analysis_with_aux_holes(x_opt_full, params);
end

%% Run Initial Design Analysis (if needed)
x0_best = all_results.runs{best_idx}.x0;
x0_full = [x0_best(1); 0; x0_best(2); x0_best(3); x0_best(4); x0_best(4)];

if ~has_initial
    fprintf('Running stress analysis for initial design...\n');
    try
        initial_result = run_analysis_with_aux_holes(x0_full, params);
        has_initial = true;
    catch ME
        fprintf('Warning: Could not analyze initial design: %s\n', ME.message);
        has_initial = false;
    end
end

%% Display Summary
fprintf('\n--- STRESS SUMMARY ---\n');
fprintf('Baseline (central hole only):\n');
fprintf('  Max von Mises stress: %.2f MPa\n', baseline.max_stress);
fprintf('  Stress concentration factor Kt: %.3f\n', baseline.Kt);

if has_initial
    fprintf('\nInitial Design:\n');
    fprintf('  x_aux = %.2f mm, a = %.2f mm, b = %.2f mm, theta = %.1f°\n', ...
        x0_best(1), x0_best(2), x0_best(3), rad2deg(x0_best(4)));
    fprintf('  Max von Mises stress: %.2f MPa\n', initial_result.max_stress);
end

fprintf('\nOptimal Design:\n');
x_opt = all_results.global_best_x;
fprintf('  x_aux = %.2f mm, a = %.2f mm, b = %.2f mm, theta = %.1f°\n', ...
    x_opt(1), x_opt(2), x_opt(3), rad2deg(x_opt(4)));
fprintf('  Max von Mises stress: %.2f MPa\n', optimal_result.max_stress);
fprintf('  Stress reduction: %.1f%%\n', all_results.summary.stress_reduction_percent);

%% Create Stress Comparison Figure
fprintf('\n--- GENERATING STRESS COMPARISON PLOT ---\n');

% Determine common color scale
max_stress_all = baseline.max_stress;  % Use baseline as reference for color scale

if has_initial
    fig = figure('Name', 'Von Mises Stress Comparison', 'Position', [50, 50, 1600, 500]);

    % --- Subplot 1: Baseline ---
    subplot(1, 3, 1);
    pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(baseline.loc_max(1), baseline.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Baseline\n\\sigma_{vM,max} = %.1f MPa (K_t = %.2f)', ...
        baseline.max_stress, baseline.Kt), 'FontSize', 12);

    % --- Subplot 2: Initial Design ---
    subplot(1, 3, 2);
    pdeplot(initial_result.mesh, 'XYData', initial_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(initial_result.loc_max(1), initial_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Initial Design\n\\sigma_{vM,max} = %.1f MPa', ...
        initial_result.max_stress), 'FontSize', 12);

    % --- Subplot 3: Optimal Design ---
    subplot(1, 3, 3);
    pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Optimal Design\n\\sigma_{vM,max} = %.1f MPa (%.1f%% reduction)', ...
        optimal_result.max_stress, all_results.summary.stress_reduction_percent), 'FontSize', 12);

    sgtitle('Von Mises Stress Comparison: Baseline → Initial → Optimal', 'FontSize', 14, 'FontWeight', 'bold');

else
    % Only baseline and optimal (no initial)
    fig = figure('Name', 'Von Mises Stress Comparison', 'Position', [50, 50, 1200, 500]);

    % --- Subplot 1: Baseline ---
    subplot(1, 2, 1);
    pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(baseline.loc_max(1), baseline.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Baseline\n\\sigma_{vM,max} = %.1f MPa (K_t = %.2f)', ...
        baseline.max_stress, baseline.Kt), 'FontSize', 12);

    % --- Subplot 2: Optimal Design ---
    subplot(1, 2, 2);
    pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal;
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Optimal Design\n\\sigma_{vM,max} = %.1f MPa (%.1f%% reduction)', ...
        optimal_result.max_stress, all_results.summary.stress_reduction_percent), 'FontSize', 12);

    sgtitle('Von Mises Stress Comparison: Baseline → Optimal', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Create Zoomed Comparison (Near Holes)
fig2 = figure('Name', 'Von Mises Stress Comparison (Zoomed)', 'Position', [50, 100, 1600, 500]);

zoom_xlim = [-100, 100];
zoom_ylim = [-60, 60];

if has_initial
    % --- Subplot 1: Baseline (Zoomed) ---
    subplot(1, 3, 1);
    pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(baseline.loc_max(1), baseline.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    % Draw central hole outline
    theta = linspace(0, 2*pi, 100);
    plot(params.centralHole.radius*cos(theta), params.centralHole.radius*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    axis equal;
    xlim(zoom_xlim); ylim(zoom_ylim);
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Baseline (Zoomed)\n\\sigma_{max} = %.1f MPa', baseline.max_stress), 'FontSize', 12);

    % --- Subplot 2: Initial Design (Zoomed) ---
    subplot(1, 3, 2);
    pdeplot(initial_result.mesh, 'XYData', initial_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(initial_result.loc_max(1), initial_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    % Draw hole outlines
    plot(params.centralHole.radius*cos(theta), params.centralHole.radius*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    axis equal;
    xlim(zoom_xlim); ylim(zoom_ylim);
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Initial Design (Zoomed)\n\\sigma_{max} = %.1f MPa', initial_result.max_stress), 'FontSize', 12);

    % --- Subplot 3: Optimal Design (Zoomed) ---
    subplot(1, 3, 3);
    pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    plot(params.centralHole.radius*cos(theta), params.centralHole.radius*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    axis equal;
    xlim(zoom_xlim); ylim(zoom_ylim);
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Optimal Design (Zoomed)\n\\sigma_{max} = %.1f MPa', optimal_result.max_stress), 'FontSize', 12);

    sgtitle('Zoomed View Near Holes', 'FontSize', 14, 'FontWeight', 'bold');

else
    % Only baseline and optimal
    subplot(1, 2, 1);
    pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(baseline.loc_max(1), baseline.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    theta = linspace(0, 2*pi, 100);
    plot(params.centralHole.radius*cos(theta), params.centralHole.radius*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    axis equal;
    xlim(zoom_xlim); ylim(zoom_ylim);
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Baseline (Zoomed)\n\\sigma_{max} = %.1f MPa', baseline.max_stress), 'FontSize', 12);

    subplot(1, 2, 2);
    pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
    hold on;
    plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    plot(params.centralHole.radius*cos(theta), params.centralHole.radius*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    axis equal;
    xlim(zoom_xlim); ylim(zoom_ylim);
    colorbar;
    caxis([0, max_stress_all]);
    xlabel('x (mm)', 'FontSize', 11);
    ylabel('y (mm)', 'FontSize', 11);
    title(sprintf('Optimal Design (Zoomed)\n\\sigma_{max} = %.1f MPa', optimal_result.max_stress), 'FontSize', 12);

    sgtitle('Zoomed View Near Holes', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Export Figures
if export_figure
    fprintf('\nExporting figures...\n');

    % Full view
    figure(fig);
    exportgraphics(fig, export_filename, 'Resolution', export_resolution);
    fprintf('  Saved: %s\n', export_filename);

    % Zoomed view
    [~, name, ext] = fileparts(export_filename);
    zoomed_filename = [name '_zoomed' ext];
    figure(fig2);
    exportgraphics(fig2, zoomed_filename, 'Resolution', export_resolution);
    fprintf('  Saved: %s\n', zoomed_filename);
end

%% Final Summary
fprintf('\n=============================================================\n');
fprintf('  STRESS COMPARISON COMPLETE\n');
fprintf('=============================================================\n');
fprintf('Figures generated:\n');
fprintf('  1. Full view stress comparison\n');
fprintf('  2. Zoomed view near holes\n');
if export_figure
    fprintf('\nExported files:\n');
    fprintf('  - %s\n', export_filename);
    fprintf('  - %s\n', zoomed_filename);
end
fprintf('\nStress Reduction Summary:\n');
fprintf('  Baseline:  %.2f MPa\n', baseline.max_stress);
fprintf('  Optimal:   %.2f MPa\n', optimal_result.max_stress);
fprintf('  Reduction: %.1f%%\n', all_results.summary.stress_reduction_percent);
fprintf('=============================================================\n');
