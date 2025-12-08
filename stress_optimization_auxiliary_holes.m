%% Shape Optimization of Plate with Stress-Relieving Auxiliary Holes
% This script optimizes the placement, size, and orientation of two
% elliptical auxiliary holes to minimize the maximum von Mises stress
% in a plate with a central circular hole.
%
% FEATURES:
% - y_aux fixed to 0 (centerline) based on optimization results
% - Symmetric rotation constraint: theta_L = -theta_R
% - Parallel multi-start optimization using parfor
% - Complete history saved for all runs (for post-processing animations)
% - Comprehensive visualization and comparison
%
% DESIGN VARIABLES (4 total):
% - x_aux: x-position of auxiliary hole center
% - a: semi-major axis of ellipse
% - b: semi-minor axis of ellipse
% - theta_R: rotation angle of right auxiliary hole (theta_L = -theta_R)
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% ========================================================================
%  SECTION 1: PROBLEM PARAMETERS
%  ========================================================================

fprintf('=============================================================\n');
fprintf('  SHAPE OPTIMIZATION: Stress-Relieving Auxiliary Holes\n');
fprintf('  (4 Variables: y_aux=0, theta_L=-theta_R)\n');
fprintf('=============================================================\n\n');

% --- Plate Geometry (Fixed) ---
params.plate.width = 50.0;                    % Half-width W (mm)
params.plate.length = 4 * params.plate.width; % Half-length L (mm)

% --- Central Hole (Fixed) ---
params.centralHole.radius = 20.0;             % Radius R_c (mm)
params.centralHole.x = 0;
params.centralHole.y = 0;

% --- Material Properties ---
params.material.E = 200e3;                    % Young's modulus (MPa)
params.material.nu = 0.25;                    % Poisson's ratio

% --- Loading ---
params.load.tension = 100;                    % Applied tensile stress sigma_0 (MPa)

% --- Constraint Parameters ---
params.constraints.min_gap_central = 3.0;     % delta_c: min gap from central hole (mm)
params.constraints.min_gap_edge = 5.0;        % delta_e: min gap from plate edges (mm)
params.constraints.min_aux_size = 6.0;        % a_min, b_min (mm)
params.constraints.max_aux_size = 15.0;       % a_max, b_max (mm)
params.constraints.max_aspect_ratio = 3.0;    % alpha_max
params.constraints.max_x_aux = 0.75 * params.plate.length;  % 150 mm

% --- Mesh Parameters ---
params.mesh.Hmax_factor = 2;                  % Quick mode (2), Fine mode (6)
params.mesh.n_ellipse_pts = 100;              % Points for ellipse polygon approximation

% --- Optimization Parameters ---
params.opt.num_starts = 20;                   % Number of multi-start points
params.opt.penalty = 1e10;                    % Penalty for infeasible/failed geometries
params.opt.use_parallel = true;               % Enable parallel processing
params.opt.save_filename = 'optimization_results.mat';  % Results file

%% ========================================================================
%  SECTION 2: DESIGN VARIABLE BOUNDS
%  ========================================================================

% Design vector: x = [x_aux, a, b, theta_R]
% Fixed: y_aux = 0 (auxiliary holes on centerline)
% Symmetric: theta_L = -theta_R

lb = [
    params.centralHole.radius + params.constraints.max_aux_size + params.constraints.min_gap_central;
    params.constraints.min_aux_size;
    params.constraints.min_aux_size;
    -pi/2
];

ub = [
    params.constraints.max_x_aux;
    params.constraints.max_aux_size;
    params.constraints.max_aux_size;
    pi/2
];

fprintf('Design Variables (4 total):\n');
fprintf('  x_aux:   [%.1f, %.1f] mm\n', lb(1), ub(1));
fprintf('  a:       [%.1f, %.1f] mm\n', lb(2), ub(2));
fprintf('  b:       [%.1f, %.1f] mm\n', lb(3), ub(3));
fprintf('  theta_R: [%.2f, %.2f] rad (theta_L = -theta_R)\n', lb(4), ub(4));
fprintf('Fixed parameters:\n');
fprintf('  y_aux:   0 mm (centerline)\n');
fprintf('\n');

%% ========================================================================
%  SECTION 3: BASELINE ANALYSIS
%  ========================================================================

fprintf('=== Running Baseline Analysis (Central Hole Only) ===\n');

tic;
baseline = run_baseline_analysis(params);
baseline.time = toc;

fprintf('Baseline Results:\n');
fprintf('  Max von Mises stress: %.2f MPa\n', baseline.max_stress);
fprintf('  Stress concentration factor (Kt): %.3f\n', baseline.Kt);
fprintf('  Computation time: %.2f s\n\n', baseline.time);

%% ========================================================================
%  SECTION 4: GENERATE INITIAL POINTS
%  ========================================================================

rng(42);  % For reproducibility
lhs_samples = lhsdesign(params.opt.num_starts, 4);
initial_points = lb' + lhs_samples .* (ub' - lb');

%% ========================================================================
%  SECTION 5: PARALLEL MULTI-START OPTIMIZATION
%  ========================================================================

fprintf('=== Starting Multi-Start Optimization ===\n');
fprintf('Number of starting points: %d\n', params.opt.num_starts);

% Initialize parallel pool if needed
if params.opt.use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool;
    end
    fprintf('Using parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Running sequentially (parallel disabled)\n');
end
fprintf('\n');

% Pre-allocate results storage
num_starts = params.opt.num_starts;
run_results = cell(num_starts, 1);

% Optimization options (no display for parallel)
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1500, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4, ...
    'FiniteDifferenceType', 'central');

% Timer
total_start_time = tic;

% Run optimization (parallel or sequential)
if params.opt.use_parallel
    parfor i = 1:num_starts
        run_results{i} = run_single_optimization(i, initial_points(i,:)', lb, ub, options, params);
    end
else
    for i = 1:num_starts
        fprintf('Running start point %d/%d...\n', i, num_starts);
        run_results{i} = run_single_optimization(i, initial_points(i,:)', lb, ub, options, params);
        fprintf('  Result: %.2f MPa (exit flag: %d)\n', run_results{i}.fval, run_results{i}.exitflag);
    end
end

total_opt_time = toc(total_start_time);

%% ========================================================================
%  SECTION 6: PROCESS RESULTS
%  ========================================================================

fprintf('\n=== Processing Results ===\n');

% Find best result
best_fval = Inf;
best_idx = 0;
total_feval = 0;

for i = 1:num_starts
    if run_results{i}.fval < best_fval
        best_fval = run_results{i}.fval;
        best_idx = i;
    end
    total_feval = total_feval + run_results{i}.funcCount;
end

% Display summary of all runs
fprintf('\nRun Summary:\n');
fprintf('%-6s | %-12s | %-10s | %-10s | %-8s\n', 'Run', 'Final Obj', 'Iterations', 'Func Evals', 'Time (s)');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:num_starts
    marker = '';
    if i == best_idx
        marker = ' *BEST*';
    end
    fprintf('%-6d | %-12.2f | %-10d | %-10d | %-8.1f%s\n', ...
        i, run_results{i}.fval, run_results{i}.iterations, ...
        run_results{i}.funcCount, run_results{i}.time, marker);
end

fprintf('\n=== Optimization Complete ===\n');
fprintf('Best result from start point %d\n', best_idx);
fprintf('Optimal objective: %.2f MPa\n', best_fval);
fprintf('Total optimization time: %.1f s\n', total_opt_time);
fprintf('Total function evaluations: %d\n\n', total_feval);

%% ========================================================================
%  SECTION 7: EXTRACT OPTIMAL SOLUTION
%  ========================================================================

x_opt = run_results{best_idx}.x_opt;
opt.x_aux = x_opt(1);
opt.y_aux = 0;  % Fixed at centerline
opt.a = x_opt(2);
opt.b = x_opt(3);
opt.theta_R = x_opt(4);
opt.theta_L = -x_opt(4);

fprintf('=== Optimal Design Variables ===\n');
fprintf('  x_aux:   %.3f mm\n', opt.x_aux);
fprintf('  y_aux:   %.3f mm (fixed at centerline)\n', opt.y_aux);
fprintf('  a:       %.3f mm\n', opt.a);
fprintf('  b:       %.3f mm\n', opt.b);
fprintf('  theta_R: %.3f rad (%.1f deg)\n', opt.theta_R, rad2deg(opt.theta_R));
fprintf('  theta_L: %.3f rad (%.1f deg) [= -theta_R]\n', opt.theta_L, rad2deg(opt.theta_L));
fprintf('\n');

% Verify constraints
[c_opt, ~] = compute_constraints(x_opt, params);
fprintf('Constraint Values at Optimum (should all be <= 0):\n');
constraint_names = {'g1: Hole-central clearance', ...
                    'g2: (same as g1, symmetric)', ...
                    'g3: Hole-right edge', ...
                    'g4: Hole-top/bottom edge', ...
                    'g5: Aspect ratio'};
for j = 1:length(c_opt)
    status = 'OK';
    if c_opt(j) > 1e-6
        status = 'VIOLATED';
    end
    fprintf('  %s: %.3f [%s]\n', constraint_names{j}, c_opt(j), status);
end
fprintf('\n');

%% ========================================================================
%  SECTION 8: FINAL DETAILED ANALYSIS
%  ========================================================================

fprintf('=== Running Final Detailed Analysis ===\n');

% Expand 4-var x_opt to 6-var format: [x_aux, y_aux=0, a, b, theta_R, theta_L]
x_opt_full = [x_opt(1); 0; x_opt(2); x_opt(3); x_opt(4); -x_opt(4)];
optimal_result = run_detailed_analysis(x_opt_full, params);

% Initial guess analysis
x0_best = run_results{best_idx}.x0;
x0_full = [x0_best(1); 0; x0_best(2); x0_best(3); x0_best(4); -x0_best(4)];
initial_result = run_analysis_with_aux_holes(x0_full, params);

%% ========================================================================
%  SECTION 9: COMPARISON TABLE
%  ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    COMPARISON TABLE\n');
fprintf('=============================================================\n');
fprintf('%-25s | %-15s | %-15s | %-15s\n', 'Metric', 'Baseline', 'Initial Guess', 'Optimal');
fprintf('%-25s | %-15s | %-15s | %-15s\n', repmat('-',1,25), repmat('-',1,15), repmat('-',1,15), repmat('-',1,15));
fprintf('%-25s | %-15.2f | %-15.2f | %-15.2f\n', 'Max vM Stress (MPa)', baseline.max_stress, initial_result.max_stress, optimal_result.max_stress);
fprintf('%-25s | %-15.3f | %-15.3f | %-15.3f\n', 'Stress Conc. Factor Kt', baseline.Kt, initial_result.Kt, optimal_result.Kt);

init_reduction = (baseline.max_stress - initial_result.max_stress)/baseline.max_stress * 100;
opt_reduction = (baseline.max_stress - optimal_result.max_stress)/baseline.max_stress * 100;
fprintf('%-25s | %-15s | %-15.1f | %-15.1f\n', 'Reduction vs Baseline (%)', 'N/A', init_reduction, opt_reduction);
fprintf('=============================================================\n\n');

%% ========================================================================
%  SECTION 10: SAVE ALL RESULTS
%  ========================================================================

fprintf('=== Saving Results ===\n');

% Compile all results
all_results = struct();
all_results.runs = run_results;
all_results.params = params;
all_results.baseline = baseline;
all_results.lb = lb;
all_results.ub = ub;
all_results.initial_points = initial_points;
all_results.global_best_x = x_opt;
all_results.global_best_fval = best_fval;
all_results.global_best_idx = best_idx;
all_results.total_time = total_opt_time;
all_results.total_feval = total_feval;
all_results.optimal_result = optimal_result;
all_results.initial_result = initial_result;
all_results.timestamp = datetime('now');

% Summary
all_results.summary.baseline_stress = baseline.max_stress;
all_results.summary.baseline_Kt = baseline.Kt;
all_results.summary.optimal_stress = optimal_result.max_stress;
all_results.summary.optimal_Kt = optimal_result.Kt;
all_results.summary.stress_reduction_percent = opt_reduction;
all_results.summary.num_starts = num_starts;
all_results.summary.total_time_seconds = total_opt_time;
all_results.summary.total_function_evals = total_feval;

% Save
save(params.opt.save_filename, 'all_results', '-v7.3');
fprintf('Results saved to: %s\n', params.opt.save_filename);
fprintf('  - Contains complete history from all %d runs\n', num_starts);
fprintf('  - Use plot_optimization_results(''%s'') to regenerate plots\n', params.opt.save_filename);
fprintf('  - Use animate_optimization_run(''%s'', run_idx) to animate any run\n\n', params.opt.save_filename);

%% ========================================================================
%  SECTION 11: GENERATE VISUALIZATIONS
%  ========================================================================

fprintf('=== Generating Visualizations ===\n');

% --- Figure 1: Geometry Comparison ---
figure('Name', 'Geometry Comparison', 'Position', [50, 50, 1400, 500]);

subplot(1, 3, 1);
[g_base, ~] = create_baseline_geometry(params);
pdegplot(g_base, 'FaceLabels', 'off');
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title('Baseline (Central Hole Only)');
grid on;

subplot(1, 3, 2);
[g_init, ~] = create_geometry(x0_full, params);
pdegplot(g_init, 'FaceLabels', 'off');
hold on;
plot(x0_best(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);  % y_aux = 0
plot(-x0_best(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);  % y_aux = 0
hold off;
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Initial Guess (\\sigma_{max} = %.1f MPa)', initial_result.max_stress));
grid on;

subplot(1, 3, 3);
[g_opt, ~] = create_geometry(x_opt_full, params);
pdegplot(g_opt, 'FaceLabels', 'off');
hold on;
plot(opt.x_aux, opt.y_aux, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(-opt.x_aux, opt.y_aux, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
quiver(opt.x_aux, opt.y_aux, opt.a*cos(opt.theta_R)*0.8, opt.a*sin(opt.theta_R)*0.8, 'r', 'LineWidth', 2, 'AutoScale', 'off');
quiver(-opt.x_aux, opt.y_aux, -opt.a*cos(opt.theta_L)*0.8, opt.a*sin(opt.theta_L)*0.8, 'b', 'LineWidth', 2, 'AutoScale', 'off');
hold off;
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Optimal (\\sigma_{max} = %.1f MPa, %.1f%% reduction)', optimal_result.max_stress, opt_reduction));
grid on;

% --- Figure 2: Stress Comparison ---
figure('Name', 'Stress Comparison', 'Position', [50, 50, 1400, 500]);

subplot(1, 3, 1);
pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
axis equal; colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Baseline: \\sigma_{vM,max} = %.1f MPa', baseline.max_stress));

subplot(1, 3, 2);
pdeplot(initial_result.mesh, 'XYData', initial_result.vonMises, 'ColorMap', 'jet');
axis equal; colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Initial: \\sigma_{vM,max} = %.1f MPa', initial_result.max_stress));

subplot(1, 3, 3);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal; colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Optimal: \\sigma_{vM,max} = %.1f MPa', optimal_result.max_stress));

% --- Figure 3: Detailed Stress Components ---
figure('Name', 'Optimal Stress Components', 'Position', [50, 50, 1400, 900]);

subplot(2, 2, 1);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.stress.sxx, 'ColorMap', 'jet', 'Contour', 'on');
axis equal; colorbar;
title('\sigma_{xx} (MPa)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2, 2, 2);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.stress.syy, 'ColorMap', 'jet', 'Contour', 'on');
axis equal; colorbar;
title('\sigma_{yy} (MPa)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2, 2, 3);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.stress.sxy, 'ColorMap', 'jet', 'Contour', 'on');
axis equal; colorbar;
title('\sigma_{xy} (MPa)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2, 2, 4);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet', 'Contour', 'on');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal; colorbar;
title(sprintf('\\sigma_{vM} (Max = %.1f MPa)', optimal_result.max_stress));
xlabel('x (mm)'); ylabel('y (mm)');

% --- Figure 4: Zoomed View ---
figure('Name', 'Stress Detail Near Holes', 'Position', [50, 50, 1000, 500]);

subplot(1, 2, 1);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
xlim([-100, 100]); ylim([-60, 60]);
colorbar;
title('von Mises Stress (Zoomed)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(1, 2, 2);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.Displacement.ux, 'ColorMap', 'jet');
axis equal;
xlim([-100, 100]); ylim([-60, 60]);
colorbar;
title('X-Displacement u_x (mm)');
xlabel('x (mm)'); ylabel('y (mm)');

% --- Figure 5: Convergence History ---
figure('Name', 'Convergence History', 'Position', [50, 50, 1200, 500]);

subplot(1, 2, 1);
hold on;
for i = 1:num_starts
    if ~isempty(run_results{i}.history.fval)
        hist_fval = run_results{i}.history.fval;
        if i == best_idx
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
best_hist = run_results{best_idx}.history;
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
title(sprintf('Best Run Convergence (Start %d)', best_idx));
legend('Objective', 'Baseline', 'Optimal', 'Location', 'northeast');
grid on;

% --- Figure 6: Boundary Stress ---
figure('Name', 'Boundary Stress Distribution', 'Position', [50, 50, 1200, 400]);

theta_pts = linspace(0, 2*pi, 200);

% Central hole
x_central = params.centralHole.radius * cos(theta_pts);
y_central = params.centralHole.radius * sin(theta_pts);
try
    stress_central = interpolateStress(optimal_result.R, [x_central; y_central]);
    vM_central = sqrt(stress_central.sxx.^2 - stress_central.sxx.*stress_central.syy + ...
                      stress_central.syy.^2 + 3*stress_central.sxy.^2);
    subplot(1, 3, 1);
    plot(rad2deg(theta_pts), vM_central, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_pts), stress_central.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)'); ylabel('Stress (MPa)');
    title('Central Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on; xlim([0, 360]);
catch
    subplot(1, 3, 1);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

% Right aux hole
x_local = opt.a * cos(theta_pts);
y_local = opt.b * sin(theta_pts);
x_right = cos(opt.theta_R)*x_local - sin(opt.theta_R)*y_local + opt.x_aux;
y_right = sin(opt.theta_R)*x_local + cos(opt.theta_R)*y_local + opt.y_aux;
try
    stress_right = interpolateStress(optimal_result.R, [x_right; y_right]);
    vM_right = sqrt(stress_right.sxx.^2 - stress_right.sxx.*stress_right.syy + ...
                    stress_right.syy.^2 + 3*stress_right.sxy.^2);
    subplot(1, 3, 2);
    plot(rad2deg(theta_pts), vM_right, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_pts), stress_right.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)'); ylabel('Stress (MPa)');
    title('Right Aux Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on; xlim([0, 360]);
catch
    subplot(1, 3, 2);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

% Left aux hole
x_left = -(cos(opt.theta_L)*x_local - sin(opt.theta_L)*y_local) - opt.x_aux;
y_left = sin(opt.theta_L)*x_local + cos(opt.theta_L)*y_local + opt.y_aux;
try
    stress_left = interpolateStress(optimal_result.R, [x_left; y_left]);
    vM_left = sqrt(stress_left.sxx.^2 - stress_left.sxx.*stress_left.syy + ...
                   stress_left.syy.^2 + 3*stress_left.sxy.^2);
    subplot(1, 3, 3);
    plot(rad2deg(theta_pts), vM_left, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_pts), stress_left.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)'); ylabel('Stress (MPa)');
    title('Left Aux Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on; xlim([0, 360]);
catch
    subplot(1, 3, 3);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

%% ========================================================================
%  SECTION 12: POST-PROCESSING ANIMATION
%  ========================================================================

fprintf('=== Generating Design Evolution Animation (Best Run) ===\n');

animate_single_run(run_results{best_idx}, params, baseline);

%% ========================================================================
%  SECTION 13: FINAL SUMMARY
%  ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    OPTIMIZATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('\n');
fprintf('OPTIMAL DESIGN (4 variables: y_aux=0, theta_L=-theta_R):\n');
fprintf('  Auxiliary hole x-position: %.2f mm from center\n', opt.x_aux);
fprintf('  Auxiliary hole y-position: %.2f mm (fixed at centerline)\n', opt.y_aux);
fprintf('  Ellipse semi-axes:         a = %.2f mm, b = %.2f mm\n', opt.a, opt.b);
fprintf('  Aspect ratio:              %.2f\n', opt.a/opt.b);
fprintf('  Right hole rotation:       %.1f degrees\n', rad2deg(opt.theta_R));
fprintf('  Left hole rotation:        %.1f degrees (symmetric)\n', rad2deg(opt.theta_L));
fprintf('\n');
fprintf('PERFORMANCE:\n');
fprintf('  Baseline max stress:      %.2f MPa (Kt = %.3f)\n', baseline.max_stress, baseline.Kt);
fprintf('  Optimal max stress:       %.2f MPa (Kt = %.3f)\n', optimal_result.max_stress, optimal_result.Kt);
fprintf('  Stress reduction:         %.1f%%\n', opt_reduction);
fprintf('\n');
fprintf('COMPUTATION:\n');
fprintf('  Number of multi-starts:   %d\n', num_starts);
fprintf('  Parallel workers:         %d\n', pool.NumWorkers);
fprintf('  Total function evals:     %d\n', total_feval);
fprintf('  Total optimization time:  %.1f seconds\n', total_opt_time);
fprintf('  Best run (start point):   %d\n', best_idx);
fprintf('\n');
fprintf('FILES:\n');
fprintf('  Results saved to: %s\n', params.opt.save_filename);
fprintf('=============================================================\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function result = run_single_optimization(run_idx, x0, lb, ub, options, params)
%RUN_SINGLE_OPTIMIZATION Run one optimization from a starting point
%   This function is designed to run independently in parallel workers

    % Initialize result structure
    result = struct();
    result.run_idx = run_idx;
    result.x0 = x0;
    result.x_opt = [];
    result.fval = Inf;
    result.exitflag = -99;
    result.iterations = 0;
    result.funcCount = 0;
    result.time = 0;
    result.history = struct('fval', [], 'x', []);

    % Check feasibility of starting point
    [c0, ~] = compute_constraints(x0, params);
    if any(c0 > 0)
        x0 = (lb + ub) / 2;
        x0(4) = 0;  % theta_R = 0 as safe default
    end

    % Create objective and constraint function handles
    history = struct('fval', [], 'x', []);

    function f = obj_with_history(x)
        % Expand 4-var to 6-var: [x_aux, y_aux=0, a, b, theta_R, theta_L=-theta_R]
        x_full = [x(1); 0; x(2); x(3); x(4); -x(4)];
        try
            res = run_analysis_with_aux_holes(x_full, params);
            f = res.max_stress;
        catch
            f = params.opt.penalty;
        end
        history.fval = [history.fval; f];
        history.x = [history.x; x'];
    end

    function [c, ceq] = nonlcon(x)
        [c, ceq] = compute_constraints(x, params);
    end

    % Run optimization
    run_start = tic;
    try
        [x_opt, fval, exitflag, output] = fmincon(@obj_with_history, x0, ...
            [], [], [], [], lb, ub, @nonlcon, options);

        result.x_opt = x_opt;
        result.fval = fval;
        result.exitflag = exitflag;
        result.iterations = output.iterations;
        result.funcCount = output.funcCount;
        result.output = output;
    catch ME
        result.error = ME.message;
    end

    result.time = toc(run_start);
    result.history = history;
end

function [c, ceq] = compute_constraints(x, params)
%COMPUTE_CONSTRAINTS Compute nonlinear inequality constraints
%   x = [x_aux, a, b, theta_R] (4 variables)
%   y_aux = 0 (fixed), theta_L = -theta_R (symmetric)

    x_aux = x(1);
    y_aux = 0;  % Fixed at centerline
    a = x(2);
    b = x(3);
    theta_R = x(4);
    theta_L = -theta_R;

    R_c = params.centralHole.radius;
    L = params.plate.length;
    W = params.plate.width;
    delta_c = params.constraints.min_gap_central;
    delta_e = params.constraints.min_gap_edge;
    alpha_max = params.constraints.max_aspect_ratio;

    ex_R = sqrt(a^2*cos(theta_R)^2 + b^2*sin(theta_R)^2);
    ey_R = sqrt(a^2*sin(theta_R)^2 + b^2*cos(theta_R)^2);
    ex_L = sqrt(a^2*cos(theta_L)^2 + b^2*sin(theta_L)^2);
    ey_L = sqrt(a^2*sin(theta_L)^2 + b^2*cos(theta_L)^2);

    e_max_R = max(ex_R, ey_R);
    e_max_L = max(ex_L, ey_L);
    d_c = x_aux;  % Since y_aux = 0, distance is simply x_aux

    c = zeros(5, 1);  % Reduced from 7 to 5 (y constraints are trivially satisfied with y_aux=0)
    c(1) = R_c + e_max_R + delta_c - d_c;  % Right hole-central clearance
    c(2) = R_c + e_max_L + delta_c - d_c;  % Left hole-central clearance (same as c(1))
    c(3) = (x_aux + ex_R) - (L - delta_e);  % Right hole-right edge
    c(4) = ey_R - (W - delta_e);  % Right hole-top/bottom edge (y_aux = 0)
    c(5) = a - alpha_max * b;  % Aspect ratio

    ceq = [];
end

function result = run_baseline_analysis(params)
%RUN_BASELINE_ANALYSIS Run FEA for plate with central hole only

    [g, ~] = create_baseline_geometry(params);

    model = femodel(AnalysisType="structuralStatic", Geometry=g);
    model.MaterialProperties = materialProperties(...
        YoungsModulus=params.material.E, PoissonsRatio=params.material.nu);
    model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
    model.EdgeBC(4) = edgeBC(XDisplacement=0);
    model.VertexBC(1) = vertexBC(YDisplacement=0);

    Hmax = params.centralHole.radius / params.mesh.Hmax_factor;
    model = generateMesh(model, Hmax=Hmax);

    R = solve(model);
    stress = R.Stress;
    vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

    result.max_stress = max(vonMises);
    result.Kt = result.max_stress / params.load.tension;
    result.vonMises = vonMises;
    result.mesh = R.Mesh;
    result.stress = stress;
end

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

function result = run_analysis_with_aux_holes(x, params)
%RUN_ANALYSIS_WITH_AUX_HOLES Run FEA with auxiliary holes

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
end

function result = run_detailed_analysis(x, params)
%RUN_DETAILED_ANALYSIS Run FEA with full output

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

    fprintf('  Max von Mises stress: %.2f MPa at (%.2f, %.2f) mm\n', ...
        max_stress, result.loc_max(1), result.loc_max(2));
    fprintf('  Kt = %.3f\n', result.Kt);
end

function [g, status] = create_geometry(x, params)
%CREATE_GEOMETRY Create plate with central and auxiliary holes

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

function animate_single_run(run_result, params, baseline)
%ANIMATE_SINGLE_RUN Create animation from a single optimization run
%   x = [x_aux, a, b, theta_R] (4 variables), y_aux = 0

    history = run_result.history;
    if isempty(history.x) || size(history.x, 1) < 2
        fprintf('Insufficient history for animation.\n');
        return;
    end

    n_frames = min(30, size(history.x, 1));
    frame_indices = round(linspace(1, size(history.x, 1), n_frames));

    fig = figure('Name', 'Design Evolution Animation', 'Position', [50, 50, 1200, 500]);

    for frame = 1:length(frame_indices)
        idx = frame_indices(frame);
        x_frame = history.x(idx, :)';
        fval_frame = history.fval(idx);

        % Expand 4-var to 6-var: [x_aux, y_aux=0, a, b, theta_R, theta_L]
        x_full = [x_frame(1); 0; x_frame(2); x_frame(3); x_frame(4); -x_frame(4)];

        clf(fig);

        subplot(1, 2, 1);
        try
            [g_frame, ~] = create_geometry(x_full, params);
            pdegplot(g_frame, 'FaceLabels', 'off');
            hold on;
            plot(x_frame(1), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);  % y_aux = 0
            plot(-x_frame(1), 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);  % y_aux = 0
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

        sgtitle(sprintf('Design Evolution - Frame %d/%d', frame, n_frames));
        drawnow;
        pause(0.15);
    end

    fprintf('Animation complete.\n');
end
