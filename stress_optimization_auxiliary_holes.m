%% Shape Optimization of Plate with Stress-Relieving Auxiliary Holes
% This script optimizes the placement, size, and orientation of two
% elliptical auxiliary holes to minimize the maximum von Mises stress
% in a plate with a central circular hole.
%
% Optimization formulation follows: Optimization_Problem_Formulation.md
% FEA methodology based on: stress_analysis_auxiliary_holes.m
%
% KEY CHANGE: Symmetric rotation constraint enforced (theta_L = -theta_R)
% This ensures symmetry for reverse loading conditions.
% Design variables reduced from 6 to 5.
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% ========================================================================
%  SECTION 1: PROBLEM PARAMETERS (from Optimization_Problem_Formulation.md)
%  ========================================================================

fprintf('=============================================================\n');
fprintf('  SHAPE OPTIMIZATION: Stress-Relieving Auxiliary Holes\n');
fprintf('  (With Symmetric Rotation Constraint: theta_L = -theta_R)\n');
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
params.opt.save_filename = 'optimization_results.mat';  % Results file

% --- Animation Parameters ---
params.anim.enabled = true;                   % Enable real-time animation
params.anim.pause_time = 0.1;                 % Pause between frames (seconds)

%% ========================================================================
%  SECTION 2: DESIGN VARIABLE BOUNDS
%  ========================================================================
%  Design vector: x = [x_aux, y_aux, a, b, theta_R]
%  Note: theta_L = -theta_R (symmetric constraint enforced internally)

% Lower bounds
lb = [
    params.centralHole.radius + params.constraints.max_aux_size + params.constraints.min_gap_central;  % x_aux >= 38 mm
    0;                                                    % y_aux >= 0
    params.constraints.min_aux_size;                      % a >= 6 mm
    params.constraints.min_aux_size;                      % b >= 6 mm
    -pi/2                                                 % theta_R >= -90 deg
];

% Upper bounds
ub = [
    params.constraints.max_x_aux;                         % x_aux <= 150 mm
    params.plate.width - params.constraints.min_gap_edge - params.constraints.max_aux_size;  % y_aux <= 30 mm
    params.constraints.max_aux_size;                      % a <= 15 mm
    params.constraints.max_aux_size;                      % b <= 15 mm
    pi/2                                                  % theta_R <= 90 deg
];

fprintf('Design Variables (5 total, with symmetric rotation):\n');
fprintf('  x_aux:   [%.1f, %.1f] mm\n', lb(1), ub(1));
fprintf('  y_aux:   [%.1f, %.1f] mm\n', lb(2), ub(2));
fprintf('  a:       [%.1f, %.1f] mm\n', lb(3), ub(3));
fprintf('  b:       [%.1f, %.1f] mm\n', lb(4), ub(4));
fprintf('  theta_R: [%.2f, %.2f] rad (theta_L = -theta_R)\n', lb(5), ub(5));
fprintf('\n');

%% ========================================================================
%  SECTION 3: BASELINE ANALYSIS (Central Hole Only)
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
%  SECTION 4: SETUP ANIMATION FIGURE
%  ========================================================================

if params.anim.enabled
    anim_fig = figure('Name', 'Optimization Progress', 'Position', [50, 50, 1400, 600]);

    % Subplot 1: Current geometry
    ax_geom = subplot(1, 3, 1);
    title('Current Design');
    xlabel('x (mm)');
    ylabel('y (mm)');
    axis equal;
    grid on;

    % Subplot 2: Stress distribution
    ax_stress = subplot(1, 3, 2);
    title('von Mises Stress');
    xlabel('x (mm)');
    ylabel('y (mm)');
    axis equal;

    % Subplot 3: Convergence
    ax_conv = subplot(1, 3, 3);
    title('Convergence History');
    xlabel('Function Evaluation');
    ylabel('Max von Mises Stress (MPa)');
    grid on;
    hold(ax_conv, 'on');
    yline(ax_conv, baseline.max_stress, 'r--', 'LineWidth', 1.5);

    drawnow;
end

%% ========================================================================
%  SECTION 5: OPTIMIZATION SETUP
%  ========================================================================

% Global storage for iteration history (accessible by nested functions)
global_history = struct();
global_history.fval = [];
global_history.x = [];
global_history.feval_count = 0;
global_history.best_fval = Inf;
global_history.best_x = [];

% --- Objective Function (with animation callback) ---
objective = @(x) objective_with_animation(x, params, anim_fig, baseline);

% --- Nonlinear Constraints ---
nonlcon = @(x) constraint_function(x, params);

% --- Optimizer Options ---
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1500, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4, ...
    'FiniteDifferenceType', 'central');

%% ========================================================================
%  SECTION 6: MULTI-START OPTIMIZATION (Sequential for Animation)
%  ========================================================================

fprintf('=== Starting Multi-Start Optimization ===\n');
fprintf('Number of starting points: %d\n', params.opt.num_starts);
fprintf('Note: Running sequentially to enable real-time animation\n\n');

% Storage for results from each start
all_results = struct();
all_results.runs = cell(params.opt.num_starts, 1);
all_results.params = params;
all_results.baseline = baseline;
all_results.lb = lb;
all_results.ub = ub;
all_results.timestamp = datetime('now');

% Generate initial points using Latin Hypercube Sampling
rng(42);  % For reproducibility
lhs_samples = lhsdesign(params.opt.num_starts, 5);
initial_points = lb' + lhs_samples .* (ub' - lb');

% Global best tracking
global_best_fval = Inf;
global_best_x = [];
global_best_idx = 0;

total_opt_time = 0;
total_feval = 0;

% Run optimization from each starting point
for i = 1:params.opt.num_starts
    fprintf('\n========== Start Point %d/%d ==========\n', i, params.opt.num_starts);

    x0 = initial_points(i, :)';

    % Reset global history for this run
    global_history.fval = [];
    global_history.x = [];
    global_history.feval_count = 0;
    global_history.best_fval = Inf;
    global_history.best_x = [];

    % Check if starting point is feasible
    [c0, ~] = nonlcon(x0);
    if any(c0 > 0)
        fprintf('  Initial point infeasible, adjusting...\n');
        x0 = (lb + ub) / 2;
        x0(5) = 0;  % theta_R = 0
    end

    fprintf('  Initial: x_aux=%.1f, y_aux=%.1f, a=%.1f, b=%.1f, theta=%.1f deg\n', ...
        x0(1), x0(2), x0(3), x0(4), rad2deg(x0(5)));

    % Run fmincon
    run_start_time = tic;
    try
        [x_opt, fval, exitflag, output] = fmincon(objective, x0, ...
            [], [], [], [], lb, ub, nonlcon, options);

        run_time = toc(run_start_time);

        % Store run results
        run_result = struct();
        run_result.x0 = x0;
        run_result.x_opt = x_opt;
        run_result.fval = fval;
        run_result.exitflag = exitflag;
        run_result.output = output;
        run_result.history = global_history;
        run_result.time = run_time;
        run_result.iterations = output.iterations;
        run_result.funcCount = output.funcCount;

        all_results.runs{i} = run_result;

        total_opt_time = total_opt_time + run_time;
        total_feval = total_feval + output.funcCount;

        fprintf('  Final objective: %.2f MPa\n', fval);
        fprintf('  Iterations: %d, Function evals: %d\n', output.iterations, output.funcCount);
        fprintf('  Exit flag: %d, Time: %.1f s\n', exitflag, run_time);

        % Update global best
        if fval < global_best_fval
            global_best_fval = fval;
            global_best_x = x_opt;
            global_best_idx = i;
            fprintf('  *** New best solution! ***\n');
        end

    catch ME
        fprintf('  Optimization failed: %s\n', ME.message);
        run_result = struct();
        run_result.x0 = x0;
        run_result.fval = Inf;
        run_result.exitflag = -99;
        run_result.error = ME.message;
        run_result.time = toc(run_start_time);
        all_results.runs{i} = run_result;
    end
end

% Store global best
all_results.global_best_x = global_best_x;
all_results.global_best_fval = global_best_fval;
all_results.global_best_idx = global_best_idx;
all_results.total_time = total_opt_time;
all_results.total_feval = total_feval;

fprintf('\n=== Optimization Complete ===\n');
fprintf('Best result from start point %d\n', global_best_idx);
fprintf('Optimal objective: %.2f MPa\n', global_best_fval);
fprintf('Total optimization time: %.1f s\n', total_opt_time);
fprintf('Total function evaluations: %d\n\n', total_feval);

%% ========================================================================
%  SECTION 7: EXTRACT OPTIMAL SOLUTION
%  ========================================================================

% Extract optimal design variables (5 variables + computed theta_L)
x_opt = global_best_x;
opt.x_aux = x_opt(1);
opt.y_aux = x_opt(2);
opt.a = x_opt(3);
opt.b = x_opt(4);
opt.theta_R = x_opt(5);
opt.theta_L = -x_opt(5);  % Symmetric constraint

fprintf('=== Optimal Design Variables ===\n');
fprintf('  x_aux:   %.3f mm\n', opt.x_aux);
fprintf('  y_aux:   %.3f mm\n', opt.y_aux);
fprintf('  a:       %.3f mm\n', opt.a);
fprintf('  b:       %.3f mm\n', opt.b);
fprintf('  theta_R: %.3f rad (%.1f deg)\n', opt.theta_R, rad2deg(opt.theta_R));
fprintf('  theta_L: %.3f rad (%.1f deg) [= -theta_R]\n', opt.theta_L, rad2deg(opt.theta_L));
fprintf('\n');

% Verify constraints at optimal point
[c_opt, ~] = nonlcon(x_opt);
fprintf('Constraint Values at Optimum (should all be <= 0):\n');
constraint_names = {'g1: Right hole-central clearance', ...
                    'g2: Left hole-central clearance', ...
                    'g3: Right hole-right edge', ...
                    'g4: Left hole-left edge', ...
                    'g5: Right hole-top/bottom edge', ...
                    'g6: Left hole-top/bottom edge', ...
                    'g7: Aspect ratio'};
all_satisfied = true;
for j = 1:length(c_opt)
    status = 'OK';
    if c_opt(j) > 1e-6
        status = 'VIOLATED';
        all_satisfied = false;
    end
    fprintf('  %s: %.3f [%s]\n', constraint_names{j}, c_opt(j), status);
end
fprintf('\n');

%% ========================================================================
%  SECTION 8: FINAL DETAILED ANALYSIS ON OPTIMAL DESIGN
%  ========================================================================

fprintf('=== Running Final Detailed Analysis on Optimal Design ===\n');

% Convert 5-variable to 6-variable format for analysis
x_opt_full = [x_opt(1:4); x_opt(5); -x_opt(5)];
optimal_result = run_detailed_analysis(x_opt_full, params);
all_results.optimal_result = optimal_result;

%% ========================================================================
%  SECTION 9: INITIAL GUESS ANALYSIS (for comparison)
%  ========================================================================

% Get initial guess from best run
x0_best = all_results.runs{global_best_idx}.x0;
x0_full = [x0_best(1:4); x0_best(5); -x0_best(5)];
initial_result = run_analysis_with_aux_holes(x0_full, params, false);
all_results.initial_result = initial_result;

%% ========================================================================
%  SECTION 10: COMPARISON TABLE
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
%  SECTION 11: SAVE RESULTS TO FILE
%  ========================================================================

fprintf('=== Saving Results to File ===\n');

% Add summary statistics
all_results.summary.baseline_stress = baseline.max_stress;
all_results.summary.baseline_Kt = baseline.Kt;
all_results.summary.optimal_stress = optimal_result.max_stress;
all_results.summary.optimal_Kt = optimal_result.Kt;
all_results.summary.stress_reduction_percent = opt_reduction;
all_results.summary.num_starts = params.opt.num_starts;
all_results.summary.total_time_seconds = total_opt_time;
all_results.summary.total_function_evals = total_feval;

% Save to file
save(params.opt.save_filename, 'all_results', '-v7.3');
fprintf('Results saved to: %s\n\n', params.opt.save_filename);

%% ========================================================================
%  SECTION 12: COMPREHENSIVE VISUALIZATION
%  ========================================================================

fprintf('=== Generating Comprehensive Visualizations ===\n');

% Close animation figure
if params.anim.enabled && ishandle(anim_fig)
    close(anim_fig);
end

% --- Figure 1: Geometry Comparison ---
figure('Name', 'Geometry Comparison', 'Position', [50, 50, 1400, 500]);

% Baseline geometry
subplot(1, 3, 1);
[g_base, ~] = create_baseline_geometry(params);
pdegplot(g_base, 'FaceLabels', 'off');
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title('Baseline (Central Hole Only)');
grid on;

% Initial guess geometry
subplot(1, 3, 2);
[g_init, ~] = create_geometry(x0_full, params);
pdegplot(g_init, 'FaceLabels', 'off');
hold on;
plot(x0_best(1), x0_best(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(-x0_best(1), x0_best(2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
hold off;
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Initial Guess (\\sigma_{max} = %.1f MPa)', initial_result.max_stress));
grid on;

% Optimal geometry
subplot(1, 3, 3);
[g_opt, ~] = create_geometry(x_opt_full, params);
pdegplot(g_opt, 'FaceLabels', 'off');
hold on;
plot(opt.x_aux, opt.y_aux, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(-opt.x_aux, opt.y_aux, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
% Draw rotation arrows
quiver(opt.x_aux, opt.y_aux, opt.a*cos(opt.theta_R)*0.8, opt.a*sin(opt.theta_R)*0.8, 'r', 'LineWidth', 2, 'AutoScale', 'off');
quiver(-opt.x_aux, opt.y_aux, -opt.a*cos(opt.theta_L)*0.8, opt.a*sin(opt.theta_L)*0.8, 'b', 'LineWidth', 2, 'AutoScale', 'off');
hold off;
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Optimal (\\sigma_{max} = %.1f MPa, %.1f%% reduction)', optimal_result.max_stress, opt_reduction));
grid on;

% --- Figure 2: Stress Contours Comparison ---
figure('Name', 'Stress Comparison', 'Position', [50, 50, 1400, 500]);

% Baseline stress
subplot(1, 3, 1);
pdeplot(baseline.mesh, 'XYData', baseline.vonMises, 'ColorMap', 'jet');
axis equal;
colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Baseline: \\sigma_{vM,max} = %.1f MPa', baseline.max_stress));

% Initial stress
subplot(1, 3, 2);
pdeplot(initial_result.mesh, 'XYData', initial_result.vonMises, 'ColorMap', 'jet');
axis equal;
colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Initial: \\sigma_{vM,max} = %.1f MPa', initial_result.max_stress));

% Optimal stress
subplot(1, 3, 3);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
colorbar;
caxis([0, max(baseline.vonMises)]);
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Optimal: \\sigma_{vM,max} = %.1f MPa', optimal_result.max_stress));

% --- Figure 3: Detailed Optimal Stress Components ---
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

% --- Figure 4: Zoomed Stress Near Holes ---
figure('Name', 'Stress Detail Near Holes', 'Position', [50, 50, 1000, 500]);

subplot(1, 2, 1);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
xlim([-100, 100]);
ylim([-60, 60]);
colorbar;
title('von Mises Stress (Zoomed)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(1, 2, 2);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.Displacement.ux, 'ColorMap', 'jet');
axis equal;
xlim([-100, 100]);
ylim([-60, 60]);
colorbar;
title('X-Displacement u_x (mm)');
xlabel('x (mm)'); ylabel('y (mm)');

% --- Figure 5: Convergence History (All Runs) ---
figure('Name', 'Convergence History', 'Position', [50, 50, 1200, 500]);

subplot(1, 2, 1);
hold on;
colors = lines(params.opt.num_starts);
for i = 1:params.opt.num_starts
    if ~isempty(all_results.runs{i}) && isfield(all_results.runs{i}, 'history') && ~isempty(all_results.runs{i}.history.fval)
        hist_fval = all_results.runs{i}.history.fval;
        if i == global_best_idx
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
legend('Best Run', 'Other Runs', 'Baseline', 'Location', 'best');
grid on;

subplot(1, 2, 2);
% Best run convergence
best_hist = all_results.runs{global_best_idx}.history;
if ~isempty(best_hist.fval)
    plot(1:length(best_hist.fval), best_hist.fval, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    hold on;
    yline(baseline.max_stress, 'r--', 'LineWidth', 1.5, 'Label', 'Baseline');
    [min_fval, min_idx] = min(best_hist.fval);
    plot(min_idx, min_fval, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    xlabel('Function Evaluation');
    ylabel('Max von Mises Stress (MPa)');
    title(sprintf('Best Run Convergence (Start %d)', global_best_idx));
    legend('Objective', 'Baseline', 'Optimal', 'Location', 'best');
    grid on;
end

% --- Figure 6: Stress Along Hole Boundaries ---
figure('Name', 'Boundary Stress Distribution', 'Position', [50, 50, 1200, 400]);

% Interpolate stress on hole boundaries
theta_pts = linspace(0, 2*pi, 200);

% Central hole boundary
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
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Central Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);
catch
    subplot(1, 3, 1);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

% Right auxiliary hole boundary
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
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Right Aux Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);
catch
    subplot(1, 3, 2);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

% Left auxiliary hole boundary
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
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Left Aux Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);
catch
    subplot(1, 3, 3);
    text(0.5, 0.5, 'Could not interpolate', 'HorizontalAlignment', 'center');
end

%% ========================================================================
%  SECTION 13: DESIGN EVOLUTION ANIMATION (Post-Processing)
%  ========================================================================

fprintf('=== Generating Design Evolution Animation ===\n');

% Create animation from best run history
best_history = all_results.runs{global_best_idx}.history;

if ~isempty(best_history.x) && size(best_history.x, 1) > 1
    % Sample frames (max 30 frames for reasonable animation)
    n_frames = min(30, size(best_history.x, 1));
    frame_indices = round(linspace(1, size(best_history.x, 1), n_frames));

    anim_evolution_fig = figure('Name', 'Design Evolution', 'Position', [50, 50, 1200, 500]);

    for frame = 1:length(frame_indices)
        idx = frame_indices(frame);
        x_frame = best_history.x(idx, :)';
        fval_frame = best_history.fval(idx);

        % Convert to 6-variable format
        x_full = [x_frame(1:4); x_frame(5); -x_frame(5)];

        clf(anim_evolution_fig);

        % Subplot 1: Geometry
        subplot(1, 2, 1);
        try
            [g_frame, ~] = create_geometry(x_full, params);
            pdegplot(g_frame, 'FaceLabels', 'off');
            hold on;
            plot(x_frame(1), x_frame(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            plot(-x_frame(1), x_frame(2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
            hold off;
        catch
            text(0.5, 0.5, 'Geometry failed', 'HorizontalAlignment', 'center');
        end
        axis equal;
        axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
        xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('Iteration %d/%d', idx, size(best_history.x, 1)));
        grid on;

        % Subplot 2: Convergence progress
        subplot(1, 2, 2);
        plot(1:idx, best_history.fval(1:idx), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        yline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
        plot(idx, fval_frame, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        xlabel('Function Evaluation');
        ylabel('Max von Mises Stress (MPa)');
        title(sprintf('\\sigma_{max} = %.1f MPa', fval_frame));
        xlim([0, size(best_history.x, 1) + 1]);
        ylim([min(best_history.fval)*0.95, max(best_history.fval)*1.05]);
        grid on;

        sgtitle(sprintf('Design Evolution - Frame %d/%d', frame, length(frame_indices)));

        drawnow;
        pause(0.2);
    end

    fprintf('Animation complete.\n\n');
else
    fprintf('Insufficient history data for animation.\n\n');
end

%% ========================================================================
%  SECTION 14: FINAL SUMMARY
%  ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    OPTIMIZATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('\n');
fprintf('OPTIMAL DESIGN (with symmetric rotation theta_L = -theta_R):\n');
fprintf('  Auxiliary hole position:  (%.2f, %.2f) mm from center\n', opt.x_aux, opt.y_aux);
fprintf('  Ellipse semi-axes:        a = %.2f mm, b = %.2f mm\n', opt.a, opt.b);
fprintf('  Aspect ratio:             %.2f\n', opt.a/opt.b);
fprintf('  Right hole rotation:      %.1f degrees\n', rad2deg(opt.theta_R));
fprintf('  Left hole rotation:       %.1f degrees (symmetric)\n', rad2deg(opt.theta_L));
fprintf('\n');
fprintf('PERFORMANCE:\n');
fprintf('  Baseline max stress:      %.2f MPa (Kt = %.3f)\n', baseline.max_stress, baseline.Kt);
fprintf('  Optimal max stress:       %.2f MPa (Kt = %.3f)\n', optimal_result.max_stress, optimal_result.Kt);
fprintf('  Stress reduction:         %.1f%%\n', opt_reduction);
fprintf('\n');
fprintf('COMPUTATION:\n');
fprintf('  Number of multi-starts:   %d\n', params.opt.num_starts);
fprintf('  Total function evals:     %d\n', total_feval);
fprintf('  Total optimization time:  %.1f seconds\n', total_opt_time);
fprintf('  Best run (start point):   %d\n', global_best_idx);
fprintf('\n');
fprintf('RESULTS SAVED TO: %s\n', params.opt.save_filename);
fprintf('=============================================================\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function f = objective_with_animation(x, params, anim_fig, baseline)
%OBJECTIVE_WITH_ANIMATION Compute objective and update animation

    global global_history

    % Convert 5-variable to 6-variable format
    x_full = [x(1:4); x(5); -x(5)];

    try
        % Run FEA analysis
        result = run_analysis_with_aux_holes(x_full, params, false);
        f = result.max_stress;

        % Update history
        global_history.feval_count = global_history.feval_count + 1;
        global_history.fval = [global_history.fval; f];
        global_history.x = [global_history.x; x'];

        if f < global_history.best_fval
            global_history.best_fval = f;
            global_history.best_x = x;
        end

        % Update animation (every 5 evaluations to reduce overhead)
        if params.anim.enabled && ishandle(anim_fig) && mod(global_history.feval_count, 5) == 0
            try
                update_animation(anim_fig, x_full, result, global_history, params, baseline);
            catch
                % Animation update failed, continue anyway
            end
        end

    catch
        % Return large penalty for failed geometries
        f = params.opt.penalty;
        global_history.feval_count = global_history.feval_count + 1;
        global_history.fval = [global_history.fval; f];
        global_history.x = [global_history.x; x'];
    end
end

function update_animation(fig, x_full, result, history, params, baseline)
%UPDATE_ANIMATION Update the animation figure

    figure(fig);

    % Subplot 1: Geometry
    subplot(1, 3, 1);
    cla;
    try
        [g, ~] = create_geometry(x_full, params);
        pdegplot(g, 'FaceLabels', 'off');
        hold on;
        plot(x_full(1), x_full(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(-x_full(1), x_full(2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        hold off;
    catch
    end
    axis equal;
    axis([-100, 100, -60, 60]);
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('Current Design (eval %d)', history.feval_count));
    grid on;

    % Subplot 2: Stress
    subplot(1, 3, 2);
    cla;
    try
        pdeplot(result.mesh, 'XYData', result.vonMises, 'ColorMap', 'jet');
        axis equal;
        xlim([-100, 100]);
        ylim([-60, 60]);
        colorbar;
    catch
    end
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('\\sigma_{vM,max} = %.1f MPa', result.max_stress));

    % Subplot 3: Convergence
    subplot(1, 3, 3);
    cla;
    if ~isempty(history.fval)
        plot(1:length(history.fval), history.fval, 'b-o', 'LineWidth', 1, 'MarkerSize', 2);
        hold on;
        yline(baseline.max_stress, 'r--', 'LineWidth', 1.5);
        hold off;
    end
    xlabel('Function Evaluation');
    ylabel('Max vM Stress (MPa)');
    title('Convergence');
    grid on;

    drawnow;
end

function [c, ceq] = constraint_function(x, params)
%CONSTRAINT_FUNCTION Compute nonlinear inequality constraints
%   x = [x_aux, y_aux, a, b, theta_R] (5 variables)
%   theta_L = -theta_R (symmetric constraint)

    % Extract design variables
    x_aux = x(1);
    y_aux = x(2);
    a = x(3);
    b = x(4);
    theta_R = x(5);
    theta_L = -theta_R;  % Symmetric constraint

    % Extract constraint parameters
    R_c = params.centralHole.radius;
    L = params.plate.length;
    W = params.plate.width;
    delta_c = params.constraints.min_gap_central;
    delta_e = params.constraints.min_gap_edge;
    alpha_max = params.constraints.max_aspect_ratio;

    % Compute AABB extents for rotated ellipses
    ex_R = sqrt(a^2*cos(theta_R)^2 + b^2*sin(theta_R)^2);
    ey_R = sqrt(a^2*sin(theta_R)^2 + b^2*cos(theta_R)^2);

    ex_L = sqrt(a^2*cos(theta_L)^2 + b^2*sin(theta_L)^2);
    ey_L = sqrt(a^2*sin(theta_L)^2 + b^2*cos(theta_L)^2);

    % Maximum extent (conservative for central hole clearance)
    e_max_R = max(ex_R, ey_R);
    e_max_L = max(ex_L, ey_L);

    % Distance from origin to auxiliary hole center
    d_c = sqrt(x_aux^2 + y_aux^2);

    % Initialize constraint vector
    c = zeros(7, 1);

    % g1: Right hole - central hole clearance
    c(1) = R_c + e_max_R + delta_c - d_c;

    % g2: Left hole - central hole clearance
    c(2) = R_c + e_max_L + delta_c - d_c;

    % g3: Right hole - right edge clearance
    c(3) = (x_aux + ex_R) - (L - delta_e);

    % g4: Left hole - left edge clearance
    c(4) = (x_aux + ex_L) - (L - delta_e);

    % g5: Right hole - top/bottom edge clearance
    c(5) = (abs(y_aux) + ey_R) - (W - delta_e);

    % g6: Left hole - top/bottom edge clearance
    c(6) = (abs(y_aux) + ey_L) - (W - delta_e);

    % g7: Aspect ratio constraint
    c(7) = a - alpha_max * b;

    % No equality constraints
    ceq = [];
end

function result = run_baseline_analysis(params)
%RUN_BASELINE_ANALYSIS Run FEA for plate with central hole only

    [g, ~] = create_baseline_geometry(params);

    % Create FE model
    model = femodel(AnalysisType="structuralStatic", Geometry=g);
    model.MaterialProperties = materialProperties(...
        YoungsModulus=params.material.E, ...
        PoissonsRatio=params.material.nu);

    % Boundary conditions
    model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
    model.EdgeBC(4) = edgeBC(XDisplacement=0);
    model.VertexBC(1) = vertexBC(YDisplacement=0);

    % Generate mesh
    Hmax = params.centralHole.radius / params.mesh.Hmax_factor;
    model = generateMesh(model, Hmax=Hmax);

    % Solve
    R = solve(model);

    % Compute von Mises stress
    stress = R.Stress;
    vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

    % Results
    result.max_stress = max(vonMises);
    result.Kt = result.max_stress / params.load.tension;
    result.vonMises = vonMises;
    result.mesh = R.Mesh;
    result.stress = stress;
end

function [g, status] = create_baseline_geometry(params)
%CREATE_BASELINE_GEOMETRY Create plate geometry with central hole only

    L = params.plate.length;
    W = params.plate.width;
    R_c = params.centralHole.radius;

    R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];
    C1 = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

    gdm = [R1, C1];
    ns = char('R1', 'C1')';
    sf = 'R1 - C1';

    [g, ~] = decsg(gdm, sf, ns);
    status = 0;
end

function result = run_analysis_with_aux_holes(x, params, verbose)
%RUN_ANALYSIS_WITH_AUX_HOLES Run FEA for plate with auxiliary holes
%   x = [x_aux, y_aux, a, b, theta_R, theta_L] (6 variables)

    if nargin < 3
        verbose = false;
    end

    % Create geometry
    [g, ~] = create_geometry(x, params);

    % Create FE model
    model = femodel(AnalysisType="structuralStatic", Geometry=g);
    model.MaterialProperties = materialProperties(...
        YoungsModulus=params.material.E, ...
        PoissonsRatio=params.material.nu);

    % Boundary conditions
    model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
    model.EdgeBC(4) = edgeBC(XDisplacement=0);
    model.VertexBC(1) = vertexBC(YDisplacement=0);

    % Generate mesh
    min_feature = min([x(3), x(4), params.constraints.min_gap_central]);
    Hmax = min_feature / params.mesh.Hmax_factor;
    model = generateMesh(model, Hmax=Hmax);

    % Solve
    R = solve(model);

    % Compute von Mises stress
    stress = R.Stress;
    vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

    % Find maximum stress location
    [max_stress, idx] = max(vonMises);
    nodes = R.Mesh.Nodes;
    loc_max = nodes(:, idx);

    % Results
    result.max_stress = max_stress;
    result.Kt = max_stress / params.load.tension;
    result.vonMises = vonMises;
    result.mesh = R.Mesh;
    result.loc_max = loc_max;
    result.stress = stress;

    if verbose
        fprintf('  Max von Mises stress: %.2f MPa at (%.2f, %.2f) mm\n', ...
            max_stress, loc_max(1), loc_max(2));
        fprintf('  Kt = %.3f\n', result.Kt);
    end
end

function result = run_detailed_analysis(x, params)
%RUN_DETAILED_ANALYSIS Run FEA with full output for final visualization

    % Create geometry
    [g, ~] = create_geometry(x, params);

    % Create FE model
    model = femodel(AnalysisType="structuralStatic", Geometry=g);
    model.MaterialProperties = materialProperties(...
        YoungsModulus=params.material.E, ...
        PoissonsRatio=params.material.nu);

    % Boundary conditions
    model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
    model.EdgeBC(4) = edgeBC(XDisplacement=0);
    model.VertexBC(1) = vertexBC(YDisplacement=0);

    % Generate mesh
    min_feature = min([x(3), x(4), params.constraints.min_gap_central]);
    Hmax = min_feature / params.mesh.Hmax_factor;
    model = generateMesh(model, Hmax=Hmax);

    % Solve
    R = solve(model);

    % Compute von Mises stress
    stress = R.Stress;
    vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

    % Find maximum stress location
    [max_stress, idx] = max(vonMises);
    nodes = R.Mesh.Nodes;
    loc_max = nodes(:, idx);

    % Results
    result.max_stress = max_stress;
    result.Kt = max_stress / params.load.tension;
    result.vonMises = vonMises;
    result.mesh = R.Mesh;
    result.loc_max = loc_max;
    result.stress = stress;
    result.Displacement = R.Displacement;
    result.R = R;  % Full result object for interpolation

    fprintf('  Max von Mises stress: %.2f MPa at (%.2f, %.2f) mm\n', ...
        max_stress, loc_max(1), loc_max(2));
    fprintf('  Kt = %.3f\n', result.Kt);
end

function [g, status] = create_geometry(x, params)
%CREATE_GEOMETRY Create plate geometry with central and auxiliary holes
%   x = [x_aux, y_aux, a, b, theta_R, theta_L] (6 variables)

    % Extract parameters
    L = params.plate.length;
    W = params.plate.width;
    R_c = params.centralHole.radius;
    n_pts = params.mesh.n_ellipse_pts;

    % Extract design variables
    x_aux = x(1);
    y_aux = x(2);
    a = x(3);
    b = x(4);
    theta_R = x(5);
    theta_L = x(6);

    % Rectangle (plate)
    R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];

    % Central circular hole
    C_central = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

    % Create elliptical auxiliary holes as polygons
    theta_pts = linspace(0, 2*pi, n_pts + 1);
    theta_pts = theta_pts(1:end-1);

    % Local ellipse coordinates
    x_local = a * cos(theta_pts);
    y_local = b * sin(theta_pts);

    % Right auxiliary hole (apply rotation and translation)
    cos_R = cos(theta_R);
    sin_R = sin(theta_R);
    x_right = cos_R * x_local - sin_R * y_local + x_aux;
    y_right = sin_R * x_local + cos_R * y_local + y_aux;

    % Left auxiliary hole (mirror x, apply rotation)
    cos_L = cos(theta_L);
    sin_L = sin(theta_L);
    x_left = -(cos_L * x_local - sin_L * y_local) - x_aux;
    y_left = sin_L * x_local + cos_L * y_local + y_aux;

    % Create polygon geometry matrices
    P_right = [2; n_pts; x_right(:); y_right(:)];
    P_left = [2; n_pts; x_left(:); y_left(:)];

    % Pad all matrices to same length
    max_len = max([length(R1), length(C_central), length(P_right), length(P_left)]);
    R1 = [R1; zeros(max_len - length(R1), 1)];
    C_central = [C_central; zeros(max_len - length(C_central), 1)];
    P_right = [P_right; zeros(max_len - length(P_right), 1)];
    P_left = [P_left; zeros(max_len - length(P_left), 1)];

    % Combine geometry
    gdm = [R1, C_central, P_right, P_left];
    ns = char('R1', 'C1', 'P_right', 'P_left')';
    sf = 'R1 - C1 - P_right - P_left';

    % Create decomposed geometry
    [g, ~] = decsg(gdm, sf, ns);
    status = 0;
end
