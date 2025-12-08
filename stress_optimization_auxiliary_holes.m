%% Shape Optimization of Plate with Stress-Relieving Auxiliary Holes
% This script optimizes the placement, size, and orientation of two
% elliptical auxiliary holes to minimize the maximum von Mises stress
% in a plate with a central circular hole.
%
% Optimization formulation follows: Optimization_Problem_Formulation.md
% FEA methodology based on: stress_analysis_auxiliary_holes.m
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% ========================================================================
%  SECTION 1: PROBLEM PARAMETERS (from Optimization_Problem_Formulation.md)
%  ========================================================================

fprintf('=============================================================\n');
fprintf('  SHAPE OPTIMIZATION: Stress-Relieving Auxiliary Holes\n');
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
params.opt.num_starts = 10;                   % Number of multi-start points
params.opt.penalty = 1e10;                    % Penalty for infeasible/failed geometries

%% ========================================================================
%  SECTION 2: DESIGN VARIABLE BOUNDS (from Optimization_Problem_Formulation.md)
%  ========================================================================

% Design vector: x = [x_aux, y_aux, a, b, theta_R, theta_L]

% Lower bounds
lb = [
    params.centralHole.radius + params.constraints.max_aux_size + params.constraints.min_gap_central;  % x_aux >= 38 mm
    0;                                                    % y_aux >= 0
    params.constraints.min_aux_size;                      % a >= 6 mm
    params.constraints.min_aux_size;                      % b >= 6 mm
    -pi/2;                                                % theta_R >= -90 deg
    -pi/2                                                 % theta_L >= -90 deg
];

% Upper bounds
ub = [
    params.constraints.max_x_aux;                         % x_aux <= 150 mm
    params.plate.width - params.constraints.min_gap_edge - params.constraints.max_aux_size;  % y_aux <= 30 mm
    params.constraints.max_aux_size;                      % a <= 15 mm
    params.constraints.max_aux_size;                      % b <= 15 mm
    pi/2;                                                 % theta_R <= 90 deg
    pi/2                                                  % theta_L <= 90 deg
];

fprintf('Design Variable Bounds:\n');
fprintf('  x_aux:   [%.1f, %.1f] mm\n', lb(1), ub(1));
fprintf('  y_aux:   [%.1f, %.1f] mm\n', lb(2), ub(2));
fprintf('  a:       [%.1f, %.1f] mm\n', lb(3), ub(3));
fprintf('  b:       [%.1f, %.1f] mm\n', lb(4), ub(4));
fprintf('  theta_R: [%.2f, %.2f] rad\n', lb(5), ub(5));
fprintf('  theta_L: [%.2f, %.2f] rad\n', lb(6), ub(6));
fprintf('\n');

%% ========================================================================
%  SECTION 3: BASELINE ANALYSIS (Central Hole Only)
%  ========================================================================

fprintf('=== Running Baseline Analysis (Central Hole Only) ===\n');

baseline = run_baseline_analysis(params);

fprintf('Baseline Results:\n');
fprintf('  Max von Mises stress: %.2f MPa\n', baseline.max_stress);
fprintf('  Stress concentration factor (Kt): %.3f\n', baseline.Kt);
fprintf('  Theoretical Kt: ~3.0\n\n');

%% ========================================================================
%  SECTION 4: OPTIMIZATION SETUP
%  ========================================================================

% --- Objective Function ---
% Wrapper that captures params and tracks iteration history
iteration_history = struct('x', [], 'fval', [], 'count', 0);

objective = @(x) objective_function(x, params);

% --- Nonlinear Constraints ---
nonlcon = @(x) constraint_function(x, params);

% --- Optimizer Options ---
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1200, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4, ...
    'FiniteDifferenceType', 'central', ...
    'OutputFcn', @(x,optimValues,state) output_function(x, optimValues, state));

%% ========================================================================
%  SECTION 5: MULTI-START OPTIMIZATION
%  ========================================================================

fprintf('=== Starting Multi-Start Optimization ===\n');
fprintf('Number of starting points: %d\n\n', params.opt.num_starts);

% Storage for results from each start
all_results = struct('x0', cell(params.opt.num_starts, 1), ...
                     'x_opt', cell(params.opt.num_starts, 1), ...
                     'fval', zeros(params.opt.num_starts, 1), ...
                     'exitflag', zeros(params.opt.num_starts, 1), ...
                     'history', cell(params.opt.num_starts, 1));

% Generate initial points using Latin Hypercube Sampling
rng(42);  % For reproducibility
lhs_samples = lhsdesign(params.opt.num_starts, 6);
initial_points = lb' + lhs_samples .* (ub' - lb');

% Global best tracking
global_best_fval = Inf;
global_best_x = [];
global_best_idx = 0;

% Run optimization from each starting point
for i = 1:params.opt.num_starts
    fprintf('--- Start Point %d/%d ---\n', i, params.opt.num_starts);

    x0 = initial_points(i, :)';

    % Check if starting point is feasible
    [c0, ~] = nonlcon(x0);
    if any(c0 > 0)
        fprintf('  Initial point infeasible, adjusting...\n');
        % Simple adjustment: move toward center of bounds
        x0 = (lb + ub) / 2;
        x0(5) = 0;  % theta_R = 0
        x0(6) = 0;  % theta_L = 0
    end

    % Reset iteration history for this run
    iteration_history.x = [];
    iteration_history.fval = [];
    iteration_history.count = 0;

    % Run fmincon
    try
        [x_opt, fval, exitflag, output] = fmincon(objective, x0, ...
            [], [], [], [], lb, ub, nonlcon, options);

        all_results(i).x0 = x0;
        all_results(i).x_opt = x_opt;
        all_results(i).fval = fval;
        all_results(i).exitflag = exitflag;
        all_results(i).history = iteration_history;

        fprintf('  Final objective: %.2f MPa\n', fval);
        fprintf('  Exit flag: %d\n\n', exitflag);

        % Update global best
        if fval < global_best_fval
            global_best_fval = fval;
            global_best_x = x_opt;
            global_best_idx = i;
        end

    catch ME
        fprintf('  Optimization failed: %s\n\n', ME.message);
        all_results(i).fval = Inf;
        all_results(i).exitflag = -99;
    end
end

fprintf('=== Optimization Complete ===\n');
fprintf('Best result from start point %d\n', global_best_idx);
fprintf('Optimal objective: %.2f MPa\n\n', global_best_fval);

%% ========================================================================
%  SECTION 6: EXTRACT AND DISPLAY OPTIMAL SOLUTION
%  ========================================================================

% Extract optimal design variables
x_opt = global_best_x;
opt.x_aux = x_opt(1);
opt.y_aux = x_opt(2);
opt.a = x_opt(3);
opt.b = x_opt(4);
opt.theta_R = x_opt(5);
opt.theta_L = x_opt(6);

fprintf('=== Optimal Design Variables ===\n');
fprintf('  x_aux:   %.3f mm\n', opt.x_aux);
fprintf('  y_aux:   %.3f mm\n', opt.y_aux);
fprintf('  a:       %.3f mm\n', opt.a);
fprintf('  b:       %.3f mm\n', opt.b);
fprintf('  theta_R: %.3f rad (%.1f deg)\n', opt.theta_R, rad2deg(opt.theta_R));
fprintf('  theta_L: %.3f rad (%.1f deg)\n', opt.theta_L, rad2deg(opt.theta_L));
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
for j = 1:length(c_opt)
    status = 'OK';
    if c_opt(j) > 0
        status = 'VIOLATED';
    end
    fprintf('  %s: %.3f [%s]\n', constraint_names{j}, c_opt(j), status);
end
fprintf('\n');

%% ========================================================================
%  SECTION 7: RUN FINAL ANALYSIS ON OPTIMAL DESIGN
%  ========================================================================

fprintf('=== Running Final Analysis on Optimal Design ===\n');

% Run detailed analysis on optimal configuration
optimal_result = run_analysis_with_aux_holes(x_opt, params, true);

%% ========================================================================
%  SECTION 8: COMPARISON TABLE
%  ========================================================================

% Initial guess results (from best starting point)
x0_best = all_results(global_best_idx).x0;
initial_result = run_analysis_with_aux_holes(x0_best, params, false);

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    COMPARISON TABLE\n');
fprintf('=============================================================\n');
fprintf('%-25s | %-15s | %-15s | %-15s\n', 'Metric', 'Baseline', 'Initial Guess', 'Optimal');
fprintf('%-25s | %-15s | %-15s | %-15s\n', repmat('-',1,25), repmat('-',1,15), repmat('-',1,15), repmat('-',1,15));
fprintf('%-25s | %-15.2f | %-15.2f | %-15.2f\n', 'Max vM Stress (MPa)', baseline.max_stress, initial_result.max_stress, optimal_result.max_stress);
fprintf('%-25s | %-15.3f | %-15.3f | %-15.3f\n', 'Stress Conc. Factor Kt', baseline.Kt, initial_result.Kt, optimal_result.Kt);
fprintf('%-25s | %-15s | %-15.1f | %-15.1f\n', 'Reduction vs Baseline (%)', 'N/A', ...
    (baseline.max_stress - initial_result.max_stress)/baseline.max_stress * 100, ...
    (baseline.max_stress - optimal_result.max_stress)/baseline.max_stress * 100);
fprintf('=============================================================\n\n');

%% ========================================================================
%  SECTION 9: CONVERGENCE PLOT
%  ========================================================================

fprintf('=== Generating Convergence Plot ===\n');

figure('Name', 'Optimization Convergence', 'Position', [100, 100, 800, 500]);

% Get history from best run
best_history = all_results(global_best_idx).history;

if ~isempty(best_history.fval)
    iterations = 1:length(best_history.fval);

    plot(iterations, best_history.fval, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    hold on;

    % Add baseline reference line
    yline(baseline.max_stress, 'r--', 'LineWidth', 1.5, 'Label', 'Baseline (no aux holes)');

    % Mark optimal point
    [min_fval, min_idx] = min(best_history.fval);
    plot(min_idx, min_fval, 'g*', 'MarkerSize', 15, 'LineWidth', 2);

    hold off;

    xlabel('Iteration');
    ylabel('Max von Mises Stress (MPa)');
    title('Optimization Convergence (Best Run)');
    legend('Objective', 'Baseline', 'Optimal', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No convergence history available', 'HorizontalAlignment', 'center');
end

%% ========================================================================
%  SECTION 10: VISUALIZATION OF OPTIMAL DESIGN
%  ========================================================================

fprintf('=== Generating Optimal Design Visualization ===\n');

% Create geometry for optimal design
[g_opt, ~] = create_geometry(x_opt, params);

figure('Name', 'Optimal Design', 'Position', [100, 100, 1200, 500]);

% Plot 1: Geometry
subplot(1, 2, 1);
pdegplot(g_opt, 'FaceLabels', 'off');
hold on;

% Mark hole centers
plot(opt.x_aux, opt.y_aux, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Right aux center');
plot(-opt.x_aux, opt.y_aux, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Left aux center');
plot(0, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Central hole center');

% Draw ellipse axes
scale = 1.0;
quiver(opt.x_aux, opt.y_aux, scale*opt.a*cos(opt.theta_R), scale*opt.a*sin(opt.theta_R), 'r', 'LineWidth', 1.5, 'AutoScale', 'off');
quiver(-opt.x_aux, opt.y_aux, -scale*opt.a*cos(opt.theta_L), scale*opt.a*sin(opt.theta_L), 'b', 'LineWidth', 1.5, 'AutoScale', 'off');

hold off;
axis equal;
axis([-1.1*params.plate.length, 1.1*params.plate.length, -1.1*params.plate.width, 1.1*params.plate.width]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Optimal Auxiliary Hole Configuration');
legend('Location', 'best');
grid on;

% Plot 2: von Mises stress
subplot(1, 2, 2);
pdeplot(optimal_result.mesh, 'XYData', optimal_result.vonMises, 'ColorMap', 'jet');
hold on;
plot(optimal_result.loc_max(1), optimal_result.loc_max(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('von Mises Stress (Max = %.1f MPa)', optimal_result.max_stress));

%% ========================================================================
%  SECTION 11: SUMMARY OUTPUT
%  ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    OPTIMIZATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('\n');
fprintf('OPTIMAL DESIGN:\n');
fprintf('  Auxiliary hole position:  (%.2f, %.2f) mm from center\n', opt.x_aux, opt.y_aux);
fprintf('  Ellipse semi-axes:        a = %.2f mm, b = %.2f mm\n', opt.a, opt.b);
fprintf('  Aspect ratio:             %.2f\n', opt.a/opt.b);
fprintf('  Right hole rotation:      %.1f degrees\n', rad2deg(opt.theta_R));
fprintf('  Left hole rotation:       %.1f degrees\n', rad2deg(opt.theta_L));
fprintf('\n');
fprintf('PERFORMANCE:\n');
fprintf('  Baseline max stress:      %.2f MPa (Kt = %.3f)\n', baseline.max_stress, baseline.Kt);
fprintf('  Optimal max stress:       %.2f MPa (Kt = %.3f)\n', optimal_result.max_stress, optimal_result.Kt);
fprintf('  Stress reduction:         %.1f%%\n', (baseline.max_stress - optimal_result.max_stress)/baseline.max_stress * 100);
fprintf('\n');
fprintf('=============================================================\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function f = objective_function(x, params)
%OBJECTIVE_FUNCTION Compute maximum von Mises stress for given design
%   f = objective_function(x, params)
%   x = [x_aux, y_aux, a, b, theta_R, theta_L]
%   Returns maximum von Mises stress or penalty if geometry fails

    try
        % Run FEA analysis
        result = run_analysis_with_aux_holes(x, params, false);
        f = result.max_stress;
    catch
        % Return large penalty for failed geometries
        f = params.opt.penalty;
    end
end

function [c, ceq] = constraint_function(x, params)
%CONSTRAINT_FUNCTION Compute nonlinear inequality constraints
%   [c, ceq] = constraint_function(x, params)
%   c(i) <= 0 for feasible design
%   Based on Optimization_Problem_Formulation.md Section 4.5

    % Extract design variables
    x_aux = x(1);
    y_aux = x(2);
    a = x(3);
    b = x(4);
    theta_R = x(5);
    theta_L = x(6);

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
%RUN_BASELINE_ANALYSIS Run FEA for plate with central hole only (no aux holes)

    % Create geometry: Rectangle minus central hole
    L = params.plate.length;
    W = params.plate.width;
    R_c = params.centralHole.radius;

    R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];
    C1 = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

    gdm = [R1, C1];
    ns = char('R1', 'C1')';
    sf = 'R1 - C1';

    g = decsg(gdm, sf, ns);

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
    Hmax = R_c / params.mesh.Hmax_factor;
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
end

function result = run_analysis_with_aux_holes(x, params, verbose)
%RUN_ANALYSIS_WITH_AUX_HOLES Run FEA for plate with central and auxiliary holes
%   x = [x_aux, y_aux, a, b, theta_R, theta_L]

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

    % Boundary conditions (edge numbers may vary with geometry)
    % Find edges at x = L and x = -L
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

function [g, status] = create_geometry(x, params)
%CREATE_GEOMETRY Create plate geometry with central and auxiliary holes
%   x = [x_aux, y_aux, a, b, theta_R, theta_L]

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

function stop = output_function(x, optimValues, state)
%OUTPUT_FUNCTION Callback to record iteration history

    persistent history_fval history_x

    stop = false;

    switch state
        case 'init'
            history_fval = [];
            history_x = [];
        case 'iter'
            history_fval = [history_fval; optimValues.fval];
            history_x = [history_x; x'];

            % Store in base workspace for access after optimization
            assignin('base', 'iteration_history', struct('fval', history_fval, 'x', history_x, 'count', length(history_fval)));
        case 'done'
            % Final update
            assignin('base', 'iteration_history', struct('fval', history_fval, 'x', history_x, 'count', length(history_fval)));
    end
end
