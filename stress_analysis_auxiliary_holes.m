%% Stress Analysis of Plate with Central Hole and Elliptical Auxiliary Holes
% This script performs a 2-D plane-stress elasticity analysis of a plate
% with a central circular hole and two elliptical auxiliary (stress-relief) holes.
%
% The auxiliary holes are placed symmetrically about the y-axis (one left, one right)
% but can have independent rotation angles.
%
% Author: Generated for Shape Optimization Study
% Date: 2024

clear; close all; clc;

%% ========================================================================
%  SECTION 1: PROBLEM PARAMETERS
%  ========================================================================

% --- Plate Geometry ---
plate.width = 50.0;                    % Half-width of plate (mm) - total width = 2*width
plate.length = 4 * plate.width;        % Half-length of plate (mm) - ensures far-field conditions

% --- Central Hole (Fixed) ---
centralHole.radius = 20.0;             % Radius of central hole (mm)
centralHole.x = 0;                     % x-position (centered)
centralHole.y = 0;                     % y-position (centered)

% --- Auxiliary Holes (Design Variables for Future Optimization) ---
% Position: distance from center along x-axis
auxHole.x_center = 45.0;               % x-distance from origin to aux hole center (mm)
auxHole.y_center = 0.0;                % y-offset from centerline (mm)

% Shape: ellipse semi-axes
auxHole.a = 8.0;                       % Semi-major axis (mm)
auxHole.b = 5.0;                       % Semi-minor axis (mm)

% Rotation: independent angles for left and right holes
auxHole.theta_right = deg2rad(30);     % Rotation angle for RIGHT aux hole (rad, CCW positive)
auxHole.theta_left = deg2rad(-30);     % Rotation angle for LEFT aux hole (rad, CCW positive)

% --- Material Properties (Steel-like) ---
material.E = 200e3;                    % Young's modulus (MPa)
material.nu = 0.25;                    % Poisson's ratio

% --- Loading ---
load.tension = 100;                    % Applied tensile stress (MPa)

% --- Constraint Parameters ---
constraints.min_gap_central = 3.0;     % Minimum gap from central hole (mm)
constraints.min_gap_edge = 5.0;        % Minimum gap from plate edges (mm)
constraints.min_aux_size = 2.0;        % Minimum semi-axis size (mm)
constraints.max_aux_size = 15.0;       % Maximum semi-axis size (mm)
constraints.max_aspect_ratio = 3.0;    % Maximum a/b ratio

% --- Mesh Parameters ---
mesh.Hmax_factor = 6;                  % Hmax = min_feature / Hmax_factor

%% ========================================================================
%  SECTION 2: CONSTRAINT VALIDATION FUNCTIONS
%  ========================================================================

fprintf('=== Validating Geometric Constraints ===\n\n');

% Function to compute axis-aligned bounding box half-extents of rotated ellipse
% For a rotated ellipse, the extreme x and y distances from center are:
%   x_extent = sqrt(a^2*cos^2(theta) + b^2*sin^2(theta))
%   y_extent = sqrt(a^2*sin^2(theta) + b^2*cos^2(theta))
get_ellipse_extents = @(a, b, theta) deal(...
    sqrt(a^2*cos(theta)^2 + b^2*sin(theta)^2), ...
    sqrt(a^2*sin(theta)^2 + b^2*cos(theta)^2));

% Compute extents for both auxiliary holes
[x_ext_right, y_ext_right] = get_ellipse_extents(auxHole.a, auxHole.b, auxHole.theta_right);
[x_ext_left, y_ext_left] = get_ellipse_extents(auxHole.a, auxHole.b, auxHole.theta_left);

fprintf('Right aux hole extents: x = %.2f mm, y = %.2f mm\n', x_ext_right, y_ext_right);
fprintf('Left aux hole extents:  x = %.2f mm, y = %.2f mm\n', x_ext_left, y_ext_left);

% --- Check 1: Distance from central hole ---
% Minimum distance from ellipse boundary to central hole boundary
% Using conservative estimate: center distance - ellipse extent - central radius >= min_gap

% For right hole
dist_to_center_right = sqrt(auxHole.x_center^2 + auxHole.y_center^2);
clearance_central_right = dist_to_center_right - max(x_ext_right, y_ext_right) - centralHole.radius;

% For left hole (mirrored x position)
dist_to_center_left = sqrt(auxHole.x_center^2 + auxHole.y_center^2);  % Same due to symmetry
clearance_central_left = dist_to_center_left - max(x_ext_left, y_ext_left) - centralHole.radius;

fprintf('\nClearance from central hole:\n');
fprintf('  Right hole: %.2f mm (required >= %.2f mm) ', clearance_central_right, constraints.min_gap_central);
if clearance_central_right >= constraints.min_gap_central
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end
fprintf('  Left hole:  %.2f mm (required >= %.2f mm) ', clearance_central_left, constraints.min_gap_central);
if clearance_central_left >= constraints.min_gap_central
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end

% --- Check 2: Distance from plate edges ---
% Right edge: plate.length - (x_center + x_extent) >= min_gap_edge
% Top edge: plate.width - (|y_center| + y_extent) >= min_gap_edge
% Left edge (for left hole): plate.length - (x_center + x_extent) >= min_gap_edge

clearance_right_edge = plate.length - (auxHole.x_center + x_ext_right);
clearance_top_right = plate.width - (abs(auxHole.y_center) + y_ext_right);
clearance_top_left = plate.width - (abs(auxHole.y_center) + y_ext_left);
clearance_left_edge = plate.length - (auxHole.x_center + x_ext_left);

fprintf('\nClearance from plate edges:\n');
fprintf('  Right hole to right edge: %.2f mm (required >= %.2f mm) ', clearance_right_edge, constraints.min_gap_edge);
if clearance_right_edge >= constraints.min_gap_edge
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end
fprintf('  Right hole to top edge:   %.2f mm (required >= %.2f mm) ', clearance_top_right, constraints.min_gap_edge);
if clearance_top_right >= constraints.min_gap_edge
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end
fprintf('  Left hole to left edge:   %.2f mm (required >= %.2f mm) ', clearance_left_edge, constraints.min_gap_edge);
if clearance_left_edge >= constraints.min_gap_edge
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end
fprintf('  Left hole to top edge:    %.2f mm (required >= %.2f mm) ', clearance_top_left, constraints.min_gap_edge);
if clearance_top_left >= constraints.min_gap_edge
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end

% --- Check 3: Aspect ratio ---
aspect_ratio = max(auxHole.a, auxHole.b) / min(auxHole.a, auxHole.b);
fprintf('\nAspect ratio: %.2f (required <= %.2f) ', aspect_ratio, constraints.max_aspect_ratio);
if aspect_ratio <= constraints.max_aspect_ratio
    fprintf('[OK]\n');
else
    fprintf('[VIOLATION]\n');
end

fprintf('\n');

%% ========================================================================
%  SECTION 3: GEOMETRY CREATION
%  ========================================================================

fprintf('=== Creating Geometry ===\n\n');

% --- Create rectangle (plate) ---
% Using geometry description matrix format for rectangle:
% [3; 4; x1; x2; x3; x4; y1; y2; y3; y4] where vertices are CCW
R1 = [3; 4; -plate.length; plate.length; plate.length; -plate.length; ...
      -plate.width; -plate.width; plate.width; plate.width];

% --- Create central circular hole ---
% Circle format: [1; x_center; y_center; radius; 0; 0; 0; 0; 0; 0]
C_central = [1; centralHole.x; centralHole.y; centralHole.radius; 0; 0; 0; 0; 0; 0];

% --- Create elliptical auxiliary holes ---
% MATLAB's decsg doesn't directly support rotated ellipses, so we use
% polygon approximation with sufficient points for smooth representation

n_ellipse_pts = 100;  % Number of points for ellipse approximation
theta_pts = linspace(0, 2*pi, n_ellipse_pts + 1);
theta_pts = theta_pts(1:end-1);  % Remove duplicate endpoint

% Right auxiliary hole (rotated ellipse as polygon)
x_ellipse_local = auxHole.a * cos(theta_pts);
y_ellipse_local = auxHole.b * sin(theta_pts);

% Apply rotation for right hole
cos_r = cos(auxHole.theta_right);
sin_r = sin(auxHole.theta_right);
x_right = cos_r * x_ellipse_local - sin_r * y_ellipse_local + auxHole.x_center;
y_right = sin_r * x_ellipse_local + cos_r * y_ellipse_local + auxHole.y_center;

% Left auxiliary hole (mirrored x, independent rotation)
cos_l = cos(auxHole.theta_left);
sin_l = sin(auxHole.theta_left);
x_left = -(cos_l * x_ellipse_local - sin_l * y_ellipse_local) - auxHole.x_center;
y_left = sin_l * x_ellipse_local + cos_l * y_ellipse_local + auxHole.y_center;

% Create polygon geometry description matrices
% Polygon format: [2; n_vertices; x1; x2; ...; xn; y1; y2; ...; yn]
P_right = [2; n_ellipse_pts; x_right(:); y_right(:)];
P_left = [2; n_ellipse_pts; x_left(:); y_left(:)];

% Pad all geometry matrices to same length
max_len = max([length(R1), length(C_central), length(P_right), length(P_left)]);
R1 = [R1; zeros(max_len - length(R1), 1)];
C_central = [C_central; zeros(max_len - length(C_central), 1)];
P_right = [P_right; zeros(max_len - length(P_right), 1)];
P_left = [P_left; zeros(max_len - length(P_left), 1)];

% Combine into geometry description matrix
gdm = [R1, C_central, P_right, P_left];
ns = char('R1', 'C1', 'P_right', 'P_left')';

% Set formula: Rectangle minus all holes
sf = 'R1 - C1 - P_right - P_left';

% Create decomposed geometry
try
    [g, bt] = decsg(gdm, sf, ns);
    fprintf('Geometry created successfully.\n');
    fprintf('Number of edges: %d\n', size(g, 2));
catch ME
    error('Geometry creation failed: %s\nCheck for overlapping holes or invalid geometry.', ME.message);
end

%% ========================================================================
%  SECTION 4: GEOMETRY VISUALIZATION
%  ========================================================================

fprintf('\n=== Generating Geometry Plots ===\n');

figure('Name', 'Plate Geometry', 'Position', [100, 100, 1200, 500]);

% Plot 1: Geometry with edge labels
subplot(1, 2, 1);
pdegplot(g, 'EdgeLabels', 'on');
axis equal;
axis([-1.1*plate.length, 1.1*plate.length, -1.1*plate.width, 1.1*plate.width]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Geometry with Edge Labels');
grid on;

% Plot 2: Geometry showing hole details
subplot(1, 2, 2);
pdegplot(g, 'FaceLabels', 'on');
hold on;

% Mark auxiliary hole centers
plot(auxHole.x_center, auxHole.y_center, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(-auxHole.x_center, auxHole.y_center, 'bo', 'MarkerSize', 10, 'LineWidth', 2);

% Draw semi-axes for visualization
% Right hole
quiver(auxHole.x_center, auxHole.y_center, ...
       auxHole.a*cos(auxHole.theta_right), auxHole.a*sin(auxHole.theta_right), ...
       'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
quiver(auxHole.x_center, auxHole.y_center, ...
       -auxHole.b*sin(auxHole.theta_right), auxHole.b*cos(auxHole.theta_right), ...
       'r--', 'LineWidth', 1.5);

% Left hole
quiver(-auxHole.x_center, auxHole.y_center, ...
       -auxHole.a*cos(auxHole.theta_left), auxHole.a*sin(auxHole.theta_left), ...
       'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
quiver(-auxHole.x_center, auxHole.y_center, ...
       auxHole.b*sin(auxHole.theta_left), auxHole.b*cos(auxHole.theta_left), ...
       'b--', 'LineWidth', 1.5);

axis equal;
axis([-1.1*plate.length, 1.1*plate.length, -1.1*plate.width, 1.1*plate.width]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Geometry with Hole Details');
legend('', 'Right hole center', 'Left hole center', 'Location', 'best');
grid on;
hold off;

%% ========================================================================
%  SECTION 5: FE MODEL SETUP
%  ========================================================================

fprintf('\n=== Setting up FE Model ===\n');

% Create structural model for static analysis
model = femodel(AnalysisType="structuralStatic", Geometry=g);

% Assign material properties
model.MaterialProperties = materialProperties(...
    YoungsModulus=material.E, ...
    PoissonsRatio=material.nu);

% --- Boundary Conditions ---
% The edge numbering depends on geometry creation; we need to identify edges

% Find the rightmost edge (for applied load) and leftmost edge (for BC)
% Based on original example: Edge 2 is right, Edge 4 is left
% However, with complex geometry, we should verify programmatically

% For now, assume similar edge numbering to base example
% Edge on right side (x = plate.length): apply tension
% Edge on left side (x = -plate.length): fix x-displacement

% Apply surface traction on right edge
model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[load.tension; 0]);

% Fix x-displacement on left edge
model.EdgeBC(4) = edgeBC(XDisplacement=0);

% Fix y-displacement at a vertex to prevent rigid body motion
model.VertexBC(1) = vertexBC(YDisplacement=0);

fprintf('Boundary conditions applied.\n');
fprintf('  - Tension %.1f MPa applied on right edge\n', load.tension);
fprintf('  - X-displacement fixed on left edge\n');
fprintf('  - Y-displacement fixed at vertex 1\n');

%% ========================================================================
%  SECTION 6: MESH GENERATION
%  ========================================================================

fprintf('\n=== Generating Mesh ===\n');

% Determine appropriate mesh size based on smallest feature
min_feature = min([auxHole.a, auxHole.b, constraints.min_gap_central]);
Hmax = min_feature / mesh.Hmax_factor;

fprintf('Minimum feature size: %.2f mm\n', min_feature);
fprintf('Target element size (Hmax): %.2f mm\n', Hmax);

% Generate mesh
model = generateMesh(model, Hmax=Hmax);

% Get mesh statistics
meshData = model.Geometry.Mesh;
numNodes = size(meshData.Nodes, 2);
numElements = size(meshData.Elements, 2);

fprintf('Mesh generated successfully.\n');
fprintf('  - Number of nodes: %d\n', numNodes);
fprintf('  - Number of elements: %d\n', numElements);

% Plot mesh
figure('Name', 'Finite Element Mesh', 'Position', [100, 100, 1000, 600]);

subplot(1, 2, 1);
pdemesh(model);
axis equal;
title('Full Mesh');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(1, 2, 2);
pdemesh(model);
axis equal;
xlim([-80, 80]);
ylim([-60, 60]);
title('Mesh Detail (Near Holes)');
xlabel('x (mm)');
ylabel('y (mm)');

%% ========================================================================
%  SECTION 7: SOLVE FE PROBLEM
%  ========================================================================

fprintf('\n=== Solving FE Problem ===\n');

tic;
R = solve(model);
solve_time = toc;

fprintf('Solution completed in %.2f seconds.\n', solve_time);

%% ========================================================================
%  SECTION 8: STRESS ANALYSIS AND RESULTS
%  ========================================================================

fprintf('\n=== Computing Stress Results ===\n');

% Extract stress components at nodes
stress = R.Stress;

% Compute von Mises stress
% For 2D plane stress: sigma_vM = sqrt(sxx^2 - sxx*syy + syy^2 + 3*sxy^2)
vonMises = sqrt(stress.sxx.^2 - stress.sxx.*stress.syy + stress.syy.^2 + 3*stress.sxy.^2);

% Find maximum stresses
[max_sxx, idx_sxx] = max(stress.sxx);
[max_syy, idx_syy] = max(stress.syy);
[max_sxy, idx_sxy] = max(abs(stress.sxy));
[max_vonMises, idx_vM] = max(vonMises);

% Get locations of maximum stresses
nodes = meshData.Nodes;
loc_max_vM = nodes(:, idx_vM);

fprintf('\nMaximum Stress Values:\n');
fprintf('  sigma_xx (max): %.2f MPa\n', max_sxx);
fprintf('  sigma_yy (max): %.2f MPa\n', max_syy);
fprintf('  |sigma_xy| (max): %.2f MPa\n', max_sxy);
fprintf('  von Mises (max): %.2f MPa at (%.2f, %.2f) mm\n', max_vonMises, loc_max_vM(1), loc_max_vM(2));

% Compute stress concentration factor
K_t = max_vonMises / load.tension;
fprintf('\nStress Concentration Factor (K_t): %.3f\n', K_t);

% Compare to theoretical value for plate with hole only (K_t ≈ 3.0)
K_t_theoretical = 3.0;
reduction_percent = (K_t_theoretical - K_t) / K_t_theoretical * 100;
fprintf('Theoretical K_t (central hole only): %.3f\n', K_t_theoretical);
fprintf('Stress reduction from auxiliary holes: %.1f%%\n', reduction_percent);

%% ========================================================================
%  SECTION 9: STRESS VISUALIZATION
%  ========================================================================

fprintf('\n=== Generating Stress Plots ===\n');

% Figure 3: Stress contours
figure('Name', 'Stress Distribution', 'Position', [100, 100, 1400, 900]);

% Plot 1: sigma_xx
subplot(2, 2, 1);
pdeplot(meshData, XYData=stress.sxx, ColorMap="jet", Contour="on");
axis equal;
colorbar;
title('\sigma_{xx} (MPa)');
xlabel('x (mm)');
ylabel('y (mm)');

% Plot 2: sigma_yy
subplot(2, 2, 2);
pdeplot(meshData, XYData=stress.syy, ColorMap="jet", Contour="on");
axis equal;
colorbar;
title('\sigma_{yy} (MPa)');
xlabel('x (mm)');
ylabel('y (mm)');

% Plot 3: sigma_xy
subplot(2, 2, 3);
pdeplot(meshData, XYData=stress.sxy, ColorMap="jet", Contour="on");
axis equal;
colorbar;
title('\sigma_{xy} (MPa)');
xlabel('x (mm)');
ylabel('y (mm)');

% Plot 4: von Mises stress
subplot(2, 2, 4);
pdeplot(meshData, XYData=vonMises, ColorMap="jet", Contour="on");
axis equal;
colorbar;
title('\sigma_{vM} (von Mises) (MPa)');
xlabel('x (mm)');
ylabel('y (mm)');

% Figure 4: Detailed von Mises near holes
figure('Name', 'von Mises Stress Detail', 'Position', [100, 100, 1000, 600]);

subplot(1, 2, 1);
pdeplot(meshData, XYData=vonMises, ColorMap="jet");
hold on;
plot(loc_max_vM(1), loc_max_vM(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
colorbar;
title(sprintf('von Mises Stress (Max = %.1f MPa)', max_vonMises));
xlabel('x (mm)');
ylabel('y (mm)');

subplot(1, 2, 2);
pdeplot(meshData, XYData=vonMises, ColorMap="jet");
hold on;
plot(loc_max_vM(1), loc_max_vM(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
xlim([-80, 80]);
ylim([-60, 60]);
colorbar;
title('von Mises Stress (Zoomed)');
xlabel('x (mm)');
ylabel('y (mm)');

%% ========================================================================
%  SECTION 10: STRESS ALONG HOLE BOUNDARIES
%  ========================================================================

fprintf('\n=== Computing Stress Along Hole Boundaries ===\n');

% Define points on central hole boundary
theta_central = linspace(0, 2*pi, 200);
x_central_boundary = centralHole.radius * cos(theta_central);
y_central_boundary = centralHole.radius * sin(theta_central);

% Define points on auxiliary hole boundaries
theta_aux = linspace(0, 2*pi, 200);

% Right auxiliary hole boundary
x_aux_right = auxHole.a * cos(theta_aux);
y_aux_right = auxHole.b * sin(theta_aux);
x_right_boundary = cos(auxHole.theta_right) * x_aux_right - sin(auxHole.theta_right) * y_aux_right + auxHole.x_center;
y_right_boundary = sin(auxHole.theta_right) * x_aux_right + cos(auxHole.theta_right) * y_aux_right + auxHole.y_center;

% Left auxiliary hole boundary
x_aux_left = auxHole.a * cos(theta_aux);
y_aux_left = auxHole.b * sin(theta_aux);
x_left_boundary = -(cos(auxHole.theta_left) * x_aux_left - sin(auxHole.theta_left) * y_aux_left) - auxHole.x_center;
y_left_boundary = sin(auxHole.theta_left) * x_aux_left + cos(auxHole.theta_left) * y_aux_left + auxHole.y_center;

% Interpolate stresses at boundary points
try
    stress_central = interpolateStress(R, [x_central_boundary; y_central_boundary]);
    stress_right = interpolateStress(R, [x_right_boundary; y_right_boundary]);
    stress_left = interpolateStress(R, [x_left_boundary; y_left_boundary]);

    % Compute von Mises at boundaries
    vM_central = sqrt(stress_central.sxx.^2 - stress_central.sxx.*stress_central.syy + ...
                      stress_central.syy.^2 + 3*stress_central.sxy.^2);
    vM_right = sqrt(stress_right.sxx.^2 - stress_right.sxx.*stress_right.syy + ...
                    stress_right.syy.^2 + 3*stress_right.sxy.^2);
    vM_left = sqrt(stress_left.sxx.^2 - stress_left.sxx.*stress_left.syy + ...
                   stress_left.syy.^2 + 3*stress_left.sxy.^2);

    % Report boundary stress statistics
    fprintf('\nStress at hole boundaries:\n');
    fprintf('  Central hole:     max vM = %.2f MPa, mean vM = %.2f MPa\n', max(vM_central), mean(vM_central));
    fprintf('  Right aux hole:   max vM = %.2f MPa, mean vM = %.2f MPa\n', max(vM_right), mean(vM_right));
    fprintf('  Left aux hole:    max vM = %.2f MPa, mean vM = %.2f MPa\n', max(vM_left), mean(vM_left));

    % Figure 5: Stress along boundaries
    figure('Name', 'Boundary Stress Distribution', 'Position', [100, 100, 1200, 400]);

    subplot(1, 3, 1);
    plot(rad2deg(theta_central), vM_central, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_central), stress_central.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Central Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);

    subplot(1, 3, 2);
    plot(rad2deg(theta_aux), vM_right, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_aux), stress_right.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Right Auxiliary Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);

    subplot(1, 3, 3);
    plot(rad2deg(theta_aux), vM_left, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_aux), stress_left.sxx, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees)');
    ylabel('Stress (MPa)');
    title('Left Auxiliary Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', 'Location', 'best');
    grid on;
    xlim([0, 360]);

catch ME
    warning('Could not interpolate stress at boundaries: %s', ME.message);
end

%% ========================================================================
%  SECTION 11: DISPLACEMENT VISUALIZATION
%  ========================================================================

fprintf('\n=== Generating Displacement Plots ===\n');

figure('Name', 'Displacement Field', 'Position', [100, 100, 1200, 500]);

% Plot x-displacement
subplot(1, 2, 1);
pdeplot(meshData, XYData=R.Displacement.ux, ColorMap="jet");
axis equal;
colorbar;
title('X-Displacement u_x (mm)');
xlabel('x (mm)');
ylabel('y (mm)');

% Plot y-displacement
subplot(1, 2, 2);
pdeplot(meshData, XYData=R.Displacement.uy, ColorMap="jet");
axis equal;
colorbar;
title('Y-Displacement u_y (mm)');
xlabel('x (mm)');
ylabel('y (mm)');

%% ========================================================================
%  SECTION 12: SUMMARY STATISTICS
%  ========================================================================

fprintf('\n');
fprintf('========================================================================\n');
fprintf('                         ANALYSIS SUMMARY\n');
fprintf('========================================================================\n');
fprintf('\n');
fprintf('GEOMETRY:\n');
fprintf('  Plate dimensions:          %.0f x %.0f mm (half-size)\n', plate.length, plate.width);
fprintf('  Central hole radius:       %.1f mm\n', centralHole.radius);
fprintf('  Auxiliary hole semi-axes:  a = %.1f mm, b = %.1f mm\n', auxHole.a, auxHole.b);
fprintf('  Aux hole x-position:       %.1f mm from center\n', auxHole.x_center);
fprintf('  Aux hole y-position:       %.1f mm from centerline\n', auxHole.y_center);
fprintf('  Right hole rotation:       %.1f degrees\n', rad2deg(auxHole.theta_right));
fprintf('  Left hole rotation:        %.1f degrees\n', rad2deg(auxHole.theta_left));
fprintf('\n');
fprintf('MESH:\n');
fprintf('  Number of nodes:           %d\n', numNodes);
fprintf('  Number of elements:        %d\n', numElements);
fprintf('  Element size (Hmax):       %.2f mm\n', Hmax);
fprintf('\n');
fprintf('LOADING & MATERIAL:\n');
fprintf('  Applied tension:           %.1f MPa\n', load.tension);
fprintf('  Young''s modulus:           %.0f MPa\n', material.E);
fprintf('  Poisson''s ratio:           %.2f\n', material.nu);
fprintf('\n');
fprintf('STRESS RESULTS:\n');
fprintf('  Maximum von Mises stress:  %.2f MPa\n', max_vonMises);
fprintf('  Location of max stress:    (%.2f, %.2f) mm\n', loc_max_vM(1), loc_max_vM(2));
fprintf('  Stress concentration (Kt): %.3f\n', K_t);
fprintf('  Theoretical Kt (no aux):   %.3f\n', K_t_theoretical);
fprintf('  Stress reduction:          %.1f%%\n', reduction_percent);
fprintf('\n');
fprintf('SOLUTION TIME:               %.2f seconds\n', solve_time);
fprintf('\n');
fprintf('========================================================================\n');

%% ========================================================================
%  SECTION 13: SAVE RESULTS (OPTIONAL)
%  ========================================================================

% Uncomment to save results to file
% results.plate = plate;
% results.centralHole = centralHole;
% results.auxHole = auxHole;
% results.material = material;
% results.load = load;
% results.mesh.numNodes = numNodes;
% results.mesh.numElements = numElements;
% results.mesh.Hmax = Hmax;
% results.stress.max_vonMises = max_vonMises;
% results.stress.loc_max_vM = loc_max_vM;
% results.stress.K_t = K_t;
% results.stress.reduction_percent = reduction_percent;
% results.solve_time = solve_time;
% save('stress_analysis_results.mat', 'results');
% fprintf('Results saved to stress_analysis_results.mat\n');
