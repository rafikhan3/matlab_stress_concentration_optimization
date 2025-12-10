%% Baseline Analysis: Plate with Central Hole Only
% This script performs FEA analysis on the baseline configuration
% (plate with central circular hole, no auxiliary holes) and generates
% visualization plots for mesh and stress distribution.
%
% Uses the same parameters as stress_optimization_auxiliary_holes.m
% to ensure consistency.
%
% Author: Generated for Shape Optimization Study
% Date: December 2024

clear; close all; clc;

%% ========================================================================
%  PROBLEM PARAMETERS (Same as optimization code)
%  ========================================================================

fprintf('=============================================================\n');
fprintf('  BASELINE ANALYSIS: Plate with Central Hole Only\n');
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

% --- Mesh Parameters ---
params.mesh.Hmax_factor = 2;                  % Same as optimization (Quick mode)

%% ========================================================================
%  DISPLAY PARAMETERS
%  ========================================================================

fprintf('Geometry Parameters:\n');
fprintf('  Plate dimensions: %.0f x %.0f mm (full size: %.0f x %.0f mm)\n', ...
    params.plate.length, params.plate.width, 2*params.plate.length, 2*params.plate.width);
fprintf('  Central hole radius: %.0f mm\n', params.centralHole.radius);
fprintf('\n');

fprintf('Material Properties:\n');
fprintf('  Young''s modulus E: %.0f MPa (%.0f GPa)\n', params.material.E, params.material.E/1000);
fprintf('  Poisson''s ratio nu: %.2f\n', params.material.nu);
fprintf('\n');

fprintf('Loading:\n');
fprintf('  Applied stress sigma_0: %.0f MPa (uniaxial tension)\n', params.load.tension);
fprintf('\n');

%% ========================================================================
%  CREATE GEOMETRY
%  ========================================================================

fprintf('Creating geometry...\n');

L = params.plate.length;
W = params.plate.width;
R_c = params.centralHole.radius;

% Rectangle (plate)
R1 = [3; 4; -L; L; L; -L; -W; -W; W; W];

% Circle (central hole)
C1 = [1; 0; 0; R_c; 0; 0; 0; 0; 0; 0];

% Combine geometries
gdm = [R1, C1];
ns = char('R1', 'C1')';
[g, ~] = decsg(gdm, 'R1 - C1', ns);

fprintf('  Geometry created successfully.\n');

%% ========================================================================
%  CREATE FE MODEL
%  ========================================================================

fprintf('Setting up FE model...\n');

model = femodel(AnalysisType="structuralStatic", Geometry=g);

% Material properties
model.MaterialProperties = materialProperties(...
    YoungsModulus=params.material.E, ...
    PoissonsRatio=params.material.nu);

% Boundary conditions
% Edge 2: Right edge (x = +L) - applied traction
% Edge 4: Left edge (x = -L) - fixed in x
model.EdgeLoad(2) = edgeLoad(SurfaceTraction=[params.load.tension; 0]);
model.EdgeBC(4) = edgeBC(XDisplacement=0);

% Fix one vertex to prevent rigid body motion in y
model.VertexBC(1) = vertexBC(YDisplacement=0);

fprintf('  Boundary conditions applied.\n');

%% ========================================================================
%  GENERATE MESH
%  ========================================================================

fprintf('Generating mesh...\n');

Hmax = params.centralHole.radius / params.mesh.Hmax_factor;
model = generateMesh(model, Hmax=Hmax);

mesh_info = model.Geometry;
num_nodes = size(model.Mesh.Nodes, 2);
num_elements = size(model.Mesh.Elements, 2);

fprintf('  Mesh generated:\n');
fprintf('    Hmax: %.2f mm\n', Hmax);
fprintf('    Nodes: %d\n', num_nodes);
fprintf('    Elements: %d\n', num_elements);

%% ========================================================================
%  SOLVE
%  ========================================================================

fprintf('Solving FEA problem...\n');
tic;
R = solve(model);
solve_time = toc;
fprintf('  Solution completed in %.2f seconds.\n', solve_time);

%% ========================================================================
%  COMPUTE STRESSES
%  ========================================================================

fprintf('Computing stresses...\n');

stress = R.Stress;
sxx = stress.sxx;
syy = stress.syy;
sxy = stress.sxy;

% Von Mises stress
vonMises = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*sxy.^2);

% Find maximum stress
[max_stress, max_idx] = max(vonMises);
nodes = R.Mesh.Nodes;
max_loc = nodes(:, max_idx);

% Stress concentration factor
Kt = max_stress / params.load.tension;

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    BASELINE RESULTS\n');
fprintf('=============================================================\n');
fprintf('  Maximum von Mises stress: %.2f MPa\n', max_stress);
fprintf('  Location of max stress: (%.2f, %.2f) mm\n', max_loc(1), max_loc(2));
fprintf('  Stress concentration factor Kt: %.3f\n', Kt);
fprintf('  Theoretical Kt (infinite plate): 3.000\n');
fprintf('  Difference from theory: %.1f%%\n', (Kt - 3.0)/3.0 * 100);
fprintf('=============================================================\n\n');

%% ========================================================================
%  FIGURE 1: GEOMETRY
%  ========================================================================

fprintf('Generating plots...\n');

figure('Name', 'Baseline Geometry', 'Position', [50, 50, 800, 400]);

pdegplot(g, 'EdgeLabels', 'on', 'FaceLabels', 'on');
axis equal;
axis([-1.1*L, 1.1*L, -1.1*W, 1.1*W]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Baseline Geometry with Edge Labels');
grid on;

% Add annotations
text(L + 10, 0, '\sigma_0 \rightarrow', 'FontSize', 12, 'Color', 'r');
text(-L - 30, 0, 'Fixed', 'FontSize', 10, 'Color', 'b');

%% ========================================================================
%  FIGURE 2: MESH
%  ========================================================================

figure('Name', 'Baseline Mesh', 'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
pdemesh(model);
axis equal;
axis([-1.1*L, 1.1*L, -1.1*W, 1.1*W]);
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Full Mesh (%d nodes, %d elements)', num_nodes, num_elements));

subplot(1, 2, 2);
pdemesh(model);
axis equal;
xlim([-50, 50]);
ylim([-40, 40]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Mesh Detail Near Central Hole');

%% ========================================================================
%  FIGURE 3: VON MISES STRESS
%  ========================================================================

figure('Name', 'Baseline Stress - Von Mises', 'Position', [150, 150, 1200, 500]);

subplot(1, 2, 1);
pdeplot(R.Mesh, 'XYData', vonMises, 'ColorMap', 'jet');
hold on;
plot(max_loc(1), max_loc(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Von Mises Stress (\\sigma_{max} = %.1f MPa, K_t = %.2f)', max_stress, Kt));

subplot(1, 2, 2);
pdeplot(R.Mesh, 'XYData', vonMises, 'ColorMap', 'jet');
hold on;
plot(max_loc(1), max_loc(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
xlim([-50, 50]);
ylim([-40, 40]);
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('Von Mises Stress - Zoomed Near Hole');

%% ========================================================================
%  FIGURE 4: STRESS COMPONENTS
%  ========================================================================

figure('Name', 'Baseline Stress Components', 'Position', [200, 200, 1400, 900]);

subplot(2, 2, 1);
pdeplot(R.Mesh, 'XYData', sxx, 'ColorMap', 'jet', 'Contour', 'on');
axis equal;
colorbar;
title('\sigma_{xx} (MPa) - Axial Stress');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(2, 2, 2);
pdeplot(R.Mesh, 'XYData', syy, 'ColorMap', 'jet', 'Contour', 'on');
axis equal;
colorbar;
title('\sigma_{yy} (MPa) - Transverse Stress');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(2, 2, 3);
pdeplot(R.Mesh, 'XYData', sxy, 'ColorMap', 'jet', 'Contour', 'on');
axis equal;
colorbar;
title('\tau_{xy} (MPa) - Shear Stress');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(2, 2, 4);
pdeplot(R.Mesh, 'XYData', vonMises, 'ColorMap', 'jet', 'Contour', 'on');
hold on;
plot(max_loc(1), max_loc(2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
axis equal;
colorbar;
title(sprintf('\\sigma_{vM} (MPa) - Von Mises (max = %.1f)', max_stress));
xlabel('x (mm)');
ylabel('y (mm)');

%% ========================================================================
%  FIGURE 5: DISPLACEMENT
%  ========================================================================

figure('Name', 'Baseline Displacement', 'Position', [250, 250, 1200, 500]);

ux = R.Displacement.ux;
uy = R.Displacement.uy;
u_mag = sqrt(ux.^2 + uy.^2);

subplot(1, 3, 1);
pdeplot(R.Mesh, 'XYData', ux, 'ColorMap', 'jet');
axis equal;
colorbar;
title('u_x (mm) - X-Displacement');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(1, 3, 2);
pdeplot(R.Mesh, 'XYData', uy, 'ColorMap', 'jet');
axis equal;
colorbar;
title('u_y (mm) - Y-Displacement');
xlabel('x (mm)');
ylabel('y (mm)');

subplot(1, 3, 3);
pdeplot(R.Mesh, 'XYData', u_mag, 'ColorMap', 'jet');
axis equal;
colorbar;
title('|u| (mm) - Displacement Magnitude');
xlabel('x (mm)');
ylabel('y (mm)');

%% ========================================================================
%  FIGURE 6: STRESS AROUND HOLE BOUNDARY
%  ========================================================================

figure('Name', 'Stress Distribution Around Hole', 'Position', [300, 300, 1000, 400]);

% Sample points around the central hole
theta_pts = linspace(0, 2*pi, 200);
x_boundary = R_c * cos(theta_pts);
y_boundary = R_c * sin(theta_pts);

% Interpolate stress at boundary points
try
    stress_boundary = interpolateStress(R, [x_boundary; y_boundary]);
    vM_boundary = sqrt(stress_boundary.sxx.^2 - stress_boundary.sxx.*stress_boundary.syy + ...
                       stress_boundary.syy.^2 + 3*stress_boundary.sxy.^2);

    subplot(1, 2, 1);
    plot(rad2deg(theta_pts), vM_boundary, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(rad2deg(theta_pts), stress_boundary.sxx, 'r--', 'LineWidth', 1);
    plot(rad2deg(theta_pts), stress_boundary.syy, 'g--', 'LineWidth', 1);
    hold off;
    xlabel('\theta (degrees from +x axis)');
    ylabel('Stress (MPa)');
    title('Stress Distribution Around Hole Boundary');
    legend('\sigma_{vM}', '\sigma_{xx}', '\sigma_{yy}', 'Location', 'best');
    grid on;
    xlim([0, 360]);

    % Highlight max stress locations
    subplot(1, 2, 2);
    polarplot(theta_pts, vM_boundary, 'b-', 'LineWidth', 1.5);
    title('Von Mises Stress (Polar Plot)');
    rlim([0, max(vM_boundary)*1.1]);
catch ME
    subplot(1, 2, 1);
    text(0.5, 0.5, 'Could not interpolate boundary stress', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
    subplot(1, 2, 2);
    text(0.5, 0.5, ME.message, ...
         'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 8);
end

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================

fprintf('Saving results...\n');

baseline_results = struct();
baseline_results.params = params;
baseline_results.max_stress = max_stress;
baseline_results.max_loc = max_loc;
baseline_results.Kt = Kt;
baseline_results.vonMises = vonMises;
baseline_results.stress = stress;
baseline_results.displacement = R.Displacement;
baseline_results.mesh = R.Mesh;
baseline_results.num_nodes = num_nodes;
baseline_results.num_elements = num_elements;
baseline_results.solve_time = solve_time;

save('baseline_results.mat', 'baseline_results');
fprintf('  Results saved to: baseline_results.mat\n');

%% ========================================================================
%  EXPORT FIGURES (Optional)
%  ========================================================================

fprintf('\nTo export figures, run:\n');
fprintf('  exportgraphics(figure(2), ''baseline_mesh.png'', ''Resolution'', 300)\n');
fprintf('  exportgraphics(figure(3), ''baseline_stress.png'', ''Resolution'', 300)\n');
fprintf('  exportgraphics(figure(4), ''baseline_stress_components.png'', ''Resolution'', 300)\n');

%% ========================================================================
%  FINAL SUMMARY
%  ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    ANALYSIS COMPLETE\n');
fprintf('=============================================================\n');
fprintf('Baseline Configuration:\n');
fprintf('  Plate: %.0f x %.0f mm\n', 2*L, 2*W);
fprintf('  Central hole radius: %.0f mm\n', R_c);
fprintf('  Applied stress: %.0f MPa\n', params.load.tension);
fprintf('\n');
fprintf('Results:\n');
fprintf('  Max von Mises stress: %.2f MPa\n', max_stress);
fprintf('  Stress concentration factor: Kt = %.3f\n', Kt);
fprintf('  Max stress location: (%.1f, %.1f) mm\n', max_loc(1), max_loc(2));
fprintf('\n');
fprintf('This baseline serves as the reference for optimization.\n');
fprintf('Target: Reduce max stress by adding auxiliary holes.\n');
fprintf('=============================================================\n');
