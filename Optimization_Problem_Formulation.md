# Optimization Problem Formulation
## Stress-Relieving Auxiliary Holes for Plate with Central Hole

**Document Version:** 1.0
**Date:** December 2024

---

## 1. Problem Statement

**Objective:** Minimize the maximum von Mises stress in a rectangular plate containing a central circular hole by optimally placing two elliptical auxiliary (stress-relief) holes.

**Physical Motivation:** A circular hole in a uniaxially loaded plate creates stress concentration with a theoretical stress concentration factor (SCF) of approximately 3.0. Strategically placed auxiliary holes can intercept and redistribute stress trajectories, reducing the peak stress and improving fatigue life.

---

## 2. Physical Problem Setup

### 2.1 Geometry

| Parameter | Symbol | Value | Units | Description |
|-----------|--------|-------|-------|-------------|
| Plate half-length | $L$ | 200 | mm | Distance from center to left/right edge |
| Plate half-width | $W$ | 50 | mm | Distance from center to top/bottom edge |
| Central hole radius | $R_c$ | 20 | mm | Fixed circular hole at origin |
| Total plate dimensions | - | 400 × 100 | mm | Full plate size |

**Coordinate System:**
- Origin at plate center (center of central hole)
- x-axis: horizontal, aligned with loading direction
- y-axis: vertical, perpendicular to loading

```
                    y
                    ^
                    |
    +---------------+---------------+
    |               |               |
    |      (L)      o      (R)      |  ---> x
    |         Auxiliary Holes       |
    |               |               |
    +---------------+---------------+
                    |
              Central Hole at (0,0)
```

### 2.2 Material Properties

| Property | Symbol | Value | Units |
|----------|--------|-------|-------|
| Young's modulus | $E$ | 200,000 | MPa |
| Poisson's ratio | $\nu$ | 0.25 | - |
| Analysis type | - | 2D Plane Stress | - |

### 2.3 Loading Conditions

| Load | Location | Value | Direction |
|------|----------|-------|-----------|
| Surface traction | Right edge ($x = +L$) | $\sigma_0 = 100$ MPa | +x (tensile) |

### 2.4 Boundary Conditions

| Boundary | Constraint | Description |
|----------|------------|-------------|
| Left edge ($x = -L$) | $u_x = 0$ | Fixed in x-direction |
| Bottom-left vertex | $u_y = 0$ | Prevents rigid body rotation |

---

## 3. Design Configuration

### 3.1 Auxiliary Hole Description

Two elliptical auxiliary holes are placed symmetrically about the y-axis:
- **Right auxiliary hole:** Center at $(x_{aux}, y_{aux})$
- **Left auxiliary hole:** Center at $(-x_{aux}, y_{aux})$ (mirrored)

Each ellipse is defined by:
- Semi-major axis: $a$
- Semi-minor axis: $b$
- Rotation angle: $\theta$ (independent for each hole)

### 3.2 Symmetry Exploitation

| Property | Symmetry Type | Benefit |
|----------|---------------|---------|
| Position | Mirrored about y-axis | Reduces 4 position variables to 2 |
| Shape (a, b) | Shared between holes | Reduces 4 shape variables to 2 |
| Rotation | **Independent** | Allows asymmetric orientations |

**Rationale for Independent Rotations:** While the stress field is symmetric about the x-axis, allowing independent rotation angles enables the optimizer to explore configurations where ellipses "lean" differently to optimally intercept stress trajectories.

---

## 4. Mathematical Formulation

### 4.1 Canonical Form

$$
\begin{aligned}
\min_{\mathbf{x}} \quad & f(\mathbf{x}) \\
\text{subject to} \quad & g_i(\mathbf{x}) \leq 0, \quad i = 1, \ldots, m \\
& h_j(\mathbf{x}) = 0, \quad j = 1, \ldots, p \\
& \mathbf{x}^L \leq \mathbf{x} \leq \mathbf{x}^U
\end{aligned}
$$

### 4.2 Design Variables

**Design Vector:**
$$\mathbf{x} = [x_{aux}, \, y_{aux}, \, a, \, b, \, \theta_R, \, \theta_L]^T \in \mathbb{R}^6$$

| Index | Variable | Symbol | Description | Units |
|-------|----------|--------|-------------|-------|
| 1 | $x_1$ | $x_{aux}$ | x-coordinate of auxiliary hole center | mm |
| 2 | $x_2$ | $y_{aux}$ | y-coordinate of auxiliary hole center | mm |
| 3 | $x_3$ | $a$ | Semi-major axis of ellipse | mm |
| 4 | $x_4$ | $b$ | Semi-minor axis of ellipse | mm |
| 5 | $x_5$ | $\theta_R$ | Rotation angle of right auxiliary hole | rad |
| 6 | $x_6$ | $\theta_L$ | Rotation angle of left auxiliary hole | rad |

### 4.3 Objective Function

$$f(\mathbf{x}) = \max_{(x,y) \in \Omega(\mathbf{x})} \sigma_{vM}(x, y; \mathbf{x})$$

Where:
- $\Omega(\mathbf{x})$ is the plate domain excluding all holes (depends on design variables)
- $\sigma_{vM}$ is the von Mises equivalent stress

**Von Mises Stress (2D Plane Stress):**
$$\sigma_{vM} = \sqrt{\sigma_{xx}^2 - \sigma_{xx}\sigma_{yy} + \sigma_{yy}^2 + 3\tau_{xy}^2}$$

**Computation:** The stress field is obtained by solving the 2D plane-stress elasticity equations using the Finite Element Method (MATLAB PDE Toolbox).

### 4.4 Variable Bounds

$$\mathbf{x}^L \leq \mathbf{x} \leq \mathbf{x}^U$$

| Variable | Lower Bound ($x^L$) | Upper Bound ($x^U$) | Rationale |
|----------|---------------------|---------------------|-----------|
| $x_{aux}$ | $R_c + a_{max} + \delta_c = 38$ mm | $0.75 \cdot L = 150$ mm | Clear of central hole; not too close to edge |
| $y_{aux}$ | $0$ mm | $W - \delta_e - b_{max} = 30$ mm | On or above centerline; within plate |
| $a$ | $6$ mm | $15$ mm | Manufacturability limits |
| $b$ | $6$ mm | $15$ mm | Manufacturability limits |
| $\theta_R$ | $-\pi/2$ rad | $+\pi/2$ rad | Full rotation range |
| $\theta_L$ | $-\pi/2$ rad | $+\pi/2$ rad | Full rotation range |

**Derived Constants:**
- $\delta_c = 3$ mm (minimum gap from central hole)
- $\delta_e = 5$ mm (minimum gap from plate edges)
- $a_{max} = 15$ mm, $b_{max} = 15$ mm

### 4.5 Inequality Constraints

**Constraint Count:** $m = 7$ inequality constraints

All constraints are formulated as $g_i(\mathbf{x}) \leq 0$:

---

#### Constraint 1-2: Central Hole Clearance

The auxiliary holes must maintain minimum clearance from the central hole.

**Auxiliary Function — Axis-Aligned Bounding Box (AABB) Extents:**

For a rotated ellipse with semi-axes $a$, $b$ and rotation $\theta$:
$$e_x(\theta) = \sqrt{a^2 \cos^2\theta + b^2 \sin^2\theta}$$
$$e_y(\theta) = \sqrt{a^2 \sin^2\theta + b^2 \cos^2\theta}$$

**Maximum extent (conservative):**
$$e_{max}(\theta) = \max(e_x(\theta), e_y(\theta))$$

**Constraints:**
$$g_1(\mathbf{x}) = R_c + e_{max}(\theta_R) + \delta_c - \sqrt{x_{aux}^2 + y_{aux}^2} \leq 0$$
$$g_2(\mathbf{x}) = R_c + e_{max}(\theta_L) + \delta_c - \sqrt{x_{aux}^2 + y_{aux}^2} \leq 0$$

*Physical meaning:* Distance from ellipse boundary to central hole boundary must exceed $\delta_c$.

---

#### Constraint 3-4: Plate Right/Left Edge Clearance

$$g_3(\mathbf{x}) = x_{aux} + e_x(\theta_R) - (L - \delta_e) \leq 0 \quad \text{(right hole to right edge)}$$
$$g_4(\mathbf{x}) = x_{aux} + e_x(\theta_L) - (L - \delta_e) \leq 0 \quad \text{(left hole to left edge)}$$

*Physical meaning:* Auxiliary holes must not extend beyond plate boundaries with clearance $\delta_e$.

---

#### Constraint 5-6: Plate Top/Bottom Edge Clearance

$$g_5(\mathbf{x}) = |y_{aux}| + e_y(\theta_R) - (W - \delta_e) \leq 0$$
$$g_6(\mathbf{x}) = |y_{aux}| + e_y(\theta_L) - (W - \delta_e) \leq 0$$

*Physical meaning:* Auxiliary holes must not extend beyond top/bottom plate boundaries.

---

#### Constraint 7: Aspect Ratio Limit

$$g_7(\mathbf{x}) = a - \alpha_{max} \cdot b \leq 0$$

Where $\alpha_{max} = 3.0$ is the maximum allowed aspect ratio.

*Physical meaning:* Prevents extremely elongated ellipses that are difficult to manufacture and may cause mesh issues.

---

### 4.6 Equality Constraints

$$p = 0 \quad \text{(no equality constraints)}$$

---

## 5. Constraint Summary Table

| Index | Constraint | Formula | Type |
|-------|------------|---------|------|
| $g_1$ | Right hole - central hole clearance | $R_c + e_{max}(\theta_R) + \delta_c - d_c \leq 0$ | Nonlinear |
| $g_2$ | Left hole - central hole clearance | $R_c + e_{max}(\theta_L) + \delta_c - d_c \leq 0$ | Nonlinear |
| $g_3$ | Right hole - right edge clearance | $x_{aux} + e_x(\theta_R) - (L - \delta_e) \leq 0$ | Nonlinear |
| $g_4$ | Left hole - left edge clearance | $x_{aux} + e_x(\theta_L) - (L - \delta_e) \leq 0$ | Nonlinear |
| $g_5$ | Right hole - top/bottom clearance | $\|y_{aux}\| + e_y(\theta_R) - (W - \delta_e) \leq 0$ | Nonlinear |
| $g_6$ | Left hole - top/bottom clearance | $\|y_{aux}\| + e_y(\theta_L) - (W - \delta_e) \leq 0$ | Nonlinear |
| $g_7$ | Aspect ratio | $a - \alpha_{max} \cdot b \leq 0$ | Linear |

Where $d_c = \sqrt{x_{aux}^2 + y_{aux}^2}$ is the distance from origin to auxiliary hole center.

---

## 6. Parameter Values Summary

### 6.1 Fixed Geometric Parameters

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| $L$ | 200 | mm | Plate half-length |
| $W$ | 50 | mm | Plate half-width |
| $R_c$ | 20 | mm | Central hole radius |

### 6.2 Constraint Parameters

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| $\delta_c$ | 3 | mm | Minimum gap from central hole |
| $\delta_e$ | 5 | mm | Minimum gap from plate edges |
| $a_{min}$ | 6 | mm | Minimum semi-axis |
| $a_{max}$ | 15 | mm | Maximum semi-axis |
| $b_{min}$ | 6 | mm | Minimum semi-axis |
| $b_{max}$ | 15 | mm | Maximum semi-axis |
| $\alpha_{max}$ | 3.0 | - | Maximum aspect ratio |
| $x_{aux,max}$ | 150 | mm | Maximum x-position (= 0.75L) |

### 6.3 Material and Loading Parameters

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| $E$ | 200,000 | MPa | Young's modulus |
| $\nu$ | 0.25 | - | Poisson's ratio |
| $\sigma_0$ | 100 | MPa | Applied tensile stress |

---

## 7. Optimization Algorithm Configuration

### 7.1 Solver Selection

**Solver:** MATLAB `fmincon`

**Rationale:**
- Handles nonlinear inequality constraints
- Supports box bounds natively
- Multiple algorithm options (SQP, interior-point, active-set)
- Parallel gradient computation capability

### 7.2 Algorithm Options

```matlab
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1200, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4, ...
    'FiniteDifferenceType', 'central', ...
    'UseParallel', true);
```

### 7.3 Multi-Start Strategy

Due to potential local minima, use multiple random starting points:

| Parameter | Value | Description |
|-----------|-------|-------------|
| Number of starts | 10-20 | For quick exploration |
| Sampling method | Latin Hypercube | Better coverage than random |
| Feasibility check | Required | Reject infeasible starting points |

### 7.4 Computational Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| Mesh Hmax factor | 2 | Hmax = min_feature / 2 (quick mode) |
| Penalty for infeasible geometry | $10^{10}$ | Large value for failed FEA |

---

## 8. Solution Output Specification

### 8.1 Primary Outputs

| Output | Symbol | Description |
|--------|--------|-------------|
| Optimal design vector | $\mathbf{x}^*$ | 6 optimal design variables |
| Optimal objective | $f^*$ | Minimum maximum von Mises stress |
| Stress concentration factor | $K_t^* = f^*/\sigma_0$ | Optimal SCF |
| Constraint values | $g_i(\mathbf{x}^*)$ | Should all be $\leq 0$ |

### 8.2 Comparison Outputs

| Configuration | Description |
|---------------|-------------|
| Baseline | Plate with central hole only (no auxiliary holes) |
| Initial guess | Starting point configuration |
| Optimal | Final optimized configuration |

**Metrics for Comparison:**
- Maximum von Mises stress
- Stress concentration factor ($K_t$)
- Stress reduction percentage: $\Delta = \frac{K_{t,baseline} - K_t^*}{K_{t,baseline}} \times 100\%$

### 8.3 Convergence Output

- Objective function vs. iteration plot
- Constraint violation vs. iteration (if applicable)

---

## 9. Implementation Notes

### 9.1 Ellipse Geometry Creation

MATLAB's `decsg` does not directly support rotated ellipses. Use polygon approximation:

```matlab
n_pts = 100;  % Number of polygon vertices
theta_pts = linspace(0, 2*pi, n_pts+1);
theta_pts = theta_pts(1:end-1);

% Local ellipse coordinates
x_local = a * cos(theta_pts);
y_local = b * sin(theta_pts);

% Apply rotation
x_rotated = cos(theta)*x_local - sin(theta)*y_local + x_center;
y_rotated = sin(theta)*x_local + cos(theta)*y_local + y_center;
```

### 9.2 Handling FEA Failures

Some design variable combinations may produce invalid geometries. Handle with try-catch and penalty:

```matlab
function f = objective(x)
    try
        % Create geometry, mesh, solve
        f = max(vonMisesStress);
    catch
        f = 1e10;  % Large penalty
    end
end
```

### 9.3 AABB Extent Functions

```matlab
function [ex, ey] = ellipse_extents(a, b, theta)
    ex = sqrt(a^2*cos(theta)^2 + b^2*sin(theta)^2);
    ey = sqrt(a^2*sin(theta)^2 + b^2*cos(theta)^2);
end
```

---

## 10. Validation Checklist

- [ ] Baseline stress analysis (central hole only) produces $K_t \approx 3.0$
- [ ] Constraint functions correctly identify feasible/infeasible designs
- [ ] Mesh convergence: stress values stable with mesh refinement
- [ ] Optimizer finds improvement over initial guess
- [ ] Optimal design satisfies all constraints ($g_i \leq 0$)
- [ ] Stress reduction is physically reasonable (15-30% expected)

---

## Appendix A: Quick Reference

### A.1 Design Variable Bounds (Numerical)

```matlab
lb = [38, 0, 6, 6, -pi/2, -pi/2];   % Lower bounds
ub = [150, 30, 15, 15, pi/2, pi/2]; % Upper bounds
```

### A.2 Constraint Function Signature

```matlab
function [c, ceq] = nonlcon(x)
    % x = [x_aux, y_aux, a, b, theta_R, theta_L]
    % c: inequality constraints (7x1), c <= 0 for feasible
    % ceq: equality constraints (empty)
end
```

### A.3 Objective Function Signature

```matlab
function f = objective(x)
    % x = [x_aux, y_aux, a, b, theta_R, theta_L]
    % f: maximum von Mises stress (scalar)
end
```
