# Shape Optimization of a Plate with Stress-Relieving Auxiliary Holes

## 1. Introduction and Problem Statement

### 1.1 Background: Stress Concentration in Plates with Holes

When a structural plate containing a geometric discontinuity (such as a hole) is subjected to loading, the stress field becomes non-uniform. The stresses near the discontinuity are significantly amplified compared to the nominal far-field stress. This phenomenon is known as **stress concentration**.

For a circular hole in an infinite plate under uniaxial tension, classical elasticity theory (Kirsch solution, 1898) predicts a **stress concentration factor (SCF) of 3.0**. This means the maximum tangential stress at the hole boundary is three times the applied far-field stress:

$$K_t = \frac{\sigma_{max}}{\sigma_{nominal}} = 3.0$$

For finite-width plates, this factor can be even higher depending on the ratio of hole diameter to plate width. These high-stress regions are critical sites for:
- **Fatigue crack initiation**
- **Fracture propagation**
- **Premature structural failure**

### 1.2 Stress-Relieving Auxiliary Holes Concept

A proven technique to mitigate stress concentration is the introduction of **auxiliary holes** (also called "defense holes" or "stress-relief holes") strategically placed near the primary discontinuity. These smaller holes:

1. **Intercept stress trajectories** - Redirect the flow of stress around the main hole
2. **Distribute peak stresses** - Spread concentrated stresses over a larger region
3. **Smooth stress gradients** - Reduce the severity of stress transitions
4. **Lower the effective SCF** - Achieve a net reduction in maximum stress

The effectiveness of auxiliary holes depends on their:
- **Shape** (circular vs. elliptical)
- **Size** (diameter or semi-axes lengths)
- **Position** (distance from main hole and angular location)
- **Orientation** (rotation angle for elliptical holes)

### 1.3 Optimization Objective

**Goal:** Determine the optimal configuration of two symmetrically placed auxiliary holes (one on each side of the central hole along the loading axis) that **minimizes the maximum von Mises stress** anywhere in the plate.

---

## 2. Geometry and Problem Setup

### 2.1 Reference Configuration

Based on the MATLAB example (`matlab_example_Stress_Concentration_Plate_with_Circular_Hole.m`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| Central hole radius | 20.0 mm | Fixed primary hole |
| Plate width | 50.0 mm | Half-width from centerline |
| Plate length | 200.0 mm | 4x width for far-field conditions |
| Applied load | 100 MPa | Uniaxial tension (x-direction) |
| Young's modulus | 200 GPa | Steel-like material |
| Poisson's ratio | 0.25 | Typical metallic material |

### 2.2 Auxiliary Hole Placement Strategy

The auxiliary holes will be placed along the **x-axis** (loading direction), symmetrically on the left and right sides of the central hole. This placement targets the regions of maximum stress concentration, which occur at the points on the hole boundary perpendicular to the loading direction (at 0 and 180 from the load axis).

```
                        Applied Load (Tension)
                              -->  -->  -->
    +----------------------------------------------------------+
    |                                                          |
    |        (o)            ( O )            (o)               |
    |     Auxiliary       Central Hole    Auxiliary            |
    |       Hole                            Hole               |
    |                                                          |
    +----------------------------------------------------------+
                              <--  <--  <--
                        (Reaction/Constraint)
```

---

## 3. Circular vs. Elliptical Auxiliary Holes: Analysis and Decision

### 3.1 Circular Auxiliary Holes

**Advantages:**
- **Simpler geometry** - Defined by only 3 parameters per hole (center x, center y, radius)
- **Easier manufacturing** - Standard drilling operations
- **Lower computational cost** - Fewer design variables in optimization
- **No orientation sensitivity** - Rotationally symmetric

**Disadvantages:**
- **Limited shape adaptability** - Cannot conform to stress flow patterns
- **Suboptimal stress redistribution** - Circular shape may not optimally intercept stress trajectories
- **Constrained optimization space** - Fewer degrees of freedom to minimize stress

**Design Variables per Hole (Circular):** 3
- `x_c`: x-coordinate of center
- `y_c`: y-coordinate of center (typically 0 for symmetric placement)
- `r`: radius

### 3.2 Elliptical Auxiliary Holes

**Advantages:**
- **Greater geometric flexibility** - Can adapt shape to stress flow
- **Superior stress redistribution** - Elongated shapes can better intercept stress trajectories
- **Larger optimization space** - More design freedom often leads to better optima
- **Physical intuition** - Stress trajectories are generally elliptical; matching this pattern can be beneficial

**Disadvantages:**
- **More complex geometry** - Defined by 5 parameters per hole
- **Manufacturing complexity** - Requires CNC machining or specialized tooling
- **Higher computational cost** - More design variables increase optimization time
- **Orientation sensitivity** - Rotation angle adds complexity

**Design Variables per Hole (Elliptical):** 5
- `x_c`: x-coordinate of center
- `y_c`: y-coordinate of center
- `a`: semi-major axis length
- `b`: semi-minor axis length
- `theta`: rotation angle from horizontal

### 3.3 Engineering Decision: Elliptical Auxiliary Holes

**Recommendation: Use Elliptical Auxiliary Holes**

**Rationale:**

1. **Generality:** Elliptical holes encompass circular holes as a special case (when a = b). If circular is truly optimal, the optimization will converge to a = b naturally. This allows the algorithm to "discover" whether circular is sufficient or if elliptical provides additional benefit.

2. **Physical Basis:** The stress field around a circular hole in a uniaxially loaded plate is not axisymmetric. The stress contours form elliptical patterns. Elliptical auxiliary holes can better conform to these natural stress trajectories.

3. **Literature Support:** Published research on defense holes consistently shows that elliptical auxiliary holes outperform circular ones in reducing stress concentration factors, typically achieving 15-25% additional reduction.

4. **Optimization Philosophy:** Starting with a richer design space and allowing the algorithm to find constraints (like a = b) is mathematically sounder than artificially restricting the space upfront.

5. **Practical Consideration:** While elliptical holes are more complex to manufacture, modern CNC machining makes this a minor concern for high-value components where fatigue life is critical.

**Counter-argument Consideration:**
Some might argue that during optimization, elliptical holes will simply converge to circular. However, this is not guaranteed because:
- The optimal ellipse orientation may align with or against the load direction
- Aspect ratios (a/b) different from 1.0 can provide superior stress redistribution
- The interaction between two elliptical holes creates complex stress fields where shape matters

---

## 4. Optimization Problem Formulation

### 4.1 Design Variables

For **two elliptical auxiliary holes** with positional symmetry but **independent rotation angles**:

| Variable | Symbol | Description | Units |
|----------|--------|-------------|-------|
| x-position | `x_aux` | Distance from origin to auxiliary hole center | mm |
| y-position | `y_aux` | Vertical offset from centerline | mm |
| Semi-major axis | `a` | Larger semi-axis of ellipse | mm |
| Semi-minor axis | `b` | Smaller semi-axis of ellipse | mm |
| Right hole rotation | `theta_right` | CCW rotation for right auxiliary hole | radians |
| Left hole rotation | `theta_left` | CCW rotation for left auxiliary hole | radians |

**Note on Symmetry:**
- **Positional symmetry** is maintained: the left hole center is at `(-x_aux, y_aux)`
- **Shape symmetry** is maintained: both holes share the same semi-axes `a` and `b`
- **Rotational independence**: each hole has its own rotation angle, allowing asymmetric orientations
- This reduces the design space from 10 variables (5 per hole) to **6 variables**

**Design Vector:**
```
X = [x_aux, y_aux, a, b, theta_right, theta_left]
```

**Rationale for Independent Rotations:**
The stress field around the central hole is symmetric about the x-axis but the optimal ellipse orientation to intercept stress trajectories may differ on each side. Allowing independent rotations enables the optimizer to find configurations where the ellipses "lean" toward or away from the stress flow differently on each side.

### 4.2 Objective Function

**Minimize:** Maximum von Mises stress anywhere in the plate

$$f(X) = \max_{(x,y) \in \Omega} \sigma_{vM}(x, y)$$

Where:
- $\Omega$ is the plate domain (excluding holes)
- $\sigma_{vM}$ is the von Mises equivalent stress

**Implementation:**
```matlab
% Objective function pseudocode
function f = objective(X)
    % 1. Create geometry with central hole and auxiliary holes
    % 2. Set up FE model with boundary conditions
    % 3. Generate mesh
    % 4. Solve for stress field
    % 5. Compute von Mises stress at all nodes/integration points
    % 6. Return maximum value
    f = max(vonMisesStress);
end
```

### 4.3 Constraints

#### Geometric Constraints (Inequality)

1. **Non-overlapping with central hole:**
   - The auxiliary hole boundary must not intersect the central hole
   - Include a minimum gap for mesh quality: `gap_min >= 3 mm`

2. **Contained within plate:**
   - Auxiliary holes must lie entirely within the plate boundaries
   - Include edge clearance: `edge_clearance >= 5 mm`

3. **Non-overlapping auxiliary holes:**
   - Left and right auxiliary holes must not intersect each other
   - Relevant if `y_aux` allows vertical overlap

4. **Minimum size:**
   - `a >= a_min` and `b >= b_min` (e.g., 2 mm) for manufacturability

5. **Maximum size:**
   - `a <= a_max` and `b <= b_max` to prevent excessive material removal

6. **Aspect ratio limit:**
   - `a/b <= AR_max` (e.g., 3.0) to prevent extremely thin ellipses

#### Variable Bounds (Box Constraints)

| Variable | Lower Bound | Upper Bound | Rationale |
|----------|-------------|-------------|-----------|
| `x_aux` | `radius + a_max + gap` | `totalLength - edge_clearance - a_max` | Keep away from central hole and plate edges |
| `y_aux` | `0` | `totalWidth - edge_clearance - b_max` | Symmetric; don't exceed plate boundary |
| `a` | `2.0 mm` | `15.0 mm` | Manufacturability and material limits |
| `b` | `2.0 mm` | `15.0 mm` | Manufacturability and material limits |
| `theta_right` | `-pi/2` | `pi/2` | Full rotation coverage |
| `theta_left` | `-pi/2` | `pi/2` | Full rotation coverage |

### 4.4 Constraint Functions for Rotated Ellipses

**Key Challenge:** Computing clearances for rotated ellipses requires finding the axis-aligned bounding box (AABB) of each ellipse. For a rotated ellipse with semi-axes `a`, `b` and rotation angle `theta`, the AABB half-extents are:

$$x_{extent} = \sqrt{a^2 \cos^2(\theta) + b^2 \sin^2(\theta)}$$
$$y_{extent} = \sqrt{a^2 \sin^2(\theta) + b^2 \cos^2(\theta)}$$

These formulas give the maximum x and y distances from the ellipse center to its boundary, accounting for rotation.

**Nonlinear Inequality Constraints:** `c(X) <= 0`

```matlab
function [c, ceq] = constraints(X)
    x_aux = X(1); y_aux = X(2); a = X(3); b = X(4);
    theta_right = X(5); theta_left = X(6);

    % Compute axis-aligned bounding box extents for each hole
    x_ext_right = sqrt(a^2*cos(theta_right)^2 + b^2*sin(theta_right)^2);
    y_ext_right = sqrt(a^2*sin(theta_right)^2 + b^2*cos(theta_right)^2);

    x_ext_left = sqrt(a^2*cos(theta_left)^2 + b^2*sin(theta_left)^2);
    y_ext_left = sqrt(a^2*sin(theta_left)^2 + b^2*cos(theta_left)^2);

    % Use max extent for conservative central hole clearance
    max_ext_right = max(x_ext_right, y_ext_right);
    max_ext_left = max(x_ext_left, y_ext_left);

    % Constraint 1 & 2: Central hole clearance
    dist_to_center = sqrt(x_aux^2 + y_aux^2);
    c(1) = (radius + max_ext_right + gap_min) - dist_to_center;  % Right hole
    c(2) = (radius + max_ext_left + gap_min) - dist_to_center;   % Left hole

    % Constraint 3: Right boundary clearance (right hole)
    c(3) = (x_aux + x_ext_right) - (totalLength - edge_clearance);

    % Constraint 4: Left boundary clearance (left hole)
    c(4) = (x_aux + x_ext_left) - (totalLength - edge_clearance);

    % Constraint 5 & 6: Top boundary clearance
    c(5) = (abs(y_aux) + y_ext_right) - (totalWidth - edge_clearance);
    c(6) = (abs(y_aux) + y_ext_left) - (totalWidth - edge_clearance);

    % Constraint 7: Aspect ratio (handled via bounds or explicit constraint)
    c(7) = a - constraints.max_aspect_ratio * b;

    % No equality constraints
    ceq = [];
end
```

**Important Implementation Note:** The bounding box approach is conservative—it may reject some feasible designs. For higher fidelity, one could sample points on the ellipse boundary and check exact distances, but the AABB method provides a good balance of accuracy and computational efficiency.

---

## 5. Optimization Algorithm: fmincon

### 5.1 Why fmincon?

MATLAB's `fmincon` is chosen for this problem because:

1. **Handles nonlinear constraints** - Essential for geometric feasibility
2. **Supports bounds** - Natural for physical design limits
3. **Multiple algorithms available** - Interior-point, SQP, active-set
4. **Gradient options** - Can use finite differences or user-supplied gradients
5. **Well-documented** - Extensive MATLAB support and examples
6. **Parallelization support** - Can use `UseParallel` option for gradient estimation

### 5.2 Recommended Algorithm Options

```matlab
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...              % Good for engineering problems
    'Display', 'iter', ...               % Show progress
    'MaxIterations', 200, ...            % Sufficient for 6 variables
    'MaxFunctionEvaluations', 1200, ...  % Allow thorough search
    'OptimalityTolerance', 1e-4, ...     % Reasonable for stress values
    'StepTolerance', 1e-4, ...           % Relative to variable scales
    'ConstraintTolerance', 1e-4, ...     % Geometric precision
    'FiniteDifferenceType', 'central', ... % Better gradient accuracy
    'UseParallel', true);                % Parallel gradient evaluation
```

### 5.3 Multi-Start Strategy

Due to potential local minima in this nonconvex problem, a **multi-start approach** is recommended:

1. Generate `N` random initial points within feasible bounds (e.g., N = 20-50)
2. Run `fmincon` from each starting point
3. Keep track of the best solution found
4. Optionally use Latin Hypercube Sampling (LHS) for better space coverage

```matlab
% Multi-start pseudocode
best_f = Inf;
best_X = [];

for i = 1:num_starts
    X0 = generate_random_feasible_point();
    [X_opt, f_opt] = fmincon(@objective, X0, [], [], [], [], lb, ub, @constraints, options);

    if f_opt < best_f
        best_f = f_opt;
        best_X = X_opt;
    end
end
```

### 5.4 Alternative: Using GlobalSearch or MultiStart

MATLAB's Global Optimization Toolbox provides:
- `GlobalSearch` - Scatter search with local solver
- `MultiStart` - Parallel multi-start with fmincon

```matlab
problem = createOptimProblem('fmincon', ...
    'objective', @objective, ...
    'x0', X0, ...
    'lb', lb, 'ub', ub, ...
    'nonlcon', @constraints, ...
    'options', options);

ms = MultiStart('UseParallel', true, 'Display', 'iter');
[X_opt, f_opt] = run(ms, problem, 50);  % 50 start points
```

---

## 6. Implementation Workflow

### 6.1 Process Overview

```
+-------------------+     +------------------+     +-------------------+
|  1. Initialize    | --> |  2. Generate     | --> |  3. Create        |
|  Design Variables |     |  Feasible X0     |     |  Geometry (decsg) |
+-------------------+     +------------------+     +-------------------+
                                                            |
                                                            v
+-------------------+     +------------------+     +-------------------+
|  6. Compute Max   | <-- |  5. Solve FE     | <-- |  4. Setup FE      |
|  von Mises Stress |     |  Problem (solve) |     |  Model (femodel)  |
+-------------------+     +------------------+     +-------------------+
        |
        v
+-------------------+     +------------------+     +-------------------+
|  7. Return to     | --> |  8. Converged?   | --> |  9. Output        |
|  Optimizer        |     |  Check criteria  |     |  Optimal Design   |
+-------------------+     +------------------+     +-------------------+
```

### 6.2 Key Implementation Steps

1. **Parameterized Geometry Function:**
   - Create a function that takes design variables and returns a geometry object
   - Handle ellipse creation using parametric equations or rotation transformations
   - Use CSG operations: `Plate - CentralHole - LeftAuxHole - RightAuxHole`

2. **Robust Meshing:**
   - Use adaptive mesh refinement near hole boundaries
   - Set `Hmax` based on smallest geometric feature
   - Include error handling for mesh failures (infeasible geometries)

3. **Stress Evaluation:**
   - Compute von Mises stress at all mesh nodes
   - Consider stress evaluation on hole boundaries specifically
   - Handle stress singularities appropriately

4. **Parallel Evaluation:**
   - Use `parfor` for multi-start initial evaluations
   - Enable parallel gradient computation in fmincon

### 6.3 Handling Infeasible Geometries

During optimization, some parameter combinations may produce invalid geometries (overlapping holes, holes outside plate, etc.). Handle this by:

```matlab
function f = objective(X)
    try
        % Create geometry and solve
        ...
        f = max(vonMisesStress);
    catch
        % Return penalty for infeasible geometry
        f = 1e10;  % Large penalty value
    end
end
```

---

## 7. Expected Outcomes and Validation

### 7.1 Expected Results

Based on literature and physical intuition:

1. **Stress Reduction:** Expect 15-30% reduction in maximum stress compared to baseline (central hole only)

2. **Optimal Position:** Auxiliary holes likely to be positioned:
   - Close to but not touching the central hole
   - Slightly above/below the x-axis (non-zero y_aux) to "catch" stress flow

3. **Optimal Shape:**
   - If elliptical is beneficial: elongated perpendicular to loading direction
   - If circular is sufficient: a approximately equal to b

4. **Optimal Size:** Auxiliary holes typically 30-60% of central hole radius

### 7.2 Validation Strategy

1. **Mesh Convergence Study:**
   - Run analysis with progressively finer meshes
   - Verify stress values converge

2. **Comparison with Literature:**
   - Compare SCF reduction percentages with published results
   - Validate stress patterns qualitatively

3. **Baseline Comparison:**
   - Compute stress for plate with central hole only (SCF approximately 3.0)
   - Quantify improvement from optimization

4. **Sensitivity Analysis:**
   - Perturb optimal design variables
   - Verify optimum is a genuine minimum

---

## 8. Summary and Recommendations

### 8.1 Problem Summary

| Aspect | Specification |
|--------|---------------|
| **Objective** | Minimize maximum von Mises stress |
| **Design Variables** | 6 (position, size, independent orientations of elliptical auxiliary holes) |
| **Constraints** | Geometric feasibility (no overlap, within plate, clearance from edges) |
| **Solver** | fmincon with SQP algorithm |
| **Strategy** | Multi-start for global exploration |

### 8.2 Key Recommendations

1. **Use elliptical auxiliary holes** - Greater design freedom leads to better solutions; circular is a special case that can emerge naturally.

2. **Exploit partial symmetry** - Reduce from 10 to 6 design variables by assuming positional and shape symmetry, while allowing independent rotation angles.

3. **Implement robust geometry handling** - Use try-catch and penalty functions for infeasible configurations.

4. **Use multi-start optimization** - Nonconvex problem requires global search strategy.

5. **Validate thoroughly** - Perform mesh convergence studies and compare with known results.

### 8.3 Future Extensions

- **Multiple auxiliary hole pairs** - Add additional hole pairs for further stress reduction
- **Non-symmetric configurations** - Allow asymmetric placement for non-uniform loading
- **Alternative hole shapes** - Consider slots, rounded rectangles, or other shapes
- **Multi-objective optimization** - Balance stress reduction against material removal
- **Topology optimization** - Let the algorithm determine optimal hole count and placement

---

## References

1. MathWorks Documentation: "Stress Concentration in Plate with Circular Hole" - Partial Differential Equation Toolbox Example

2. Peterson, R.E. (1974). *Stress Concentration Factors*. John Wiley & Sons.

3. Meguid, S.A. (1986). "Finite element analysis of defence hole systems for the reduction of stress concentration in a uniaxially-loaded plate with two coaxial holes." *Engineering Fracture Mechanics*, 25(4), 403-413.

4. Providakis, C.P., Sotiropoulos, D.A. (1997). "A Finite Element Approach on the Stress Concentration Factor of Plates with Multiple Circular Holes." *Archive of Applied Mechanics*, 67, 594-602.

5. Rajaiah, K., Naik, N.K. (1984). "Hole-shape optimization in a finite plate in the presence of auxiliary holes." *Experimental Mechanics*, 24(2), 157-161.
