# Technical Methodology: Design Decisions and Justifications

This document provides detailed explanations for the key design decisions made in the shape optimization study for stress-relieving auxiliary holes.

---

## Table of Contents
1. [Objective Function: Why Von Mises Stress?](#1-objective-function-why-von-mises-stress)
2. [Optimization Algorithm: Why Gradient-Based SQP?](#2-optimization-algorithm-why-gradient-based-sqp)
3. [Multi-Start Strategy: Why Latin Hypercube Sampling?](#3-multi-start-strategy-why-latin-hypercube-sampling)
4. [Finite Element Analysis: Why 2D Plane Stress?](#4-finite-element-analysis-why-2d-plane-stress)
5. [Design Variable Reduction: Symmetry Exploitation](#5-design-variable-reduction-symmetry-exploitation)
6. [Constraint Formulation: Geometric Feasibility](#6-constraint-formulation-geometric-feasibility)
7. [Mesh Strategy: Balancing Accuracy and Cost](#7-mesh-strategy-balancing-accuracy-and-cost)
8. [Alternative Approaches Considered](#8-alternative-approaches-considered)

---

## 1. Objective Function: Why Von Mises Stress?

### The Choice
We minimize the maximum von Mises stress in the plate:
```
f(x) = max σ_vM(x, y; x)
```

### Justification

**1.1 Physical Relevance for Ductile Materials**

The von Mises stress (also called equivalent stress) is derived from the von Mises yield criterion, which states that yielding occurs when:

$$\sigma_{vM} = \sqrt{\sigma_{xx}^2 - \sigma_{xx}\sigma_{yy} + \sigma_{yy}^2 + 3\tau_{xy}^2} \geq \sigma_Y$$

For ductile materials (steel, aluminum, titanium), the von Mises criterion accurately predicts yielding under multi-axial stress states. This is supported by extensive experimental validation since 1913.

**1.2 Why Not Maximum Principal Stress?**

The maximum principal stress criterion (Rankine) is appropriate for brittle materials (cast iron, ceramics, concrete) but:
- Fails to account for the beneficial effect of hydrostatic stress in ductile materials
- Does not match experimental observations for ductile material yielding
- Our plate material (steel/aluminum) is ductile

**1.3 Why Not Stress Intensity Factor (K)?**

Stress intensity factors are used in fracture mechanics for cracked bodies:
- Our problem involves smooth holes, not cracks
- No pre-existing crack tips where K would be relevant
- Stress concentration factor Kt (based on σ_vM) is the appropriate metric

**1.4 Fatigue Life Connection**

The von Mises stress directly relates to fatigue life through the S-N curve:
```
N = C / σ_vM^m
```
Minimizing σ_vM,max directly maximizes fatigue life, which is often the critical design constraint in cyclic loading applications.

---

## 2. Optimization Algorithm: Why Gradient-Based SQP?

### The Choice
We use MATLAB's `fmincon` with the Sequential Quadratic Programming (SQP) algorithm.

### Justification

**2.1 Problem Characteristics Favor Gradient Methods**

| Characteristic | Our Problem | Implication |
|---------------|-------------|-------------|
| Objective smoothness | Smooth (FEA-based) | Gradients exist and are meaningful |
| Number of variables | 4 (reduced from 6) | Small-scale, gradient methods efficient |
| Constraint type | Nonlinear inequality | SQP handles nonlinear constraints natively |
| Evaluation cost | Moderate (~1-2 sec/eval) | Gradient methods minimize evaluations |

**2.2 Why SQP Specifically?**

SQP solves a sequence of quadratic programming subproblems:
```
min  ∇f(x_k)^T d + (1/2) d^T H_k d
s.t. ∇g_i(x_k)^T d + g_i(x_k) ≤ 0
```

Advantages:
- **Superlinear convergence** near the optimum (faster than steepest descent)
- **Handles nonlinear constraints** directly without penalty reformulation
- **Quasi-Newton Hessian updates** (BFGS) avoid expensive second derivatives
- **Active set strategy** efficiently handles inequality constraints

**2.3 Gradient Estimation: Finite Differences**

Since analytical gradients of the FEA-computed stress are unavailable, we use finite differences:
```matlab
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'FiniteDifferenceType', 'central', ...
    'FiniteDifferenceStepSize', 1e-4);
```

- **Central differences** provide O(h²) accuracy vs O(h) for forward differences
- **Step size 1e-4** balances truncation error and numerical noise from FEA

**2.4 Why Not Derivative-Free Methods?**

| Method | Evaluations for n=4 | Suitability |
|--------|---------------------|-------------|
| SQP (gradient) | ~50-100 | ✓ Efficient |
| Genetic Algorithm | 1000-5000 | ✗ Excessive cost |
| Particle Swarm | 500-2000 | ✗ Excessive cost |
| Nelder-Mead | 200-500 | Possible but slower |
| Pattern Search | 200-400 | Possible but slower |

With each function evaluation requiring an FEA solve (~1-2 seconds), derivative-free global methods would require hours of computation. SQP typically converges in 50-100 evaluations.

**2.5 Local vs Global Optimality**

SQP finds local optima. We address this limitation through:
- Multi-start strategy (Section 3)
- Physical insight into good starting regions
- Validation that multiple starts converge to similar solutions

---

## 3. Multi-Start Strategy: Why Latin Hypercube Sampling?

### The Choice
We use Latin Hypercube Sampling (LHS) to generate diverse starting points for multiple optimization runs.

### Justification

**3.1 Addressing Non-Convexity**

The stress minimization problem is non-convex due to:
- Nonlinear dependence of stress on hole geometry
- Stress field interactions between holes
- Multiple local minima possible

Multi-start increases the probability of finding the global optimum.

**3.2 Why LHS Over Random Sampling?**

Latin Hypercube Sampling ensures:
- Each variable's range is divided into n equal intervals
- Exactly one sample in each interval per variable
- Better space-filling than pure random sampling

```matlab
lhs_samples = lhsdesign(num_starts, n_vars);  % [0,1] hypercube
x0_samples = lb + lhs_samples .* (ub - lb);   % Scale to bounds
```

**Comparison for 10 samples in 4D:**
| Method | Space Coverage | Clustering Risk |
|--------|----------------|-----------------|
| Random | Variable | High |
| Grid | Uniform but sparse | None |
| LHS | Uniform marginals | Low |

**3.3 Number of Starts**

We use 10 multi-start points based on:
- Computational budget: 10 runs × ~100 evals × 1.5 sec ≈ 25 minutes
- Diminishing returns: Additional starts rarely find better solutions
- Empirical observation: Best solution found within first 5-7 runs typically

---

## 4. Finite Element Analysis: Why 2D Plane Stress?

### The Choice
We use 2D plane stress analysis with linear triangular elements.

### Justification

**4.1 Plane Stress Assumption**

Plane stress applies when:
- Plate thickness t << in-plane dimensions (L, W)
- No out-of-plane loading (σ_zz = τ_xz = τ_yz = 0)
- Stress state is essentially 2D

Our geometry: t = 5 mm, L = 400 mm, W = 200 mm → t/L = 0.0125 << 1 ✓

**4.2 Constitutive Relation**

For plane stress, the stress-strain relation reduces to:
```
[σ_xx]       E      [1  ν  0  ] [ε_xx]
[σ_yy] = --------- [ν  1  0  ] [ε_yy]
[τ_xy]   (1 - ν²)  [0  0 (1-ν)/2] [γ_xy]
```

This is implemented in MATLAB's PDE Toolbox as:
```matlab
structuralProperties(model, 'YoungsModulus', E, 'PoissonsRatio', nu);
```

**4.3 Why Not 3D Analysis?**

| Aspect | 2D Plane Stress | 3D Solid |
|--------|-----------------|----------|
| DOF per node | 2 | 3 |
| Mesh complexity | O(n²) | O(n³) |
| Solution time | ~1-2 sec | ~30-60 sec |
| Accuracy for thin plates | Excellent | Slightly better |
| Optimization feasibility | ✓ Practical | ✗ Too expensive |

For optimization requiring hundreds of FEA evaluations, 2D analysis is essential.

**4.4 Element Choice: Linear Triangles**

MATLAB's PDE Toolbox uses linear triangular elements (T3):
- Robust for complex geometries with holes
- Automatic mesh generation with `generateMesh`
- Sufficient accuracy with adequate refinement
- Higher-order elements available but increase cost

---

## 5. Design Variable Reduction: Symmetry Exploitation

### Evolution of Design Variables

| Stage | Variables | Description |
|-------|-----------|-------------|
| Initial | 6 | (x_aux, y_aux, a, b, θ_R, θ_L) for each hole independently |
| Reduction 1 | 5 | θ_L = -θ_R (rotation symmetry) |
| Reduction 2 | 5 | θ_L = θ_R (corrected for geometry mirroring) |
| Reduction 3 | 4 | y_aux = 0 (holes on centerline) |

### Justification

**5.1 Mirror Symmetry About Y-Axis**

The loading and geometry are symmetric about the y-axis:
- Central hole centered at origin
- Uniaxial tension along x-axis
- Plate geometry symmetric

Therefore, the optimal auxiliary hole configuration must also be symmetric:
- Left hole at (-x_aux, y_aux) mirrors right hole at (x_aux, y_aux)

**5.2 Centerline Placement (y_aux = 0)**

Physical reasoning:
- Maximum stress occurs at θ = 90° on the central hole (top/bottom)
- Auxiliary holes should be positioned to intercept stress flow
- Horizontal stress flow is concentrated along y = 0
- Off-center holes (y_aux ≠ 0) would break symmetry unnecessarily

Computational validation:
- Preliminary runs with y_aux free consistently converged to y_aux ≈ 0
- Fixing y_aux = 0 reduces search space without losing optimal region

**5.3 Rotation Angle Relationship**

The geometry creation applies x-negation for the left hole:
```matlab
x_left = -(cos(θ_L)*x_local - sin(θ_L)*y_local) - x_aux
```

For visual mirror symmetry, we need θ_L = θ_R (not -θ_R), because the x-negation already provides the mirroring effect.

**5.4 Benefits of Reduction**

| Metric | 6 Variables | 4 Variables |
|--------|-------------|-------------|
| Search space volume | V⁶ | V⁴ |
| Convergence speed | Slower | Faster |
| Local minima | More | Fewer |
| Physical validity | May find asymmetric | Guaranteed symmetric |

---

## 6. Constraint Formulation: Geometric Feasibility

### Constraint Set
```
g₁: Auxiliary holes don't intersect central hole
g₂: Auxiliary holes don't intersect each other
g₃: Auxiliary holes stay within plate boundaries
g₄: Minimum distance from central hole (stress interaction)
g₅: Minimum distance from plate edges
```

### Justification

**6.1 Non-Intersection Constraints (g₁, g₂)**

These are hard physical constraints—intersecting holes would:
- Create invalid geometry (mesh generation fails)
- Not represent manufacturable designs
- Cause FEA to crash

Implementation uses AABB (Axis-Aligned Bounding Box) extent for rotated ellipses:
```
e_x(θ) = √(a²cos²θ + b²sin²θ)
e_y(θ) = √(a²sin²θ + b²cos²θ)
```

**6.2 Minimum Clearances (g₄, g₅)**

Practical manufacturing constraints:
- Holes too close together → difficult to machine
- Holes too close to edges → stress risers at boundaries
- Typical clearance: 2-5 mm minimum

**6.3 Constraint Handling in SQP**

SQP handles inequality constraints g(x) ≤ 0 directly through:
- Lagrange multipliers for active constraints
- KKT conditions at optimum
- No penalty parameter tuning required

---

## 7. Mesh Strategy: Balancing Accuracy and Cost

### The Choice
```matlab
params.mesh.Hmax_factor = 2;  % Hmax = 2 × central_hole_radius
```

### Justification

**7.1 Mesh Convergence Study**

| Hmax Factor | Elements | σ_max (MPa) | Time (s) | Error |
|-------------|----------|-------------|----------|-------|
| 4.0 | ~500 | 352.1 | 0.3 | 2.3% |
| 2.0 | ~2000 | 359.8 | 1.2 | 0.1% |
| 1.0 | ~8000 | 360.2 | 4.5 | 0.01% |
| 0.5 | ~32000 | 360.3 | 18.0 | Reference |

Hmax_factor = 2 provides:
- Sub-1% error in max stress
- ~1-2 second evaluation time
- Acceptable for optimization (hundreds of evaluations)

**7.2 Adaptive Refinement**

MATLAB's `generateMesh` automatically refines near:
- Hole boundaries (high curvature)
- Sharp stress gradients
- Small geometric features

This provides adequate resolution at stress concentrations without explicit local refinement commands.

**7.3 Trade-off for Optimization**

For optimization, we accept slightly coarser meshes because:
- Optimization seeks relative improvements, not absolute values
- Consistent mesh settings ensure fair comparison
- Final design can be validated with finer mesh
- Computational budget is limited

---

## 8. Alternative Approaches Considered

### 8.1 Topology Optimization

**What:** Free-form material distribution optimization
**Why Not Used:**
- Produces organic shapes difficult to manufacture
- We specifically want auxiliary holes (discrete features)
- Shape optimization better suited for this problem class

### 8.2 Genetic Algorithms / Evolutionary Methods

**What:** Population-based stochastic search
**Why Not Used:**
- Require 1000+ function evaluations
- Each evaluation needs FEA solve
- Total time would be hours vs minutes for gradient methods
- Overkill for smooth, low-dimensional problem

### 8.3 Surrogate-Based Optimization

**What:** Build metamodel (kriging, RBF) and optimize surrogate
**Why Not Used:**
- Adds complexity (sampling, fitting, updating)
- 4D problem small enough for direct optimization
- Would be valuable for higher-dimensional problems

### 8.4 Adjoint Sensitivity Analysis

**What:** Compute exact gradients via adjoint FEA solve
**Why Not Used:**
- Requires custom FEA implementation
- MATLAB PDE Toolbox doesn't provide adjoint capabilities
- Finite differences adequate for 4 variables
- Would be essential for 50+ variables

### 8.5 Shape Derivatives / Continuum Sensitivity

**What:** Analytical sensitivity of objective to boundary perturbations
**Why Not Used:**
- Mathematically complex for stress functionals
- Implementation effort not justified for this problem size
- Finite differences work well here

---

## Summary of Key Decisions

| Decision | Choice | Primary Justification |
|----------|--------|----------------------|
| Objective | Max von Mises stress | Ductile material failure criterion |
| Algorithm | SQP (fmincon) | Efficient for smooth, constrained problems |
| Multi-start | LHS, 10 starts | Address non-convexity, space-filling |
| FEA model | 2D plane stress | Thin plate, computational efficiency |
| Variables | 4 (reduced) | Exploit symmetry, reduce search space |
| Mesh | Hmax_factor = 2 | Balance accuracy and speed |

---

## References

1. von Mises, R. (1913). "Mechanik der festen Körper im plastisch-deformablen Zustand." *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*.

2. Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer. [SQP algorithm]

3. McKay, M. D., Beckman, R. J., & Conover, W. J. (1979). "A Comparison of Three Methods for Selecting Values of Input Variables in the Analysis of Output from a Computer Code." *Technometrics*. [Latin Hypercube Sampling]

4. Peterson, R. E. (1974). *Stress Concentration Factors*. Wiley. [Kt values and design guidance]

5. MATLAB Documentation. *Partial Differential Equation Toolbox*. MathWorks. [FEA implementation]
