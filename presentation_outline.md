# Presentation Outline: Shape Optimization of Stress-Relieving Auxiliary Holes

**Duration:** ~15 minutes
**Audience:** Graduate students in structural optimization / computational mechanics

---

## Slide 1: Title Slide
**Title:** Shape Optimization of Stress-Relieving Auxiliary Holes in a Plate with Central Hole

**Subtitle:** A Gradient-Based Multi-Start Optimization Approach

**Content:**
- Your name / Institution
- Course name
- Date

---

## Slide 2: Motivation & Problem Context (1.5 min)

**Title:** Why Stress Concentration Matters

**Key Points:**
- Holes, notches, and geometric discontinuities cause stress concentration
- Stress concentration factor (Kt) for circular hole in infinite plate: Kt ≈ 3.0
- Fatigue failures initiate at stress concentration sites
- **Question:** Can we reduce stress concentration by adding MORE holes?

**[FIGURE: Simple schematic showing stress flow lines around a circular hole]**
*Annotation: Show stress trajectories bunching at the hole boundary, illustrating why σ_max occurs at the hole edge perpendicular to loading*

**Speaker Notes:**
- "Counter-intuitive idea: adding holes to reduce stress"
- Reference aerospace and automotive applications where weight/stress optimization is critical

---

## Slide 3: Physical Concept (1.5 min)

**Title:** Stress-Relief Through Auxiliary Holes

**Key Points:**
- Auxiliary holes intercept and redistribute stress flow
- Strategic placement "shields" the critical central hole region
- Trade-off: auxiliary holes create their own (smaller) stress concentrations

**[FIGURE: Geometry Comparison - 3 panels showing baseline, initial, optimal]**
*Annotation: Use the geometry comparison figure from the optimization output. Highlight the symmetric placement of auxiliary holes.*

**Speaker Notes:**
- Explain that this is a well-known technique in fracture mechanics
- Mention that historically this was done by trial-and-error; now we optimize

---

## Slide 4: Problem Setup - Geometry (1 min)

**Title:** Problem Geometry

**Content - Table:**
| Parameter | Symbol | Value |
|-----------|--------|-------|
| Plate dimensions | 2L × 2W | 400 × 100 mm |
| Central hole radius | R_c | 20 mm |
| Material | Steel | E = 200 GPa, ν = 0.25 |
| Loading | σ₀ | 100 MPa (uniaxial tension) |

**[FIGURE: Annotated schematic of plate geometry with dimensions]**
*Annotation: Create a clean schematic showing:
- Rectangular plate with coordinate system at center
- Central circular hole at origin
- Two auxiliary holes (left/right) with design variables labeled
- Arrows showing tensile loading on right edge
- Fixed BC symbols on left edge*

---

## Slide 5: Design Variables (1.5 min)

**Title:** Design Variables & Parameterization

**Content:**
Design vector: **x** = [x_aux, a, b, θ]ᵀ ∈ ℝ⁴

| Variable | Description | Bounds |
|----------|-------------|--------|
| x_aux | Distance from center | [38, 150] mm |
| a | Semi-major axis | [6, 15] mm |
| b | Semi-minor axis | [6, 15] mm |
| θ | Rotation angle | [-90°, 90°] |

**Fixed/Derived:**
- y_aux = 0 (centerline, based on preliminary results)
- θ_L = θ_R (mirror symmetry via geometry)

**[FIGURE: Close-up of auxiliary hole showing a, b, θ parameters]**
*Annotation: Show an ellipse with labeled semi-axes a, b and rotation angle θ from horizontal*

---

## Slide 6: Mathematical Formulation (2 min)

**Title:** Optimization Problem Formulation

**Content - Canonical Form (Symbolic):**
$$
\begin{aligned}
\min_{\mathbf{x}} \quad & f(\mathbf{x}) = \max_{(x,y) \in \Omega} \sigma_{vM}(x, y; \mathbf{x}) \\
\text{subject to} \quad & g_i(\mathbf{x}) \leq 0, \quad i = 1, \ldots, 5 \\
& \mathbf{x}^L \leq \mathbf{x} \leq \mathbf{x}^U
\end{aligned}
$$

**Constraints in Detail:**
- g₁: R_c + e_max(θ) + δ_c - x_aux ≤ 0  (central hole clearance)
- g₂: x_aux + e_x(θ) - (L - δ_e) ≤ 0  (right edge clearance)
- g₃: e_y(θ) - (W - δ_e) ≤ 0  (top/bottom edge clearance)
- g₄: a - α_max · b ≤ 0  (aspect ratio limit)

**Box Bounds (Side Constraints):**
$$\mathbf{x}^L \leq \mathbf{x} \leq \mathbf{x}^U$$

This notation means each design variable is bounded:
- x^L = lower bounds = [38, 6, 6, -π/2]
- x^U = upper bounds = [150, 15, 15, +π/2]

**Von Mises Stress (2D Plane Stress):**
$$\sigma_{vM} = \sqrt{\sigma_{xx}^2 - \sigma_{xx}\sigma_{yy} + \sigma_{yy}^2 + 3\tau_{xy}^2}$$

**Speaker Notes:**
- Emphasize this is a "min-max" problem (minimax) - we minimize the maximum stress
- Box bounds are simple inequality constraints on individual variables
- Von Mises stress is computed at every node in the FE mesh; we take the maximum

---

## Slide 7: Constraint Handling (1 min)

**Title:** Geometric Constraints

**Content - AABB Extent Formulas:**
For rotated ellipse:
- e_x(θ) = √(a²cos²θ + b²sin²θ)
- e_y(θ) = √(a²sin²θ + b²cos²θ)

**Constraints:**
1. **Central clearance:** R_c + e_max + δ_c ≤ x_aux
2. **Right edge:** x_aux + e_x ≤ L - δ_e
3. **Top/bottom edge:** e_y ≤ W - δ_e
4. **Aspect ratio:** a ≤ 3·b

**[FIGURE: Diagram showing AABB around rotated ellipse]**
*Annotation: Show a rotated ellipse with its bounding box, labeling e_x and e_y extents*

---

## Slide 8: Solution Methodology (1.5 min)

**Title:** Solution Approach

**Content - Algorithm:**
1. **Analysis:** 2D Plane-Stress FEA (MATLAB PDE Toolbox)
2. **Optimizer:** fmincon with SQP algorithm
3. **Gradients:** Finite differences (central)
4. **Multi-start:** Latin Hypercube Sampling (248 points)
5. **Parallelization:** parfor across starting points

**Why Multi-Start?**
- Non-convex objective (stress field is complex)
- Multiple local minima expected
- Global optimum not guaranteed, but robust exploration

**Speaker Notes:**
- SQP (Sequential Quadratic Programming) is efficient for smooth nonlinear problems
- Central differences more accurate than forward differences (but 2× cost)
- LHS provides better coverage of design space than random sampling

---

## Slide 9: FEA & Baseline Analysis (1.5 min)

**Title:** Finite Element Analysis & Baseline Reference

**FEA Implementation:**
| Parameter | Value |
|-----------|-------|
| Mesh type | Adaptive triangular |
| Elements | ~38,000 |
| Nodes | ~77,000 |
| Solve time | ~0.3-0.5 s per evaluation |

**Baseline Analysis (Central Hole Only):**
| Metric | Value |
|--------|-------|
| Max σ_vM | **360.28 MPa** |
| Stress Concentration Factor | K_t = 3.60 |
| Location | Top of central hole (0, 20) mm |

*This baseline serves as the reference for measuring optimization improvement.*

**[FIGURE: Baseline stress contour showing central hole only]**
*Annotation: Show the stress distribution for the baseline case. Highlight the max stress location at the top/bottom of the central hole where σ_max ≈ 3σ₀*

**Speaker Notes:**
- Theoretical K_t for circular hole in infinite plate is 3.0
- Our K_t = 3.6 is higher due to finite plate width effects
- Every optimization evaluation is compared against this baseline

---

## Slide 10: Convergence History (1 min)

**Title:** Multi-Start Optimization Convergence

**[FIGURE: Convergence history showing ALL runs - dedicated full slide]**
*Annotation: Use the convergence history plot with:
- Gray lines: All 248 optimization runs
- Blue bold line: Best run (Run #199)
- Red dashed line: Baseline stress (360.28 MPa)
- Highlight the spread of final objectives indicating multiple local minima*

**Key Observations:**
- Wide spread in final objectives (320 - 413 MPa)
- Many runs converge to local minima
- Best run achieves 11.1% reduction below baseline
- Justifies the multi-start approach

**Statistics:**
- Total function evaluations: 19,864
- Average per run: ~80 evaluations
- Total computation time: 2.6 hours

---

## Slide 11: Results - Optimal Design (1.5 min)

**Title:** Optimization Results

**Content - Comparison Table:**
| Metric | Baseline | Optimized | Change |
|--------|----------|-----------|--------|
| Max σ_vM | 360.3 MPa | 320.5 MPa | **-11.1%** |
| Kt | 3.60 | 3.20 | -11.0% |

**Optimal Parameters:**
- x_aux = 46.0 mm (2.3× central hole radius)
- a = b = 15.0 mm (**circular!**)
- θ = -35.8° (irrelevant for circle)

**Key Finding:** Optimizer converged to **circular** auxiliary holes at maximum allowed size.

**[FIGURE: Stress Comparison - 3 panels: baseline, initial, optimal]**
*Annotation: Use the stress comparison figure. Add arrows pointing to max stress locations.*

---

## Slide 12: Stress Distribution Analysis (1 min)

**Title:** Detailed Stress Analysis

**[FIGURE: Stress Detail Near Holes - zoomed view]**
*Annotation: Use the zoomed stress contour plot. Highlight:
- Max stress still at central hole (top/bottom)
- Auxiliary holes create local stress concentrations but lower than central
- Stress "shielding" effect visible*

**Key Observations:**
- Maximum stress remains at central hole boundary
- Auxiliary holes redistribute load away from critical region
- 11% reduction achieved with simple circular holes

---

## Slide 13: Why Circular? (1 min)

**Title:** Physical Interpretation

**Why did the optimizer choose circular holes?**

1. **No preferred orientation** - robust to load direction uncertainty
2. **Uniform stress distribution** around perimeter
3. **No sharp features** - ellipse tips create higher local Kt
4. **Maximum size** - hit upper bound (15 mm), suggesting larger = better

**[FIGURE: Final Objective Distribution histogram]**
*Annotation: Use the histogram showing spread of final objectives. Highlight that best solutions cluster around 320 MPa.*

**Speaker Notes:**
- Discuss that for more complex loading (biaxial, shear), elliptical shapes might emerge
- The "circular is optimal" result is specific to uniaxial tension

---

## Slide 14: Fatigue Life Implications (1.5 min)

**Title:** Fatigue Life Extension

**Content - S-N Curve Approach:**
Life improvement: N_opt/N_base = (σ_base/σ_opt)^m

| Fatigue Exponent (m) | Application | Life Extension |
|---------------------|-------------|----------------|
| 3 | Welded structures | **1.4×** |
| 5 | Machined parts | **1.8×** |
| 6 | Typical steel | **2.0×** |
| 8 | High-quality steel | **2.6×** |

**Practical Impact:**
- For m ≈ 5-6: Fatigue life approximately **doubles**
- Component lasting 100k cycles → **180k-200k cycles**

**Speaker Notes:**
- Emphasize the nonlinear relationship between stress reduction and life extension
- 11% stress reduction → 40-155% life increase depending on material

---

## Slide 15: Computational Cost (0.5 min)

**Title:** Computational Summary

| Metric | Value |
|--------|-------|
| Multi-start runs | 248 |
| Total FEA evaluations | 19,864 |
| Total time | 2.6 hours |
| Average per run | ~80 evaluations |
| Best run found | Run #199 |

**Observations:**
- Wide spread in final objectives (std = 19.7 MPa)
- Confirms multiple local minima
- Multi-start essential for this problem

---

## Slide 16: Conclusions (1 min)

**Title:** Conclusions

**Key Takeaways:**
1. ✓ Auxiliary holes can reduce stress concentration by **11%**
2. ✓ Optimal shape is **circular** for uniaxial tension
3. ✓ Fatigue life improvement of **1.8-2× for typical steel**
4. ✓ Multi-start optimization essential due to multiple local minima

**Limitations & Future Work:**
- 2D analysis (3D effects near surfaces)
- Single load case (multi-objective for varying loads)
- Manufacturing constraints not fully considered
- Could explore more auxiliary holes

---

## Slide 17: Questions?

**Title:** Questions & Discussion

**Content:**
- Contact information
- Code available at: [repository link]

**[FIGURE: Animation frame or attractive stress plot]**
*Annotation: Use a visually appealing frame from the optimization animation showing the final optimal design*

---

# Backup Slides

## Backup Slide A: AABB Derivation

**Title:** AABB Extent Derivation

Show the mathematical derivation of the axis-aligned bounding box formulas for a rotated ellipse.

## Backup Slide B: Algorithm Parameters

**Title:** fmincon Settings

```matlab
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1500, ...
    'OptimalityTolerance', 1e-4, ...
    'FiniteDifferenceType', 'central');
```

## Backup Slide C: Boundary Stress Distribution

**[FIGURE: Boundary stress plot around central hole]**
*Annotation: Use the boundary stress distribution plot showing σ_vM vs θ around the central hole perimeter*

---

# Figure Checklist

Figures to export from MATLAB (use `exportgraphics` or `saveas`):

1. **geometry_comparison.png** - 3-panel geometry comparison
2. **baseline_stress.png** - Baseline stress contour (central hole only)
3. **stress_comparison.png** - 3-panel stress field comparison
4. **stress_detail_zoomed.png** - Zoomed view near holes
5. **convergence_history.png** - All runs + best run convergence (FULL SIZE for dedicated slide)
6. **objective_distribution.png** - Histogram of final objectives
7. **optimization_animation.gif** - Animation of best run

**MATLAB export commands:**
```matlab
% After regenerating plots
exportgraphics(gcf, 'figure_name.png', 'Resolution', 300)

% For specific figures
fig = figure(1);
exportgraphics(fig, 'geometry_comparison.png', 'Resolution', 300)
```

---

# Timing Guide

| Section | Slides | Time |
|---------|--------|------|
| Introduction & Motivation | 1-3 | 3.0 min |
| Problem Formulation | 4-7 | 4.5 min |
| Methodology & Baseline | 8-10 | 4.0 min |
| Results & Analysis | 11-14 | 4.5 min |
| Conclusions | 15-17 | 1.5 min |
| **Total** | **17** | **~17.5 min** |

*Note: Slightly over 15 min - can trim fatigue slide or computational cost slide if needed*

---

# Presentation Tips

1. **Start with the counter-intuitive hook:** "We're going to reduce stress by adding holes"

2. **Use animations sparingly:** The GIF is great for showing optimization progress, but limit to one slide

3. **Emphasize practical impact:** Graduate students appreciate real-world relevance (fatigue life doubling)

4. **Be prepared for questions about:**
   - Why not use gradient-based global optimization?
   - How sensitive is the result to mesh refinement?
   - What about manufacturing tolerances?
   - Would topology optimization give better results?

5. **Key insight to emphasize:** The circular result is elegant - it suggests that for symmetric loading, the simplest auxiliary hole shape is optimal
