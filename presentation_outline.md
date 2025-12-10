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

**Content - Canonical Form:**
```
minimize    f(x) = max σ_vM(x,y; x)
   x        over domain Ω

subject to  g₁: Central hole clearance ≤ 0
            g₂: Edge clearance (x) ≤ 0
            g₃: Edge clearance (y) ≤ 0
            g₄: Aspect ratio ≤ 0
            x^L ≤ x ≤ x^U
```

**Key Points:**
- Objective: Minimize maximum von Mises stress
- 5 nonlinear inequality constraints
- Box bounds on all variables
- **Note:** Objective requires FEA solve at each evaluation

**Speaker Notes:**
- Emphasize this is a "min-max" problem (minimax)
- Constraints use AABB (Axis-Aligned Bounding Box) for rotated ellipse clearance

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

**[FIGURE: Convergence history showing all runs]**
*Annotation: Use the convergence history plot showing gray lines (other runs) and blue line (best run). Highlight the spread indicating multiple local minima.*

---

## Slide 9: FEA Implementation (1 min)

**Title:** Finite Element Analysis

**Content:**
- **Mesh:** Adaptive with Hmax = min_feature / 2
- **Elements:** ~38,000 triangular elements
- **Nodes:** ~77,000
- **Solve time:** ~0.3-0.5 seconds per evaluation

**Von Mises Stress:**
σ_vM = √(σ_xx² - σ_xx·σ_yy + σ_yy² + 3τ_xy²)

**[FIGURE: Mesh visualization or stress contour plot]**
*Annotation: Show the mesh near the holes, or a stress contour plot from the optimal design*

---

## Slide 10: Results - Optimal Design (1.5 min)

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

## Slide 11: Stress Distribution Analysis (1 min)

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

## Slide 12: Why Circular? (1 min)

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

## Slide 13: Fatigue Life Implications (1.5 min)

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

## Slide 14: Computational Cost (0.5 min)

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

## Slide 15: Conclusions (1 min)

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

## Slide 16: Questions?

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
2. **stress_comparison.png** - 3-panel stress field comparison
3. **stress_detail_zoomed.png** - Zoomed view near holes
4. **convergence_history.png** - All runs + best run convergence
5. **objective_distribution.png** - Histogram of final objectives
6. **optimization_animation.gif** - Animation of best run

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
| Methodology | 8-9 | 2.5 min |
| Results & Analysis | 10-13 | 4.0 min |
| Conclusions | 14-16 | 1.0 min |
| **Total** | **16** | **15.0 min** |

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
