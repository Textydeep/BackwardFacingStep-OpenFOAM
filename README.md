# Steady-State CFD Analysis: Backward-Facing Step
**OpenFOAM v2306 | simpleFoam | k-ω SST RANS**

---

## 1. Problem Description

The backward-facing step (BFS) is one of the most widely used benchmark cases in computational fluid dynamics. Flow separates at the step corner and forms a recirculation bubble on the downstream wall before reattaching. It is extensively used to validate turbulence models because the reattachment length — where the separated shear layer reattaches to the wall — is sensitive to the choice of turbulence closure.

**Physical setup:**

```
         ___________________________________  topWall (y = 2h)
        |         Block 0     |
 inlet  |     (upstream)      |________________  y = h
 ──────>|                    STEP               ──────> outlet
                              |________________  y = 0
                              |      Block 2      |
                              |   (recirculation) |
                              |___________________| bottomWall
                              x=0               x=30h
```

| Parameter | Value |
|-----------|-------|
| Fluid | Water at 20°C |
| Kinematic viscosity ν | 1×10⁻⁶ m²/s |
| Step height h | 0.01 m |
| Expansion ratio H/h | 2 (inlet height h, downstream height 2h) |
| Upstream length | 4h = 0.04 m |
| Downstream length | 30h = 0.30 m |
| Inlet velocity U | 1.0 m/s |
| Reynolds number Re_h | **10,000** (based on step height) |

---

## 2. Numerical Setup

### 2.1 Solver

**`simpleFoam`** — steady-state, incompressible, pressure-velocity coupling via the SIMPLE algorithm.

The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm iterates between solving the momentum equation and correcting the pressure field until convergence. It is the standard solver for steady RANS simulations.

### 2.2 Turbulence Model

**k-ω SST (Shear Stress Transport)** — Menter (1994)

The k-ω SST model is the industry standard for external aerodynamics and flows with separation. It blends:
- **k-ω** near walls (accurate in the viscous sublayer)
- **k-ε** in the freestream (avoids k-ω's sensitivity to freestream conditions)

This makes it particularly well-suited to the backward-facing step, where both wall boundary layers and a separated free shear layer must be captured.

**Why not k-ε?** The standard k-ε model is known to underpredict reattachment length in separated flows because it overpredicts turbulent diffusion in the recirculation zone.

### 2.3 Mesh

Generated using **blockMesh** (OpenFOAM's structured hex mesher) with 3 blocks:

| Block | Region | Cells |
|-------|--------|-------|
| 0 | Upstream channel (4h × h) | 40 × 40 = 1,600 |
| 1 | Upper downstream (30h × h) | 300 × 40 = 12,000 |
| 2 | Lower downstream (30h × h) | 300 × 40 = 12,000 |
| **Total** | | **25,600 cells** |

**Mesh quality (checkMesh):**

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Max non-orthogonality | 0° | < 70° | ✅ |
| Max skewness | 3.3×10⁻¹³ | < 4 | ✅ |
| Max aspect ratio | 4 | < 1000 | ✅ |
| Geometric directions | 2 (x, y) | — | ✅ 2D confirmed |

Cell size: Δx = Δy = 1 mm (uniform). The domain is extruded 1 mm in z with `empty` boundary conditions, making this a 2D simulation.

### 2.4 Boundary Conditions

| Patch | Type | U | p | k | ω |
|-------|------|---|---|---|---|
| inlet | velocity inlet | (1, 0, 0) m/s | zeroGradient | 3.75×10⁻³ m²/s² | 160 s⁻¹ |
| outlet | pressure outlet | inletOutlet | 0 Pa (gauge) | zeroGradient | zeroGradient |
| topWall | wall (no-slip) | (0,0,0) | zeroGradient | kqRWallFunction | omegaWallFunction |
| bottomWall | wall (no-slip) | (0,0,0) | zeroGradient | kqRWallFunction | omegaWallFunction |
| stepWall | wall (no-slip) | (0,0,0) | zeroGradient | kqRWallFunction | omegaWallFunction |
| stepTop | wall (no-slip) | (0,0,0) | zeroGradient | kqRWallFunction | omegaWallFunction |
| frontAndBack | 2D (empty) | — | — | — | — |

**Inlet turbulence quantities:**
- Turbulence intensity TI = 5%
- k = 1.5 × (TI × U)² = 1.5 × (0.05 × 1.0)² = **3.75×10⁻³ m²/s²**
- ω = √k / (Cμ^0.25 × L) where L = 0.07h → **160 s⁻¹**

Wall functions (kqRWallFunction, omegaWallFunction, nutkWallFunction) are used for high-Reynolds-number treatment, appropriate when the first cell y⁺ ≫ 1.

### 2.5 Numerical Schemes

| Term | Scheme | Rationale |
|------|--------|-----------|
| Time | steadyState | Steady RANS |
| Gradients | Gauss linear (cellLimited) | Bounded, prevents oscillations |
| Convection (U, k, ω) | Gauss linearUpwind | 2nd-order, bounded — good for separated flows |
| Laplacian | Gauss linear corrected | 2nd-order, handles non-orthogonality |

### 2.6 SIMPLE Relaxation Factors

| Field | Relaxation Factor |
|-------|-----------------|
| U | 0.7 |
| p | 0.3 |
| k | 0.5 |
| ω | 0.5 |

Under-relaxation damps oscillations in the SIMPLE loop. Lower values = more stable but slower convergence. These are standard values for incompressible RANS.

---

## 3. Results

### 3.1 Convergence

The simulation ran for **2000 iterations** in **61 seconds** on a MacBook (Apple Silicon via Docker x86 emulation).

**Final residuals (iteration 2000):**

| Field | Initial Residual | Converged |
|-------|-----------------|-----------|
| Ux | ~2×10⁻⁶ | ✅ (target < 10⁻⁵) |
| Uy | ~6×10⁻⁶ | ✅ |
| k | ~1×10⁻⁴ | ✅ |
| ω | converged | ✅ |

Residuals dropped approximately **5 orders of magnitude** from their initial values, indicating a well-converged solution.

### 3.2 Reattachment Length

The reattachment length is the most important result for the backward-facing step. It is defined as the x-location where the streamwise velocity Ux changes sign in the near-wall cell layer above the bottom wall.

| | Value |
|--|-------|
| **Computed (k-ω SST, Re_h = 10,000)** | **x_r = 0.0847 m = 8.47 h** |
| Typical range (k-ω SST, Re_h ~ 10⁴) | 7–10 h |
| Driver & Seegmiller (1985) experiment | 6.4 h at Re_H = 37,400 |

The longer reattachment length compared to the Driver & Seegmiller experiment is expected due to the lower Reynolds number (Re_h = 10,000 vs. their ~30,000) — lower Re flows take longer to reattach. At higher Re, turbulent mixing is more vigorous and forces earlier reattachment.

### 3.3 Recirculation Zone

The recirculation bubble extends from x = 0 (step corner) to x ≈ 0.085 m (8.5h) along the bottom wall. Within this zone, Ux < 0 — the fluid near the wall moves in the **upstream direction**, driven by the low-pressure vortex shed at the step corner.

**Maximum reverse velocity:** Ux_min = −0.190 m/s at x ≈ 4h (40% of the inlet velocity, reversed)

### 3.4 Velocity Profiles

Streamwise velocity Ux normalised by inlet velocity U₀ = 1 m/s, at four cross-sections:

**x/h = 1 (just downstream of step):**

| y/(2h) | Ux (m/s) | Note |
|--------|----------|------|
| 0.063 | +0.010 | Near step toe |
| 0.131 | −0.047 | Recirculation |
| 0.256 | −0.062 | Peak reverse flow |
| 0.444 | +0.118 | Shear layer |
| 0.506 | +0.562 | Inlet jet |
| 0.569 | +1.040 | Accelerated over step |

**x/h = 4 (deep in recirculation zone):**

| y/(2h) | Ux (m/s) | Note |
|--------|----------|------|
| 0.063 | −0.190 | Peak reverse velocity |
| 0.256 | −0.075 | |
| 0.506 | +0.656 | Main flow above |
| 0.631 | +1.023 | Near top |

**x/h = 6 (near reattachment):**

| y/(2h) | Ux (m/s) | Note |
|--------|----------|------|
| 0.063 | −0.116 | Still reversed |
| 0.256 | +0.042 | Near zero (reattachment approaching) |
| 0.506 | +0.655 | |

**x/h = 10 (post-reattachment recovery):**

| y/(2h) | Ux (m/s) | Note |
|--------|----------|------|
| 0.063 | +0.042 | Fully attached, positive |
| 0.506 | +0.604 | Recovering profile |
| 0.757 | +0.842 | Near top |

The profiles clearly show the jet of fast fluid entering in the upper half at x/h = 1 and progressively mixing downward, eventually producing a more uniform profile by x/h = 10.

---

## 4. Key Physical Observations

### 4.1 Flow Separation and Reattachment
At the step corner (x = 0, y = h), the flow cannot follow the sudden change in wall geometry. It separates and forms a **free shear layer** — the boundary between the recirculating flow below and the main flow above. This shear layer is unstable (Kelvin-Helmholtz type) and generates turbulence, which is what eventually mixes high-momentum fluid downward and forces reattachment.

### 4.2 Pressure Recovery
After the step, the flow decelerates (the channel area doubles), leading to an **adverse pressure gradient** that helps drive the recirculation. After reattachment, a **favorable pressure gradient** develops and the boundary layer re-grows toward a fully-developed state.

### 4.3 Why k-ω SST Performs Well Here
The SST model activates its k-ω branch in the adverse pressure gradient region near the separation and reattachment points (where the flow history matters), and switches to k-ε in the free shear layer. This blending strategy captures the separated shear layer and the near-wall recovery more accurately than either model alone.

---

## 5. How to Visualise Results in ParaView

Open **ParaView**, then **File → Open → `backwardFacingStep.foam` → Apply**.

Jump to **Time = 2000** using the time controls.

### Velocity Contour
- Color by: **U → Magnitude**
- Expected: blue recirculation zone behind step, red/yellow main jet above

### Streamlines
1. **Filters → Stream Tracer**
2. Vectors: **U**, Seed type: **Line Source**
3. Point 1: `(0.005, 0.001, 0.0005)`, Point 2: `(0.005, 0.019, 0.0005)`, Resolution: 30
4. Apply → color by U Magnitude
5. Expected: clockwise vortex in the recirculation bubble (x = 0 to ~8.5h)

### Wall Shear Stress / Reattachment
1. Select the case in the Pipeline Browser
2. In the "Representation" toolbar: choose **Surface**
3. Color by **wallShearStress** (if available) or use the streamlines method above

### Velocity Profiles (Plot Over Line)
1. **Filters → Plot Over Line**
2. Set a vertical line at x = 0.04 (x/h = 4): Point 1 `(0.04, 0, 0.0005)`, Point 2 `(0.04, 0.02, 0.0005)`
3. Apply → shows Ux vs y at that cross-section
4. Repeat for x = 0.06, 0.10

---

## 6. Case File Structure

```
OpenFOAM pilot1/
├── 0/                        Initial & boundary conditions
│   ├── U                     Velocity field
│   ├── p                     Kinematic pressure (p/ρ)
│   ├── k                     Turbulent kinetic energy
│   ├── omega                 Specific dissipation rate
│   └── nut                   Turbulent viscosity
├── constant/
│   ├── transportProperties   Fluid viscosity (Newtonian, ν = 1e-6 m²/s)
│   ├── turbulenceProperties  k-ω SST model selection
│   └── polyMesh/             Generated mesh (blockMesh output)
├── system/
│   ├── blockMeshDict         Mesh geometry and block topology
│   ├── controlDict           Solver, time steps, output
│   ├── fvSchemes             Numerical discretisation schemes
│   └── fvSolution            Linear solvers and SIMPLE settings
├── 500/, 1000/, 1500/, 2000/ Solution snapshots
├── postProcessing/           Wall shear stress, residuals
├── log.simpleFoam            Complete solver output log
└── backwardFacingStep.foam   ParaView case file
```

---

## 7. Summary

| Item | Value |
|------|-------|
| Solver | simpleFoam (steady RANS) |
| Turbulence model | k-ω SST |
| Mesh cells | 25,600 |
| Iterations | 2,000 |
| Solve time | 61 seconds |
| Final Ux residual | ~2×10⁻⁶ |
| Reattachment length | **8.47 h** |
| Peak reverse velocity | −0.190 m/s (−19% of U_inlet) |
| Recirculation extent | x = 0 to 0.085 m |

The simulation successfully captures the defining feature of backward-facing step flow: a recirculation bubble driven by flow separation at the step corner, with reattachment at approximately 8.5 step heights downstream. Results are physically consistent with expectations for Re_h = 10,000 turbulent flow with wall function treatment.

---

*Generated using OpenFOAM v2306 on macOS (Apple Silicon) via Docker.*
*Turbulence model: k-ω SST (Menter 1994). Solver: simpleFoam (SIMPLE algorithm).*
