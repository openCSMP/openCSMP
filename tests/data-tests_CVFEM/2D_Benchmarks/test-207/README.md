<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Convection onset: growth or decay rate of a convective mode (CVFEM_tests_main, 2D tests test-207a–e)

Transient, analytical check of the buoyancy (gravity) term together with conduction and the bulk heat capacity, on the
Horton-Rogers-Lapwood problem (Horton & Rogers 1945; Lapwood 1948): a closed porous layer heated from below. Weatherill,
Simmons, Voss & Robinson (2004, *Advances in Water Resources* 27, 547–562) propose this problem as a benchmark for
density-driven flow codes because its stability criteria are analytical (onset at Ra_c = 4π²). This test measures the
growth or decay RATE of the first mode below and above onset, which linear stability theory gives exactly:

| test | Ra | σ = Ra/2 − 2π² | behaviour |
|---|---|---|---|
| test-207a | 20 | −9.70 | decays |
| test-207b | 30 | −4.68 | decays |
| test-207c | 60 | +10.39 | grows |
| test-207d | 100 | +30.47 | grows |
| test-207e | 30, 10 % perturbation | −4.68 | must decay (any amplitude) |

with τ = λ t / (C_bulk H²), C_bulk = (1 − φ) ρ_r c_r + φ ρ_f c_f, Ra = ρ_f β_T ΔT g H k (ρ c)_f / (μ λ).

- Box 200 m × 100 m (L = 2 H), impermeable on all sides; base 150 °C, top 140 °C; insulated side walls; fluid pressure
  fixed at the top-left corner node only (a pressure reference, not an open boundary), as in test-206.
- Initial state: conductive profile ± A0 cos(πx/H) sin(πz/H), hydrostatic pressure. The mode satisfies every boundary
  condition and is an exact eigenmode of the linearised Darcy–Boussinesq equations.
- Each case runs twice, with the mode added and subtracted. Half the difference of the two mode amplitudes is the
  linear response; half the sum collects anything that does not depend on the sign (e.g. spurious buoyancy flow of the
  discretisation, quadratic nonlinearity) and is reported.
- The run stops at |σ| τ = 3; A0 is 1 % of ΔT for a decaying mode and 1 % × e⁻³ for a growing one (linear regime).
  Time step capped at |σ| Δτ = 0.005. ln A is fitted against τ over the last 75 % of the run.
- Pass: |Ra_σ − Ra| < 3 % of Ra_c, with Ra_σ = 2 (σ + 2π²) the Rayleigh number implied by the fitted rate (|σ − σ_theory|
  < 0.59). The tolerance is relative to Ra_c, not Ra, because near onset σ is a small difference of two large terms.
- Ra uses the thermal density contrast at constant pressure from the model's own equation of state (9.07 kg m⁻³); the
  hydrostatic compression with depth (top-to-base difference 8.56 kg m⁻³) does not drive flow.
- Real water (not Boussinesq): 140–150 °C at 5 MPa keeps it close to that limit, as in test-206. Rock heat capacity
  100 J kg⁻¹ K⁻¹ (thermal time unit about 100 yr); fluid included in the conduction heat capacity.
- Mesh `test-207`: 40 × 20 squares split into triangles with alternating diagonals (structured). On 80 × 40 the
  default PETSc settings do not converge; with an algebraic multigrid preconditioner the result is the same.
- test-207e: below onset the Darcy–Boussinesq layer is also globally stable: the energy method gives the same threshold
  4π² for all initial data (Straughan 2004, *The Energy Method, Stability, and Nonlinear Convection*, 2nd ed.,
  Springer, Applied Mathematical Sciences 91, sec. 4.2.1, pp. 64–68, "R_L = R_E", shown for a periodic layer; the
  impermeable, insulated side walls of the box leave the same boundary terms zero), so a perturbation of ANY size must
  decay. It starts at 10 % of ΔT (1 K, both signs) and passes if both runs
  have decayed by at least 5 when |σ| τ = 3 (linear theory: e³ ≈ 20). A scheme can pass the small-amplitude cases
  and still sustain convection below onset once the perturbation is large enough; this case catches that.
- Fields of the first run are written at 40 checkpoints (`test-207?_Model_(REGION)_<n>.vtu`).
