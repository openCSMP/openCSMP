<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Square-cavity thermal convection (CVFEM_tests_main, 2D tests test-206a/b/c)

Benchmark of **Geiger, Driesner, Heinrich & Matthäi (2006)**, *Multiphase thermohaline convection in the Earth's
crust: II. Benchmarking and application of a finite element – finite volume solution technique with a NaCl–H2O
equation of state*, Transport in Porous Media 63, 435–461, section 3.2.2 and Table II. It was used to verify the
predecessor of this code; here it assesses how the current CVFEM PHX scheme, and in particular its buoyancy
(gravity) term, compares.

| Ra_th | Caltagirone (1975) | Geiger et al. (2006) | test limit |
|---|---|---|---|
| 100 | 2.65 | 2.560 | ±5 % of Caltagirone |
| 200 | 3.81 | 3.818 | ±5 % |
| 300 | 4.52 | 4.503 | ±5 % |

- Closed square 100 m × 100 m, impermeable on all sides; base 150 °C, top 140 °C; insulated side walls.
- Closed top: fluid pressure fixed at the top-left corner node only (a pressure reference, not an open boundary).
- Nu = −dT̄/dz̄ at the top (eq. 11); Ra = (ρ_top − ρ_base) g H k / (μ κ), κ = λ / (ρ c_p)_water (eq. 10).
- Rock heat capacity 100 J kg⁻¹ K⁻¹ (small, to shorten the transient; steady Nu does not depend on it).
- Geiger et al. solved the Boussinesq problem (linear ρ(T), constant μ) and give no temperatures. openCSMP uses the
  full water equation of state, so 140–150 °C at 5 MPa is chosen to stay near the Boussinesq limit (μ varies ±3.6 %,
  Δρ = 9.07 kg m⁻³). Permeability per case from IAPWS water at 145 °C, 5 MPa, scaled (×1.056) so the Ra realised by the model (its Δρ includes the hydrostatic pressure difference) is the target; the test reports it.
- Mesh `test-206`: 40 × 40 squares split into triangles with alternating diagonals (structured, as required for
  the top-boundary gradient).
