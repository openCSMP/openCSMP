<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Hydrostatic column with a conductive temperature gradient (CVFEM_tests_main, 2D tests test-205a-d)

A pure-water column 500 m × 1000 m (top 10 °C and 1 bar fixed, base 110 °C fixed, insulated no-flow sides,
permeability 10⁻¹⁵ m²), initialised with the linear conductive temperature profile and the hydrostatic pressure.
Its Rayleigh number (~3) is far below the onset of convection (4π², Horton & Rogers 1945, Lapwood 1948), so the exact
solution is no flow at all times; any flow is discretisation error of the buoyancy (gravity) term.

The test passes if the mass transported between interior control volumes (the transport scheme's
"fluid flux integral", as a mean Darcy speed) stays below 10⁻² of the column buoyancy velocity scale
k g (ρ_max − ρ_min) / μ (bulk fluid density, smallest viscosity of the phases present). The stored element velocity ("velocity liquid", a visualisation quantity) and the transport
at Dirichlet nodes are reported for information. Each mesh is run with the donor-node gravity density (default) and
with CVFEM_PHX_Scheme::UseFacetMeanGravityDensity(true).

On this column at rest the transported mass with either density is round-off (10⁻¹³–10⁻¹² m/s): the upwind
direction at facets with zero velocity is decided by round-off, so the ratio changes with memory layout (for example
with the number of model properties) between 0 and 2 × 10⁻³ of the velocity scale. The limit, 10⁻², leaves room above
that noise; a faulty gravity term gives spurious flows of about 0.1 of the scale. The donor density's real errors show
in the near-critical column (test-205g/h) and in the convection-onset test (test-207).

- `test-205-structured`: 20 × 40 squares split into triangles with alternating diagonals (25 m).
- `test-205-unstructured`: gmsh Frontal-Delaunay triangles, target size 25 m.

## Near-critical column (test-205g/h)

A column 100 m × 200 m at 25 MPa (above the critical pressure of water), 420 °C at the top and 360 °C at the base:
hot, light fluid above cold, dense fluid, so it is stable at any Rayleigh number and the exact solution is again no
flow. It crosses the pseudo-critical temperature at 25 MPa (~385 °C), where the density falls from 590 to 134 kg/m³,
much of it within 10 K. The equation of state labels the low-density part of this single supercritical fluid vapour,
so the column also carries a liquid/vapour label boundary that is not a phase transition. Meshes in `test-205-nc/`
(the test-205 meshes scaled by 0.2: 5 m). Run with the facet-mean density only: with the donor density the
spurious velocities are 85–120 times larger and the time step collapses (17–42 times more steps), which makes it too
slow for the test suite.
