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
"fluid flux integral", as a mean Darcy speed) stays below 10⁻³ of the column buoyancy velocity scale
k g (ρ_max − ρ_min) / μ. The stored element velocity ("velocity liquid", a visualisation quantity) and the transport
at Dirichlet nodes are reported for information. Each mesh is run with the donor-node gravity density (default) and
with CVFEM_PHX_Scheme::UseFacetMeanGravityDensity(true).

- `test-205-structured`: 20 × 40 squares split into triangles with alternating diagonals (25 m).
- `test-205-unstructured`: gmsh Frontal-Delaunay triangles, target size 25 m.
