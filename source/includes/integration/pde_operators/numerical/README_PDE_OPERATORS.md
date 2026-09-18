<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# PDE Operator Naming Convention and Catalogue

## Overview

This directory contains numerically integrated finite element PDE operators
used to assemble the global stiffness matrix, mass matrix, and right-hand
side vectors for flow, transport and mechanics problems.

All operators follow a strict naming convention that encodes the mathematical
form of the integral directly in the class name.

---

## Naming Convention

### Prefix

| Prefix | Meaning |
|--------|---------|
| `Num`  | Numerically integrated (Gauss quadrature) |

### Tokens

| Token   | Mathematical symbol | Meaning |
|---------|---------------------|---------|
| `N`     | N                   | Scalar interpolation (test/trial) function |
| `NT`    | Nᵀ                  | Transposed scalar interpolation function |
| `dN`    | ∇N                  | Gradient matrix of scalar solution variable |
| `dNT`   | (∇N)ᵀ               | Transposed gradient matrix of scalar solution variable |
| `B`     | B                   | Strain-displacement matrix for **vector** solution variable (mechanics) |
| `BT`    | Bᵀ                  | Transposed strain-displacement matrix |
| `P`     | P                   | Interpolation matrix for **vector** solution variable (mechanics) |
| `PT`    | Pᵀ                  | Transposed interpolation matrix for vector solution variable |
| `op`    | [σ]                 | Material operand (scalar, vector or tensor property) |
| `lhsop` | [σ]                 | Material operand contributing to the left-hand side matrix |
| `rhsop` | [σ]                 | Material operand contributing to the right-hand side vector |
| `D`     | D                   | Constitutive (elasticity) matrix, function of E and ν |
| `v`     | **v**               | Advection velocity vector |
| `dV`    | dV                  | Volume integral over element domain Ω |
| `dS`    | dS                  | Surface (face) integral over boundary Γ |

### Key distinction: `B`/`P` vs `dN`/`N`

| Symbol | DOF type        | Typical application |
|--------|-----------------|---------------------|
| `B`    | **Vector** DOFs | Mechanics — strain-displacement matrix relating nodal displacements to strains |
| `P`    | **Vector** DOFs | Mechanics — interpolation matrix relating nodal displacements to displacement field |
| `dN`   | **Scalar** DOFs | Flow and transport — gradient of scalar field (pressure, concentration) |
| `N`    | **Scalar** DOFs | Flow and transport — interpolation of scalar field |

This distinction follows standard mechanics textbook notation where `B` is
the strain-displacement matrix and `P` is the displacement interpolation
matrix, while `dN` and `N` are used for scalar gradient and interpolation
operators in flow and transport problems.

### Side suffix

| Suffix   | Meaning |
|----------|---------|
| _(none)_ | Left-hand side matrix operator (`MathOperatorLHS`) |
| `_rhs`   | Right-hand side vector operator (`MathOperatorRHS`) |

In practice the LHS/RHS distinction is encoded in the base class rather
than the name, but the mathematical form of the integral makes it clear.

---

## Numerical Integration

All operators prefixed with `Num` use **Gauss quadrature** in parametric
reference space. They require isoparametric finite elements where the
geometry and the solution field are interpolated using the same shape
functions:

```cpp
assert( e.FE()->Isoparametric() == true );  // enforced in ComputeContribution
