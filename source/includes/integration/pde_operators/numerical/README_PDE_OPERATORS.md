/# PDE Operator Naming Convention and Catalogue

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

| Token  | Mathematical symbol | Meaning |
|--------|--------------------|---------| 
| `N`    | $$N$$              | Scalar interpolation (test/trial) function |
| `NT`   | $$N^T$$            | Transposed scalar interpolation function |
| `DN`   | $$\nabla N$$       | Gradient matrix of scalar solution variable |
| `DNT`  | $$(\nabla N)^T$$   | Transposed gradient matrix of scalar solution variable |
| `B`    | $$B$$              | Strain-displacement matrix for **vector** solution variable (mechanics) |
| `BT`   | $$B^T$$            | Transposed strain-displacement matrix |
| `op`   | $$[\sigma]$$       | Material operand (scalar, vector or tensor property) |
| `lhsop`| $$[\sigma]$$       | Material operand on the left-hand side matrix |
| `rhsop`| $$[\sigma]$$       | Material operand on the right-hand side vector |
| `v`    | $$\mathbf{v}$$     | Advection velocity vector |
| `dV`   | $$dV$$             | Volume integral |
| `dS`   | $$dS$$             | Surface (face) integral |

### Key distinction: `B` vs `DN`

| Symbol | DOF type        | Typical application |
|--------|-----------------|---------------------|
| `B`    | **Vector** DOFs | Mechanics (displacement, strain) |
| `DN`   | **Scalar** DOFs | Flow and transport (pressure, concentration) |

This distinction follows standard mechanics textbook notation where `B` is
the strain-displacement matrix, while `DN` is used for scalar gradient
operators in flow and transport problems.

### Side suffix

| Suffix | Meaning |
|--------|---------|
| _(none)_ | Left-hand side matrix operator (`MathOperatorLHS`) |
| `_rhs`   | Right-hand side vector operator (`MathOperatorRHS`) |

In practice the LHS/RHS distinction is encoded in the base class rather
than the name, but the mathematical form of the integral makes it clear.

---

## Numerical vs Analytical Integration

CSMP PDE operators are divided into two families distinguished by their
name prefix: b

### `NumIntegral` — Numerical Integration (Isoparametric Elements)

Operators prefixed with `Num` use **Gauss quadrature** in a parametric
reference space. They require isoparametric finite elements where the
geometry and the solution field are interpolated using the same shape
functions.

```cpp
assert( e.FE()->Isoparametric() == true );  // enforced in ComputeContribution


## Operator Catalogue

### Left-Hand Side (Matrix) Operators

These operators inherit from `MathOperatorLHS<dim,CELL>` and assemble
contributions into the global stiffness or mass matrix.

---

#### `NumIntegral_dNT_lhsop_dN_dV`

**Known as:** Stiffness matrix / conductance matrix

**Integral:**
$$K_{jk} = \int_{\Omega^e} (\nabla N_j)^T \, [\sigma] \, \nabla N_k \, dV$$

**Operand `op`:** Scalar, vector (diagonal anisotropy) or full tensor
diffusivity/conductivity/permeability — element-placed.

**Test variable:** Scalar, node-placed (e.g. `"fluid pressure"`).

**Application:** Pressure diffusion, heat conduction, species transport.

**Constructor:**
```cpp
NumIntegral_dNT_lhsop_dN_dV( const PropertyDatabase<dim>&,
                           const char* oper,   // e.g. "conductivity"
                           const char* basic,  // e.g. "fluid pressure"
                           const char* test ); // e.g. "fluid pressure"

---

#### `NumIntegral_NT_lhsop_N_dV`

**Known as:** Capacitance matrix / mass matrix (LHS)

**Integral:**
$$C_{jk} = \int_{\Omega^e} N_j \, \sigma \, N_k \, dV$$

**Operand `op`:** Scalar — element, integration point or node-placed.

**Test variable:** Scalar, node-placed.

**Lumped formulation:** Supported — produces diagonal matrix whose
diagonal equals the row sums of the consistent mass matrix.

**Application:** Transient storage term in pressure diffusion,
compressibility matrix in poromechanics.

**Constructor:**
```cpp
NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>&,
                            const char* oper,   // e.g. "compressibility"
                            const char* basic,  // e.g. "fluid pressure"
                            const char* test ); // e.g. "fluid pressure"
