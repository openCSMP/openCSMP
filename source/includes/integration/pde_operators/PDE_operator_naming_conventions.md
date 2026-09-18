<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# CSMP PDE Operator Naming Convention — Summary

## Token Reference

| Token   | Mathematical symbol          | Meaning                                              | DOF type       |
|---------|------------------------------|------------------------------------------------------|----------------|
| `N`     | $$N$$                        | Scalar interpolation (test/trial) function           | Scalar         |
| `NT`    | $$N^T$$                      | Transposed scalar interpolation function             | Scalar         |
| `dN`    | $$\nabla N$$                 | Gradient matrix of scalar solution variable          | Scalar         |
| `dNT`   | $$(\nabla N)^T$$             | Transposed gradient matrix of scalar solution variable | Scalar       |
| `B`     | $$B$$                        | Strain-displacement matrix                           | Vector (mechanics) |
| `BT`    | $$B^T$$                      | Transposed strain-displacement matrix                | Vector (mechanics) |
| `op`    | $$[\sigma]$$                 | Material operand (scalar, vector or tensor)          | —              |
| `lhsop` | $$[\sigma]$$                 | Material operand — used only when both LHS and RHS exist for the same integral form | — |
| `rhsop` | $$[\sigma]$$                 | Material operand — used only when both LHS and RHS exist for the same integral form | — |
| `v`     | $$\mathbf{v}$$               | Advection velocity vector                            | —              |
| `dV`    | $$dV$$                       | Volume integral                                      | —              |
| `dS`    | $$dS$$                       | Surface (face) integral                              | —              |

---

## Key Distinctions

| Symbol pair | Rule |
|---|---|
| `B` / `BT` vs `dN` / `dNT` | `B`/`BT` reserved for **vector DOF** mechanics operators only. `dN`/`dNT` used for all **scalar DOF** flow and transport operators. |
| `op` vs `lhsop` / `rhsop` | Use plain `op` when only one side (LHS or RHS) exists for a given integral form. Use `lhsop`/`rhsop` as a **pair** only when both the matrix (LHS) and vector (RHS) versions of the same integral exist. |
| `Num` prefix vs no prefix | `Num` = numerical Gauss quadrature in parametric space (isoparametric elements). No prefix = exact analytical integration in global coordinates (straight-sided elements only). |
| `var` prefix | Indicates a **variable nodal prefactor** multiplying the integral — e.g. `Integral_var_NT_lhsop_N_dV` has a nodal scalar prefactor in addition to the element material operand. |

---

## LHS / RHS Tagging Rule

| Situation | LHS operator name | RHS operator name |
|---|---|---|
| Only LHS matrix exists | `..._op_...` | — |
| Only RHS vector exists | — | `..._op_...` |
| Both LHS matrix and RHS vector exist for the same integral | `..._lhsop_...` | `..._rhsop_...` |

---

## Operator Catalogue with Final Names

### Left-Hand Side (Matrix) Operators — `MathOperatorLHS`

| Operator | Integral | DOF type |
|---|---|---|
| `Integral_dNT_op_dN_dV` | $$\int (\nabla N)^T [\sigma] \nabla N \, dV$$ | Scalar |
| `Integral_NT_lhsop_N_dV` | $$\int N^T \sigma N \, dV$$ | Scalar |
| `Integral_var_NT_lhsop_N_dV` | $$\int N^T \sigma \lambda N \, dV$$ ($$\lambda$$ = nodal prefactor) | Scalar |
| `Integral_DNT_op_DN_NT_v_DN_dV` | $$\int \left[(\nabla N)^T [D] \nabla N + N^T \mathbf{v} \cdot \nabla N\right] dV$$ | Scalar |
| `NumIntegral_dNT_op_dN_dV` | $$\int (\nabla N)^T [\sigma] \nabla N \, dV$$ | Scalar |
| `NumIntegral_NT_lhsop_N_dV` | $$\int N^T \sigma N \, dV$$ | Scalar |
| `NumIntegral_var_NT_lhsop_N_dV` | $$\int N^T \sigma \lambda N \, dV$$ | Scalar |
| `NumIntegral_dNT_op_dN_NT_v_dN_dV` | $$\int \left[(\nabla N)^T [D] \nabla N + N^T \mathbf{v} \cdot \nabla N\right] dV$$ | Scalar |
| `NumIntegral_BT_D_B_dV` | $$\int B^T [D] B \, dV$$ | Vector (mechanics) |

### Right-Hand Side (Vector) Operators — `MathOperatorRHS`

| Operator | Integral | DOF type |
|---|---|---|
| `Integral_NT_op_N_dV` | $$\int N^T \sigma N \, dV \cdot \{1\}$$ (row-summed) | Scalar |
| `Integral_var_NT_rhsop_N_dV` | $$\int N^T \sigma \lambda N \, dV \cdot \{u\}$$ | Scalar |
| `Integral_dNT_op_dV` | $$\int (\nabla N)^T \mathbf{v} \, dV$$ | Scalar |
| `NumIntegral_NT_op_N_dV` | $$\int N^T \sigma N \, dV \cdot \{1\}$$ (row-summed) | Scalar |
| `NumIntegral_var_NT_rhsop_N_dV` | $$\int N^T \sigma \lambda N \, dV \cdot \{u\}$$ | Scalar |
| `NumIntegral_dNT_rhsop_dN_dV` | $$\int (\nabla N)^T [\sigma] \nabla N_k u_k \, dV$$ | Scalar |
| `NumIntegral_dNT_op_dV` | $$\int (\nabla N)^T \mathbf{v} \, dV$$ | Scalar |
| `NumIntegral_BT_D_op_dV` | $$\int B^T [D] \{\varepsilon\} \, dV$$ | Vector (mechanics) |

---

## Integration Method Reference

| Prefix | Integration method | Element requirement | Performance |
|---|---|---|---|
| `Num` | Gauss quadrature in parametric space | Isoparametric elements required | General — handles curved elements |
| _(none)_ | Exact analytical formula in global coordinates | Straight-sided elements only | Highest — no quadrature loop |


### Missing

| Matrix | Missing? | Priority |
| --- | --- | --- |
| Mass / capacitance | ✅ Have | — |
| Stiffness / conductance | ✅ Have | — |
| Elasticity stiffness | ✅ Have | — |
| Advection-dispersion | ✅ Have | — |
| Gradient-vector RHS | ✅ Have | — |
| Thermal/pressure expansion RHS | ✅ Have | — |
| Divergence / gradient coupling | ❌ Missing | High — mixed formulations |
| Robin / convection boundary | ❌ Missing | High — boundary conditions |
| Consistent source vector | ❌ Missing | High — body forces |
| Neumann boundary vector | ❌ Missing | High — flux BCs |
| Geometric stiffness | ❌ Missing | Medium — buckling |
| Coupling matrix (poromechanics) | ❌ Missing | Medium — Biot |
| Damping matrix | ✅ Covered by mass matrix | — |
