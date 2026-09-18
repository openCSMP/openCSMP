# CVFEM Hydrothermal System Modelling

This document describes the layout of `hydrothermal_system_modelling`, the CVFEM (Control Volume Finite Element Method) benchmarks that validate it, and examples to demonstrate how to run a model.

## Folder Structure

```
hydrothermal_system_modelling/
├── includes/
│   ├── cvfe_hydrothermal/              # Core CVFEM discretization for hydrothermal physics
│   └── hydrothermal_system_simulator/  # Simulator-level model components
│       ├── helper_classes/             # Model construction / time handling utilities
│       ├── magmatic_fluid_production/  # Magmatic fluid source models
│       ├── permeability_porosity/      # Permeability / porosity 
│       └── tracer_models/              # Tracer models
└── modules/                            # .cpp implementations, mirroring includes/
    ├── cvfe_hydrothermal/
    └── hydrothermal_system_simulator/
        ├── helper_classes/
        ├── magmatic_fluid_production/
        ├── permeability_porosity/
        ├── tracer_models/
```

**`includes/cvfe_hydrothermal` & `modules/cvfe_hydrothermal`**
The CVFEM transport scheme as described in [Weis et al., 2014](https://doi.org/10.1111/gfl.12080). 
The transport scheme got extended to also support unsaturated flow (i.e., an air phase is introduced). 

**`includes/hydrothermal_system_simulator` & `modules/hydrothermal_system_simulator`**
Simulator-level features built on top of the CVFEM layer, organized by process:
- `helper_classes/` — model building and time-step/integer conversion utilities.
- `magmatic_fluid_production/` — magma chamber and fluid-production models (`MagmaModel`, `ParmigianiVisitor`, `RatioVisitor*`).
- `permeability_porosity/` — different permeability models, porosity change, and a simplified permeability-porosity coupling.
- `tracer_models/` — different tracer models (`TracerModel`, `TracerVisitor`, `LithiumModel`, `LithiumVisitor`, `QuartzModel`).


## Benchmarks

```
tests/CVFEM/
├── CVFEM_1D_VVCase.cpp / .h
└── CVFEM_2D_VVCase.cpp / .h

data-tests_CVFEM/                  # Reference data compared against
├── 1D_Benchmarks/                 # tests 101–104 (each with a -wg variant; i.e., with-gravity)
└── 2D_Benchmarks/                 # tests 201a/b, 202a/b, 203a/b
    └── test-XXX/description.txt   # brief description of the test setup
```

The benchmarks reproduce the hydrothermal test cases presented in [Weis et al., 2014](https://doi.org/10.1111/gfl.12080). Please refer to the publication for the full model setup (geometry, boundary conditions, and physical parameters) of each case.

For every test case, three data snapshots are provided, each split into `boundaries`, `regions`, and `variables` `.dat` files plus a `.vset` file:
- **Initial** – the starting state of the model.
- **First-iteration** – the state after the first solver iteration.
- **Final** – the final solution at the end of the simulation.

Each case has a `test-XXX-configuration.txt` defining the setup (i.e., initial values and boundary conditions).

**How to use them:** build and run `CVFEM_tests_main` (requires `-DCSMP_WITH_CVFEM_HYDROTHERMAL:BOOL=ON` to build with CVFEM). The executable will prompt you interactively for:
1. **Benchmark data path** — the path to `data-tests_CVFEM/` (the folder containing `1D_Benchmarks/` and `2D_Benchmarks/`). This is located in `/tests/data-tests_CVFEM/` per default but could be moved to a different directory.
2. **Suite** — `0` for 1D tests, `1` for 2D tests.
3. **Test type** — `0` (QUICK: checks model initialization and the first-iteration output only) or `1` (EXTENDED: runs the full simulation and compares the Initial, First-iteration, and Final states; slower but gives a more detailed comparison against the benchmark data).
4. **Which tests to run** — a comma-separated list of test names (e.g., `test-201a, test-202b`), or `all` to run every test in the selected suite.
Each requested test is added to a `TestSuite`, run, and checked against the corresponding reference data in `1D_Benchmarks`/`2D_Benchmarks`; the total number of failures is printed at the end. A result of zero failures confirms the code reproduces the corresponding benchmark result from [Weis et al., 2014](https://doi.org/10.1111/gfl.12080). In case of differences between the calculated values and the benchmark data, the tests fail and differences are printed to the screen.
 
 
 
## Examples

```
examples/source/CVFEM_geothermal_cooling_magma_chamber_example_.cpp / .h
examples/source/CVFEM_fault_well_lithium_example.cpp / .h
examples/source/CVFEM_topo_magma_air_example.cpp / .h
```

`CVFEM_geothermal_cooling_magma_chamber_example` demonstrates how to configure and run a simulation with a magmatic intrusion that drives hydrothermal convection. It reads the following input files:
- `examples/example_inputs/CVFEM_examples/CVFEM_geothermal-configuration.txt` — configuration of initial values and boundary conditions.
- `examples/example_inputs/CVFEM_examples/2D_intrusion.dat / .asc` — unstructured grid.
- `examples/example_inputs/CVFEM_examples/2D_intrusion-regions.txt` — list of regions within the model.

`CVFEM_fault_well_lithium_example` demonstrates how to configure and run a simulation with a fault combined with injection and production wells. The fault is represented as a split-boundary object and a LithiumModel is applied to simulate metal transport via the brine. It reads the following input files:
- `examples/example_inputs/CVFEM_examples/CVFEM_fault_well_lithium_example-wells.txt` — configuration of wells.
- `examples/example_inputs/CVFEM_examples/CVFEM_fault_well_lithium_example_mesh.dat / .asc` — unstructured grid.
- `examples/example_inputs/CVFEM_examples/CVFEM_fault_well_lithium_example-regions.txt` — list of regions within the model.

`CVFEM_topo_magma_air_example` demonstrates how to configure and run a simulation with unsaturated flow using a mesh with topography and magmatic intrusion. It reads the following input files:
- `examples/example_inputs/CVFEM_examples/CVFEM_topo_magma_air_example_mesh.dat / .asc` — unstructured grid.
- `examples/example_inputs/CVFEM_examples/CVFEM_topo_magma_air_example-regions.txt` — list of regions within the model.

These examples are built and run as part of the example suite (`examples/source/ExamplesMain.cpp`), and serve as a starting point for setting up hydrothermal system simulations with this codebase. CVFEM examples can be found in the CVFEM Examples categorie when running `ExamplesMain`. Please copy the required input files listed above into your working directory. 
