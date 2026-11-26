# openCSMP Examples

The examples are the quickest way to learn openCSMP. They are compiled into one
program, `openCSMP-examples`, which shows a menu and runs the example you select.
Every example is a short, documented C++ class in `source/` that you can read, run
and modify.

This page explains what is in this directory and how the examples work. For
installing and building openCSMP see the
[Quick Start](../README.md#quick-start) in the top-level README.

## What is in this directory

| Directory         | Content |
|-------------------|---------|
| `source/`         | The examples, one class per `<Name>_Example.h/.cpp`, and `ExamplesMain.cpp`, which registers them in the menu |
| `example_inputs/` | The input files the examples read: meshes, variable and configuration files. They are copied into your workspace (see below); do not run the examples here. |
| `data/`           | Further models (ANSYS, Eclipse, gOcad, Rhino, Triangle), bitmaps, tutorial data and test fixtures |

## Running the examples

After building openCSMP, from the top-level `openCSMP/` directory:

```bash
cmake --build build --target setup_examples   # once: creates the workspace
build/openCSMP-examples                       # shows the menu of examples
```

**First run `Interfaces_Example`** (in the group of software interfaces). It converts
selected input files into the openCSMP native format, which most other examples
start from. If an example reports a missing model, this step was skipped.

Then work through the "CSMP for beginners" examples from top to bottom.

## The workspace: where inputs are read and results are written

The examples do not run inside the repository. They use a *workspace* directory,
by default `~/openCSMP/examples`, so that results are not committed by accident
and survive when the build directory is deleted.

```
~/openCSMP/examples/
├── example_inputs/                  a copy of examples/example_inputs/
│   ├── input_meshes/                meshes for the openCSMP interfaces
│   ├── variables_and_configuration_files/
│   ├── bitmaps_and_other/
│   ├── CVFEM_examples/
│   └── csmp_native_format_models/   created by Interfaces_Example
└── example_outputs/
    └── <example name>/              one directory per example
```

When you run an example it

1. creates the directory `example_outputs/<example name>/` (an existing directory
   of that name is deleted first, so copy results you want to keep),
2. copies the input files it needs from `example_inputs/` into that directory,
3. writes all of its results there.

**Viewing results.** Open the files in `example_outputs/<example name>/` with
[ParaView](https://www.paraview.org). `openCSMP-examples` prints the full path of
the workspace when it starts. A ParaView how-to is in `../doc/howto/`.

**Using your own inputs.** Put your files into the matching subdirectory of
`example_inputs/` in the workspace, not into the repository. `setup_examples`
never overwrites files that already exist in the workspace, so your changes are
safe when you run it again, for instance to pick up new inputs after updating
openCSMP.

**Using a different workspace.** The program looks for the workspace in this
order:

1. `build/openCSMP-examples --workspace <dir>`
2. the environment variable `CSMP_EXAMPLES_DIR`
3. the current working directory, if it is named `example_inputs`
4. the default, `~/openCSMP/examples`, or what was set at configure time with
   `-DCSMP_EXAMPLES_WORKSPACE=<dir>` (this is also where `setup_examples` copies
   the inputs)

**In an IDE.** Add `--workspace <dir>` to the program arguments of the
`openCSMP-examples` target, or set its working directory to
`<workspace>/example_inputs`.

## The groups of examples

The groups below follow the order in which the examples are registered in
`source/ExamplesMain.cpp`; the menu may word the headings slightly differently.

| Menu group | What you learn |
|------------|----------------|
| CSMP for beginners | C++ as used in openCSMP, building a model, physical variables, a steady-state heat conduction, pressure diffusion, mesh convergence. Start here, in the listed order. |
| C++ examples | Language techniques the code relies on: the curiously recurring template pattern, policy-based design, test-driven development |
| Applied math / numerical methods | Dense matrices, error metrics, linear solvers (including PETSc), time-stepping approaches |
| Software interfaces | Reading and converting models: ANSYS, Triangle, Eclipse, Rhino, SKUA and others. `Interfaces_Example` belongs here. |
| Software functionality | Variables, regions, visitors, averaging, permeability tensors, monitoring and statistics |
| Simulation of physical processes | Pressure diffusion, tracer advection, two-phase flow, elasticity and poro-elasticity, thermal convection, geothermal systems and more |
| CVFEM | Hydrothermal systems with faults, wells and magma bodies |
| Tutorials | Tutorials 1 to 4 |
| Write your own examples | `Experimental_Example`, an empty template for your own code |

## Reading and writing an example

Each example class has two functions worth reading first:

- `Specifications()` describes what the example does and which inputs it uses,
- `Run()` contains the code.

To try out your own ideas, edit `source/Experimental_Example.cpp`; it is already in
the menu. To add a new example of your own:

1. create `source/MyTopic_Example.h` and `.cpp`, using an existing example as a
   pattern,
2. include the header in `source/Examples.h`,
3. register the class in `source/ExamplesMain.cpp` with
   `examplesSuite.RegisterExample( new MyTopic_Example() );`,
4. configure and build again (`cmake -S . -B build`, then `cmake --build build`);
   the configure step is what finds new files in `source/`.

Put input files for your example into the workspace. Add them to
`example_inputs/` in the repository only if the example is meant to be shared.

## If something does not work

| Symptom | Cause and remedy |
|---------|------------------|
| "no example inputs found in the examples workspace" | The workspace does not exist yet. Run `cmake --build build --target setup_examples`, or give the right directory with `--workspace`. |
| An example cannot find a `.vset` or other model file | `Interfaces_Example` has not been run yet in this workspace. |
| "'input_meshes/' directory is missing" | The program was started with a working directory that is not an `example_inputs/` directory. Use `--workspace`, or correct the working directory in your IDE. |
| Results of an earlier run are gone | Each run of an example re-creates its output directory. Copy results you want to keep elsewhere. |

---

Based on the original README by Qi Shao (2023).
