# openCSMP

## Background

openCSMP (Complex Systems Modelling Platform) supports the coupled simulation of
thermal–hydrological–mechanical–chemical (THMC) processes in complex subsurface
systems on a range of length scales.

This GitHub repository contains the C++23 ANSI / ISO standard-compliant
application programming interface (API) of openCSMP with examples and applications
to Discrete Fracture & Matrix (DFM) modelling (`DFM_and_fault_modelling/`),
reservoir simulation in 2D (`FECFVM_simulator/`), and hydrothermal fluid flow
including magmatic intrusions, faults, and wells (`hydrothermal_system_modelling/`).

The core of this repository is the geologic Model, discretised with a collocated
finite-element (FEM) and finite-volume (FVM) mesh. Combining FEM with FVM enables
locally and globally conservative simulation.

This data model supports computations on regions, boundaries and internal
discontinuities. Material properties and dependent variables are stored using tree
structures representing the mesh. Computational domains can evolve in shape and
properties over time, and multi-domain simulations, where different physics are
modelled in different regions, are supported using multi-models.

In addition to solution methods for elliptic-parabolic and hyperbolic partial
differential equations, openCSMP provides libraries of:

- Polytype linear and quadratic finite elements, with exact integration as well as
  with quadrature rules
- Discrete fracture and matrix (DFM) and geologic fault modelling provisions
- Constitutive relationships for THMC processes
- Equations of state
- Property modelling and simulation analysis functionality
- Local (asynchronous) time-stepping methods and support for dynamic mesh
  adaptivity, facilitating goal-based simulation

The openCSMP library forms the basis of a wide variety of simulation programmes,
including discrete fracture and matrix modelling, upscaling and reservoir
simulation. 50+ examples illustrate these capabilities and the input and output
interfaces, and also provide a structured introduction to openCSMP.

---

## Quick Start

The seven steps below take you from an empty machine to running examples. The
commands are the same on Linux, Windows (through WSL) and macOS, except for the
installation of the prerequisites in step 2.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/openCSMP/openCSMP.git
   cd openCSMP
   ```

2. **Install the prerequisites** (details: [Prerequisites](#prerequisites)):

   *Ubuntu / Debian, also inside the Windows Subsystem for Linux (WSL):*
   ```bash
   sudo apt install build-essential cmake git pkg-config \
                    libopenmpi-dev openmpi-bin \
                    petsc-dev libjpeg-dev libcgns-dev
   ```
   openCSMP needs CMake 3.29 or later. Check with `cmake --version`; Ubuntu 24.04
   ships CMake 3.28, in which case a current version can be installed with
   `pip install cmake` (or from <https://cmake.org/download/>).

   *macOS (Homebrew):*
   ```bash
   brew install cmake pkg-config open-mpi petsc jpeg cgns eigen
   ```

   The header-only libraries plf_colony and RapidXML are downloaded automatically
   when you configure the build; so is Eigen if it is not installed.

3. **Configure with CMake** (default: PETSc solver and JPEG output enabled):
   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   ```
   This creates the directory `build/` and reports which libraries were found.

4. **Build the openCSMP library, the examples and the tests:**
   ```bash
   cmake --build build --parallel
   ```

5. **Set up the examples workspace:**

   The examples are run in a *workspace* outside of the repository, so that
   simulation results are neither committed by accident nor lost when `build/` is
   deleted. Create it with
   ```bash
   cmake --build build --target setup_examples
   ```
   This copies the example inputs to `~/openCSMP/examples/example_inputs/` and
   prints the location. Files that already exist there are never overwritten. See
   [Running the Examples](#running-the-examples) to use a different directory.

6. **Explore openCSMP via the examples:**
   ```bash
   build/openCSMP-examples
   ```
   The program reports where it reads inputs and writes results, and then shows a
   menu of examples:

   - Start with `Interfaces_Example`, after studying its `Specifications()` in the
     source (`.cpp`) file. It converts selected input files into the "csmp native
     format" that most other examples start from.
   - Then work through the "CSMP for beginners" examples in the order in which
     they are listed. This gives you a working knowledge of the essential
     functionality of the code.
   - Use `Experimental_Example` to write your own program, for instance by
     combining functionality from the other examples.

   Each example writes its results to
   `~/openCSMP/examples/example_outputs/<example name>/`; open the files there
   with ParaView.

7. **Run the tests:**
   ```bash
   ctest --test-dir build --output-on-failure
   ```
   The large suite of unit and integration tests also serves as a collection of
   code usage examples.

To work with the code and write your own programs you will want an integrated
development environment (IDE). CLion, for example, can open the project directly
from its `CMakeLists.txt`.

Much success! Please provide feedback so that we can improve the code and its
documentation.

---

## Running the Examples

A description of the examples themselves, of how to write your own and of what to
do when something does not work is in [`examples/README.md`](examples/README.md).

### The examples workspace

```
~/openCSMP/examples/                 the workspace (default location)
├── example_inputs/                  copied from examples/example_inputs/ by 'setup_examples'
│   ├── input_meshes/                mesh files that can be read by the openCSMP interfaces
│   ├── variables_and_configuration_files/
│   ├── bitmaps_and_other/
│   ├── CVFEM_examples/
│   └── csmp_native_format_models/   created when you run Interfaces_Example
└── example_outputs/
    └── <example name>/              created by each example; deleted and re-created
                                     when the example is run again
```

Each example copies the input files it needs from `example_inputs/` into its own
directory `example_outputs/<example name>/` and writes all of its results there.
To use your own inputs, place them into the corresponding subdirectory of
`example_inputs/`.

### Choosing a different workspace

`openCSMP-examples` determines the workspace in this order:

1. the command-line option `--workspace <dir>`,
2. the environment variable `CSMP_EXAMPLES_DIR`,
3. the current working directory, if it is a directory named `example_inputs`
   (its parent is then the workspace),
4. the location set when openCSMP was configured: `~/openCSMP/examples`, or the
   value of `-DCSMP_EXAMPLES_WORKSPACE=<dir>`.

The `setup_examples` target always uses the location of point 4. For example, to
keep everything in `~/csmp/sandbox`:

```bash
cmake -S . -B build -DCSMP_EXAMPLES_WORKSPACE=~/csmp/sandbox
cmake --build build --target setup_examples
build/openCSMP-examples
```

The workspace must be located outside of the repository; `setup_examples` refuses
a directory inside of it.

### Running from an IDE

Either add `--workspace <dir>` to the program arguments of the `openCSMP-examples`
target, or set its working directory to `<workspace>/example_inputs` (Xcode:
Product → Scheme → Edit Scheme → Options → custom working directory; CLion: Run →
Edit Configurations → Working directory).

### Setting up the workspace by hand

Without CMake, copy the directory itself (not only its contents):

```bash
mkdir -p ~/openCSMP/examples
cp -r examples/example_inputs ~/openCSMP/examples/
```

---

## Prerequisites

The CMake build locates these automatically using `find_package`, `pkg-config`
or Homebrew (macOS). Nothing needs to be placed next to the repository.

### Required for all builds

| Dependency     | Purpose                        | Ubuntu / Debian / WSL                         | macOS (Homebrew)          |
|----------------|--------------------------------|-----------------------------------------------|---------------------------|
| C++23 compiler | GCC ≥ 13.3, Clang, or Intel    | `sudo apt install build-essential`            | Xcode Command Line Tools  |
| CMake ≥ 3.29   | Build system generator         | `sudo apt install cmake` or `pip install cmake` | `brew install cmake`    |
| Git            | Downloading dependencies       | `sudo apt install git`                        | `brew install git`        |
| `pkg-config`   | Library discovery              | `sudo apt install pkg-config`                 | `brew install pkg-config` |
| MPI            | Parallel execution             | `sudo apt install libopenmpi-dev openmpi-bin` | `brew install open-mpi`   |

### Downloaded automatically at configure time

These header-only libraries are fetched by CMake into `build/_deps/`; an internet
connection is needed the first time you configure.

| Library | Purpose | Note |
|---|---|---|
| [plf::colony](https://github.com/mattreecebentley/plf_colony) | Fast, memory-efficient non-contiguous container (to become `std::hive` in C++26) | pinned to a fixed commit |
| [RapidXML](https://github.com/discord/rapidxml) | XML parsing | pinned to a fixed commit |
| [Eigen](https://gitlab.com/libeigen/eigen) 3.4.0 | Linear algebra | only downloaded if no installed Eigen is found (`brew install eigen`, `sudo apt install libeigen3-dev`); force the download with `-DCSMP_FETCH_EIGEN=ON` |

### Enabled by default, can be switched off

| Dependency | CMake option                   | Purpose               | Ubuntu / Debian / WSL          | macOS (Homebrew)     |
|------------|--------------------------------|-----------------------|--------------------------------|----------------------|
| PETSc      | `-DCSMP_WITH_PETSC_SOLVER=ON`  | Default linear solver | `sudo apt install petsc-dev`   | `brew install petsc` |
| libjpeg    | `-DCSMP_WITH_IMAGE_OUTPUT=ON`  | Image output support  | `sudo apt install libjpeg-dev` | `brew install jpeg`  |

### Further dependencies

| Dependency   | CMake option         | Purpose                                            | Notes |
|--------------|----------------------|----------------------------------------------------|-------|
| CGNS         | *(auto-detected)*    | CFD General Notation System I/O                    | `sudo apt install libcgns-dev` / `brew install cgns`. Without the library, openCSMP compiles against a bundled header, but programs that use the CGNS interface do not link. Treat it as required for the examples and tests. |
| IAPS H2O EOS | *(auto-detected)*    | H₂O equation of state                              | Built from the source in `thirdparty/iaps_h2o_eos/` |
| OpenNURBS    | `-DOPENNURBS_ENABLE=ON` | NURBS geometry kernel (McNeel & Associates)     | Uses the checkout given with `-DOPENNURBS_ROOT=<dir>`; otherwise cloned into `build/_deps/` (branch `OPENNURBS_GIT_BRANCH`, default `8.x`). Needs zlib. |
| SAMG         | `-DSAMG_ENABLE=ON`   | Fraunhofer algebraic multigrid solver (commercial) | See [SAMG Solver](#samg-solver-commercial) |

### Compilers

Avoid Visual Studio's native compiler; on Windows use WSL, or Intel's free oneAPI
plug-in (see `doc/CSMP_on_Windows2026.docx`). Please share your experience if you
build successfully with other compilers, for example NVIDIA's NVC++.

---

## How the Build Finds Its Libraries on Each Platform

`CMakeLists.txt` takes a different path depending on the platform and on whether
SAMG is enabled. The configure output and the summary at its end show what was
chosen.

| | Linux and WSL | macOS, Apple Silicon (default) | macOS with SAMG |
|---|---|---|---|
| Architecture | native | arm64 | x86_64 under Rosetta 2 (set automatically when `SAMG_ENABLE=ON`) |
| Package source | `apt` | Homebrew in `/opt/homebrew` | a second, x86_64 Homebrew in `/usr/local` |
| Eigen | installed copy, else downloaded | `brew --prefix eigen`, else downloaded | as on the left, using the `brew` of the Rosetta shell |
| MPI | `find_package(MPI)` | compilers from `brew --prefix open-mpi` | `/usr/local/opt/open-mpi` |
| PETSc | `pkg-config` | `pkg-config`, pointed at `brew --prefix petsc` | `pkg-config`, pointed at `/usr/local/opt/petsc` |
| libjpeg | `find_package(JPEG)` | `brew --prefix jpeg` | `/usr/local/opt/jpeg` |
| CGNS | `find_package`, `pkg-config`, then system search | as on the left, then `brew --prefix cgns` | as on the left |
| Linear solver | PETSc | PETSc | SAMG (PETSc can stay enabled) |

**Linux.** Everything comes from the package manager, apart from the downloaded
header-only libraries.

**Windows (WSL).** Inside WSL the build is the Linux build. Clone the repository
into the Linux file system (for example `~/openCSMP`), not under `/mnt/c/`, where
compilation is much slower.

**macOS, Apple Silicon.** The native arm64 build uses the Homebrew packages listed
above and the PETSc solver. SAMG is not available in this build because its
binaries are x86_64 only.

**macOS with SAMG.** See the next section.

---

## SAMG Solver (Commercial)

SAMG (algebraic multigrid) is a high-performance commercial solver from Fraunhofer
SCAPOS. It is **not required** for a standard build; openCSMP defaults to the
PETSc solver.

To enable SAMG you need:

1. a valid SAMG licence (contact [Fraunhofer SCAPOS](https://www.scapos.de)),
2. the prebuilt SAMG binaries in a directory of your choice, with the libraries
   in its subdirectory `linux/` and/or `macosx/`.

Then tell CMake where that directory is:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DSAMG_ENABLE=ON -DSAMG_ROOT=/path/to/samg
```

If `SAMG_ROOT` is not given, the former location `../external_support_libs/samg`
(next to the repository) is still used when it exists. If the libraries are not
found, CMake warns and continues without SAMG. Add `-DSAMG_LEGACY=ON` for the
original SAMG integration with multiple instances; the default is the modern,
single-instance integration.

### SAMG on macOS (Apple Silicon)

The SAMG binaries are x86_64, so the whole build, including every library
openCSMP links against, must be x86_64 and run under Rosetta 2. The extra
requirements are:

1. **Rosetta 2** (`softwareupdate --install-rosetta`).
2. **A second Homebrew for x86_64**, installed from within a Rosetta shell. It
   lives in `/usr/local`, next to the arm64 Homebrew in `/opt/homebrew`:
   ```bash
   arch -x86_64 /bin/zsh          # enter a Rosetta shell; leave it again with 'exit'
   /usr/local/bin/brew install cmake pkg-config open-mpi petsc jpeg cgns eigen
   ```
3. **The SAMG libraries** `libamg_mult.dylib` and `libiomp5.dylib` in
   `<SAMG_ROOT>/macosx/`.

Configure and build from within the Rosetta shell:

```bash
arch -x86_64 /bin/zsh
cmake -S . -B build-samg -DCMAKE_BUILD_TYPE=Release \
      -DSAMG_ENABLE=ON -DSAMG_ROOT=/path/to/samg
cmake --build build-samg --parallel
```

With `SAMG_ENABLE=ON` on an Apple Silicon Mac, CMake sets the architecture to
x86_64 itself and takes MPI, PETSc and libjpeg from `/usr/local/opt/`. Use a
separate build directory (here `build-samg`) so that the arm64 and the x86_64
build do not get mixed up.

### SAMG runtime environment

To run executables linked against SAMG, set:

```bash
export DYLD_LIBRARY_PATH=/path/to/samg/macosx   # macOS
export LD_LIBRARY_PATH=/path/to/samg/linux      # Linux
export SVD_LICENSE_FILE=/path/to/samg/license.dat
```

Add these to your login shell (`~/.bashrc` or `~/.zshrc`). For IDEs such as Xcode
or CLion, add them to the run environment of your project. `ctest` sets the two
library paths itself.

---

## More on Building

### Build options

| Option                    | Default               | Description                                              |
|---------------------------|-----------------------|----------------------------------------------------------|
| `CMAKE_BUILD_TYPE`        | *(none)*              | `Release`, `Debug`, or `RelWithDebInfo`                  |
| `CSMP_WITH_PETSC_SOLVER`  | `ON`                  | Enable the PETSc linear solver                           |
| `CSMP_WITH_IMAGE_OUTPUT`  | `ON`                  | Enable JPEG image output via libjpeg                     |
| `SAMG_ENABLE`             | `OFF`                 | Enable the commercial Fraunhofer SAMG solver             |
| `SAMG_ROOT`               | *(empty)*             | Directory with the SAMG binaries                         |
| `SAMG_LEGACY`             | `OFF`                 | Use the legacy SAMG integration (multiple instances)     |
| `OPENNURBS_ENABLE`        | `OFF`                 | Enable the OpenNURBS geometry kernel                     |
| `OPENNURBS_ROOT`          | *(empty)*             | Existing OpenNURBS checkout to use                       |
| `CSMP_FETCH_EIGEN`        | `OFF`                 | Download Eigen even if an installed copy is found        |
| `CSMP_EXAMPLES_WORKSPACE` | `~/openCSMP/examples` | Where `setup_examples` puts the examples workspace       |
| `BUILD_EXAMPLES`          | `ON`                  | Build the examples executable                            |
| `BUILD_TESTING`           | `ON`                  | Build the test executable                                |
| `BUILD_CVFEM_TESTS`       | `OFF`                 | Build the CVFEM verification & validation tests          |
| `BUILD_SHARED_LIBS`       | `OFF`                 | Build shared libraries (default is static)               |
| `CSMP_WITH_ASAN`          | `OFF`                 | Build with AddressSanitizer                              |
| `CSMP_WITH_LTO`           | `OFF`                 | Enable link-time optimisation                            |

Examples:

```bash
# Debug build in its own directory
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug

# Without PETSc and JPEG
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCSMP_WITH_PETSC_SOLVER=OFF -DCSMP_WITH_IMAGE_OUTPUT=OFF

# Ninja instead of Make (faster on many systems); use a new build directory
cmake -S . -B build-ninja -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-ninja
```

CMake is not a build system itself; it *generates* build files (Makefiles or Ninja
files) for your operating system and compiler from `CMakeLists.txt`. The workflow
is always: configure (`cmake -S . -B build`), build (`cmake --build build`), and
optionally install. Building in a separate directory keeps the source tree clean.

### Tests

```bash
ctest --test-dir build                         # run all tests
ctest --test-dir build --output-on-failure     # show output for failing tests
ctest --test-dir build -V                      # verbose output
ctest --test-dir build -R Suite_Fundamentals   # run a specific suite
```

The test suites run in strict dependency order:

```
Refactoring → Fundamentals → Interdependent1 → Interdependent2 → Interfaces → Composite
```

### Cleaning the build

```bash
cmake --build build --target clean   # remove the compiled files only
rm -rf build                         # full clean; configure again afterwards
```

Your examples workspace is not affected by either.

---

## Installation

Installing is optional; openCSMP can be used directly from the build directory.

### User-specific installation (recommended for development)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~/local/csmp
cmake --build build --parallel
cmake --install build
```

### System-wide installation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build --parallel
sudo cmake --install build
```

This installs the library to `<prefix>/lib`, the headers to
`<prefix>/include/openCSMP`, the executables `openCSMP-examples` and
`openCSMP-tests` to `<prefix>/bin`, and the example inputs to
`<prefix>/example_inputs`.

### Setting PATH after the installation

```bash
echo 'export PATH=~/local/csmp/bin:$PATH' >> ~/.zshrc    # zsh; use ~/.bashrc for bash
source ~/.zshrc
```

---

## Repository Structure

```
.
├── CMakeLists.txt                   build description for library, examples and tests
├── cmake                            CMake helper scripts
├── source
│   ├── includes                     headers, by module: analysis, discretisation, eos,
│   │                                integration, interfaces, math, model,
│   │                                property_manipulation, regular_grids, run_control,
│   │                                thmc, utilities
│   └── modules                      sources, same module structure
├── non-GIT-files                    variable-set headers
├── examples
│   ├── source                       the examples (ExamplesMain.cpp starts the menu)
│   ├── example_inputs               inputs copied into the examples workspace
│   └── data                         models: ANSYS, Eclipse, gOcad, Rhino, Triangle,
│                                    bitmaps, tutorials, fixtures
├── tests
│   ├── unit                         unit tests
│   ├── integration                  integration tests
│   ├── dataset                      input data for the tests
│   ├── CVFEM, data-tests_CVFEM      CVFEM verification & validation cases and data
│   └── data, doc
├── DFM_and_fault_modelling          discrete fracture & matrix and fault modelling
├── FECFVM_simulator                 2D reservoir simulator
├── hydrothermal_system_modelling    hydrothermal CVFEM modules and simulator
├── thirdparty
│   └── iaps_h2o_eos                 H₂O equation of state
├── doc                              user guide, interface guide, conventions, how-tos
├── latex
└── LICENSES
```

---

## Usage and Documentation

Learn how to use openCSMP with the examples executable (`examples/`), the models
in `examples/data`, and the configuration files. Modify and combine these to
create new functionality, and place new code into the example template
`Experimental_Example.cpp`.

Do **not** create builds of your own programs or simulation output inside the
openCSMP repository. Work in the examples workspace (see
[Running the Examples](#running-the-examples)) or in another directory outside of
the repository. This way you do not risk committing unwieldy files once you have
developer rights.

Documentation:

- `doc/CSMP_UserGuide.docx`: in-depth description of openCSMP's functionality
- `doc/CSMP_InterfaceGuide.docx` and `doc/CSMP_ModelBuildingAndInterfacesGuide.docx`:
  input and output interfaces and model building
- `doc/Open_CSMP_detailed_overview.pptx`: overview presentation
- `doc/howto/`: how-to guides, including one on ParaView
- Doxygen comments throughout the code: most functions have an abbreviated
  triple-slash (`///`) explanation in the header file, and classes are documented
  in the header and source files. Your IDE can show these together with command
  completion. The syntax is explained in
  `doc/CSMP_doxygen_doc_guidelines_basic.pdf` and its advanced version.

### More on the Complex Systems Modelling Platform

Multi-physics simulation of complex geologic systems needs to express dynamic
behaviours manifesting as plumes, convection cells and/or fracturing and faulting.
The "playing field" is the geomodel. Instabilities emerge across multiple length
scales and create unique complex patterns and statistical distributions. The
Complex Systems Modelling Platform (CSMP++) employs a space-time adaptive
collocated finite element – finite volume discretisation of a geomodel to simulate
subsurface behaviour, with a focus on geo-energy systems. This discretisation is
implemented with dynamic polymorphic data structures. The code uses policy-based
class design, with elements, nodes and other data structures implemented as class
templates.

---

## Support

Contact us via email (stephan.matthai@icloud.com), seek help via the GitHub
repository site (also consulting its issue tracker), or contact other members of
the openCSMP community.

## Roadmap

This is the first release, which has space to grow and to improve. The goal is
now to improve and extend the code base and the documentation step by step.

## Contributing

[Contributions](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project)
to openCSMP are most welcome. A strict code-review process is in place. To qualify
for review, source code must be standard-compliant, documented using Doxygen, and
accompanied by passing unit or integration tests as applicable.

---

## Licence

openCSMP is licensed under LGPL-3.0-only.
