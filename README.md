# openCSMP
## Background
openCSMP (Complex Systems Modeling Platform) supports  
the coupled simulation of thermal–hydrological–mechanical–chemical (THMC) processes, 
in subterranean complex earth's systems on a range of length scales. 

This GitHub repository contains the C++23 ANSI / ISO standard-compliant 
application programmer interface (API) of openCSMP with examples and applications
to Discrete Fracture & Matrix (DFM) modelling (`DFM_and_fault_modelling/`), reservoir simulation in 2D (`FECFVM_simulator/`),
and hydrothermal fluid flow including magmatic intrusions, faults, and wells (`hydrothermal_system_modelling/`). 

The core of this repository is the geologic Model discretised with a collocated finite-element(FEM) and
finite volume (FVM)  mesh enabling locally and globally conservative simulation combining FEM with FVM. 

This data model supports computations on regions, boundaries and
internal discontinuities. Material properties and dependent variables are stored
using tree structures representing the mesh. Computational domains can evolve in
shape and properties over time, and multi-domain simulations, where different
physics are modelled in different regions are supported using multi-models.
  
In addition to solution methods for elliptic-parabolic and hyperbolic partial differential equations, 
openCSMP provides libraries of:

- polytype linear and quadratic finite elements with exact integration as well as one with quadrature rules

- discrete fracture and matrix (DFM) and geologic fault modelling provisions

- Constitutive relationships for THMC processes

- Equations of state

- Property modelling and simulation analysis functionality

- Local (asynchronous) time-stepping methods and support for dynamic mesh
  adaptivity, facilitating goal-based simulation

The openCSMP library forms the basis of a wide variety of simulation programmes, including
discrete fracture and matrix modelling, upscaling and reservoir simulation. 50+ examples illustrate these capabilities,
input and output interfaces, and a also a structured introduction to openCSMP.

---

## Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/openCSMP/openCSMP.git
   cd openCSMP
   ```

2. **Install Prerequisites** (see [Prerequisites](#prerequisites) for full
   details):

   *WindowsSubsytemLinux/Ubuntu/Debian:*
   ```bash
   sudo apt install build-essential cmake git pkg-config \
                    libopenmpi-dev openmpi-bin \
                    eigen \
                    petsc-dev \
                    libjpeg-dev \
                    libcgns-dev \
                    nlohmann-json
   ```

You also need an integrated development environment (IDE) to work with the code and write your own programs.
Try Clion for this purpose as it can use the project CMake file for setup.

   *macOS (Homebrew):*
   On arm64 (Apple Silicon) build with PETSc or use Rosetta 2 to build with SAMG.
   
   ```bash
   brew install cmake open-mpi petsc jpeg cgns eigen
   ```
see further details below under 'Prerequisites', and place the required header only
libraries <plf_colony.h> (that will become <hive> in C++26), as well as those for 
'rapid-xml' into a directory where they can be readily accessed. 

IMPORTANT: for openCSMP with Fraunhofer's SAMG, install support libs and perform built from within
ROSETTA 2 and zsh:
   ```bash
   arch -x86_64 /bin/zsh
   ```
use 'exit' to leave this shell again. Rosetta 2 will be deprecated in 2028.


3. **Create and enter build Directory:**
   ```bash
   mkdir build && cd build
   ```

4. **Configure with CMake** (default: PETSc solver enabled, JPEG output enabled):
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCSMP_WITH_PETSC_SOLVER=ON \
            -DCSMP_WITH_IMAGE_OUTPUT=ON
   ```

5. **Build the openCSMP Library, Examples and Tests:**
   ```bash
   make -j$(nproc)              # Linux
   make -j$(sysctl -n hw.ncpu)  # macOS
   ```
use Ninja for much improved build speed.

   ```bash
   cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DSAMG_ENABLE=ON
   cmake --build .
   ```

6. **Explore openCSMP via Examples**

Copy input files from repository into a suitable working directory outside of the repository:
   ```bash
   cp -r openCSMP/examples/example_inputs/*  working_dir/
   ```
Set your project's working directory to this new example_inputs/ directory.

After studying its "Specifications() in the source (.cpp) file,
run Interfaces_Example to convert selected input files into "csmp native format".
Selecting the interface example from the menu that will be shown 
when you launch 'openCSMP-examples'.
Note that there is extensive DOxygen documentation for most CSMP classes in the code.
Your IDE can parse that and show it to you in conjunction with command completions.

There are different categories of examples to choose from.
Please start with the "CSMP for beginners" examples in the order in which they 
are listed. This will guarantee that you get a working knowledge of the essential functionality of the code.
Refer to the CSMP_UserGuide in openCSMP/doc/ to get more in-depth information on openCSMP's functionality.

Use Experimental_Example to write your own program; try this by combining functionality 
from the other examples.

Much success! - and please provide feedback so we can improve the code and its documentation.

7. **Run ctest:**

And examine the rather large suite of unit and integration tests which also serve as code
usage examples.
   ```bash
   ctest --output-on-failure
   ```

The following is a more detailed description of the openCSMP build process and the 
libraries that it builds on. Of course, like Eigen, these all have their own comprehensive 
documentation are will be known to AI tools for coding.

---

## 8. **More Details**

### What is CMake?

CMake is not a build system itself, but rather a tool that *generates* build
system files (such as Makefiles or Ninja build files) for your specific operating
system and compiler. You write a `CMakeLists.txt` file once, and CMake can
generate build files for many different environments.

### The CMake Workflow

1. **Configure:** Run `cmake` in a build directory to process `CMakeLists.txt`
   and generate the build system.
2. **Build:** Run the generated build system (e.g., `make`) to compile the source
   code.
3. **Install:** Optionally install the built software to a specified location.

### Step-by-Step Build Guide

It is strongly recommended to perform an **out-of-source build** — building in a
separate directory from the source keeps things clean and makes it easy to start
over.

1. **Open a terminal** and navigate to the root directory of the project.

2. **Create a build directory:**
   ```bash
   mkdir build && cd build
   ```

3. **Configure the project:**

   *Default build (PETSc solver + JPEG output enabled):*
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```

   *Explicit default options (same result):*
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCSMP_WITH_PETSC_SOLVER=ON \
            -DCSMP_WITH_IMAGE_OUTPUT=ON
   ```

   *Debug build:*
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Debug
   ```

   *With SAMG solver enabled (requires licence and prebuilt binaries — see
   below):*
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release -DSAMG_ENABLE=ON
   ```

   *Disable PETSc or JPEG if not needed:*
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCSMP_WITH_PETSC_SOLVER=OFF \
            -DCSMP_WITH_IMAGE_OUTPUT=OFF
   ```

   *Using Ninja instead of Make (faster on many systems):*
   ```bash
   cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
   cmake --build .
   ```

4. **Build the project:**
   ```bash
   make -j$(nproc)              # Linux — use all available cores
   make -j$(sysctl -n hw.ncpu)  # macOS
   ```

5. **Run tests (optional):**
   ```bash
   ctest                              # run all tests
   ctest --output-on-failure          # show output for failing tests
   ctest -V                           # verbose output
   ctest -R Suite_Fundamentals        # run a specific suite
   ```

   Test suites run in strict dependency order:
   ```
   Fundamentals → Interdependent1 → Interdependent2 → Interfaces → Composite
   ```

6. **Install (optional):**
   ```bash
   make install
   ```

### CMake Build Options Summary

| Option                   | Default  | Description                                        |
|--------------------------|----------|----------------------------------------------------|
| `CMAKE_BUILD_TYPE`       | *(none)* | `Release`, `Debug`, or `RelWithDebInfo`            |
| `CSMP_WITH_PETSC_SOLVER` | `ON`     | Enable PETSc linear solver                         |
| `CSMP_WITH_IMAGE_OUTPUT` | `ON`     | Enable JPEG image output via libjpeg               |
| `SAMG_ENABLE`            | `OFF`    | Enable commercial Fraunhofer SAMG solver           |
| `SAMG_LEGACY`            | `ON`     | Use legacy SAMG integration (when SAMG is enabled) |
| `BUILD_TESTING`          | `ON`     | Build the test executable                          |
| `BUILD_EXAMPLES`         | `ON`     | Build the examples executable                      |
| `BUILD_SHARED_LIBS`      | `OFF`    | Build shared libraries (default is static)         |

### Cleaning the Build

```bash
# Option 1 — clean targets only
make clean

# Option 2 — full clean (safest)
cd ..
rm -rf build
mkdir build && cd build
```

---

Before configuring the build, ensure the following are available on your system.
The CMake build system will locate them automatically using `find_package`,
`pkg-config`, or Homebrew (macOS).

### Required for All Builds

C++23 and the associated standard libraries.

| Dependency     | Purpose                        | Ubuntu/Debian                              | macOS (Homebrew)            |
|----------------|--------------------------------|--------------------------------------------|-----------------------------|
| CMake ≥ 3.29   | Build system generator         | `sudo apt install cmake`                   | `brew install cmake`        |
| C++23 compiler | GCC ≥ 13.3, Clang, or Intel    | `sudo apt install build-essential`         | Xcode Command Line Tools    |
| Git            | Source fetching (FetchContent) | `sudo apt install git`                     | `brew install git`          |
| `pkg-config`   | Library discovery              | `sudo apt install pkg-config`              | `brew install pkg-config`   |
| MPI            | Parallel execution             | `sudo apt install libopenmpi-dev openmpi-bin` | `brew install open-mpi`  |
| Eigen 3        | Linear algebra (header-only)   | fetched automatically on Linux             | `brew install eigen`        |

> **Note — Eigen on Linux:** Eigen is fetched automatically at configure time via
> CMake `FetchContent`. No manual installation is required on Linux. On macOS it
> is located via Homebrew.

Avoid Visual Studio's native compiler. Use Intel's free OneAPI plugin instead.
Share with us if you successfully exploare other options like, for example, using Nvidea's NVC++.

### Fetched Automatically at Configure Time

The following header-only libraries are downloaded automatically by CMake during
the configure step and require no manual installation:

| Library | Purpose |
|---|---|
| [plf::colony](https://github.com/mattreecebentley/plf_colony) | Fast, memory-efficient non-contiguous container |
| [RapidXML](https://github.com/discord/rapidxml) | XML parsing (header-only) |

### Default Optional Dependencies (Enabled by Default)

| Dependency | CMake Flag                    | Purpose              | Ubuntu/Debian                    | macOS (Homebrew)    |
|------------|-------------------------------|----------------------|----------------------------------|---------------------|
| PETSc      | `-DCSMP_WITH_PETSC=ON`        | Default linear solver | `sudo apt install petsc-dev`    | `brew install petsc` |
| libjpeg    | `-DCSMP_WITH_IMAGE_OUTPUT=ON` | Image output support  | `sudo apt install libjpeg-dev`  | `brew install jpeg`  |

### Further Optional Dependencies

| Dependency   | CMake Flag          | Purpose                                    | Notes |
|--------------|---------------------|--------------------------------------------|-------|
| CGNS         | *(auto-detected)*   | CFD General Notation System I/O            | `sudo apt install libcgns-dev` / `brew install cgns`. Falls back to bundled headers if not found. |
| IAPS H2O EOS | *(auto-detected)*   | H₂O equation of state                     | Built from source if present in `thirdparty/iaps_h2o_eos/` |
| SAMG         | `-DSAMG_ENABLE=ON`  | Fraunhofer algebraic multigrid solver (commercial) | See [SAMG Solver](#samg-solver-commercial) |

---

On OSX on the Mac, openCSMP can be built for Apple Silicon (Eigen & PETSC) or ROSETTA 2 needed for SAMG.
If you want both, you must install 2 version of homebrew. The X86_64 will live in /usr/local/bin and accessed with zsh.
Use the following instructions:

| lib name | install command | include path | source-dir | linklib name |
| :--- | :--- | :--- | :--- | :--- |
| **Eigen** | `brew install eigen` | `/opt/homebrew/opt/eigen/include/eigen3` | `/opt/homebrew/Cellar/eigen/<version>` | *(Header-only)* |
| **CGNS** | `brew install cgns` | `/opt/homebrew/opt/cgns/include` | `/opt/homebrew/Cellar/cgns/<version>` | `-lcgns` *(also requires `-lhdf5`)* |
| **JPEG** | `brew install jpeg` *(or `jpeg-turbo`)* | `/opt/homebrew/opt/jpeg/include` | `/opt/homebrew/Cellar/jpeg/<version>` | `-ljpeg` |
| **MPI** | `brew install open-mpi` *(or `mpich`)* | `/opt/homebrew/opt/open-mpi/include` | `/opt/homebrew/Cellar/open-mpi/<version>` | `-lmpi` |
| **open-mpi** | `brew install open-mpi` | `/opt/homebrew/opt/open-mpi/include` | `/opt/homebrew/Cellar/open-mpi/<version>` | `-lmpi` *(or `-lmpi_cxx`)* |
| **rapid-xml** | *Manual / `vcpkg install rapidxml`* | `./include/rapidxml` *(or `/opt/homebrew/include/rapidxml` if copied)* | Custom project vendor folder | *(Header-only)* |
| **json** | `brew install nlohmann-json` | `/opt/homebrew/opt/nlohmann-json/include` | `/opt/homebrew/Cellar/nlohmann-json/<version>` | *(Header-only)* |
| **open-nurbs** |  |  |  |  header & source from RhinocerosTM |
| :--- | :--- | :--- | :--- | :--- |

## Installation

openCSMP can be installed to a custom location using CMake's installation feature.

### User-Specific Installation (Recommended for Development)

1. **Create an installation directory:**
   ```bash
   mkdir -p ~/local/csmp
   ```

2. **Configure with the installation prefix:**
   ```bash
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=~/local/csmp \
            -DCMAKE_BUILD_TYPE=Release \
            -DCSMP_WITH_PETSC=ON \
            -DCSMP_WITH_IMAGE_OUTPUT=ON
   ```

3. **Build and install:**
   ```bash
   make -j$(nproc)              # Linux
   make -j$(sysctl -n hw.ncpu)  # macOS
   make install
   ```

### System-Wide Installation

For system-wide installation (requires `sudo` privileges):

```bash
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local \
         -DCMAKE_BUILD_TYPE=Release \
         -DCSMP_WITH_PETSC=ON \
         -DCSMP_WITH_IMAGE_OUTPUT=ON
make -j$(nproc)
sudo make install
```

### openCSMP Repository Structure

```
.
├── cmake
├── DFM_and_fault_modelling
│   ├── doc
│   │   ├── FFSA-Liem-matlab
│   │   └── Liem-ForwardSimilator
│   ├── main
│   ├── source
│   ├── test_data
│   │   ├── aperture_20_frac
│   │   ├── aperture20_frac (broken model)
│   │   ├── box_hrz_fault_model
│   │   ├── box_hrz_fault_model2
│   │   ├── fracs19
│   │   ├── Irregular2D-2-frac-model
│   │   ├── multi_dim_region_model
│   │   ├── test_3_frac
│   │   └── test_fractures_with_shared_neighbors
│   └── tests
│       └── aperture_20_frac
├── doc
│   ├── figures
│   ├── guiding_publications_and_manuals
│   ├── howto
│   ├── models
│   └── thirdparty_libraries
│       ├── JSON
│       └── SAMG
├── examples
│   ├── data
│   │   ├── ansys_models
│   │   ├── bitmaps
│   │   ├── Eclipse
│   │   ├── fixtures
│   │   ├── gOcad_models
│   │   ├── physical_variables
│   │   ├── Rhino
│   │   ├── triangle_models
│   │   └── tutorials
│   ├── example_inputs
│   │   ├── bitmaps_and_other
│   │   ├── input_meshes
│   │   └── variables_and_configuration_files
│   └── source
│       ├── deprecated
│       ├── interrelations
│       └── visitors
├── FECFVM_simulator
│   ├── cmake
│   ├── doc
│   ├── includes
│   └── modules
├── latex
├── non-GIT-files
├── source
│   ├── includes
│   │   ├── analysis
│   │   ├── discretisation
│   │   ├── eos
│   │   ├── integration
│   │   ├── interfaces
│   │   ├── math
│   │   ├── model
│   │   ├── property_manipulation
│   │   ├── regular_grids
│   │   ├── run_control
│   │   ├── thmc
│   │   └── utilities
│   └── modules
│       ├── analysis
│       ├── discretisation
│       ├── eos
│       ├── integration
│       ├── interfaces
│       ├── math
│       ├── model
│       ├── property_manipulation
│       ├── regular_grids
│       ├── run_control
│       ├── thmc
│       └── utilities
├── tests
│   ├── data
│   │   ├── fixtures
│   │   └── meshes
│   ├── dataset
│   ├── doc
│   ├── integration
│   └── unit
└── thirdparty
    ├── iaps_h2o_eos
```

### Post-Installation: Setting PATH

```bash
# User-specific installation
echo 'export PATH=~/local/csmp/bin:$PATH' >> ~/.bashrc   # bash
echo 'export PATH=~/local/csmp/bin:$PATH' >> ~/.zshrc    # zsh

# System-wide installation
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bashrc     # bash
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.zshrc      # zsh
```

After editing, reload your shell:
```bash
source ~/.bashrc   # bash
source ~/.zshrc    # zsh
```

---

## Building openCSMP (Library, Examples and Tests) with CMake



## SAMG Solver (Commercial)

SAMG (Algebraic Multigrid) is a high-performance commercial solver from
Fraunhofer SCAPOS. It is **not required** for a standard build — openCSMP
defaults to the PETSc solver.

To enable SAMG you need:

1. A valid SAMG licence (contact
   [Fraunhofer SCAPOS](https://www.scapos.de))
2. Prebuilt SAMG binaries placed in
   `../external_support_libs/samg/linux` (Linux) or
   `../external_support_libs/samg/macosx` (macOS)

Then configure with:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DSAMG_ENABLE=ON
```

### SAMG on macOS (Apple Silicon)

SAMG requires x86_64 binaries. On Apple Silicon Macs, you must build under
Rosetta 2:

```bash
# Enter a Rosetta shell first
arch -x86_64 /bin/zsh

# Install x86_64 Homebrew dependencies under /usr/local if not already present
# Then configure:
cmake .. \
  -DCMAKE_OSX_ARCHITECTURES=x86_64 \
  -DCMAKE_BUILD_TYPE=Release \
  -DSAMG_ENABLE=ON \
  -DMPI_CXX_COMPILER=/usr/local/opt/open-mpi/bin/mpicxx \
  -DMPI_C_COMPILER=/usr/local/opt/open-mpi/bin/mpicc \
  -DMPI_Fortran_COMPILER=/usr/local/opt/open-mpi/bin/mpifort
```

### SAMG Runtime Environment

To run executables linked against SAMG, the following environment variables must
be set:

```bash
export PATH="$HOME/local/csmp/bin:$PATH"
export DYLD_LIBRARY_PATH=/path/to/external_support_libs/samg/macosx   # macOS
export LD_LIBRARY_PATH=/path/to/external_support_libs/samg/linux      # Linux
export SVD_LICENSE_FILE=/path/to/samg/license.dat
```

Add these to your login shell (`~/.bashrc` or `~/.zshrc`) so they are set
automatically. For IDEs such as Xcode, CLion, or Visual Studio, add these paths
to your project's run environment settings.

---

## Adding Your Own New Files

- We welcome [contributions](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) to the project.

---

9. **More on the Complex Systems Modelling Platform (openCSMP)**

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

openCSMP is a C++23 ANSI / ISO standard-compliant application programmer interface
(API). It is a library of numerical methods for the combined simulation of
thermal–hydrological–mechanical–chemical (THMC) processes, designed for complex
models of the subsurface spanning a range of length scales.

openCSMP is documented using Doxygen, using this tool's specific syntax explained in
`doc/howto/CSMP_doxygen_doc_guidelines_basic.pdf` and its advanced version. Most
functions have an abbreviated triple-slash (`///`) explanation in the header file,
while class documentation is in the header and source files.

Additional documentation is in `doc/CSMP_UserGuide.docx`, the Interface Guide,
and the `Open_CSMP_detailed_overview.ppt` presentation. Some documentation of the
third-party libraries is provided in `doc/thirdparty_libraries`.

## Usage

Learn how to use openCSMP with the examples executable (`examples/`), models
(`examples/data`), and configurations. Modify and combine these to create new
functionality. Place new code into the example template
`Experimental_Example.cpp`. Do **not** create related builds and simulation output
inside the openCSMP repository; instead, copy the directory `examples/example_inputs/`
and the files contained therein to a new location (for example `~/csmp/sandbox/`).
This way you do not risk committing unwieldy files to the repository once you have
developer rights.

## Support

Contact via email (stephan.matthai@icloud.com), seek help via the GutHub
repository site (also consulting its issue tracker), or contact other members of
the openCSMP community.

## Roadmap

This is only the first release which will has space to grow and to improve.
The goal is now improve and extend the code-base and documentation step-by-step.

## Contributing

Contributions to openCSMP are most welcome. A strict code-review process is in place.
To qualify for review, source code must be standard-compliant, documented using
Doxygen, and accompanied by passing unit or integration tests as applicable.

---

## Licence

openCSMP is licensed under LGPL-3.0-only.



