# Open-CSMP++ README

Some background. CSMP (Complex Systems Modeling Platform) is a C++17 ANSI / ISO standard-compliant application programmer interface (API).
It is a library of numerical methods for the combined simulation of thermal–hydrological–mechanical–chemical (THMC) processes, and has been designed for complex models of the subsurface, spanning a range of length scales.

This CSMP repository contains implementations of the finite-element- (FEM) and finite volume (FVM) methods and combinations thereof.
In addition to solution methods for elliptic-parabolic and hyberbolic partial differential equations, it provides implementations of:

CSMP's' unique data model supporting computations on regions, boundaries and internal discontinuities. Material properties and dependent variables  stored using a tree structures representing the mesh. Computational domains can evolve in shape and properties over time, and multi-domain simulations, where different physics are modelled in different regions, can be implemented.

Constitutive relationships for THMC processes
Equations of state
Property modeling and simulation analysis functionality
Local (asynchronous) time-stepping methods and support for dynamic mesh adaptivity, facilitating goal-based simulation.

CSMP forms the basis of a wide variety of simulation programs, including discrete fracture and matrix modelling, upscaling and reservoir simulation. Publications on these applications are furnished upon request.

## Quick Start

1. **Clone repository and associated third-party libraries Eigen, Catch2, and plf::colony:**
Get Eigen from https://eigen.tuxfamily.org, Catch2 (https://github.com/catchorg), and plf colony from https://github.com/mattreecebentley/plf_colony/blob/master/plf_colony.h
   ```bash
   git clone https://gitlab.com/csmp/open-csmp-2024.git
   cd open-csmp-2024
   cd thirdparty/
   git clone https://gitlab.com/libeigen/eigen.git
   git clone https://github.com/catchorg/Catch2.git  
   cd thirdparty/plf_colony
    ```
  place plf_colony.h in there or replace existing file with the latest version

2. **Create and enter build directory:**
   ```bash
   mkdir build && cd build
   ```

3. **Configure with CMake:**
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```

4. **Build the project:**
   ```bash
   make -j$(nproc)  # For Linux
   make -j$(sysctl -n hw.ncpu)  # For macOS
   ```

5. **Run tests (optional):**
   ```bash
   ctest
   ```

6. **Start with examples:**
   - Check the `examples/` directory for sample code
   - Copy `examples/example_inputs/` to a new location for your own experiments
   - Use `examples/Experimental_Example.cpp` as a template for new development

For more detailed instructions, see the [Building CSMP](#building-csmp-library-examples-and-tests-with-cmake) section below.

## Installation

CSMP can be installed to a custom location using CMake's installation feature. This is useful for system-wide or user-specific installations.

### User-specific Installation (Recommended for Development)

1. **Create an installation directory:**
   ```bash
   mkdir -p ~/local/csmp
   ```

2. **Configure with the installation prefix:**
   ```bash
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=~/local/csmp -DCMAKE_BUILD_TYPE=Release
   ```

3. **Build and install:**
   ```bash
   make -j$(nproc)  # For Linux
   make -j$(sysctl -n hw.ncpu)  # For macOS
   make install
   ```

### System-wide Installation

For system-wide installation (requires sudo privileges):

```bash
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)  # For Linux
make -j$(sysctl -n hw.ncpu)  # For macOS
sudo make install
```

### Installation Structure

The installation creates the following directory structure under the chosen prefix:

```
prefix/
├── bin/          # Executables
├── include/      # Header files
├── lib/          # Libraries
└── share/        # Documentation and other resources
```

### Post-Installation

To use the installed executables, add the installation's `bin` directory to your PATH:

```bash
# For user-specific installation
echo 'export PATH=~/local/csmp/bin:$PATH' >> ~/.bashrc  # For bash
echo 'export PATH=~/local/csmp/bin:$PATH' >> ~/.zshrc   # For zsh

# For system-wide installation
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bashrc    # For bash
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.zshrc     # For zsh
```

After adding to PATH, either restart your terminal or run:
```bash
source ~/.bashrc  # For bash
source ~/.zshrc   # For zsh
```

# Building CSMP (Library, Examples and Tests) with CMake

This guide explains how to build the CSMP project using CMake, a tool that helps manage the build process across different platforms. This is especially useful if you are new to building software on Linux or other Unix-like systems.

## What is CMake?

CMake is not a build system itself, but rather a tool that *generates* build system files (like Makefiles) for your specific operating system and compiler. This means you write a `CMakeLists.txt` file once, and CMake can generate build files for many different environments (like Make on Linux, Ninja, or even Xcode projects on macOS).

## The CMake Workflow

The typical CMake workflow involves three main steps:

1.  **Configure:** Run `cmake` in a build directory to process the `CMakeLists.txt` files and generate the build system (e.g., Makefiles).
2.  **Build:** Run the generated build system (e.g., `make`) to compile the source code and create executables and libraries. There are multiple options here.
3.  **Install:** Install the built software to a specified location.

## Prerequisites

Before you begin, make sure you have the following installed on your system:

*   **CMake:** You can download it from [cmake.org](https://cmake.org/) or install it using your system's package manager (e.g., `sudo apt-get install cmake` on Ubuntu, `brew install cmake` on macOS).
*   **A C++ Compiler:** Such as GCC or Clang. Install using your system's package manager (e.g., `sudo apt-get install build-essential` on Ubuntu).
*   **Make:** This is usually included with build tools (like `build-essential`).

*   **Eigen:**  download this linear algebra library from https://eigen.tuxfamily.org, 
*   **Catch2:** download this unit-testing framework from https://github.com/catchorg
*   **Colony:** download this container class from https://github.com/mattreecebentley/plf_colony/blob/master/plf_colony.h

Make sure also that you have access to the directories where CSMP++ will be installed. You may need administrator rights on your computer.

## Step-by-Step Build Guide

It is highly recommended to perform an **out-of-source build**. This means you build the project in a separate directory from your source code. This keeps your source directory clean and makes it easy to start over if something goes wrong.

1.  **Open a Terminal:** Navigate to the root directory of the project in your terminal.

2.  **Create a Build Directory:** Create a new directory where the build files will be generated. A common name is `build`.
    ```bash
    mkdir build
    ```

3.  **Navigate into the Build Directory:**
    ```bash
    cd build
    ```

4.  **Configure the Project:** Run `cmake` from inside the `build` directory, pointing it to the root of your source code (which is the directory above your current location, represented by `..`).

    ```bash
    cmake ..
    ```

    By default, this will configure a `Debug` build. You can specify the build type using the `-DCMAKE_BUILD_TYPE` option:

    *   For a release build (optimized):
        ```bash
        cmake .. -DCMAKE_BUILD_TYPE=Release
        ```
    *   For a debug build (with debugging symbols):
        ```bash
        cmake .. -DCMAKE_BUILD_TYPE=Debug
        ```

    You can also enable specific options defined in the `CMakeLists.txt` files using the `-D` flag. For example:

    *   To enable the SAMG solver (if you a license to this commercial linear algebra package and configured it with CMake):
        ```bash
        cmake .. -DCMAKE_BUILD_TYPE=Release -Dsamg=ON
        ```
    *   To enable the SKM unit tests:
        ```bash
        cmake .. -DBUILD_SKM_TESTS=ON
        ```
    *   You can combine options:
        ```bash
        cmake .. -DCMAKE_BUILD_TYPE=Release -Dsamg=ON -DBUCCESS_SKM_TESTS=ON
        ```

    After running the configure step, CMake will generate the build files in the `build` directory. You should see output indicating that the configuration was successful.

5.  **Build the Project:** Now that the build files are generated, you can use the `make` command (or the command for your chosen build system) to compile the code.

    ```bash
    make
    ```

    To speed up the build process by using multiple processor cores, you can use the `-j` flag followed by the number of cores you want to use. A common practice is to use all available cores:

    ```bash
    make -j$(nproc) # For Linux
    make -j$(sysctl -n hw.ncpu) # For macOS
    ```

    The build process will compile the source files and create the libraries and executables.

6.  **Run Tests (Optional):** If testing was enabled during configuration (`-DBUILD_TESTING=ON`), you can run the tests using CTest from the `build` directory.

    ```bash
    ctest
    ```

    To see more detailed output during testing, use the `-V` flag:

    ```bash
    ctest -V
    ```

7.  **Install the Project (Optional):** If an install target is defined in the `CMakeLists.txt` file, you can install the built project.

    ```bash
    make install
    ```

    You can specify an installation prefix during the CMake configure step using `-DCMAKE_INSTALL_PREFIX=/path/to/install`.
    One normally installs the executables, include files and library into your home directory (~user/local/bin, ~user/local/include, ~user/local/lib) or globally for all users on the machine (/usr/local/bin, /usr/local/include, /usr/local/lib).
    You must have access rights to the install directory. 
    One some systems, you need to use the sudo command for this.

## Cleaning the Build

If you need to clean the build directory (e.g., to start a fresh build), you can either:

1.  Run the clean target from the build directory:
    ```bash
    make clean
    ```
2.  Simply delete the `build` directory and create it again. This is often the safest way to ensure a completely clean build.
    ```bash
    cd ..
    rm -rf build
    mkdir build
    cd build
    ```

This covers the basic steps to configure and build the project using CMake. If you encounter specific errors, they often provide clues about missing dependencies or configuration issues. Much success!

## Search paths

To run CSMP++ from any directory on your computer you need to set paths to the installation directory. For the csmp executables and the (commercial) SAMG libraries you need
  ```bash
  export PATH="$HOME/local/csmp/bin:$PATH"
  export DYLD_LIBRARY_PATH=/usr/local/thirdparty/samg/linux
  export SVD_LICENSE_FILE=/usr/local/sofware_licensinfg/samg/license.dat
  echo "DYLD_LIBRARY_PATH: $DYLD_LIBRARY_PATH" && echo "SVD_LICENSE_FILE: $SVD_LICENSE_FILE" 
   ```

you can add these commands to your login shell so that you do not have to repeat them everytime you run CSMP++ in a new window.
For XCode, CLion, Visual Studio or other integrated development environments, you must add these paths to your project settings. 

## Adding your own new files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.com/csmp/open-csmp-2024.git
git branch -M main
git push -uf origin main
```

## Integrate with other tools

- [ ] [Set up project integrations](https://gitlab.com/csmp/open-csmp-2024/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing (SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Details

## Name
The Complex Systems Modelling Platform (Open-CSMP++)


## Description
Multi-physics simulation of complex geologic systems needs to express dynamic behaviours manifesting as plumes, convection cells and-or fracturing and faulting. The "playing field" is the geomodel. Instabilities emerge across multiple length scales and create unique complex patterns and statistical distributions. The Complex Systems Modelling Platform (CSMP++) employs a space-time adaptive collocated finite element – finite volume discretization of geomodel to simulate subsurface behaviour, with a focus on geo-energy systems. This discretisation is implemented with dynamic polymorphic data structures. The code uses policy-based class design, and elements, nodes and other data structures implemented as class templates. 

CSMP is a C++17 ANSI / ISO standard-compliant application programmer interface (API).
It is a library of numerical methods for the combined simulation of thermal–hydrological–mechanical–chemical (THMC) processes, and has been designed for complex models of the subsurface, spanning a range of length scales.

This CSMP repository contains implementations of the finite-element- (FEM) and finite volume (FVM) methods and combinations thereof.
In addition to solution methods for elliptic-parabolic and hyberbolic partial differential equations, it provides implementations of:

CSMP's' unique data model supporting computations on regions, boundaries and internal discontinuities. Material properties and dependent variables  stored using a tree structures representing the mesh. Computational domains can evolve in shape and properties over time, and multi-domain simulations, where different physics are modelled in different regions, can be implemented.

Constitutive relationships for THMC processes
Equations of state
Property modeling and simulation analysis functionality
Local (asynchronous) time-stepping methods and support for dynamic mesh adaptivity, facilitating goal-based simulation.

CSMP forms the basis of a wide variety of simulation programs, including discrete fracture and matrix modelling, upscaling and reservoir simulation. Publications on these applications are furnished upon request.

CSMP is documented using DOxgen, using this tool's specific syntax that is explained in doc/howto/CSMP_doxygen_doc_guidelines_basic.pdf and its advanced version.
Thus, most functions have an abbreviated triple slash (///) explanation in the header file, while class documentation is in the header and source files.

Additional documentation is in doc/CSMP_UserGuide.docx, the Interface_Guide and the Open_CSMP_detailed_overview.ppt presentation. Some documentation of the third-party libraries is provided in doc/thirdparty_libraries.

## Usage
Try to learn how to use CSMP with the examples executable (examples/), models (examples/data), and configurations. Modify and combine these to create new functionality. Place such new code into the example template "Experimental_Example.cpp." Do not create related builds and simulation output inside of the CSMP repository, but, instead, copy the directory examples/example_inputs/ and the files contained therein into a new location (for example C:/csmp/sandbox/). This way you do not run the risk of committing unwieldy files to the repository once you will get developer rights.

## Support
Contact me via email (stephan.matthai@icloud.com), seek help via the Gitlab repository site (also consulting its issue tracker), or contact other members of the CSMP community.

## Roadmap
The communities goal is to release Open-CSMP++ in 2025, following the successful dissolution of the CSMP Originator's Group of Universities developers agreement.

## Contributing
Of course we are open to- and most grateful for contributions to CSMP and we are practicing a strict code-review process for this purpose. To qualify for this review, the source code must be standard-compliant, documented using DOxygen, and accompanied by passing unit or integration tests, depending on what is applicable.

## Authors and acknowledgment
CSMP would not have become what it is today through the inspiration and hard work of many individuals.

Conceived in 1995, by the author at Stanford University, CSMP's' application programmer interface (API) follows an object-oriented design using C++ templates. In collaboration with S. G. Roberts (ANU, Australia) its prototype was finalised and applied to fluid flow in faults [36-38] and fractures [67]. 
From 1996-99 as a member of C. A. Heinrich's Fluids and Ore Deposits Group at the ETHZ, Switzerland, SKM extended CSP, applying it to fracture flow [35] and reactive transport simulation of the formation of the giant Mount Isa copper deposit, Australia [34]. These simulations support the hypothesis that the Cu-orebody formed by the mixing of a reduced fluid with an oxidised one transporting the copper.
From 2001-08 at Imperial College, London (UK), CSP adopted C++11, and inter-faced with the geomodelling software gOcad and the finite-element meshing tool ICEM Tetra (ICEM Technologies, now ANSYS) and the fracture-geometry generators FracMan and FRED (Golder Associates). To employ hybrid-element meshes, A. Me-zentsev implemented missing finite-element types in the library. Advised by C. Pain, Mezentsev and the author added higher-order accurate implicit immiscible two-phase flow and transport modules to CSP [32]. These were applied to discrete-fracture and matrix (DFM) models [5,33] and in water flooding simulations of naturally fractured and structurally complex reservoirs [30-33]. Palusny et al. [49] subsequently refac-tored CSP's collocated finite element – finite volume framework and applied it to fracture-propagation modelling [47-48]. Paluszny [46] further extended this code to 3D fragmentation, including contact detection, transfer of tractions, and a J-integral based formulation to evaluate crack-tip stresses, combined with dynamic remeshing. 
In the same decade, at ETHZ, T. Driesner developed an H2O-NaCl equation of state [13] and S. Geiger a thermohaline convection module [18-19], IMPES multi-phase fluid-flow module [20-21], and a black-oil model for fracture flow [17]. Cou-mou, Matthai, Geiger, and Driesner [6] parallelised CSP using MPI-CH [10] and ran simulations of thermohaline convection at mid-ocean ridges [7,9]. This simulation study is an excellent illustration how emergent behaviour explains the frequency and variation of seafloor hydrothermal venting as well as characteristic black- and white-smoker vent temperatures that naturally arise from the nonlinear pressure-temperature dependent density and viscosity of brine, obviating the need to invoke sudden perme-ability changes to explain episodicity of high-Raleigh number convection.
From 2009-14, at Montanuniversitaet Leoben, P. Lang, R. Manasipov, J. Mindel and SKM refactored CSP using policy-based class design, introducing the distributed variable storage system that is studied in this paper. Manasipov wrote a Petrel/Eclipse interface, complementing the unstructured grid gOcad - SKUA-ANSYS-CSP fracture upscaling [29] and reservoir simulation workflows [39-40]. Inspired by C. Cordes & W. Kinzelbach [11], S. Bazr-Afkan et al. [3] devised a new method to simulate satu-ration-jump discontinuities in fractured media, developing a finite element-centred FECFVM transport scheme, applying it to Gas-Oil Gravity Drainage (GOGD) of naturally fractured reservoirs [28]. The ETHZ and Berlin teams made major progress modelling high-enthalpy magmatic hydrothermal systems, e.g., [71-72], including the formation of ore shells in porphyry copper systems [59-60]. A new strategic alli-ance with Fiona Whitaker (Bristol, UK) was formed on reactive transport modelling and D. Kulik, A. Leal, and G. Kosakovski (now D. Miron) at the Paul-Scherer Insti-tute (PSI), Villigen, Switzerland, interfacing CSP (now CSMP++) with their open-source C++ code GEMS. First applications of GEMS-CSMP focused on carbonate dissolution and dolomitization [76-77], but also on other reactive transport processes [75].
Over the last decade, at Melbourne University, SKM refactored the core of CSMP to draw on C++17 language features such as rvalue (&&) references / move semantics, lambda functions, ranges, constexpr and more. Mesh management was improved to better support the modelling of discontinuities such as faults (Fig. 1) and mesh evo-lution with time. For crack propagation modelling [53-54] and the representation of saturation and capillary pressure discontinuities [68], an earlier embedded discontinui-ty discretisation [43] was improved. Based on the CSMP API, the Australian Carbon Geo-Sequestration compositional (CO2-H2O-NaCl) simulator (ACGSS) was devel-oped and benchmarked [62] and applied to understand plume spreading in highly heterogeneous porous media [4,61]. ACGSS is space-time adaptive, implementing Discrete Event Simulation (DES, [24]) and pre-emptive event processing (PEP, [44]), in CSMP's unstructured collocated finite element – finite volume discretisation [63-64].  Most recently, this framework has been applied successfully to numerical up-scaling of relative permeability and capillary pressure using periodic flux [27] and double-periodic flux and pressure boundary conditions [78]. In the FluidFlower benchmarking study [15], CSMP simulations ("Melbourne") were the only ones con-ducted with a hybrid method as opposed to a FV-only discretisation, delivering rea-sonable forecasts with only about one third of the model degrees of freedom used by the other participants.
In 2022, A. Tertois (AspenTech) developed a direct interface between SKUA (former-ly Gocad) and CSMP facilitating lossless transfer of geomodels.  Youssef et al. [79] demonstrate its application to multiphase fluid flow simulation through cross-bedded sandstone, cf., [69].
      In the meantime, the ETHZ - Melbourne University collaboration produced CSMP applications to thermal fluid convection in fault zones [50-51], the simulation of high-enthalpy geothermal systems [12, 58-60,73], and the propagation of magmat-ic-volatile driven fractures in the carapace of magmatic intrusions [54]. Liem et al. [25] used CSMP's discrete fracture and matrix modelling capabilities in an ensemble-Kalman filter-based (ESMDA) calibration of fracture aperture in a field-data based DFM model, using far-field stress-based apertures as a prior [40, submitted]. At Im-perial, Salimzadeh et al. [57] presented the first thermo-hydro-mechanical CSP model for deformable fractured geothermal systems. Salimzadeh and Nick [56] also consid-ered reactive transport processes in these. At the Danish Institute of Technology (DTU), Peters et al. [52] studied the performance of multilateral geothermal wells. The topic was revisited in [74].
 As CSP developers moved institutions, in 2010, the community spanned ETHZ, Imperial (A. Palusnzny), HWU (S. Geiger, R. Annewandtner), and TU Berlin (D. Coumou, P. Weiss). In 2009 the Montanuniversitaet Leoben joined (L. Mosser, P. Lang, J. Mindel, R. Manasipov, A Yapparova), then DTU, Denmark in 2012 (H.M. Nick and students), and Melbourne University in 2014 (Q. Shao, H. Agheshliu, L.K. Tran). In 2016, ENSG, France, joined the community (M. Raguenel, F. Bonneau), then UQ (L. Gross), Australia. In 2017, S. Salimzadeh joined CSIRO, Clayton, Aus-tralia, leading a team developing CSMP hydraulic fracture propagation models and a new Gmsh interface. To foster inter-institutional collaboration many workshops were organised, starting in 2004 at Seealpsee, Switzerland. At the 2010 Bad Gams work-shop, Austria, in search for a unique name, CSP was renamed CSMP.  J. Mindel introduced a github-based continuous integration system, replacing the previous SVN-based system, also adding routine unit and integration testing.
 
 ## Cited publications
 
 1.  Alexandrescu, A.: Modern C++ Design: Generic Programming and Design Patterns Applied. Addison-Wesley (2001).
2.  Arabas, S., Jarecka, D., Jaruga, A., Fijalkowski, M.: Formula translation in Blitz++, Numpy and modern Fortran: a case study of the language choice tradeoffs. Scientific Programming 22, pp. 201-22. DOI 10.3233/SPR-140379 (2014).
3.  Bazrafkan, S., Matthai, S.K., Mindel, J.E.: The finite-element-centered finite-volume discretization method (FECFVM) for multiphase transport in porous media with sharp material discontinuities. In: ECMOR XIV—14th European Conference on the Mathematics of Oil Recovery, vol. 2014, no. 1, pp. 1–22. European Association of Geoscientists & Engineers (2014).
4.  Boon, M., Matthäi, S.K., Shao, Q., Youssef, A.A., Mishra, A., Benson, S.M.: Aniso-tropic rate-dependent saturation functions for compositional simulation of sandstone composites. Journal of Petroleum Science and Engineering, 209, p.109934 (2022).
5.  Bourbiaux, B.: Fractured reservoir simulation: A challenging and rewarding issue. Oil & Gas Science and Technology—Revue de l'Institut Français du Pétrole, 65(2), pp. 227–238 (2010).
6.  Coumou, D., Matthäi, S., Geiger, S. and Driesner, T.: A parallel FE–FV scheme to solve fluid flow in complex geologic media. Computers & Geosciences, 34(12), pp.1697-1707 (2008).
7.  Coumou, D., Driesner, T. and Heinrich, C.A.: The structure and dynamics of mid-ocean ridge hydrothermal systems. Science, 321(5897), pp.1825-1828 (2008).
8.  Chauvin, B.P., Lovely, P.J., Stockmeyer, J.M., Plesch, A., Caumon, G., Shaw, J.H.: Val-idating novel boundary conditions for three-dimensional mechanics-based restora-tion: An extensional sandbox model example. AAPG Bulletin, 102(2), pp. 245–266 (2018).
9.  Coumou, D., Driesner, T., Heinrich, C.A.: The structure and dynamics of mid-ocean ridge hydrothermal systems. Science, 321(5897), pp. 1825–1828 (2008).
10.  Coumou, D., Matthäi, S., Geiger, S., Driesner, T.: A parallel FE–FV scheme to solve fluid flow in complex geologic media. Computers & Geosciences, 34(12), pp. 1697–1707 (2008).
11.  Cordes, C., Kinzelbach, W.: Continuous groundwater velocity fields and path lines in linear, bilinear, and trilinear finite elements. Water Resources Research, 28(11), pp. 2903–2911 (1992).
12.  Driesner, T., Weis, P., Scott, S.: A new generation of numerical simulation tools for studying the hydrology of geothermal systems to "supercritical" and magmatic condi-tions. In: World Geothermal Congress (2015).
13.  Driesner, T., Geiger, S., Heinrich, C.A., Matthäi, S.K.: Modeling multiphase flow of H2O–NaCl fluids by combining CSP5.0 with SoWat2.0. Geochimica et Cosmo-chimica Acta, 70(18/Supplement), p. A147 (2006).
14.  Durlofsky, L.J., Aziz, K.: Advanced techniques for reservoir simulation and modeling of nonconventional wells. Stanford University (US) (2004).
15.  Flemisch, B., Nordbotten, J.M., Fernø, M., Juanes, R., Both, J.W., Class, H., Delshad, M., Doster, F., Ennis-King, J., Franc, J., Geiger, S.: The FluidFlower Validation Benchmark Study for the Storage of CO₂. Transport in Porous Media, pp. 1–48 (2023).
16.  Gamma, E. Helm, R., Johnson, R., Vlissides, J.: Design Patterns: Elements of Reusable
Object-Oriented Software. Boston: Addison-Wesley (1995).
17.  Geiger, S., Matthäi, S., Niessner, J., Helmig, R.: Black-oil simulations for three-component, three-phase flow in fractured porous media. SPE Journal, 14(02), pp. 338–354 (2009).
18.  Geiger, S., Driesner, T., Heinrich, C.A., Matthäi, S.K.: Multiphase thermohaline con-vection in the earth's crust: I. A new finite element–finite volume solution technique combined with a new equation of state for NaCl–H2O. Transport in Porous Media, 63, pp. 399–434 (2006).
19.  Geiger, S., Driesner, T., Heinrich, C.A., Matthäi, S.K.: On the dynamics of NaCl–H₂O fluid convection in the Earth's crust. Journal of Geophysical Research: Solid Earth, 110(B7) (2005).
20.  Geiger, S., Roberts, S., Matthäi, S.K., Zoppou, C., Burri, A.: Combining finite element and finite volume methods for efficient multiphase flow simulations in highly heter-ogeneous and structurally complex geologic media. Geofluids 4(4), 284–299 (2004).
21.  Geiger, S., Roberts, S., Matth, S.K. and Zoppou, C.: Combining finite volume and finite element methods to simulate fluid flow in geologic media. Anziam Journal, 44, pp.C180-C201 (2002).
22.  Gries, S., Metsch, B., Terekhov, K.M., Tomin, P.: System-AMG for fully coupled reser-voir simulation with geomechanics. In: SPE Reservoir Simulation Conference, p. D021S011R003. SPE (2019).
23.  Juanes, R., Samper, J., Molinero, J.: A general and efficient formulation of fractures and boundary conditions in the finite element method. International Journal for Numerical Methods in Engineering 54(12), 1751–1774 (2002).
24.  Karimabadi, H., Driscoll, J., Omelchenko, Y.A., Omidi, N.: A new asynchronous methodology for modeling of physical systems: breaking the curse of Courant condi-tion. Journal of Computational Physics 205(2), 755–775 (2005).
25.  Liem, M., Matthai, S.K., Jenny, P.: Estimation of fracture aperture in naturally frac-tured reservoirs using an ensemble smoother with multiple data assimilation. In: ECMOR 2022, vol. 2022, pp. 1–18. European Association of Geoscientists & Engi-neers (2022).
26.  Lovely, P.J., Jayr, S.N., Medwedeff, D.A.: Practical and efficient three-dimensional structural restoration using an adaptation of the GeoChron model. AAPG Bulletin 102(10), 1985–2016 (2018).
27.  Matthai, S. K., & Tran, L. K. (2023). Numeric determination of relative permeability of heterogeneous porous media with capillary discontinuities. Advances in Water Re-sources, 175, 104430.
28.  Matthai, S., Bazrafkan, S.: Simulation of gas oil gravity drainage: comparison of the dual continuum with the discrete-fracture and matrix approach. In: Second EAGE Workshop on Naturally Fractured Reservoirs, pp. cp-371. European Association of Geoscientists & Engineers (2013).
29.  Matthai, S. K., & Nick, H. M. (2009). Upscaling two-phase flow in naturally fractured reservoirs. AAPG bulletin, 93(11), 1621-1632.
30.  Matthai, S.K., Mezentsev, A., Belayneh, M.: Finite element–node-centered finite-volume two-phase-flow experiments with fractured rock represented by unstructured hybrid-element meshes. SPE Reservoir Evaluation & Engineering 10(06), 740–756 (2007a).
31.  Matthai, S.K., Geiger, S., Roberts, S.G., Paluszny, A., Belayneh, M., Burri, A., Me-zentsev, A., Lu, H., Coumou, D., Driesner, T., Heinrich, C.A.: Numerical simulation of multi-phase fluid flow in structurally complex reservoirs. Geological Society, Lon-don, Special Publications 292(1), 405–429 (2007b).
32.  Matthai, S.K., Mezentsev, A., Belayneh, M.: Control-volume finite-element two-phase flow experiments with fractured rock represented by unstructured 3D hybrid meshes. In: SPE Reservoir Simulation Conference, pp. SPE-93341. SPE (2005).
33.  Matthai, S.K., Belayneh, M.: Fluid flow partitioning between fractures and a permea-ble rock matrix. Geophysical Research Letters 31(7) (2004).
34.  Matthai, S.K., Heinrich, C.A., Driesner, T.: Is the Mount Isa copper deposit the product of forced brine convection in the footwall of a major reverse fault? Geology 32(4), 357–360 (2004).
35.  Matthäi, S. K. (2003). Fluid flow and (reactive) transport in fractured and faulted rock. Journal of Geochemical Exploration, 78, 179-182.
36.  Matthäi, S.K., Geiger, S. and Roberts, S.G.: Complex Systems Platform: CSP3D3. 0: user's guide. ETH Zurich, https://doi.org/10.3929/ethz-a-004432279 (2001).
37.  Matthai, S.K., Aydin, A., Pollard, D.D., Roberts, S.G.: Simulation of transient well-test signatures for geologically realistic faults in sandstone reservoirs. SPE Journal 3(01), 62–76 (1998).
38.  Matthai, S.K., Roberts, S.G.: Transient versus continuous fluid flow in seismically ac-tive faults: an investigation by electric analogue and numerical modelling. In: Fluid Flow and Transport in Rocks: Mechanisms and Effects, pp. 263–295. Springer, Dor-drecht (1997).
39.  Matthai, S.K., Roberts, S.G.: The influence of fault permeability on single-phase fluid flow near fault-sand intersections: results from steady-state high-resolution models of pressure-driven fluid flow. AAPG Bulletin 80(11), 1763–1779 (1996).
40.  Milliotte, C., Jonoud, S., Wennberg, O.P., Matthäi, S.K., Jurkiw, A., Mosser, L.: Well-data-based discrete fracture and matrix modelling and flow-based upscaling of multi-layer carbonate reservoir horizons (2018).
41.  Milliotte, C., Matthai, S.: From seismic interpretation to reservoir model: an integrated study accounting for the structural complexity of the Vienna Basin using an unstruc-tured reservoir grid. First Break 32(5) (2014).
42.  Mindel, J.E., Alt-Epping, P., Les Landes, A.A., et al.: Benchmark study of simulators for thermo-hydraulic modelling of low enthalpy geothermal processes. Geothermics 96, Article 102130 (2021).
43.  Nick, H.M., Matthai, S.K.: A hybrid finite‐element finite‐volume method with em-bedded discontinuities for solute transport in heterogeneous media. Vadose Zone Journal 10(1), 299–312 (2011).
44.  Omelchenko, Y.A., Karimabadi, H.: A time-accurate explicit multi-scale technique for gas dynamics. Journal of Computational Physics 226(1), 282–300 (2007).
45.  Palsberg, J., Jay, C.B.: The essence of the visitor pattern. In: Proceedings. The Twenty-Second Annual International Computer Software and Applications Conference (Compsac'98), pp. 9–15. IEEE (1998).
46.  Paluszny, A. & Zimmerman, R. W. (2011). Numerical simulation of multiple 3D frac-ture propagation using arbitrary meshes. Computer Methods in Applied Mechanics and Engineering, 200(9-12), 953-966.
47.  Paluszny, A., Matthai, S.K.: Impact of fracture development on the effective permea-bility of porous rocks as determined by 2‐D discrete fracture growth modeling. Jour-nal of Geophysical Research: Solid Earth 115(B2) (2010).
48.  Paluszny, A., Matthai, S.K.: Numerical modeling of discrete multi-crack growth ap-plied to pattern formation in geological brittle media. International Journal of Solids and Structures 46(18–19), 3383–3397 (2009).
49.  Paluszny, A., Matthai, S.K., Hohmeyer, M.: Hybrid finite element–finite volume dis-cretization of complex geologic structures and a new simulation workflow demon-strated on fractured rocks. Geofluids 7(2), 186–208 (2007).
50.  Patterson, J.W., Driesner, T., Matthai, S., Tomlinson, R.: Heat and fluid transport in-duced by convective fluid circulation within a fracture or fault. Journal of Geophysi-cal Research: Solid Earth 123(4), 2658–2673 (2018a).
51.  Patterson, J.W., Driesner, T., Matthai, S.K.: Self‐organizing fluid convection patterns in an en echelon fault array. Geophysical Research Letters 45(10), 4799–4808 (2018b).
52.  Peters, E., Blöcher, G., Salimzadeh, S., Egberts, P.J., Cacace, M.: Modelling of multi-lateral well geometries for geothermal applications. Advances in Geosciences 45, 209–215 (2018).
53.  Pezzulli, E., Nejati, M., Salimzadeh, S., Matthäi, S.K., Driesner, T.: An enhanced J‐integral for hydraulic fracture mechanics. International Journal for Numerical and Analytical Methods in Geomechanics 46(11), 2163–2190 (2022).
54.  Pezzulli, E.: Simulating hydraulic fracture propagation in crustal processes. Doctoral Dissertation, ETH Zurich (2022).
55.  Ponting, D.K.: Corner point geometry in reservoir simulation. In: ECMOR I-1st Euro-pean Conference on the Mathematics of Oil Recovery, pp. cp-234. European Associa-tion of Geoscientists & Engineers (1989).
56.  Salimzadeh, S., Nick, H.M.: A coupled model for reactive flow through deformable fractures in enhanced geothermal systems. Geothermics 81, 88–100 (2019).
57.  Salimzadeh, S., Paluszny, A., Nick, H.M., Zimmerman, R.W.: A three-dimensional coupled thermo-hydro-mechanical model for deformable fractured geothermal sys-tems. Geothermics 71, 212–224 (2018).
58.  Scott, S.W., Driesner, T.: Permeability changes resulting from quartz precipitation and dissolution around upper crustal intrusions. Geofluids 2018(1), Article 6957306 (2018).
59.  Scott, S., Driesner, T., Weis, P.: The thermal structure and temporal evolution of high-enthalpy geothermal systems. Geothermics 62, 33–47 (2016).
60.  Scott, S., Driesner, T., Weis, P.: Geologic controls on supercritical geothermal re-sources above magmatic intrusions. Nature Communications 6(1), Article 7837 (2015).
61.  Shao, Q., Boon, M., Youssef, A., Kurtev, K., Benson, S.M., Matthai, S.K.: Modelling CO2 plume spreading in highly heterogeneous rocks with anisotropic, rate-dependent saturation functions: A field-data based numeric simulation study of Otway. Interna-tional Journal of Greenhouse Gas Control 119, Article 103699 (2022).
62.  Shao, Q., Matthai, S., Driesner, T., Gross, L.: Predicting plume spreading during CO2 geo-sequestration: benchmarking a new hybrid finite element–finite volume composi-tional simulator with asynchronous time marching. Computational Geosciences 25, pp. 299–323 (2021).
63.  Shao, Q., Matthai, S.: Numerical modelling of CO2 migration through faulted storage strata with a new asynchronous FE-FV compositional simulator. In: ECMOR XVII, vol. 2020, pp. 1–16. European Association of Geoscientists & Engineers (2020).
64.  Shao, Q., Matthai, S.K., Gross, L.: Efficient modelling of solute transport in heteroge-neous media with discrete event simulation. Journal of Computational Physics 384, pp. 134–150 (2019).
65.  Shao, Q., Matthai, S., Gross, L.: Efficient modelling of CO2 injection and plume spreading with discrete event simulation (DES). In: 14th Greenhouse Gas Control Technologies Conference Melbourne, pp. 21–26 (2018).
66.  Stüben, K., Ruge, J.W., Clees, T., Gries, S.: Algebraic multigrid: from academia to in-dustry. In: Scientific Computing and Algorithms in Industrial Simulations: Projects and Products of Fraunhofer SCAI, pp. 83–119 (2017).
67.  Taylor, W.L., Pollard, D.D., Aydin, A.: Fluid flow in discrete joint sets: field observa-tions and numerical simulations. Journal of Geophysical Research: Solid Earth 104(B12), 28983–29006 (1999).
68.  Tran, L.K., Kim, J.C., Matthai, S.K.: Simulation of two-phase flow in porous media with sharp material discontinuities. Advances in Water Resources 142, Article 103636 (2020).
69.  Tertois, A.: Cross-bedding, bedforms, and tetrahedral meshes. In: 85th EAGE Annual Conference & Exhibition (including the Workshop Programme), vol. 2024, pp. 1–5. European Association of Geoscientists & Engineers (2024).
70.  Vandevoorde, D., Josuttis, N.M.: C++ templates: the complete guide. Addison-Wesley Professional (2002).
71.  Weis, P., Driesner, T., Coumou, D., Geiger, S.: Hydrothermal, multiphase convection of H2O‐NaCl fluids from ambient to magmatic temperatures: A new numerical scheme and benchmarks for code comparison. Geofluids 14(3), 347–371 (2014).
72.  Weis, P., Driesner, T., Heinrich, C.A.: Porphyry-copper ore shells form at stable pres-sure-temperature fronts within dynamic fluid plumes. Science 338(6114), 1613–1616 (2012).
73.  Yapparova, A., Lamy-Chappuis, B., Scott, S.W., Gunnarsson, G., Driesner, T.: Cold water injection near the magmatic heat source can enhance production from high-enthalpy geothermal fields. Geothermics 112, Article 102744 (2023).
74.  Yapparova, A., Lamy-Chappuis, B., Scott, S.W., Driesner, T.: A Peaceman-type well model for the 3D Control Volume Finite Element Method and numerical simulations of supercritical geothermal resource utilization. Geothermics 105, Article 102516 (2022).
75.  Yapparova, A., Miron, G.D., Kulik, D.A., Kosakowski, G., Driesner, T.: An advanced reactive transport simulation scheme for hydrothermal systems modelling. Geother-mics 78, 138–153 (2019).
76.  Yapparova, A., Gabellone, T., Whitaker, F., Kulik, D.A., Matthäi, S.K.: Reactive transport modelling of dolomitisation using the new CSMP++ GEM coupled code: governing equations, solution method and benchmarking results. Transport in Po-rous Media 117(3), 385–413 (2017a).
77.  Yapparova, A., Gabellone, T., Whitaker, F., Kulik, D.A., Matthäi, S.K.: Reactive transport modelling of hydrothermal dolomitisation using the CSMP++ GEM coupled code: effects of temperature and geological heterogeneity. Chemical Geology 466, 562–574 (2017b).
78.  Youssef, A. A., Shao, Q., & Matthäi, S. K. (2024). Computing Relative Permeability and Capillary Pressure of Heterogeneous Rocks Using Realistic Boundary Condi-tions. Transport in Porous Media, 1-26.
79.  Youssef, A., Tertois, A., Shao, Q., Matthai, S.: Simulating unsteady CO2 flow through brine-saturated cross-bedded sandstones: towards relative permeability curves for unstable displacement. In: ECMOR 2022, earthdoc.
80.  Knuth, D.E., 2007. Computer programming as an art. In ACM Turing award lec-tures (p. 1974).

 

## License
Open-CSMP++ is going to be licensed under L-GPL.

## Project status
CSMP's development is ongoin at a number of international institutions.
This library is the stable branch of CSMP, supporting a range of applications developed on top of it such as the Australian Carbon Geo-Sequestration Simulator (ACGSS).
Another branch is a work in progress, refactoring the finite element and finite volume libraries for massively parallel computations and introducing Eigen's expression templates for matrix algebra into the code, replacing csmp::DenseMatrix and its associated functionality.
