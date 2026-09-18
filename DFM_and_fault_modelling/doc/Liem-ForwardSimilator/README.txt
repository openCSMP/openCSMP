Forward simulators developed on the basis of OpenCSMP.
- Load CSMP native model
- Read aperture values from file
- Calculate steady-state pressure and velocity fields
- Advect tracer along the steady-state velocity field
- Write measurement data to file

run_well_fracture.cpp: Run the forward simulator of the history matching scenario, where fluid is injected into the domain through the well fracture.

run_left_right.cpp: Run the forward simulator of the forecast scenario, where fluid enters the domain from the left and exits at the right.