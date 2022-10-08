
To run the examples, please copy the 'example_inputs/' directory to a location outside of the open-csmp source diectory, and set it as your working directory.

This input directory contains all default input files required by the examples, in its 3 sub directories:

	- 'input_meshes/' contains mesh files that can be read by Open-CSMP interfaces;

	- 'csmp_native_format_models/' contains Open-CSMP native format binary files, that most examples start with;

	- 'variables_and_configuration_files/' contains the variables and configuration files.


Each example will 

	1} create an example directory from its name (i.e., 'example_outputs/example_name/');

	2) find its default input files from 'example_inputs/' directory, and copy them into the example directory (to use your own inputs, please copy the input files into corresponding sub directories of 'example_inputs/').

	3) produce all outputs in the created example directory (note: if an example directory already exists, it will be overwritten). 
