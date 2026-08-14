*** Files in this directory ***

The files contained in this directory are SKM's current unit_testing draft protocol, and  variable sets (header files) generated with source/variable_set_preprocessor.cpp (main file) from .csv (comma-delimited) text files output from Excel containing the specifications of the variables used in specific programs (name indicates which).

Once compiled, the variable_set_preprocessor.exe becomes a stand-alone command-line program that converts the csv files into CSMP header files.  Run -h options to see how the program is used:

prompt: variable_set_preprocessor VariableStructName variables.csv variables.h

The minimum viable specifications that must be included into the input file (csv), include variable name, notation, units, type (scalar..array), meaningful min-max values, and placement (Element, Node..). Beyond these extra categories will be read if they are indicated by column headers.

