#include <iostream>
#include <fstream>
#include <cassert>
#include <deque>
#include <string>

/* IMPORTANT NOTE:
 *
 * This should not depend on any CSMP library for build reasons.
 * Header files are generally okay.
 */

// The CSV processor is distributed under the BSD 3-clause licence.
// This means that documentation will need to be updated if it is
// ever used as part of the distribution.
//
// This is a build tool, and as such is not distributed.
//
#include "csv.h"


using namespace std;


struct variable {
    string name;
    string notation;
    string units;
    size_t dimension;
    double minval;
    double maxval;
    string placement;
};


string headerGuard(const char* fname)
{
    string guard(fname);
    for (size_t i = 0; i < guard.size(); ++i) {
      guard[i] = toupper(guard[i]);
        if (guard[i] == '.') {
            guard[i] = '_';
        }
    }
    return string("CSMP_") + guard;
}


const char* varType(size_t dimension)
{
    if (dimension == 1) {
        return "SCALAR";
    }
    if (dimension == 2) {
        return "VECTOR";
    }
    if (dimension == 3) {
        return "TENSOR";
    }
    return "ARRAY";
}


int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " VariableStructName variables.csv variables.h\n";
        return 1;
    }


    std::deque<variable> vars;
    try
    {
        io::CSVReader<7> csv(argv[2]);
        csv.read_header(0, "name", "notation", "units", "dimension", "minval", "maxval", "placement");
        variable var;
        while (csv.read_row(var.name, var.notation, var.units, var.dimension,
                    var.minval, var.maxval, var.placement)) {
            vars.emplace_back(std::move(var));
        }
    }
    catch (std::exception e)
    {
        std::cerr << "Error in processing input file: " << e.what() << '\n';
        return 1;
    }

    ofstream ofs(argv[3]);

    string header_guard = headerGuard(argv[3]);

    ofs <<"#ifndef " << header_guard << '\n';
    ofs <<"#define " << header_guard << "\n\n";
    ofs <<"/**\n";
    ofs <<"@file " << argv[3] << '\n';
    ofs <<"Automatically generated from " << argv[2] << '\n';
    ofs <<"DO NOT EDIT!\n";
    ofs <<"*/\n\n";
    ofs <<"#include \"Index.h\"\n";
    ofs <<"#include \"Exception.h\"\n";
    ofs <<"#include \"PropertyDatabase.h\"\n\n";
    ofs <<"namespace csmp { namespace variables {\n\n";

    ofs << "struct " << argv[1] << " {\n";
    for (auto& var : vars) {
        ofs << "  csmp::INDEX<" << varType(var.dimension) << ','
            << var.placement << "> key_" << var.notation << "; // " << var.name << "\n";
    }

    ofs << "\n  template<size_t dim>\n";
    ofs << "  explicit " << argv[1] << "( const PropertyDatabase<dim>& db )\n";
    bool first = true;
    for (auto& var : vars) {
        ofs << (first ? "    : " : "    , ")
            << "key_" << var.notation << "( INDEX<"
            << varType(var.dimension) << ',' << var.placement
            << ">( db.StorageKey(\"" << var.name << "\") ))\n";
        first = false;
    }
    ofs << "  {\n";
    for (auto& var : vars) {
        ofs << "    if ( key_" << var.notation << ".place != " << var.placement
            << " || key_" << var.notation << ".type != " << varType(var.dimension)
            << " )\n      throw csmp::Exception( FATAL_ERROR, \""
            << argv[3] << "::" << argv[3] << ":\",\n        "
            << "\"The '" << var.name << "' variable must be " << varType(var.dimension)
            << " and placed on " << var.placement << "\"  );\n";
    }
    ofs << "  }\n";
    ofs << "};\n";

    ofs <<"\n} } // end namespace csmp\n\n";

    ofs <<"#endif // " << header_guard << "\n";

    cout <<"\nCSMP header file '"<< argv[3] <<"' written successfully.\n";

    return 0;

} // end main
