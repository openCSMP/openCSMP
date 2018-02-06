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
  string type;
  double minval;
  double maxval;
  string placement;
  string usage;
  string explanation;
  string reference;
  string notes;
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


int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " VariableStructName variables.csv variables.h\n";
        return 1;
    }


    std::deque<variable> vars;
    try
    {
        io::CSVReader<11> csv(argv[2]);
        csv.read_header(0, "variable name", "notation", "units", "type", "minval", "maxval", "placement", "usage", "explanation", "reference", "notes");
        variable var;
        while (csv.read_row(var.name, var.notation, var.units, var.type,
                    var.minval, var.maxval, var.placement, var.usage, var.explanation,
                    var.reference, var.notes)) {
            vars.emplace_back(std::move(var));
        }
    }
    catch (io::error::base& e)
    {
        std::cerr << "Error in processing input file: " << e.what() << '\n';
        return 1;
    }
  
    std::sort(vars.begin(), vars.end(), [](auto& lhs, auto& rhs) { return lhs.explanation < rhs.explanation; });

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
    string explanation("");
    for (auto& var : vars) {
      if (var.explanation != explanation) {
        explanation = var.explanation;
        ofs << "  // " << explanation << '\n';
      }
        ofs << "  csmp::INDEX<" << var.type << ','
            << var.placement << "> key_" << var.notation << "; // " << var.name << "\n";
    }

    ofs << "\n  template<size_t dim>\n";
    ofs << "  explicit " << argv[1] << "( const PropertyDatabase<dim>& db )\n";
    bool first = true;
    for (auto& var : vars) {
        ofs << (first ? "    : " : "    , ")
            << "key_" << var.notation << "( INDEX<"
            << var.type << ',' << var.placement
            << ">( db.StorageKey(\"" << var.name << "\") ))\n";
        first = false;
    }
    ofs << "  {\n";
    for (auto& var : vars) {
        ofs << "    if ( key_" << var.notation << ".place != " << var.placement
            << " || key_" << var.notation << ".type != " << var.type
            << " )\n      throw csmp::Exception( FATAL_ERROR, \""
            << argv[1] << "::" << argv[1] << ":\",\n        "
            << "\"The '" << var.name << "' variable must be " << var.type
            << " and placed on " << var.placement << "\"  );\n";
    }
    ofs << "  }\n";
    ofs << "};\n";

    ofs <<"\n} } // end namespace csmp\n\n";

    ofs <<"#endif // " << header_guard << "\n";

    cout <<"\nCSMP header file '"<< argv[3] <<"' written successfully.\n";

    return 0;

} // end main
