#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include <iostream>
#include <stdexcept>
#include <exception>
#include <string>
#include <map>
#include <cassert>

using std::string;

namespace csmp {

  class CommandLineParser {
  public:
	CommandLineParser(int argc, const char** argv);

	template <typename T>
	T get(string shortOption, string longOption);
	template <typename T>
	T get(string shortOption, string longOption, const T defaultValue);
	template <typename T>
	T get(string option);
	template <typename T>
	T get(string option, const T defaultValue);

	bool operands() const;
	string getOperands() const;
	bool unused() const;
	string getUnused() const;
	
	void print(std::ostream& os) const;
	
  private:
	typedef std::map<string, string>::const_iterator optIterator;
	
	// disallowed
	CommandLineParser(const CommandLineParser&);
	CommandLineParser& operator=(const CommandLineParser&);
	
	bool containsOption(string opt) const;
	bool containsShortOption(string opt) const;
	
	bool isOption(string opt) const;
	bool isShortOption(string opt) const;
	bool isLongOption(string opt) const;
	
	bool isEndMarker(string str) const;

	void store(string opt, string arg);
	optIterator findArgument(string opt) const;
	optIterator findArgument(string shortOpt, string longOpt) const;

  template <typename T>
  T processArgument(optIterator it);
	template <typename T>
	T convertArgument(string arg);
	
	const char** argv_;
	const int argc_;
	string operands_;
	std::map<string, string> options_;
};

std::ostream& operator<< (std::ostream& os, const CommandLineParser& parser);
/**
 
@class CommandLineParser CommandLineParser "utilities/CommandLineParser.h"
@author A. Burri
@date 2003

@section motivation Motivation
 
Command line parser which allows for type-safe Unix-like command line
parsing. 

 
@section design Design Intent
 
The command line parser makes command line parsing safer and easier.
Beside alleviating the need to do error-prone string parsing it offers
sophisticated error handling capabilities through exception handling. 

Most Unix-like command line parsing feature are supported, like:
- long and/or short options (like "-d" and "--double")
- type-safe retrieval
- optional arguments
- listing of unused options
- listing of non-arguments 

For a more complete overview of the functionality see the unit-test
class of the CommandLineParser class. 

One specialty is that boolean options (i.e. command line flags) aren't
allowed to take arguments, since they wouldn't be used anyway. So you
are forced to use the end of options marker if you want to pass non-
option arguments after a flag. 

 
@section applicability Applicability
 
Parsing of command lines in the main function of a program. 

@section implementation Implementation
 
Stand-alone class. Member template methods used. The command line is
parsed in the constructor and the items can be retrieved using get
functions. 

 
@section examples Application Examples

@code
int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  
  // Retrieving arguments
  // ! Note that you have to specify the template parameter (return
  // value parametrization) !
  try {
  	double d = parser.get<double>("-d"); // mandatory argument
  	string s = parser.get<const char*>("-s", "Default");
  											// optional argument
  	int i = parser.get<int>("-i", "--int");
  											// choice between long and
  											// short option
  	bool b = parser.get<bool>("-b");
  } catch (std::runtime_error& e) {
  	cout << e.what() << endl;
    exit(1);
  }
  
  // Display if there are arguments which don't belong to an option
  if (parser.operands()) {
    cout << "Non-arguments: " << parser.getOperands() << endl;
  }
  
  // Display if there are unused arguments
  if (parser.unused()) {
  	cout << "Unused: " << parser.getUnused() << endl;
  }
  
  // Output command line to screen
  // (could be redirected to a file stream as well, as long as the file
  // stream is opened before)
  cout << parser;

}  
@endcode
 */
} // end namespace csp

#endif
