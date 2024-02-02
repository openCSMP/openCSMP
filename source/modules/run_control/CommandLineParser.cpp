#include <cstdlib>
#include <cstring>
#include "CommandLineParser.h"

namespace csmp {

CommandLineParser::CommandLineParser(int argc, const char** argv) :
  argv_(argv),
  argc_(argc),
  operands_()
{ 
  int i = 1;
  while (i < argc && containsOption(argv[i])) {
    string str = string(argv[i]);
    ++i;
    if (isOption(str)) {
      if (!(i < argc) || containsOption(argv[i]) || isEndMarker(argv[i])) {
	    store(str, "");
      } else {
	    store(str, argv[i]);
	    ++i;
      }
    } else {
      assert(containsShortOption(str));
      string opt = str.substr(0,2);
      string arg = str.substr(2);
      store(opt, arg);
    }
  } // end while
  
  if (i < argc && isEndMarker(argv[i])) {
  	++i;
  }
  for (; i < argc; ++i) {
    operands_ += string(argv[i]);
  }
}

template <typename T>
T CommandLineParser::get(string shortOption, string longOption) {
  if (!isShortOption(shortOption) || !isLongOption(longOption)) {
    throw std::runtime_error("Option " + shortOption + " or " + longOption +
			     " is not a valid option");
  }
  
  CommandLineParser::optIterator it = findArgument(shortOption, longOption);
  
  if (it == options_.end()) {
    throw std::runtime_error("Option " + shortOption + " and " + longOption 
			     + " not found");
  }
  
  return processArgument<T>(it);
}

template <typename T>
T CommandLineParser::get(string shortOption, string longOption, const T defaultValue) {
  if (!isShortOption(shortOption) || !isLongOption(longOption)) {
    throw std::runtime_error("Option " + shortOption + " or " + longOption +
			     " is not a valid option");
  }

  CommandLineParser::optIterator it = findArgument(shortOption, longOption);
  
  if (it == options_.end()) {
    return defaultValue;
  }

  return processArgument<T>(it);
}

template <typename T>
T CommandLineParser::get(string option) {
  if (!isOption(option)) {
    throw std::runtime_error("Option " + option + " is not a valid option");
  }
  
  CommandLineParser::optIterator it = findArgument(option);

  if (it == options_.end()) {
    throw std::runtime_error("Option " + option + " not found");
  }

  return processArgument<T>(it);
}

template <typename T>
T CommandLineParser::get(string option, const T defaultValue) {
  if (!isOption(option)) {
    throw std::runtime_error("Option " + option + " is not a valid option");
  }
  
  CommandLineParser::optIterator it = findArgument(option);

  if (it == options_.end()) {
    return defaultValue;
  }

  return processArgument<T>(it);
}

bool CommandLineParser::operands() const {
  return !operands_.empty();
}

string CommandLineParser::getOperands() const {
  return operands_;
}

bool CommandLineParser::unused() const {
  return !options_.empty();
}

string CommandLineParser::getUnused() const {
  string result;
  CommandLineParser::optIterator it = options_.begin();
  if (it != options_.end()) {
  	result += it->first + " " + it->second;
  	++it;
  	while (it != options_.end()) {
  		result += " " + it->first + " " + it->second;
  		++it;
  	}
  }
  return result;
}

void CommandLineParser::print(std::ostream& os) const {
  // * Need to add stuff here
  if (os.fail()) {
    throw std::runtime_error("Error on output stream");
  }
  assert(argc_ > 0);
  os << argv_[0];
  for (int i = 1; i < argc_; ++i) {
    os << " " << argv_[i];
  }
}

bool CommandLineParser::containsOption(string opt) const {
  return (containsShortOption(opt) || isLongOption(opt));
}

bool CommandLineParser::containsShortOption(string opt) const {
  return (opt.find_first_not_of("-") == 1);
}

bool CommandLineParser::isOption(string opt) const {
  return (isShortOption(opt) || isLongOption(opt));
}

bool CommandLineParser::isShortOption(string opt) const {
  return (containsShortOption(opt) && opt.length() == 2);
}

bool CommandLineParser::isLongOption(string opt) const {
  return (opt.find_first_not_of("-") == 2);
}

bool CommandLineParser::isEndMarker(string str) const {
  return (str == string("--"));
}

void CommandLineParser::store(string opt, string arg) {
  if (options_.find(opt) != options_.end()) {
  	throw std::runtime_error("Option " + opt + " specified twice");
  }
  options_.insert(std::pair<string, string>(opt, arg));
}

CommandLineParser::optIterator CommandLineParser::findArgument(string shortOption, string longOption) const {
  CommandLineParser::optIterator itShort = findArgument(shortOption);
  CommandLineParser::optIterator itLong = findArgument(longOption);

  if (itShort != options_.end() && itLong != options_.end()) {
    throw std::runtime_error("Related short option (" + shortOption +
			     ") and long option (" + longOption + 
			     ") simultaneously present");
  }
  
  return (itShort == options_.end()) ? itLong : itShort;
}

CommandLineParser::optIterator CommandLineParser::findArgument(string opt) const {
  return options_.find(opt);
}

template <typename T>
T CommandLineParser::processArgument(CommandLineParser::optIterator it) {
  T result = convertArgument<T>(it->second);
  options_.erase(it->first);
  return result;
}

// Specialisations
template <>
int CommandLineParser::convertArgument(string arg) {
  return std::atoi(arg.c_str());
}

template <>
double CommandLineParser::convertArgument(string arg) {
  return std::atof(arg.c_str());
}

template <>
const char* CommandLineParser::convertArgument(string arg) {
  char* result = new char[arg.size()+1]; // null terminator of string
  std::strcpy(result, arg.c_str());
  return result;
}

template <>
bool CommandLineParser::convertArgument(string arg) {
  if (!arg.empty()) {
    throw std::runtime_error("Boolean option with argument");
  }
  return true;
}

// int specialisation
template int CommandLineParser::get(string shortOption, string longOption);
template int CommandLineParser::get(string shortOption, string longOption, const int defaultValue);
template int CommandLineParser::get(string option);
template int CommandLineParser::get(string option, const int defaultValue);

// double specialisation
template double CommandLineParser::get(string shortOption, string longOption);
template double CommandLineParser::get(string shortOption, string longOption, const double defaultValue);
template double CommandLineParser::get(string option);
template double CommandLineParser::get(string option, const double defaultValue);

// bool specialisation
template <>
bool CommandLineParser::get<bool>(string shortOption, string longOption) {
	return get<bool>(shortOption, longOption, false);
}	
template bool CommandLineParser::get(string shortOption, string longOption, const bool defaultValue);
template <>
bool CommandLineParser::get<bool>(string option) {
	return get<bool>(option, false);
}
template bool CommandLineParser::get(string option, const bool defaultValue);

// string specialisation
template const char* CommandLineParser::get(string shortOption, string longOption);
template const char* CommandLineParser::get(string shortOption, string longOption, const char* defaultValue);
template const char* CommandLineParser::get(string option);
template const char* CommandLineParser::get(string option, const char* defaultValue);

// Output function
std::ostream& operator<< (std::ostream& os, const CommandLineParser& parser) {
  parser.print(os);
  return os;
}

} // end namespace csmp
