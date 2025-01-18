#ifndef INTREPID_EXAMPLE_H
#define INTREPID_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<uint32_t> class Model;

class  Intrepid_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();

private:
  void CopyInputFiles(std::string& model_name, std::string& variable_file);
};

} // csmp

#endif // INTREPID_EXAMPLE_H
