// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "LimitVisitor.h"
#include "Exception.h"
#include "ModelTime.h"

using namespace std;

namespace csmp {

  /** custom constructor */
  template<uint32_t dim>
  LimitVisitor<dim>::LimitVisitor(Model<dim>& model,
	  const char* var_name, double min_value, double max_value)
	  : min_val(min_value), max_val(max_value)
  { 

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

	var_key_ = model.Database().StorageKey(var_name);
    T_key_ = model.Database().StorageKey("temperature");

	if (var_key_.type != SCALAR || var_key_.place != NODE)
      throw csmp::Exception( ERROR, "LimitVisitor::(constructor)",
	  var_name, " must be a nodal scalar property.");

  }

  /** deconstructor */
  template<uint32_t dim>
  LimitVisitor<dim>::~LimitVisitor() 
  {}

  /** visit function for Node */
  template<uint32_t dim>
  void LimitVisitor<dim>::Visit(Node<dim>* n) 
  { 

	  n->Read(var_key_, var);

      n->Read(T_key_, T);
      if(var()<min_val || var()>max_val)
      {
          double model_time = ModelTime::Instance().modelTime;

          // open file in read and write mode
          std::fstream file("PressureLimiter.log", std::ios::in | std::ios::out);

          // move the write pointer to the end of the file to append new entries
          file.seekp(0, std::ios::end);

          if (dim == 2U)
              file << "\t" << model_time << "\t" << n->x() << "\t" << n->y() << "\t" << "0." << "\t" <<
                  var() << "\t" << T() << std::endl;
          else if (dim == 3U)
              file << "\t" << model_time << "\t" << n->x() << "\t" << n->y() << "\t" << n->z() << "\t"
                   << var() << "\t" << T() << std::endl;
          file.close();

          cerr << n->Coordinate() << endl;
          cerr << "fluid pressure = " << var() << endl;
          cerr << "fluid temperature = " << T() << endl;
      } // end write to LimiterLogFile

	  var() = std::max(var(), min_val);
	  var() = std::min(var(), max_val);
	  n->Store(var_key_, var);

  }

  /** visit function for Region */
  template<uint32_t dim>
  void LimitVisitor<dim>::Visit(Model<dim>* n)
  {
      // no calcution for the model
  }

  /** visit function for Region */
  template<uint32_t dim>
  void LimitVisitor<dim>::Visit(Region<dim>* n)
  {
	  // no calcution for the region
  }

  template class LimitVisitor<1U>;
  template class LimitVisitor<2U>;
  template class LimitVisitor<3U>;

} // csmp
