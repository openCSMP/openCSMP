#ifndef INTERFACE_FRACTURE_VISITOR_H
#define INTERFACE_FRACTURE_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"


namespace csmp {

  template<uint32_t> class PropertyDatabase;
  class ScalarVariable;
  template<uint32_t> class Model;

  template<uint32_t dim>
  class InterFaceFractureVisitor : public Visitor<dim> {
  public:
    InterFaceFractureVisitor( const PropertyDatabase<dim>& p_ref,
                              const char* displacement,
                              const char* aperture,
                              INTERFACE_SIDE aperture_side);


    InterFaceFractureVisitor( const PropertyDatabase<dim>& p_ref,
                              const char* displacement,
                              const char* aperture,
                              INTERFACE_SIDE aperture_side,
                              const char* conductivity,
                              const char* viscosity,
                              double min_aperture = 0.0);



    virtual ~InterFaceFractureVisitor();

    virtual void Visit(InterFace<dim>* );

  private:
    InterFaceFractureVisitor();

  private:

    Index                displacement_key_;
    Index                aperture_key_;
    const INTERFACE_SIDE aperture_side_;

    const bool           calculate_conductivity_;
    Index                viscosity_key_;
    Index                conductivity_key_;
    double               min_aperture_;

  };



} // csmp

#endif
