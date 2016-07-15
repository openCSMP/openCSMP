#include "FractureCapillaryPressureVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"

namespace csmp{


  template<size_t dim>
  FractureCapillaryPressureVisitor<dim>::FractureCapillaryPressureVisitor( Model<dim>& model, const char* fracApTag,
                                                                         const char* fracPcTag, double pd, double ift, double ca )
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ), fracApKey_( model.Database().StorageKey(fracApTag) ),
      fracPcKey_( model.Database().StorageKey(fracPcTag) ), Pd_matrix_(pd), ift_(ift), ca_(ca)
    {
    }


  template<size_t dim>
  void FractureCapillaryPressureVisitor<dim>::Visit( Element<dim>* element )
    {
      const double fracAp( element->Read( fracApKey_ ) );
      fracPc_ = std::min ( 2*cos(ca_)*ift_ / (fracAp/2), Pd_matrix_);
      element->Store( fracPcKey_, makeScalar(PLAIN, fracPc_) );
    }

  template class FractureCapillaryPressureVisitor<1>;
  template class FractureCapillaryPressureVisitor<2>;
  template class FractureCapillaryPressureVisitor<3>;

}
