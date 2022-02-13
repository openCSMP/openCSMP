#include "SourceVisitor.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
SourceVisitor<dim>::SourceVisitor( Model<dim>& model, std::vector<std::string>* to_initialize_keys):
    Visitor<dim>( MODEL, ELEMENT ),
    model_ (model),
    to_initialize_keys_(to_initialize_keys),
    timestepKey_(csmp::Index()),
    thicknessKey_(csmp::Index())
{
    if (!this->to_initialize_keys_){
        porosityKey_    =           model.Database().StorageKey( "porosity"                            ) ;
        densityDiffKey_ =           model.Database().StorageKey( "density difference"                  ) ;
        nfvsKey_        =           model.Database().StorageKey( "nodal fluid volume source"           ) ;
        thicknessKey_   =           model.Database().StorageKey( "thickness"                           ) ;
    }
    else if (this->to_initialize_keys_->size()==5)
    {
        porosityKey_ =           model.Database().StorageKey( to_initialize_keys_->operator [](0).c_str());
                densityDiffKey_ =           model.Database().StorageKey( to_initialize_keys_->operator [](1).c_str()) ;
                nfvsKey_        =           model.Database().StorageKey( to_initialize_keys_->operator [](2).c_str()) ;
                timestepKey_    =           model.Database().StorageKey( to_initialize_keys_->operator [](3).c_str()) ;
                thicknessKey_   =           model.Database().StorageKey( to_initialize_keys_->operator [](4).c_str()) ;
    }
    else
    {
        throw csmp::Exception(FATAL_ERROR,"SourceVisitor(constructor)"," Incorrect number of variable names provided"," Please use the correct number.");
    }
}

template<uint32_t dim>
SourceVisitor<dim>::SourceVisitor(Model<dim>& model,
                                  Index porosityKey,
                                  Index densityDiffKey,
                                  Index nfvsKey,
                                  Index thicknessKey)
    : Visitor<dim>( MODEL, ELEMENT ),
      model_ (model),
      porosityKey_    (porosityKey),
      densityDiffKey_ (densityDiffKey),
      nfvsKey_        (nfvsKey),
      thicknessKey_   (thicknessKey),
      timestepKey_    (csmp::Index())
{

}

template<uint32_t dim>
SourceVisitor<dim>::~SourceVisitor()
{
}

template<uint32_t dim>
void SourceVisitor<dim>::Visit( Model<dim>* model )
{

    timestep()=1.0; // this is equivalent to not dividing the source term by the timestep, as used by the Geothermal Example.
    if (timestepKey_!=csmp::Index())
        model_.Read(timestepKey_,timestep);
    else
        timestep()=1.0; // this is equivalent to not dividing the source term by the timestep, as used by the Geothermal Example.
}

template<uint32_t dim>
void SourceVisitor<dim>::Visit( Element<dim>* e )
{
    e->Read( porosityKey_, phi);
    e->Read( thicknessKey_,thickness);
    //! if lower dimension elements are included in transport calculation thickness needs to be taken into account
    /// @todo need to find a more efficient way to perform this calculation. Julian 19.02.2016

    for (auto i=0;i< e->Nodes (); i++)
    {
        e->N(i)->Read(nfvsKey_, nfvs);
        nfvs() += e->N(i)->Read(densityDiffKey_) * e->SectorVolume(i) * phi() * thickness() / timestep();
        e->N(i)->Store(nfvsKey_, nfvs);
    }
    
}


template class SourceVisitor<1>;
template class SourceVisitor<2>;
template class SourceVisitor<3>;

}// end csmp

