#ifndef CSMP_BO_LOOKUP_VISITOR_H
#define CSMP_BO_LOOKUP_VISITOR_H

#include "Visitor.h"
#include "Index.h"
#include "ScalarVariable.h"

namespace csmp {

class BO_PropertyCalculatorInterface;
class ReservoirWaterPropertyCalculatorInterface;

  /**
   @class BO_LookupVisitor BOLookupVisitor.h
   
   @author L.J. Mosser
   @section motivation Motivation
   
   This Visitor populates Nodes of a given csmp mesh
   with temperature and pressure dependent oil and 
   reservoir brine properties.
   
   @section design Design Intent
   Modular design takes as an input to its constructor
   a pointer to an object implementing the BOPropertyCalculatorInterface
   as well as an object pointer implementing the ReservoirWaterPropertyCalculator.
   
   @attention Attention
   Initializing objects outside of the scope of this using normal pointers requires
   the implementer to handle garbage collection manually. Use shared pointers instead
   leads to safe execution and garbage collection.
   */
  
template<uint32_t dim>
class BO_LookupVisitor : public Visitor<dim>
  {
    public:
    
      /// Visitor is endowed with oil and water properties upon its construction
      BO_LookupVisitor( const Model<dim>&,
                        BO_PropertyCalculatorInterface&,
                        ReservoirWaterPropertyCalculatorInterface& );

      void Visit( Node<dim>* node ) override final;
      void Visit( Model<dim>* model ) override final;
      
    private:
      Index pressureFluidKey_,
            tempFluidKey_,
            viscosityOilKey_,
            compressibilityOilKey_,
            densityOilKey_,
            viscosityWaterKey_,
            compressibilityWaterKey_,
            densityWaterKey_;
    
      ScalarVariable presFluid_,
                     tempFluid_,
                     visOil_,
                     comOil_,
                     denOil_,
                     visWater_,
                     comWater_,
                     denWater_;

      BO_PropertyCalculatorInterface& bo_;
      ReservoirWaterPropertyCalculatorInterface& water_;    
  };
  
} //csmp

#endif // CSMP_BO_LOOKUP_VISITOR_H
