#ifndef CSMP_TRANSPORT_VISITOR_2D_H
#define CSMP_TRANSPORT_VISITOR_2D_H

#include "Visitor.h"
#include "ElementGrid.h"
#include "FiniteDifferenceGrid.h"
#include "AP_BoolVector.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 2002
*/

class TransportVisitor2D : public Visitor<2U> {
  public:
    TransportVisitor2D( Model<2U>&, 
                        const char* advected_prop, 
                        const char* advecting_prop );
                          
    ~TransportVisitor2D();
    
    virtual void Visit(Element<2U>* );   
    virtual void Visit(Model<2U>* ); 
    
    bool AdvectUntil( Model<2U>&, double final_time ); 

    size_t    MaximumIncrements() const;
    void      MaximumIncrements( size_t maxi );
    
    void      InterpolateOnlyWithinGrid( bool interpolate_just_within_grid );
    
  private:
     const PropertyDatabase<2U>&       pref;
     FiniteDifferenceGrid              grid;
     std::vector<ElementGrid>          egrids;
     BoolVector                        visited;
     double                          resolution;
     DenseMatrix<DM_MIN>               XY, NN;
     std::vector<VectorVariable<2U> >  P;
     std::vector<double>             xy;
     double                          time_increment;
     csmp::Index                       v_key, prop_key;
     size_t                            max_increments;
     bool                              interpolate_only_within_grid;   
     char                              transp_prop[200], adv_prop[200];

     void StoreResultsInGrid();
     
     bool IsInsideTriangle( double x, double y, bool update=true );
     
     void MinMaxCoordinates( double& min_x, double& max_x, 
                             double& min_y, double& max_y );
                             
     void InitializeElementGrid( size_t idx );
     
     void InputPropertyFromGrid( Model<2U>&, const char* prop ) const;
     void WritePropertyToGrid( Model<2U>&, const char* prop );
};

} // csmp

/* copyright (c) 2002 by S. K. Matthaei & S. G. Roberts & S. Geiger */

#endif



















