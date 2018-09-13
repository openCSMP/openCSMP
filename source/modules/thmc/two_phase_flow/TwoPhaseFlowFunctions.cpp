#include "TwoPhaseFlowFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::TangentOfFractionalFlowFunction( const TARGET_PLACEMENT& p, double64 S ) const
  {
     const double64 srH2O = p.Obtain(User()->key_srH2O);
     return dfds_at(p, S) - (f_at(p, 0U, S) - f_at(p, 0U, srH2O)) / (S - srH2O);
  }
 
template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<1U,NODE>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<2U,NODE>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<3U,NODE>&, double64 ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentOfFractionalFlowFunction( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const ;




/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::FindRootSecantMethod( const TARGET_PLACEMENT& p, double64 S1, double64 S2 ) const
  {
    double64 F1 = 10000.;

    while ( abs(F1)>1e-10 ) {
      
        F1 = TangentOfFractionalFlowFunction(p, S1);
        double64 F2 = TangentOfFractionalFlowFunction(p, S2);
        
        double64 NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
        
        S2 = S1;
        S1 = NewPoint;
      }
    
    return  S1;
}

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<1U,NODE>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<2U,NODE>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<3U,NODE>&, double64 , double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<1U,ELEMENT>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<2U,ELEMENT>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<3U,ELEMENT>&, double64 , double64 ) const;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 , double64 ) const;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::FindRootSecantMethod( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 , double64 ) const;




/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::InflectionPointSaturation( const TARGET_PLACEMENT& p ) const
  {
    double64 S = 1.-p.Obtain(User()->key_srCO2) ;
    
    double64 Swmin = p.Obtain(User()->key_srH2O) ;
    double64 Swmax = 1.0-p.Obtain(User()->key_srCO2);
    
    double64 DS(0.001);
    double64 Fold(-10000.);
    
    int Maxiter(4) ;
    int it(1) ;
    
    while (it < Maxiter){
      
      double64 F1 = dfds_at(p, S);
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = dfds_at(p, S);
        
        if((S<Swmin)||(S>Swmax)) return S=0;
        
      }
      
      S = S + 2*DS ;
      DS = DS/10. ;
      it++ ;
      
      if((S<Swmin)||(S>Swmax)) return S=0;
      
    }
    
    S = S-10*DS ;
    
    if((S<Swmin)||(S>Swmax)) return S=0;
    
    
    return S ;
}
  
template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<1U,NODE>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<2U,NODE>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<3U,NODE>& ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<3U,ELEMENT>& ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::InflectionPointSaturation( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const ;



/**
 
  The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.

 */

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::TangentPointSaturation( const TARGET_PLACEMENT& p ) const
  {
    const double64 Si = InflectionPointSaturation(p);
    const double64 Sf = 1-p.Obtain(User()->key_srCO2) ;
    
    double64 St = FindRootSecantMethod(p,Si,Sf) ;
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
}

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<1U,NODE>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<2U,NODE>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<3U,NODE>& ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<3U,ELEMENT>& ) const ;

template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::TangentPointSaturation( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const ;






template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockHeight( const TARGET_PLACEMENT& p ) const
{
  return TangentPointSaturation(p);
}
  
template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::ShockHeight(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::ShockHeight(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::ShockHeight(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;
  


/**

  The Shock front wave calculated after estimation of tangent Saturation point.
  
*/

template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockFrontVelocity( const TARGET_PLACEMENT& p ) const
  {
    double64 S = TangentPointSaturation(p);
    S = std::min( std::max( S, 0. ), 1. );
    
    return dfds_at(p, S);
}
  
template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::ShockFrontVelocity( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::ShockFrontVelocity( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::ShockFrontVelocity( const FiniteElementPlacement<3U,ELEMENT>& ) const ;






/// shock speed base on Buckley Leverett theory
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 TwoPhaseFlowFunctions<dim,USER>::ShockSpeed( const  TARGET_PLACEMENT& p ) const
{
  return ShockFrontVelocity(p) ;
}
  
template double64 TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::ShockSpeed(  const  FiniteElementPlacement<1U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::ShockSpeed(  const  FiniteElementPlacement<2U,ELEMENT>& ) const;
template double64 TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::ShockSpeed(  const  FiniteElementPlacement<3U,ELEMENT>& ) const;




  
  
  
  
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
void TwoPhaseFlowFunctions<dim,USER>::ShockSpeedAndHeight( const TARGET_PLACEMENT& p, double64& speed, double64& height) const
{
  height = ShockHeight(p) ;
  speed  = ShockSpeed(p)  ;
}

template void TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>::ShockSpeedAndHeight(  const  FiniteElementPlacement<1U,ELEMENT>& , double64& , double64& ) const;
template void TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>::ShockSpeedAndHeight(  const  FiniteElementPlacement<2U,ELEMENT>& , double64& , double64& ) const;
template void TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>::ShockSpeedAndHeight(  const  FiniteElementPlacement<3U,ELEMENT>& , double64& , double64& ) const;





template class TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>;


} // end namespace csmp





