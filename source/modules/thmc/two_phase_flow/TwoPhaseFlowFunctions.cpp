#include "TwoPhaseFlowFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::TangentOfFractionalFlowFunction( const Element<dim>* const e, double64 S ) const
  {
     const double64 srH2O = e->Read(User()->key_srH2O);
     return dfds_at(e, S) - (f_at(e, 0U, S) - f_at(e, 0U, srH2O)) / (S - srH2O);
  }




/**
 
 Using the the SecantMethod to find the root of The Buckley- Leverett function See Eq. 1.86 in page 44 from Guinot book.
 
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::FindRootSecantMethod( const Element<dim>* const e, double64 S1, double64 S2 ) const
  {
    double64 F1 = 10000.;

    while ( abs(F1)>1e-10 ) {
      
        F1 = TangentOfFractionalFlowFunction(e, S1);
        double64 F2 = TangentOfFractionalFlowFunction(e, S2);
        
        double64 NewPoint = S1 - F1*(S1-S2)/(F1-F2) ;
        
        S2 = S1;
        S1 = NewPoint;
      }
    
    return  S1;
}




/**
   
  The Inflection Saturation Point, calculated from the maxima of 1st derivative fractional flow function See page 144 from Helmig book.
  This can be more accurate by puting in the loop of more and more finer maxima serach algorithm.
 
*/
template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::InflectionPointSaturation( const Element<dim>* const e ) const
  {
    double64 S     = 1.-e->Read(User()->key_srCO2) ;
    double64 Swmin = e->Read(User()->key_srH2O) ;
    double64 Swmax = 1. - S;
    
    double64 DS(0.001);
    double64 Fold(-10000.);
    
    int Maxiter(4) ;
    int it(1) ;
    
    while (it < Maxiter){
      
      double64 F1 = dfds_at(e, S);
      while (F1 > Fold) {
        
        S = S - DS ;
        Fold = F1 ;
        F1 = dfds_at(e, S);
        
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




/**
 
  The Tanget Saturation Point, calculated from the Buckley-Leverett problem See Eq. 1.86 in page 44 from Guinot book. To find root of this nonlinear function, I use The Secant Algorithm.

 */

template<size_t dim, template<size_t> class USER>
double64 TwoPhaseFlowFunctions<dim,USER>::TangentPointSaturation( const Element<dim>* const e ) const
  {
    const double64 Si = InflectionPointSaturation(e);
    const double64 Sf = 1-e->Read(User()->key_srCO2) ;
    
    double64 St = FindRootSecantMethod(e,Si,Sf) ;
    St = std::min( std::max( St, 0. ), 1. );
    
    return  St;
    
}









template class TwoPhaseFlowFunctions<1U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<2U,CO2H2O_FunctionsModule2>;
template class TwoPhaseFlowFunctions<3U,CO2H2O_FunctionsModule2>;


} // end namespace csmp





