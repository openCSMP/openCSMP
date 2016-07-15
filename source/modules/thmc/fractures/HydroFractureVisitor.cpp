#include "HydroFractureVisitor.h"
#include "InterFace.h"
#include "Region.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {


template<size_t dim>
HydroFractureVisitor<dim>::HydroFractureVisitor( Model<dim>& sg )
  :  pref(sg.Database()),
     // property keys
     K_key(pref.StorageKey("conductivity")),
     Pe_key(pref.StorageKey("excess pressure")),
     S_key(pref.StorageKey("storativity")),
     V_key(pref.StorageKey("velocity")),
     fractured( sg.Region("Model").Elements() ),
     GRAD_LIMIT(101325.), // one bar m-1
     K_LIMIT(1.),         // maximum hydraulic conductivity that can be induced by 
     pres(3)              // hydrofracture process
  { 
     
     fractured.SetAll( false );
  }
         


      
template<size_t dim>
HydroFractureVisitor<dim>::~HydroFractureVisitor() 
 {  
 }




template<size_t dim>
void HydroFractureVisitor<dim>::Visit(Model<dim>* n)   
     { 
        char name[NAME_STRING];
        // Model has no id's therefore give element id's
        // get into the target element
        cout <<"\nHydroFractureVisitor::VisitModel: please enter name of ";
        cout <<" region from which you would like to spark streamlines\n";
        cin >> name;
        n->Region( name ).Accept( *this );
     }   
    



template<size_t dim>
void  HydroFractureVisitor<dim>::HydroFracturedElements( vector<size_t>& elmts ) const
 {
     elmts.clear();
     elmts.reserve(fractured.Size());

     for ( size_t n=0U; n<fractured.Size(); n++ )
       if ( fractured.GetBit(n) == true ) elmts.push_back( n );
       
     if ( elmts.empty() ) {
          cout <<"\nHydroFractureVisitor::HydroFracturedElements: Info: ";
          cout <<"No hydrofractured elements exist."<< endl; 
       }
     vector<size_t>( elmts ).swap( elmts );
 }
    
    
    
    
/**

@section implementation Implementation

1. test whether at least one element node is overpressured. 

2. only where an overpressure exists conductivity shall be diffused. 
     
3. back calculating dP/dxy from transport velocity. 

4. the conductivity is enhanced propertional to pressure gradient over the 
storativity, since  = dQ / bc. 
              
5. calculating desired flux in element with bc as the proportionality 
constant between desired pressure change and flux. 

6. calculating conductivity needed to accomodate this flux. This makes 
sense only where there is an appreciable pressure gradient. 
Otherwise an increased conductivity will not foster the flow.  

7. only if the new conductivity is greater than the old one it is 
assigned. 
 
*/
template<size_t dim>
void HydroFractureVisitor<dim>::Visit( Element<dim>* n )   
  { 
     // 1. test whether at least one element node is overpressured
     for ( over_pressured=false, i=0; i<n->Nodes(); i++ ) 
       {
          (n->N(i))->Read( Pe_key, Pe );
          pres[i] = Pe();
          if ( Pe() > 0.0 ) over_pressured = true;
       }
     // 2. only where an overpressure exists
     //    conductivity shall be diffused
     if ( over_pressured )
       { 
          // 0. finding the highest node overpressure
          Pe() = *max_element( pres.begin(), pres.end() );
          
          // 1. recording element which will be hydrofractured
          fractured.SetBit( n->Idx(), true );

          // 2. getting the conductivity
          n->Read( K_key, K );

          // only if the conductivity is not already very high
          if ( K() < K_LIMIT )
           {
              // 3. back calculating dP/dxy from transport velocity
              n->Read( V_key, V );
              dPdxy  = V;
              dPdxy /= K(); // get back to the pressure gradient

              // 4. the conductivity is enhanced propertional to 
              //    pressure gradient over the storativity, since
              //    dP = dQ / bc
              n->Read( S_key, S );
              dp = dPdxy.Length();
              
              // 5. calculating desired flux in element with bc
              //    as the proportionality constant between
              //    desired pressure change and flux
              q = S() * Pe();
              
              // 6. calculating conductivity needed to accomodate
              //    this flux
              //    NOTE: this makes sense only where there is an
              //    appreciable pressure gradient. Otherwise an
              //    increased conductivity will not foster the
              //    flow 
              if ( dp >= GRAD_LIMIT )
                {
                    Kf  = q / dp;
                    Kf += K;
                    // 7. only if the new conductivity is greater than the
                    //    old one it is assigned
                    if ( Kf > K ) n->Store( K_key, Kf );
                }
           }
       } // end if overpressured
    
  } // end VisitElement
    
   
template class HydroFractureVisitor<1U>;
template class HydroFractureVisitor<2U>;
template class HydroFractureVisitor<3U>;

} // csmp
