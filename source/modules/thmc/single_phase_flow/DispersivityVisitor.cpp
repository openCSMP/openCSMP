#include "DispersivityVisitor.h"
#include "Exception.h"
#include "InterFace.h"
#include "Element.h"
#include "Model.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
DispersivityVisitor< dim>::DispersivityVisitor( Model< dim>& sg, 
                                                const char* dispersivity, const char* pore_velocity, 
                                                double diffusivity, double dispersion_long, 
                                                double dispersion_trans )
  : pref(sg.Database()),
    vel_key(pref.StorageKey(pore_velocity)),
    disp_key(pref.StorageKey(dispersivity)),
    read_values(false)
  { 
    dp       = diffusivity;
    alpha_t  = dispersion_trans;
    alpha_l  = dispersion_long;
    
    if ( vel_key.place != ELEMENT || vel_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
          "The pore velocity must be a vector property placed on the element" ); 		       

    if ( disp_key.place != ELEMENT || disp_key.type != TENSOR )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
          "The dispersivity must be a tensor property placed on the element" ); 
    
    if ( alpha_t() > alpha_l() ) {
              cout << "\nalpha trans: " << alpha_t() << ", alpha long: " << alpha_l() << endl;
		throw csmp::Exception( INFO, "DispersivityVisitor< dim>::(constructor)", 
          "The transversal dispersivity is larger than the longitudinal dispersivity, is this really intended?" ); 

          
      }
    
    		       
  }

template<uint32_t dim>
DispersivityVisitor< dim>::DispersivityVisitor( Model< dim>& sg, 
                                                 const char* dispersivity, const char* pore_velocity, 
                                                 const char* diffusivity,  const char* dispersion_long, 
                                                 const char* dispersion_trans )
  : pref(sg.Database()),
    vel_key(pref.StorageKey(pore_velocity)),
    disp_key(pref.StorageKey(dispersivity)),
    dp_key(pref.StorageKey(diffusivity)),
    at_key(pref.StorageKey(dispersion_trans)),
    al_key(pref.StorageKey(dispersion_long)),
    read_values(true)
  { 

    if ( al_key.place != ELEMENT || al_key.type != SCALAR  )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 	       
          "The longitudinal dispersion coefficient must be a scalar property placed on the element" ); 	

    if ( at_key.place != ELEMENT || at_key.type != SCALAR  )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
 	      "The transveral dispersion coefficient must be a scalar property placed on the element" );    

	if ( vel_key.place != ELEMENT || vel_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
          "The pore velocity must be a vector property placed on the element" ); 		       
    
    if ( disp_key.place != ELEMENT || disp_key.type != TENSOR )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
          "The dispersivity must be a tensor property placed on the element" ); 		       
    
    if ( dp_key.place != ELEMENT || dp_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "DispersivityVisitor< dim>::(constructor)", 
          "The pore diffusion coefficient must be a scalar property placed on the element" ); 		       
   
   

  }


template<uint32_t dim>
DispersivityVisitor< dim>::~DispersivityVisitor() 
 {}


template<uint32_t dim>
void DispersivityVisitor< dim>::Visit(Element< dim>* n)   
 {     
   
   // rock properties
   n->Read( vel_key, vel );     
   if ( read_values ) { 
       n->Read( dp_key, dp );      
       n->Read( at_key, alpha_t );     
       n->Read( al_key, alpha_l );   
     }  


   v_abs = vel.Length();   
   vx = vel[0] * vel[0];
   if ( dim != 1U ) vy = vel[1] * vel[1];
   if ( dim == 3U ) vz = vel[2] * vel[2];


   // equations 3.14a-c from Ingebritsen's book Groundwater in Geologic Porcesses
   disp = 0.0;

   double cross_term(0.0);

   // diffusion
   for ( auto i=0; i<dim; i++ ) disp(i,i) = dp();

   if ( v_abs != 0.0 ) {
       // xx
       disp(0,0) += alpha_l() * vx / v_abs;
       disp(0,0) += alpha_t() * vy / v_abs;
       if ( dim == 3 ) disp(0,0) += alpha_t() * vz / v_abs;

       // yy
       disp(1,1) += alpha_t() * vx / v_abs;
       disp(1,1) += alpha_l() * vy / v_abs;
       if ( dim == 3 ) disp(1,1) += alpha_t() * vz / v_abs;

       // zz
       if ( dim == 3 ) {
           disp(2,2) += alpha_t() * vx / v_abs;
           disp(2,2) += alpha_t() * vy / v_abs;
           disp(2,2) += alpha_l() * vz / v_abs;
         }

       // off diagonal terms
       off_diag  = alpha_l();
       off_diag -= alpha_t();

       // 2D case
       if ( dim > 1 ) {
           //xy
           cross_term = vel[0] * vel[1];
           cross_term *= off_diag;
           cross_term /= v_abs;
           disp(0,1) = disp(1,0) = cross_term;
         }
       // 3D case
       if ( dim > 2 ) {
           // xz
           cross_term = vel[0] * vel[2];
           cross_term *= off_diag;
           cross_term /= v_abs;
           disp(0,2) = disp(2,0) = cross_term;
           // yz
           cross_term = vel[1] * vel[2];
           cross_term *= off_diag;
           cross_term /= v_abs;
           disp(1,2) = disp(2,1) = cross_term;
         }
     }

   n->Store( disp_key, disp );         
    
  }


template class DispersivityVisitor<1U>;
template class DispersivityVisitor<2U>;
template class DispersivityVisitor<3U>;

}
