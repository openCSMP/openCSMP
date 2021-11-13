#include "GenericNodePropertyGradientLimiter.h"
#include "Model.h"
#include "Region.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

template<size_t dim>
GenericNodePropertyGradientLimiter<dim>::GenericNodePropertyGradientLimiter( Model<dim>& sg, const char* region, const char* prop_name )
  : gref_(sg.Region(region)),
    tolerance(1.0e-15),
    //gradient( sg_, prop, VECTOR, NODE ), // tbd NOT to be called with  prop, but with name of gradient!
    u_key(sg.Database().StorageKey(prop_name)),
    grad_key(sg.Database().StorageKey( (string(prop_name) + " gradient").c_str() )),
    lim_key(sg.Database().StorageKey( (string(prop_name) + " limiter").c_str() )),
    node_prop_grad( sg, prop_name ),
    limiter( 1U, NULL )
  {
    // initialize gradient PropHandle and set prop key:
    // ------------------------------------------------  
    string   lim_name;
    lim_name  =  prop_name;
    lim_name += " limiter";
    limiter[0] = new PropertyHandle<dim>( sg, region, lim_name.c_str(), SCALAR, NODE );
  }




template<size_t dim>
GenericNodePropertyGradientLimiter<dim>::GenericNodePropertyGradientLimiter( Model<dim>& sg, const char* region, vector<char*> prop_names )
  : gref_(sg.Region(region)),
    tolerance(1.0e-15),
    //gradient( sg, prop, VECTOR, NODE ), // tbd NOT to be called with  prop, but with name of gradient!
    node_prop_grad( sg, region, prop_names ),
    limiter( prop_names.size(), NULL )
  {   
    vector<string >  limiter_strings(prop_names.size());
    string   str2;
    str2 = " limiter";

    for( size_t i = 0; i < prop_names.size(); i++ ){
      // construct the name of the property's gradient:
      limiter_strings[i] = prop_names[i]; 
      limiter_strings[i] = limiter_strings[i] + str2;
      // property handle
      limiter[i] = new PropertyHandle<dim>( sg, limiter_strings[i].c_str() , SCALAR, NODE );
    }     
  }


  
template<size_t dim>
GenericNodePropertyGradientLimiter<dim>::~GenericNodePropertyGradientLimiter()
  {
  }


template<size_t dim>
void GenericNodePropertyGradientLimiter<dim>::CalculateGenericNodalGradient()
{
    node_prop_grad.CalculateGenericNodalGradient();
}



/**

Computes the 2d slope limiter for all generic finite volumes if a higher order
finite volume method is used. The slope limiting factor is stored in an
automatically at the property database.

The method is implemented from the various FiniteVolume<fT, dim>Visitors

*/
template<size_t dim>
void GenericNodePropertyGradientLimiter<dim>::CalculateSlopeLimiter( const vector<pair<double64,double64> >& MINMAX, int counter )
  {
    ScalarVariable                                  val1, val2, phi;
    VectorVariable<dim>                             grad, dist;
    typename std::vector<size_t>::iterator          sit;
    Node<dim>                                       fv;
    typename std::vector<size_t>::iterator          nit;
    std::pair<size_t, size_t>                       ids;
    double64                                        max, min, val_left, phi_temp;
    Node<dim>*                                      cvit, neighbor_node;
    size_t                                          nloc_id, global_el_id, local_facet_id,node_id;
    
    // loop over finite volumes
    // -------------------------
    for( typename std::vector<Node<dim>*>::const_iterator
         it=gref_.NodesBegin(); it!=gref_.NodesEnd(); it++){
         
        cvit = (*it); // node pointer
        node_id =  (*(*it)).Idx();
        // if finite volume is at a node with DIRICHLET or NEUMANN boundary conditions,
        // its phi value is zero (first order upwind scheme)
        if ( (*cvit).Status(u_key) == PLAIN || (*cvit).Status(u_key) == ANY ) {
        
            // read the gradient at the node
            (*cvit).Read( grad_key,  grad );
            
            //double64  node_x =  (*cvit).x();
            //double64  node_y =  (*cvit).y();

            // read concentration at the fv at its center of mass
            (*cvit).Read( u_key, val1 );

            //string name = sg_.Database().Name( u_key );
                    
            // set phi to 1.0
            phi = 1.0;

            // get the min/max values:
            ///*
            min = val1();
            max = val1();

            for ( nit = node_prop_grad.neighbors_[node_id].begin(); nit != node_prop_grad.neighbors_[node_id].end(); nit++ ) {
                  gref_.N(*nit)->Read( u_key, val2 );
                if ( val2() < min ) { min = val2(); }
                if ( val2() > max ) { max = val2(); }
            }
            //*/

            // get the min/max values:
            //min = MINMAX[ node_id ].first;
            //max = MINMAX[ node_id ].second;
              
            // sitting at a node, loop over parent el's
            // ------------------------------------------           
            for(  size_t p = 0; p< (*cvit).Parents() ; p++ ){

                // get the global parent id:
                global_el_id = (*cvit).Parent( p )->Idx();

                // get the corresponding element:
                if(global_el_id<gref_.Elements()){

                    Element<dim>* eptr = gref_.E( global_el_id);

                    // get local node number
                    nloc_id = (*cvit).ParentNodeNumber( p );

                    // at that parent element, loop over all facets that belong to the current node/fv
                    // -------------------------------------------------------------------------------
                    for ( size_t i=0U; i < eptr->FV()->FacetsPerSector(nloc_id); i++ ) {

                         // at that facet, get local facet_id used for determining distance:
                         local_facet_id = eptr->FV()->FacetSurroundingSector( nloc_id,i );

                         // get distance barycenter - facetcenter for that facet:
                         node_prop_grad.GenericDistanceFacetFVBarycenter( global_el_id, local_facet_id, nloc_id, dist );

                         // construct linear interpolant:
                         val_left  = val1();

                         for (size_t j=0;j<dim;j++)
                             val_left += grad[j]*dist[j];

                         /*
                         // calculate r for that facet:
                         if ( val_left > val1() ) {
                             r  = max;
                             r -= val1();
                             r /= ( val_left - val1() );
                         }
                         else if ( val_left < val1() ) {
                             r  = min;
                             r -= val1();
                             r /= ( val_left - val1() );
                         }else
                             r=1.0;

                         // ...and store that r for that node
                         // find the actual value for phi
                         if ( r < 1.0 ) phi_temp = r;
                         else   phi_temp = 1.0;
                         */
                         phi_temp = limitProperty(val_left, val1(),min,max);

                         // if necessary, set new value to phi
                         if ( phi_temp < phi() ) phi=phi_temp;
                         if ( fabs( phi() ) < tolerance ) phi = 0.0;
                         if ( essentiallyEqual( phi(), 0.0 ) ) break;
                    } // end loop over facet
                }
            } // end loop over parents -> this node is done            
        } // end if statement
        // if finite volume is at a node with DIRICHLET or NEUMANN boundary conditions, 
        // its phi value is zero (first order upwind scheme)  
        else phi = 0.0;
        
        // store the gradient limiter
        // ----------------------------
        (*cvit).Store( lim_key, phi );

    } // end loop over finite volumes
    
} // end function
  



/// provide the csmp::Index for the gradient variable
template<size_t dim>
void GenericNodePropertyGradientLimiter<dim>::SetPropertyKey( csmp::Index& key )  { u_key = key; }





/**

@return double64 Returns the storage in bytes required by the GenericNodePropertyGradientLimiter object

*/    
template<size_t dim>
double64 GenericNodePropertyGradientLimiter<dim>::SizeOf() const
  {
    return sizeof( *this );
  }

template class GenericNodePropertyGradientLimiter<1U>;
template class GenericNodePropertyGradientLimiter<2U>;
template class GenericNodePropertyGradientLimiter<3U>;


}
