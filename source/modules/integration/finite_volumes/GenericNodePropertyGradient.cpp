#include "GenericNodePropertyGradient.h"
#include "FiniteVolumeTraits.h"

using namespace std;

namespace csmp {

/**
    Sets the calculator up for computation on the entire model
    (a throughgoing node numbering is assumed.
*/
template<size_t dim>
GenericNodePropertyGradient<dim>::GenericNodePropertyGradient( Model<dim>& sg,
                                                               const FiniteVolumeStencilManager<dim>& fv_connectivity,
                                                               const char* prop )
  : sg_(sg),
  fcv( fv_connectivity ),
  gradient( 1U, NULL ),
  sum_x2( sg.Region("Model").Nodes(), 0.0 ),
  sum_y2( sg.Region("Model").Nodes(), 0.0 ),
  sum_z2( sg.Region("Model").Nodes(), 0.0 ),
  sum_xy( sg.Region("Model").Nodes(), 0.0 ),
  sum_xz( sg.Region("Model").Nodes(), 0.0 ),
  sum_yz( sg.Region("Model").Nodes(), 0.0 ),
  det( sg.Region("Model").Nodes() ),
  zero_grad_index_( sg.Region("Model").Nodes(), 6U ),
  distance( sg.Region("Model").Nodes() ),
  center_of_mass_( sg.Region("Model").Nodes() ),
  tolerance(1.0e-12),
  inner_(1U, VectorVariable<dim>(PLAIN,0.0) ),
  middle_( 1U, inner_ ),
  distance_facet_FVBarycenter_( sg.Region("Model").Elements() ),
  neighbors_( sg.Region("Model").Nodes() ),
  mass_center( sg, "mass center", VECTOR, NODE ) //PropertyHandle for mass center
  {

    // initialize gradient PropHandle and set prop key:
    // ------------------------------------------------  
    string   gradient_name;
    gradient_name  =  prop;
    gradient_name += " gradient";
    gradient[0] = new PropertyHandle<dim>( sg_, gradient_name.c_str(), VECTOR, NODE );
               
    u_key = sg_.Database().StorageKey( prop );
     
    // --------------------------------------------------
    // initialization and calculation of mass centers...
    // --------------------------------------------------
    std::pair<VectorVariable<dim>, double64>   init_pair;
    init_pair.first = VectorVariable<dim>( PLAIN, 0.0);
    init_pair.second = 0.0;

    fill(  center_of_mass_.begin(), center_of_mass_.end(), init_pair );

    CalculateCenterOfMass( fcv );

    // ---------------------------------------------------------------------------
    // initialization and calculation of distances facet centers-barycenters
    // ---------------------------------------------------------------------------

    // init of vector
    // ====================================================

    Element<dim>  e;

    distance_facet_FVBarycenter_.resize( sg_.Region("Model").Elements() ); // resize outer vector
    size_t length_;
    length_ =  distance_facet_FVBarycenter_.size();

    for ( typename vector<Element<dim>*>::const_iterator
          eit=sg_.Region("Model").ElementsBegin();
          eit!=sg_.Region("Model").ElementsEnd(); eit++ ){

          e = *(*eit);

          //resize middle vector
          distance_facet_FVBarycenter_[ e.Idx() ].resize( e.FV_Stencil()->Facets() );

          // loop over facets
          // -------------------
          for( size_t fi=0U; fi < e.FV_Stencil()->Facets(); fi++ ){

                 // resize inner vector
                distance_facet_FVBarycenter_[ e.Idx() ][ fi ].resize( e.Nodes() );
                length_ = distance_facet_FVBarycenter_[ e.Idx()][ fi ].size();

                fill( distance_facet_FVBarycenter_[ e.Idx() ][ fi ].begin(),
                      distance_facet_FVBarycenter_[ e.Idx() ][ fi ].end(),
                      VectorVariable<dim>( PLAIN, 0.0)       );
          }
    }

    // ====================================================

    // ... and distances:
    CalculateDistanceFacetFVBary( fcv, distance_facet_FVBarycenter_ );

    // -----------------------------
    // neighbours and least squares
    // -----------------------------


    // initialize vector with global node id's of neighbors:
    // =====================================================
    Element<dim>     current_el;
    size_t           global_neighb_el_id, current_n_id;
    std::vector<size_t> ids;

    // loop over all nodes/fv's
    // --------------------------
    for( typename std::vector<Node<dim>*>::const_iterator
         cvit=sg_.Region("Model").NodesBegin();
         cvit!=sg_.Region("Model").NodesEnd(); cvit++){

         current_n_id = (*(*cvit)).Idx();
         // loop over parents
         for(  size_t p = 0; p< (*cvit)->Parents() ; p++ ){
             // get the global parent id:
             global_neighb_el_id = (*cvit)->Parent( p)->Idx();

             if(global_neighb_el_id<sg_.Region("Model").Elements()){
                 // get the corresponding element:
                 current_el = (*sg_.Region("Model").E( global_neighb_el_id ));
                 // ...and the global node_i's of that element
                 ids.clear();
                 for(size_t i=0;i<current_el.Nodes();i++)
                     if(current_n_id!=current_el.N(i)->Idx()) ids.push_back(current_el.N(i)->Idx());

                 // insert into vector without duplicates:
                 PushBackAvoidDuplicate( neighbors_[ current_n_id ], ids );
             }
       }
    }
    // =====================================================

    CalculateGenericLeastSquareSums();
}




/**

Constructor for several variables. This is needed if, e.g. several components
and fluid phases need to be limited.
The name of the properties to be limited needs to be passed on in the vector of characters
vector<char* > properties.

*/
template<size_t dim>
GenericNodePropertyGradient<dim>::GenericNodePropertyGradient(  Model<dim>& sg,
                                                                const FiniteVolumeStencilManager<dim>& fv_connectivity,
                                                                std::vector<char* > properties )
: sg_(sg),
  fcv(fv_connectivity), 
  gradient( 1U, NULL ),
  sum_x2( sg.Region("Model").Nodes(), 0.0 ),
  sum_y2( sg.Region("Model").Nodes(), 0.0 ),
  sum_z2( sg.Region("Model").Nodes(), 0.0 ),
  sum_xy( sg.Region("Model").Nodes(), 0.0 ),
  sum_xz( sg.Region("Model").Nodes(), 0.0 ),
  sum_yz( sg.Region("Model").Nodes(), 0.0 ),
  det( sg.Region("Model").Nodes() ),
  zero_grad_index_( sg.Region("Model").Nodes(), 6U ),
  distance( sg.Region("Model").Nodes() ),
  center_of_mass_( sg.Region("Model").Nodes() ),
  tolerance(1.0e-15),
  inner_(1U, VectorVariable<dim>() ),
  middle_( 1U, inner_ ),
  distance_facet_FVBarycenter_( sg.Region("Model").Elements() ),
  neighbors_( sg.Region("Model").Nodes() ),
  mass_center( sg, "mass center", VECTOR, NODE ) //PropertyHandle for mass center
  { 
    // --------------------------------------------------
    // initialization of gradients 
    // --------------------------------------------------
    vector<string >  gradient_strings(properties.size()); 
    string   str2;
    str2 = " gradient";

    for( size_t i = 0; i < properties.size(); i++ ){
          // construct the name of the property's gradient:
          gradient_strings[i] = properties[i];
          gradient_strings[i] = gradient_strings[i] + str2;

          // property handle
          gradient[i] = new PropertyHandle<dim>( sg_, gradient_strings[i].c_str() , VECTOR, NODE );
    }
     
    // --------------------------------------------------
    // initialization and calculation of mass centers... 
    // --------------------------------------------------
    std::pair<VectorVariable<dim>, double64>   init_pair;
    init_pair.first = VectorVariable<dim>( PLAIN, 0.0);
    init_pair.second = 0.0;
   
    fill(  center_of_mass_.begin(), center_of_mass_.end(), init_pair );
       
    CalculateCenterOfMass( fcv );  
    
    // ---------------------------------------------------------------------------
    // initialization and calculation of distances facet centers-barycenters
    // ---------------------------------------------------------------------------
    
    // init of vector 
    // ====================================================
    
    Element<dim>  e;
    
    distance_facet_FVBarycenter_.resize( sg_.Region("Model").Elements() ); // resize outer vector
    size_t length_;
    length_ =  distance_facet_FVBarycenter_.size();
     
    for ( typename vector<Element<dim>*>::const_iterator
          eit=sg_.Region("Model").ElementsBegin();
          eit!=sg_.Region("Model").ElementsEnd(); eit++ ){
          
          e = *(*eit);
          cout<< e.Idx() <<endl;

          //resize middle vector
          distance_facet_FVBarycenter_[ e.Idx() ].resize( e.FV_Stencil()->Facets() );

          // loop over facets
          // -------------------
          for( size_t fi=0U; fi < e.FV_Stencil()->Facets(); fi++ ){

                 // resize inner vector
                distance_facet_FVBarycenter_[ e.Idx() ][ fi ].resize( e.Nodes() );
                length_ = distance_facet_FVBarycenter_[ e.Idx()][ fi ].size();

                fill( distance_facet_FVBarycenter_[ e.Idx() ][ fi ].begin(),
                      distance_facet_FVBarycenter_[ e.Idx() ][ fi ].end(),
                      VectorVariable<dim>( PLAIN, 0.)       );
          }
    }
    
    // ====================================================   
   
    // ... and distances:
    CalculateDistanceFacetFVBary( fcv, distance_facet_FVBarycenter_ ); 
    
    // -----------------------------
    // neighbours and least squares
    // -----------------------------
     
    
    // initialize vector with global node id's of neighbors:
    // =====================================================
    Element<dim>     current_el;
    size_t           global_neighb_el_id, current_n_id;
    std::vector<size_t> ids;
    
    // loop over all nodes/fv's
    // --------------------------
    for( typename std::vector<Node<dim>*>::const_iterator
         cvit=sg_.Region("Model").NodesBegin();
         cvit!=sg_.Region("Model").NodesEnd(); cvit++){
         
         current_n_id = (*(*cvit)).Idx();
         // loop over parents
         for(  size_t p = 0; p< (*cvit)->Parents() ; p++ ){
             // get the global parent id:
             global_neighb_el_id = (*cvit)->Parent( p)->Idx();
             if(global_neighb_el_id<sg_.Region("Model").Elements()){
                 // get the corresponding element:
                 current_el = (*sg_.Region("Model").E( global_neighb_el_id ));
                 // ...and the global node_i's of that element
                 ids.clear();
                 for(size_t i=0;i<current_el.Nodes();i++)
                     if(current_n_id!=current_el.N(i)->Idx()) ids.push_back(current_el.N(i)->Idx());

                 // insert into vector without duplicates:
                 PushBackAvoidDuplicate( neighbors_[ current_n_id ], ids );
             }
         }
    }
    // =====================================================
     
    CalculateGenericLeastSquareSums(); 
}



template<size_t dim>
void GenericNodePropertyGradient<dim>::SetPropertyKey( csmp::Index& key )  { u_key = key; }

  
template<size_t dim>
GenericNodePropertyGradient<dim>::~GenericNodePropertyGradient()
{}
  
  
  
  

/**

Calculates the center of mass of generic finite volumes in 2d.

*/
template<size_t dim>
void GenericNodePropertyGradient<dim>::CalculateCenterOfMass( const FiniteVolumeStencilManager<dim>& fcv )
{

 Element<dim>  e;
 size_t  glob_n_id;
 double64 volume_;
 std::vector<double64>   current_bc; // stencil barycenter in global coord's
 std::vector<size_t> ids;
 VectorVariable<dim> temp_;
 Point<dim> tmp_p;
 // loop over stencils
 // ------------------
 for ( typename vector<Element<dim>*>::const_iterator
       eit=sg_.Region("Model").ElementsBegin();
       eit!=sg_.Region("Model").ElementsEnd(); eit++ ){

        e = *(*eit);
        ids.resize( e.Nodes() );
        // get vector of global id's
        for(size_t i=0;i<e.Nodes();i++)
            ids[i]=e.N(i)->Idx();


        // loop over sectors
        // -------------------
        for( size_t i=0U; i < e.FV_Stencil()->Sectors(); i++ ){

              //get the global node id for the current segment
              glob_n_id = ids[ i ];

              if( sg_.Region("Model").N( glob_n_id)->AtBoundary() ){

                  tmp_p = sg_.Region("Model").N( glob_n_id)->Coordinate();
                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for(size_t k=0;k<dim;k++)
                      temp_.Component(k,tmp_p.Coordinates()[k]);
                  center_of_mass_[ glob_n_id ].first  = temp_;
                  center_of_mass_[ glob_n_id ].second = 1.;

              }else{

                  //get the barycenter of current segment
                  ConvertToGlobalCoordinates( e, e.FV_Stencil()->SectorIntegrationPoint( i, 0U), current_bc );

                  //get volume of current sector
                  volume_ = e.SectorVolume( i );

                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for(size_t k=0;k<dim;k++)
                      temp_.Component(k,current_bc[k] * volume_);

                  center_of_mass_[ glob_n_id ].first += temp_;
                  center_of_mass_[ glob_n_id ].second += volume_;

             }
        }
     }

     for( size_t i = 0U;  i< sg_.Region("Model").Nodes(); i++ )
     {
        //  center_of_mass_[ i ].first.Out();
        //  calculate center: sum_i(x_i * A_i) / sum_i(A_i)
        center_of_mass_[ i ].first *=  1. / center_of_mass_[ i ].second;

        // write back to node:
        sg_.Region("Model").N( i )->Store( sg_.Database().StorageKey("mass center"), center_of_mass_[ i ].first );
     }

}


template<size_t dim>
void GenericNodePropertyGradient<dim>::CalculateDistanceFacetFVBary(  const FiniteVolumeStencilManager<dim>& fcv,
                                     std::vector<std::vector<std::vector<VectorVariable<dim> > > >& d_facet_FVBarycenter  )
{
 Element<dim>  e;
 // local node id's:
 size_t inside_node_, outside_node_;
 std::vector<double64>  facet_bc(dim);  // global coordinates of facet barycenter
 VectorVariable<dim>  mass_center;


 // loop over elements
 // ------------------
 for ( typename vector<Element<dim>*>::const_iterator
       eit=sg_.Region("Model").ElementsBegin();
       eit!=sg_.Region("Model").ElementsEnd(); eit++ ){

        e = *(*eit);

        // loop over facets
        // -------------------
        for( size_t fi=0U; fi < e.FV_Stencil()->Facets(); fi++ ){

              // get facet barycenter in global coordinates:
              ConvertToGlobalCoordinates( e, e.FV_Stencil()->FacetIntegrationPoint( fi, 0U ),  facet_bc);

              // get local node id's
              e.FV_Stencil()->FacetEdgeNodes( fi, inside_node_, outside_node_ );


              // get mass center for FV of inside_node_:
              GenericCenterOfMass(  e.N(inside_node_)->Idx(), mass_center );

              VectorVariable<dim> face_dist_inside_node( PLAIN,  0.0);
              for(size_t k=0;k<dim;k++)
                  face_dist_inside_node.Component(k,facet_bc[k] - mass_center[k]);

              d_facet_FVBarycenter[ e.Idx() ][ fi ][ inside_node_ ] =face_dist_inside_node;


              // get mass center for FV of outside_node_:
              GenericCenterOfMass(  e.N(outside_node_)->Idx(), mass_center );

              // calculate distance vector:
              VectorVariable<dim> face_dist_outside_node( PLAIN,  0.0);
              for(size_t k=0;k<dim;k++)
                  face_dist_outside_node.Component(k,facet_bc[k] - mass_center[k]);

              d_facet_FVBarycenter[ e.Idx() ][ fi ][ outside_node_ ] = face_dist_outside_node;

        }
 }
// bool hallo;

//  return_d_facet_FVBarycenter =  d_facet_FVBarycenter;
}




/**

Inserts all components of possible_new_entries into old_vector if they
are not there already.

*/
template<size_t dim>
void GenericNodePropertyGradient<dim>::PushBackAvoidDuplicate( std::vector<size_t>& old_vector,
                                                                 std::vector<size_t>& new_vector )
{
    bool already_exists(false);

    for( typename std::vector<size_t>::const_iterator
            new_ = new_vector.begin();
            new_ != new_vector.end(); new_++ ){
        for(  typename std::vector<size_t>::const_iterator
              old_ = old_vector.begin(); old_ != old_vector.end(); old_++   ){
             if( *new_ == *old_ ) already_exists = true;
        }
        // current element does not exist yet, insert:
        if( !already_exists ) old_vector.push_back( *new_ );
        already_exists = false;
    }
}





template<size_t dim>
void GenericNodePropertyGradient<dim>::ConvertToGlobalCoordinates( Element<dim>& el, const Point<dim>& local_c_point, std::vector<double64>& global_c )
{


  std::vector<double64> temp(el.Nodes()); //has the local interp. function values
  std::vector<double64> local_c(local_c_point.Coordinates());

  if( el.IsLineElement() ){
      el.FE()->Nr( local_c[0], temp );
  }else if( el.IsSurfaceElement()){
      el.FE()->Nrs( local_c[0], local_c[1], temp );
  }else{
      el.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
  }

  el.CoordinateMatrix();

  global_c.assign( dim, 0.0);

  // transform local c's to global c's
  for (size_t i = 0; i<el.Nodes(); i++)
      for (size_t j = 0; j<dim; j++){
        global_c[j] += el.FE()->XY(i,j) * temp[i];
  }
}



/**

Returns the calculated centers of mass of generic finite volumes in 2d.

@param size_t_global_id:       global_node_id: the node id of the corresponding FV

@param VectorVariable:            the Vector onto which the xyz-coordinates of the mass center
                         of the FV that belongs to the the node are written

*/
template<size_t dim>
void GenericNodePropertyGradient<dim>::GenericCenterOfMass( size_t global_node_id, VectorVariable<dim>& mass_center_out ) const
{
 mass_center_out =  center_of_mass_[global_node_id].first;
}

template<size_t dim>
void GenericNodePropertyGradient<dim>::GenericDistanceFacetFVBarycenter( size_t global_el_id,
                                                                           size_t local_facet_id,
                                                                           size_t local_node_id,
                                                                           VectorVariable<dim>& distance ) const
{
 distance =   distance_facet_FVBarycenter_[ global_el_id ][ local_facet_id ][ local_node_id ] ;
}







/**

Returns the storage required by the GenericNodePropertyGradient object

*/
template<size_t dim>
double64 GenericNodePropertyGradient<dim>::SizeOf() const
{
  double64 storage;

  storage  = sizeof( *this );
  storage += static_cast<double64>(sum_xy.size()) * sizeof(double64);
  storage += static_cast<double64>(sum_x2.size()) * sizeof(double64);
  storage += static_cast<double64>(sum_y2.size()) * sizeof(double64);
  storage += static_cast<double64>(det.size())    * sizeof(double64);

  for ( size_t i=0; i<distance.size(); i++ )
      storage += static_cast<double64>(distance[i].size()) * 2.0 * sizeof(double64);

  return storage;

}






/**

A private method that computes the least squares sums, i.e., the sum of the distances between
the center of mass of a finite volume to the centers of mass of all 
neighboring finite volumes. The sums are needed for the subseqent finite 
volume computations and are stored in STL vectors.

If the mesh doesn't change, this only needs to be calculated once for this class.

This computation is done automatically during the construction of the object.

A message states the successful computation of the least squares sums.

*/
template<>
void GenericNodePropertyGradient<1U>::CalculateGenericLeastSquareSums()
  {
    enum{DIM=1U};
    std::vector<Node<DIM>*>::const_iterator  cvit;
    std::vector<size_t>::iterator            nit;
    VectorVariable<DIM>                      xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);
    size_t                                   id, node_id;

    for( cvit=sg_.Region("Model").NodesBegin();
         cvit!=sg_.Region("Model").NodesEnd(); cvit++){

        node_id = (*(*cvit)).Idx() ;

        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );
        
        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient 
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );
        id = 0;           
        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        for ( nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {
        
            GenericCenterOfMass( (*nit), xyz2 );
            
            // subtract x,y,z-coordinate_neighbour_center - x,y,z-coordinate_current_fv
            dxyz.Component(0, xyz2[0U]-xyz1[0]);

            // save the distance to the neighboring FVs; Used in CalculateNodalGradient
            // for the RHS
            distance[ node_id ][id] = dxyz;

            sum_x2[ node_id ] += dxyz[0] * dxyz[0];                        
        }
 
        // calculate determinante
        det[  node_id ]  = sum_x2[ node_id ];

    }
    
    cout << "\nGenericNodePropertyGradient<dim>::CalculateLeastSquareSums:";
    cout << "\nLeast squares sums successfully calculated... ";
    cout << endl;
      
  }  // end CalculateLeastSquareSums
  
  
  

/** GenericNodePropertyGradient<dim>::CalculateNodalGradient()

@section Description

Computes the gradient of a property (concentration, saturation, etc.),
by summing up the difference in the property value in a single finite
volume and all its neighbors times the according least square sum.
Those gradients are needed for the higher order finite volume methods
and are stored in STL vectors.

Note: For both the case where just one gradient is needed as for the one where 
several ones are needed: The user has to set the Property Key to the respective value!

@section Implementation

The method is implemented from the various FiniteVolume<fT, dim>Visitors and the STL
vectors that store the x and y components of the gradients can be accessed
through the interfaces GradientComponentX() and GradientComponentY()

@test tested OK */
template<>
void GenericNodePropertyGradient<1U>::CalculateGenericNodalGradient()
{
  
    enum{DIM=1U};
    std::vector<Node<DIM>*>::const_iterator     cvit;
    std::vector<size_t>::iterator               nit;
    double64                                    val1, val2, dc;
    VectorVariable<DIM>                         grad(PLAIN,0.0);
    VectorVariable<DIM>                         dcxyz(PLAIN,0.0);
    size_t                                      id, node_id;
    VARIABLE_FLAG                               status;
    string                                      grad_name;
    
    for( std::vector<Node<DIM>*>::const_iterator
           cvit=sg_.Region("Model").NodesBegin();
           cvit!=sg_.Region("Model").NodesEnd(); cvit++){

           status = (*(*cvit)).Status( u_key );
           node_id =  (*(*cvit)).Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY) {
        
                // get neighbor ids at every cv and read concentration
                val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.0;
                id  = 0;
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for ( nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    val2 = (sg_.Region("Model").N( (*nit) ))->Read (u_key);

                    dc   = (val2 -val1);

                    // computing temporary variables for least squares calculation; needed for RHS of LGS
                    // NOTE: id's are the same as in CalculateLeastSquares since list is gone thorugh in same order
                    dcxyz.Component(0, dcxyz[0]+dc*distance[ node_id ][id][0]);

                }

                if(det[ node_id]!=0.0){

                    grad(0)  = dcxyz[0]/det[ node_id ];

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;

                }else grad=0.0;

        }
        // gradient is zero if nodes are flagged DIRICH or NEUMANN
        else grad = 0.0;
        
        // store the gradient
        // ------------------

        grad_name = sg_.Database().Name( u_key );
        grad_name = grad_name + " gradient";
        (*(*cvit)).Store( sg_.Database().StorageKey( grad_name.c_str() ), grad );
    }
} // end CalculateNodalGradient
  
  
  
  
  
  

/**

A private method that computes the least squares sums, i.e., the sum of the distances between
the center of mass of a finite volume to the centers of mass of all
neighboring finite volumes. The sums are needed for the subseqent finite
volume computations and are stored in STL vectors.

If the mesh doesn't change, this only needs to be calculated once for this class.

This computation is done automatically during the construction of the object.

A message states the successful computation of the least squares sums

*/
template<>
void GenericNodePropertyGradient<2U>::CalculateGenericLeastSquareSums()
  {
    enum{DIM=2U};
    std::vector<Node<DIM>*>::const_iterator  cvit;
    std::vector<size_t>::iterator            nit;
    VectorVariable<DIM>                      xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);
    size_t                                   id, node_id;

    for( cvit=sg_.Region("Model").NodesBegin();
         cvit!=sg_.Region("Model").NodesEnd(); cvit++){

        node_id = (*(*cvit)).Idx() ;

        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );

        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );
        id = 0;
        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        for ( nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {

            GenericCenterOfMass( (*nit), xyz2 );

            // subtract x,y,z-coordinate_neighbour_center - x,y,z-coordinate_current_fv

            dxyz.Component(0, xyz2[0U]-xyz1[0]);
            dxyz.Component(1, xyz2[1U]-xyz1[1]);

            // save the distance to the neighboring FVs; Used in CalculateNodalGradient
            // for the RHS
            distance[ node_id ][id] = dxyz;

            // Accumulate LHS matrix ([sum_x2,sum_xy],[sum_xy,sum_y2])
            sum_x2[ node_id ] += dxyz[0] * dxyz[0];
            sum_y2[ node_id ] += dxyz[1] * dxyz[1];
            sum_xy[ node_id ] += dxyz[0] * dxyz[1];
        }

        // calculate determinante
        det[  node_id ]  = sum_x2[ node_id ]*sum_y2[ node_id ]-sum_xy[ node_id ] * sum_xy[ node_id ];

        if(det[ node_id ] == 0.0){

            if(sum_x2[ node_id ] == 0.0)
                zero_grad_index_[ node_id]=0U;
            else if(sum_y2[ node_id ] == 0.0)
                zero_grad_index_[ node_id]=1U;
            else
                zero_grad_index_[ node_id ]=6U;
        }

    }

    cout << "\nGenericNodePropertyGradient<dim>::CalculateLeastSquareSums:";
    cout << "\nLeast squares sums successfully calculated... ";
    cout << endl;

  }  // end CalculateLeastSquareSums






/**

Computes the gradient of a property (concentration, saturation, etc.),
by summing up the difference in the property value in a single finite
volume and all its neighbors times the according least square sum.
Those gradients are needed for the higher order finite volume methods
and are stored in STL vectors.

Note: For both the case where just one gradient is needed as for the one where
several ones are needed: The user has to set the Property Key to the respective value!

The method is implemented from the various FiniteVolume<fT, dim>Visitors and the STL
vectors that store the x and y components of the gradients can be accessed
through the interfaces GradientComponentX() and GradientComponentY()

*/
template<>
void GenericNodePropertyGradient<2U>::CalculateGenericNodalGradient()
{
    enum{DIM=2U};
    std::vector<Node<DIM>*>::const_iterator     cvit;
    std::vector<size_t>::iterator               nit;
    double64                                    val1, val2, dc;
    VectorVariable<DIM>                         grad(PLAIN,0.0);
    VectorVariable<DIM>                         dcxyz(PLAIN,0.0);
    size_t                                      id, node_id;
    VARIABLE_FLAG                               status;
    string                                      grad_name;

    for( std::vector<Node<DIM>*>::const_iterator
           cvit=sg_.Region("Model").NodesBegin();
           cvit!=sg_.Region("Model").NodesEnd(); cvit++){

           status = (*(*cvit)).Status( u_key );
           node_id =  (*(*cvit)).Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY) {

                // get neighbor ids at every cv and read concentration
                val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.0;
                id  = 0;
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for ( nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    val2 = (sg_.Region("Model").N( (*nit) ))->Read (u_key);

                    dc   = (val2 -val1);

                    // computing temporary variables for least squares calculation; needed for RHS of LGS
                    // NOTE: id's are the same as in CalculateLeastSquares since list is gone thorugh in same order
                    dcxyz.Component(0, dcxyz[0]+dc*distance[ node_id ][id][0]);
                    dcxyz.Component(1, dcxyz[1]+dc*distance[ node_id ][id][1]);
                }

                if(det[ node_id]!=0.0){

                    grad(0)  = (dcxyz[0] * sum_y2[  node_id ] - dcxyz[1] * sum_xy[  node_id ])/det[ node_id ];
                    grad(1)  = (dcxyz[1] * sum_x2[  node_id ] - dcxyz[0] * sum_xy[  node_id ])/det[ node_id ];

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 0U){

                    grad(0) = 0.0;
                    grad(1) = dcxyz[1]/ sum_y2[ node_id ];

                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 1U){

                    grad(0) = dcxyz[0]/ sum_x2[ node_id ];
                    grad(1) = 0.0;

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;

                }else
                    grad = 0.0;

        }
        // gradient is zero if nodes are flagged DIRICH or NEUMANN
        else grad = 0.0;

        // store the gradient
        // ------------------

        grad_name = sg_.Database().Name( u_key );
        grad_name = grad_name + " gradient";
        (*(*cvit)).Store( sg_.Database().StorageKey( grad_name.c_str() ), grad );
    }
} // end CalculateNodalGradient




/**

void GenericNodePropertyGradient<dim>::CalculateLeastSquareSums()

@section Description

A private method that computes the least squares sums, i.e., the sum of the distances between
the center of mass of a finite volume to the centers of mass of all
neighboring finite volumes. The sums are needed for the subseqent finite
volume computations and are stored in STL vectors.

If the mesh doesn't change, this only needs to be calculated once for this class.

@section Implementation:

This computation is done automatically during the construction of the object.

@section Messages

A message states the successful computation of the least squares sums

*/
template<>
void GenericNodePropertyGradient<3U>::CalculateGenericLeastSquareSums()
  {
    enum{DIM=3U};
    std::vector<Node<DIM>*>::const_iterator  cvit;
    std::vector<size_t>::iterator            nit;
    VectorVariable<DIM>                      xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);
    size_t                                   id, node_id;

    for( cvit=sg_.Region("Model").NodesBegin();
         cvit!=sg_.Region("Model").NodesEnd(); cvit++){

        node_id = (*(*cvit)).Idx() ;

        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );

        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );
        id = 0;
        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        for ( nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {

            GenericCenterOfMass( (*nit), xyz2 );

            // subtract x,y,z-coordinate_neighbour_center - x,y,z-coordinate_current_fv
            for(size_t i=0U;i<DIM;i++)
                dxyz.Component(i, xyz2[i]-xyz1[i]);

            // save the distance to the neighboring FVs; Used in CalculateNodalGradient
            // for the RHS
            distance[ node_id ][id] = dxyz;

            // Accumulate LHS matrix ([sum_x2,sum_xy,sum_xz],[sum_xy,sum_y2,sum_yz],[sum_xz,sum_yz,sum_z2])
            sum_x2[ node_id ] += dxyz[0] * dxyz[0];
            sum_y2[ node_id ] += dxyz[1] * dxyz[1];
            sum_z2[ node_id ] += dxyz[2] * dxyz[2];
            sum_xy[ node_id ] += dxyz[0] * dxyz[1];
            sum_xz[ node_id ] += dxyz[0] * dxyz[2];
            sum_yz[ node_id ] += dxyz[1] * dxyz[2];
        }

        // calculate determinante
        det[  node_id ]   = ( sum_x2[ node_id ]*sum_y2[ node_id ]-sum_xy[ node_id ] * sum_xy[ node_id ] )*sum_z2 [ node_id ];
        det[  node_id ]  -= ( sum_x2[ node_id ]*sum_yz[ node_id ]-sum_xy[ node_id ] * sum_xz[ node_id ] )*sum_yz [ node_id ];
        det[  node_id ]  += ( sum_xy[ node_id ]*sum_yz[ node_id ]-sum_y2[ node_id ] * sum_xz[ node_id ] )*sum_xz [ node_id ];

        // the indexes of zero gradient components
        //0 - 0, 1- 1, 2 - 2, 0&1 - 3, 0&2 - 4, 1&2 - 5, else 6
        if(det[ node_id ] == 0.0){

            if(sum_x2[ node_id ] == 0.0){

                if(sum_y2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=3U;
                else if(sum_z2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=4U;
                else
                    zero_grad_index_[ node_id]=0U;
            }
            else if(sum_y2[ node_id ] == 0.0){

                if(sum_x2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=3U;
                else if(sum_z2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=5U;
                else
                    zero_grad_index_[ node_id]=1U;
            }
            else if(sum_z2[ node_id ] == 0.0){

                if(sum_x2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=4U;
                else if(sum_y2 [ node_id ] == 0.0)
                    zero_grad_index_[ node_id]=5U;
                else
                    zero_grad_index_[ node_id]=2U;
            }else
                zero_grad_index_[ node_id ]=6U;
        }

    }

    cout << "\nGenericNodePropertyGradient<dim>::CalculateLeastSquareSums:";
    cout << "\nLeast squares sums successfully calculated... ";
    cout << endl;

  }  // end CalculateLeastSquareSums






/**

Computes the gradient of a property (concentration, saturation, etc.),
by summing up the difference in the property value in a single finite
volume and all its neighbors times the according least square sum.
Those gradients are needed for the higher order finite volume methods
and are stored in STL vectors.

Note: For both the case where just one gradient is needed as for the one where
several ones are needed: The user has to set the Property Key to the respective value!

The method is implemented from the various FiniteVolume<fT, dim>Visitors and the STL
vectors that store the x and y components of the gradients can be accessed
through the interfaces GradientComponentX() and GradientComponentY()

*/
template<>
void GenericNodePropertyGradient<3U>::CalculateGenericNodalGradient()
{

    enum{DIM=3U};
    std::vector<Node<DIM>*>::const_iterator     cvit;
    std::vector<size_t>::iterator               nit;
    double64                                    val1, val2, dc;
    VectorVariable<DIM>                         grad(PLAIN,0.0);
    VectorVariable<DIM>                         dcxyz(PLAIN,0.0);
    size_t                                      id, node_id;
    VARIABLE_FLAG                               status;
    string                                      grad_name;

    for( std::vector<Node<DIM>*>::const_iterator
           cvit=sg_.Region("Model").NodesBegin();
           cvit!=sg_.Region("Model").NodesEnd(); cvit++){

           status = (*(*cvit)).Status( u_key );
           node_id =  (*(*cvit)).Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY ) {

                // get neighbor ids at every cv and read concentration
                val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.0;
                id  = 0;
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for ( nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    val2 = (sg_.Region("Model").N( (*nit) ))->Read (u_key);

                    dc   = (val2 -val1);

                    // computing temporary variables for least squares calculation; needed for RHS of LGS
                    // NOTE: id's are the same as in CalculateLeastSquares since list is gone thorugh in same order
                    for(size_t i=0U;i<DIM;i++)
                        dcxyz.Component(i, dcxyz[i]+dc*distance[ node_id ][id][i]);

                }

                if(det[ node_id]!=0.0){

                    grad(0)   = ( dcxyz[0]*sum_y2[ node_id ]-dcxyz[1]*sum_xy[ node_id ] )*sum_z2 [ node_id ];
                    grad(0)  -= ( dcxyz[0]*sum_yz[ node_id ]-dcxyz[1]*sum_xz[ node_id ] )*sum_yz [ node_id ];
                    grad(0)  += ( sum_xy[ node_id ]*sum_yz[ node_id ]-sum_y2[ node_id ]*sum_xz [ node_id ] )*dcxyz[2];
                    grad(0)  /= det[ node_id ];

                    grad(1)   = ( sum_x2[ node_id ]*dcxyz[1]-sum_xy[ node_id ]*dcxyz[0] )*sum_z2 [ node_id ];
                    grad(1)  -= ( sum_x2[ node_id ]*sum_yz[ node_id ]-sum_xy[ node_id ]*sum_xz [ node_id ] )*dcxyz[2];
                    grad(1)  += ( dcxyz[0]*sum_yz[ node_id ]-dcxyz[1]*sum_xz[ node_id ] )*sum_xz [ node_id ];
                    grad(1)  /= det[ node_id ];

                    grad(2)   = ( sum_x2[ node_id ]*sum_y2[ node_id ]-sum_xy[ node_id ]*sum_xy [ node_id ] )*dcxyz[2];
                    grad(2)  -= ( sum_x2[ node_id ]*dcxyz[1]-sum_xy[ node_id ]*dcxyz[0] )*sum_yz [ node_id ];
                    grad(2)  += ( sum_xy[ node_id ]*dcxyz[1]-sum_y2[ node_id ]*dcxyz[0] )*sum_xz [ node_id ];
                    grad(2)  /= det[ node_id ];


                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
                    if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 0U){

                    grad(0)   = 0.0;

                    grad(1)   = ( dcxyz[1]*sum_z2[ node_id ]-dcxyz[2]*sum_yz[ node_id ] );
                    grad(1)  /= ( sum_y2[ node_id ]*sum_z2[ node_id ]-sum_yz[ node_id ]*sum_yz[ node_id ] );

                    grad(2)   = ( sum_y2[ node_id ]*dcxyz[2]-sum_yz[ node_id ]*dcxyz[1] );
                    grad(2)  /= ( sum_y2[ node_id ]*sum_z2[ node_id ]-sum_yz[ node_id ]*sum_yz[ node_id ] );


                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;
                    if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 1U){

                    grad(0)   = ( dcxyz[0]*sum_z2[ node_id ]-dcxyz[2]*sum_xz[ node_id ] );
                    grad(0)  /= ( sum_x2[ node_id ]*sum_z2[ node_id ]-sum_xz[ node_id ]*sum_xz[ node_id ] );

                    grad(1)   = 0.0;

                    grad(2)   = ( sum_x2[ node_id ]*dcxyz[2]-sum_xz[ node_id ]*dcxyz[0] );
                    grad(2)  /= ( sum_x2[ node_id ]*sum_z2[ node_id ]-sum_xz[ node_id ]*sum_xz[ node_id ] );

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                    if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 2U){

                    grad(0)   = ( dcxyz[0]*sum_y2[ node_id ]-dcxyz[1]*sum_xy[ node_id ] );
                    grad(0)  /= ( sum_x2[ node_id ]*sum_y2[ node_id ]-sum_xy[ node_id ]*sum_xy[ node_id ] );

                    grad(1)   = ( sum_x2[ node_id ]*dcxyz[1]-sum_xy[ node_id ]*dcxyz[0] );
                    grad(1)  /= ( sum_x2[ node_id ]*sum_y2[ node_id ]-sum_xy[ node_id ]*sum_xy[ node_id ] );

                    grad(2)   = 0.0;

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 3U){

                    grad(0)   = 0.0;

                    grad(1)   = 0.0;

                    grad(2)   = dcxyz[2]/sum_z2[ node_id ];

                    if ( fabs( grad(2) ) < tolerance ) grad(2) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 4U){

                    grad(0)   = 0.0;

                    grad(1)   = dcxyz[1]/sum_y2[ node_id ];

                    grad(2)   = 0.0;

                    if ( fabs( grad(1) ) < tolerance ) grad(1) = 0.0;

                }else if(zero_grad_index_[ node_id ] == 5U){

                    grad(0)   = dcxyz[0]/sum_x2[ node_id ];

                    grad(1)   = 0.0;

                    grad(2)   = 0.0;

                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;

                }else
                    grad = 0.0;
        }
        // gradient is zero if nodes are flagged DIRICH or NEUMANN
        else grad = 0.0;

        // store the gradient
        // ------------------

        grad_name = sg_.Database().Name( u_key );
        grad_name = grad_name + " gradient";
        (*(*cvit)).Store( sg_.Database().StorageKey( grad_name.c_str() ), grad );
    }
} // end CalculateNodalGradient


template class GenericNodePropertyGradient<1U>;
template class GenericNodePropertyGradient<2U>;
template class GenericNodePropertyGradient<3U>;

} // end namespace













































