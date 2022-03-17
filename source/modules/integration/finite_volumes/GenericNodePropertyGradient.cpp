#include "GenericNodePropertyGradient.h"
#include "Region.h"

using namespace std;

namespace csmp {

/**
    Sets the calculator up for computation on the entire model
    (a throughgoing node numbering is assumed.
*/
template<uint32_t dim>
GenericNodePropertyGradient<dim>::GenericNodePropertyGradient( Model<dim>& sg, const char* region, const char* prop )
  :  gref_(sg.Region(region)),
     gradient( 1U, NULL ),
     sum_x2( gref_.Nodes(), 0.0 ),
     sum_y2( gref_.Nodes(), 0.0 ),
     sum_z2( gref_.Nodes(), 0.0 ),
     sum_xy( gref_.Nodes(), 0.0 ),
     sum_xz( gref_.Nodes(), 0.0 ),
     sum_yz( gref_.Nodes(), 0.0 ),
     det( gref_.Nodes() ),
     zero_grad_index_( gref_.Nodes(), 6U ),
     distance( gref_.Nodes() ),
     center_of_mass_( gref_.Nodes() ),
     tolerance(1.0e-12),
     inner_(1U, VectorVariable<dim>(PLAIN,0.0) ),
     middle_( 1U, inner_ ),
     distance_facet_FVBarycenter_( gref_.Elements() ),
     neighbors_( gref_.Nodes() ),
     mass_center( sg, region, "mass center", VECTOR, NODE ), //PropertyHandle for mass center
     u_key(sg.Database().StorageKey( prop )),
     grad_key(sg.Database().StorageKey( (string(prop) + " gradient").c_str() )),
     mctr_key(sg.Database().StorageKey("mass center"))
  {
    // initialize gradient PropHandle and set prop key:
    // ------------------------------------------------  
    string   gradient_name;
    gradient_name  =  prop;
    gradient_name += " gradient";
    gradient[0] = new PropertyHandle<dim>( sg, gradient_name.c_str(), VECTOR, NODE );
               
     
    // --------------------------------------------------
    // initialization and calculation of mass centers...
    // --------------------------------------------------
    std::pair<VectorVariable<dim>, double>   init_pair;
    init_pair.first = VectorVariable<dim>( PLAIN, 0.0);
    init_pair.second = 0.0;

    fill(  center_of_mass_.begin(), center_of_mass_.end(), init_pair );

    CalculateCenterOfMass();

    // ---------------------------------------------------------------------------
    // initialization and calculation of distances facet centers-barycenters
    // ---------------------------------------------------------------------------

    // init of vector
    // ====================================================
    distance_facet_FVBarycenter_.resize( gref_.Elements() ); // resize outer vector
    size_t length = distance_facet_FVBarycenter_.size();

    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.ElementsBegin();
          eit!=gref_.ElementsEnd(); eit++ ){
          //resize middle vector
          distance_facet_FVBarycenter_[ (*eit)->Idx() ].resize( (*eit)->FV()->Facets() );

          // loop over facets
          // -------------------
          for( size_t fi=0U; fi < (*eit)->FV()->Facets(); fi++ ){

                 // resize inner vector
                distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].resize( (*eit)->Nodes() );
                length = distance_facet_FVBarycenter_[ (*eit)->Idx()][ fi ].size();

                fill( distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].begin(),
                      distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].end(),
                      VectorVariable<dim>( PLAIN, 0.0)       );
          }
    }

    // ====================================================

    // ... and distances:
    CalculateDistanceFacetFVBary( distance_facet_FVBarycenter_ );

    // -----------------------------
    // neighbours and least squares
    // -----------------------------


    // initialize vector with global node id's of neighbors:
    // =====================================================
    size_t               global_neighb_el_id, current_n_id;
    std::vector<size_t>  ids;

    // loop over all nodes/fv's
    // --------------------------
    for( auto cvit=gref_.NodesBegin(); cvit!=gref_.NodesEnd(); cvit++){

         current_n_id = (*cvit)->Idx();
         // loop over parents
         for(  auto p = 0; p< (*cvit)->Parents(); p++ ){
             // get the global parent id:
             global_neighb_el_id = (*cvit)->Parent( p)->Idx();

             if(global_neighb_el_id<gref_.Elements()){
                 // get the corresponding element:
                 const Element<dim>* current_el = gref_.E( global_neighb_el_id );
                 // ...and the global node_i's of that element
                 ids.clear();
                 for(auto i=0; i<current_el->Nodes(); i++)
                   if(current_n_id!=current_el->N(i)->Idx())
                     ids.push_back(current_el->N(i)->Idx());

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
template<uint32_t dim>
GenericNodePropertyGradient<dim>::GenericNodePropertyGradient( Model<dim>& sg, const char* region, vector<char*> properties )
  : gref_(sg.Region(region)),
    gradient( 1U, NULL ),
    sum_x2( gref_.Nodes(), 0.0 ),
    sum_y2( gref_.Nodes(), 0.0 ),
    sum_z2( gref_.Nodes(), 0.0 ),
    sum_xy( gref_.Nodes(), 0.0 ),
    sum_xz( gref_.Nodes(), 0.0 ),
    sum_yz( gref_.Nodes(), 0.0 ),
    det( gref_.Nodes() ),
    zero_grad_index_( gref_.Nodes(), 6U ),
    distance( gref_.Nodes() ),
    center_of_mass_( gref_.Nodes() ),
    tolerance(1.0e-15),
    inner_(1U, VectorVariable<dim>() ),
    middle_( 1U, inner_ ),
    distance_facet_FVBarycenter_( gref_.Elements() ),
    neighbors_( gref_.Nodes() ),
    mass_center( sg, region, "mass center", VECTOR, NODE ), //PropertyHandle for mass center
     u_key(sg.Database().StorageKey( properties[0] )),
     grad_key(sg.Database().StorageKey( (string( properties[0] ) + " gradient").c_str() )),
     mctr_key(sg.Database().StorageKey("mass center"))

  {
    // --------------------------------------------------
    // initialization of gradients 
    // --------------------------------------------------
    vector<string >  gradient_strings(properties.size()); 
    string   str2;
    str2 = " gradient";

    for( auto i = 0; i < properties.size(); i++ ){
          // construct the name of the property's gradient:
          gradient_strings[i] = properties[i];
          gradient_strings[i] = gradient_strings[i] + str2;

          // property handle
          gradient[i] = new PropertyHandle<dim>( sg, gradient_strings[i].c_str() , VECTOR, NODE );
    }
     
    // --------------------------------------------------
    // initialization and calculation of mass centers... 
    // --------------------------------------------------
    std::pair<VectorVariable<dim>, double>   init_pair;
    init_pair.first = VectorVariable<dim>( PLAIN, 0.0);
    init_pair.second = 0.0;
   
    fill(  center_of_mass_.begin(), center_of_mass_.end(), init_pair );
       
    CalculateCenterOfMass();
    
    // ---------------------------------------------------------------------------
    // initialization and calculation of distances facet centers-barycenters
    // ---------------------------------------------------------------------------
    
    // init of vector 
    // ====================================================
    
    distance_facet_FVBarycenter_.resize( gref_.Elements() ); // resize outer vector
    size_t length_;
    length_ =  distance_facet_FVBarycenter_.size();
     
    for ( typename vector<Element<dim>*>::const_iterator
          eit=gref_.ElementsBegin();
          eit!=gref_.ElementsEnd(); eit++ ){
          
          cout<< (*eit)->Idx() <<endl;

          //resize middle vector
          distance_facet_FVBarycenter_[ (*eit)->Idx() ].resize( (*eit)->FV()->Facets() );

          // loop over facets
          // -------------------
          for( size_t fi=0U; fi < (*eit)->FV()->Facets(); fi++ ){

                 // resize inner vector
                distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].resize( (*eit)->Nodes() );
                length_ = distance_facet_FVBarycenter_[ (*eit)->Idx()][ fi ].size();

                fill( distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].begin(),
                      distance_facet_FVBarycenter_[ (*eit)->Idx() ][ fi ].end(),
                      VectorVariable<dim>( PLAIN, 0.)       );
          }
    }
    
    // ====================================================   
   
    // ... and distances:
    CalculateDistanceFacetFVBary( distance_facet_FVBarycenter_ );
    
    // -----------------------------
    // neighbours and least squares
    // -----------------------------
     
    
    // initialize vector with global node id's of neighbors:
    // =====================================================
    std::vector<size_t>  ids;
    
    // loop over all nodes/fv's
    // --------------------------
    for( auto cvit=gref_.NodesBegin(); cvit!=gref_.NodesEnd(); cvit++ )
      {
         auto current_n_id = (*cvit)->Idx();
         // loop over parents
         for( auto p = 0; p< (*cvit)->Parents(); p++ ){
             // get the global parent id:
             auto global_neighb_el_id = (*cvit)->Parent(p)->Idx();
             if(global_neighb_el_id<gref_.Elements()){
                 // get the corresponding element:
                 const Element<dim>* current_el = gref_.E( global_neighb_el_id );
                 // ...and the global node_i's of that element
                 ids.clear();
                 for( auto i=0;i<current_el->Nodes(); i++ )
                     if(current_n_id!=current_el->N(i)->Idx()) ids.push_back(current_el->N(i)->Idx());

                 // insert into vector without duplicates:
                 PushBackAvoidDuplicate( neighbors_[ current_n_id ], ids );
             }
         }
    }
    // =====================================================
     
    CalculateGenericLeastSquareSums();
}



template<uint32_t dim>
void GenericNodePropertyGradient<dim>::SetPropertyKey( csmp::Index& key )
  {
    u_key = key;
  
  }

  
template<uint32_t dim>
GenericNodePropertyGradient<dim>::~GenericNodePropertyGradient()
{     // getting rid of the property handles
     for ( auto& it : gradient ) delete it;
}
  
  
  
  

/**

Calculates the center of mass of generic finite volumes in 2d.

*/
template<uint32_t dim>
void GenericNodePropertyGradient<dim>::CalculateCenterOfMass()
{
 size_t  glob_n_id;
 double volume_;
 std::vector<double>   current_bc; // stencil barycenter in global coord's
 std::vector<size_t> ids;
 VectorVariable<dim> temp_;
 Point<dim> tmp_p;

 // loop over stencils
 // ------------------
 for ( auto eit=gref_.ElementsBegin();
       eit!=gref_.ElementsEnd(); eit++ ){

        ids.resize( (*eit)->Nodes() );
        // get vector of global id's
        for( uint32_t i{0}; i<(*eit)->Nodes(); i++ )
            ids[i]=(*eit)->N(i)->Idx();

        // loop over sectors
        // -------------------
        for( auto i{0}; i < (*eit)->FV()->Sectors(); i++ ){

              //get the global node id for the current segment
              glob_n_id = ids[ i ];

              if( gref_.N( glob_n_id)->AtBoundary() ){

                  tmp_p = gref_.N( glob_n_id)->Coordinate();
                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for( auto k=0;k<dim;k++)
                      temp_.Component(k,tmp_p.Coordinates()[k]);
                  center_of_mass_[ glob_n_id ].first  = temp_;
                  center_of_mass_[ glob_n_id ].second = 1.;

              }else{

                  //get the barycenter of current segment
                  ConvertToGlobalCoordinates( *(*eit), (*eit)->FV()->SectorIntegrationPoint( i, 0U), current_bc );

                  //get volume of current sector
                  volume_ = (*eit)->SectorVolume( i );

                  temp_ = VectorVariable<dim>( PLAIN,  0.0);
                  for( uint32_t k=0; k<dim; k++ )
                      temp_.Component(k,current_bc[k] * volume_);

                  center_of_mass_[ glob_n_id ].first += temp_;
                  center_of_mass_[ glob_n_id ].second += volume_;

             }
        }
     }

     for( auto i = 0U;  i< gref_.Nodes(); i++ )
     {
        //  center_of_mass_[ i ].first.Out();
        //  calculate center: sum_i(x_i * A_i) / sum_i(A_i)
        center_of_mass_[ i ].first *=  1. / center_of_mass_[ i ].second;

        // write back to node:
        gref_.N( i )->Store( mctr_key, center_of_mass_[ i ].first );
     }

}





template<uint32_t dim>
void GenericNodePropertyGradient<dim>::CalculateDistanceFacetFVBary( vector<vector<vector<VectorVariable<dim> > > >& d_facet_FVBarycenter )
{
 for ( auto eit=gref_.ElementsBegin();
       eit!=gref_.ElementsEnd(); eit++ )
   {
        // loop over facets
        // -------------------
        for( auto fi=0U; fi < (*eit)->FV()->Facets(); fi++ ){

              // get facet barycenter in global coordinates:
              vector<double>  facet_bc(dim);  // global coordinates of facet barycenter
              ConvertToGlobalCoordinates( *(*eit), (*eit)->FV()->FacetIntegrationPoint( fi, 0U ),  facet_bc);

              // get local node id's
              uint32_t inside_node_, outside_node_;
              (*eit)->FV()->FacetEdgeNodes( fi, inside_node_, outside_node_ );


              // get mass center for FV of inside_node_:
              VectorVariable<dim> mass_ctr;
              GenericCenterOfMass(  (*eit)->N(inside_node_)->Idx(), mass_ctr );

              VectorVariable<dim> face_dist_inside_node( PLAIN,  0.0);
              for(auto k=0;k<dim;k++)
                  face_dist_inside_node.Component(k,facet_bc[k] - mass_ctr[k] );

              d_facet_FVBarycenter[ (*eit)->Idx() ][ fi ][ inside_node_ ] =face_dist_inside_node;


              // get mass center for FV of outside_node_:
              GenericCenterOfMass(  (*eit)->N(outside_node_)->Idx(), mass_ctr );

              // calculate distance vector:
              VectorVariable<dim> face_dist_outside_node( PLAIN,  0.0);
              for(auto k=0;k<dim;k++)
                  face_dist_outside_node.Component(k,facet_bc[k] - mass_ctr[k] );

              d_facet_FVBarycenter[ (*eit)->Idx() ][ fi ][ outside_node_ ] = face_dist_outside_node;

        }
 }
// bool hallo;

//  return_d_facet_FVBarycenter =  d_facet_FVBarycenter;
}




/**

Inserts all components of possible_new_entries into old_vector if they
are not there already.

@todo logic?

*/
template<uint32_t dim>
void GenericNodePropertyGradient<dim>::PushBackAvoidDuplicate( vector<size_t>& old_vector,
                                                               const vector<size_t>& new_vector )
{
    bool already_exists(false);

    for( auto new_ = new_vector.begin();
              new_ != new_vector.end(); new_++ ){
        for( auto old_ = old_vector.begin(); old_ != old_vector.end(); old_++ ){
             if( *new_ == *old_ ) already_exists = true;
        }
        // current element does not exist yet, insert:
        if( !already_exists ) old_vector.push_back( *new_ );
        already_exists = false;
    }
}





template<uint32_t dim>
void GenericNodePropertyGradient<dim>::ConvertToGlobalCoordinates( Element<dim>& el, const Point<dim>& local_c_point, std::vector<double>& global_c )
{


  std::vector<double> temp(el.Nodes()); //has the local interp. function values
  std::vector<double> local_c(local_c_point.Coordinates());

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
  for (auto i = 0; i<el.Nodes(); i++)
      for (auto j = 0; j<dim; j++){
        global_c[j] += el.FE()->XY(i,j) * temp[i];
  }
}



/**

Returns the calculated centers of mass of generic finite volumes in 2d.

@param global_node_id  the node id of the corresponding FV

@param mass_center_out        the Vector onto which the xyz-coordinates of the mass center
                         of the FV that belongs to the the node are written
*/
template<uint32_t dim>
void GenericNodePropertyGradient<dim>::GenericCenterOfMass( size_t global_node_id, VectorVariable<dim>& mass_center_out ) const
{
 mass_center_out =  center_of_mass_[global_node_id].first;
}

template<uint32_t dim>
void GenericNodePropertyGradient<dim>::GenericDistanceFacetFVBarycenter( size_t global_el_id,
                                                                         uint32_t local_facet_id,
                                                                         uint32_t local_node_id,
                                                                         VectorVariable<dim>& fdistance ) const
{
   fdistance =   distance_facet_FVBarycenter_[ global_el_id ][ local_facet_id ][ local_node_id ] ;
}







/**

Returns the storage required by the GenericNodePropertyGradient object

*/
template<uint32_t dim>
double GenericNodePropertyGradient<dim>::SizeOf() const
{
  double storage;

  storage  = sizeof( *this );
  storage += static_cast<double>(sum_xy.size()) * sizeof(double);
  storage += static_cast<double>(sum_x2.size()) * sizeof(double);
  storage += static_cast<double>(sum_y2.size()) * sizeof(double);
  storage += static_cast<double>(det.size())    * sizeof(double);

  for ( size_t i=0; i<distance.size(); i++ )
      storage += static_cast<double>(distance[i].size()) * 2.0 * sizeof(double);

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
    VectorVariable<DIM> xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);

    for( auto cvit=gref_.NodesBegin(); cvit!=gref_.NodesEnd(); cvit++){
        const size_t node_id = (*(*cvit)).Idx() ;
        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );
        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient 
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );

        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        size_t id{0};
        for ( auto nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {
        
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
    VectorVariable<DIM>  grad(PLAIN,0.0);
    VectorVariable<DIM>  dcxyz(PLAIN,0.0);

    for( std::vector<Node<DIM>*>::const_iterator
           cvit=gref_.NodesBegin();
           cvit!=gref_.NodesEnd(); cvit++){

           VARIABLE_FLAG status  = (*(*cvit)).Status( u_key );
           size_t        node_id = (*cvit)->Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY) {
        
                // get neighbor ids at every cv and read concentration
                double val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.;
                size_t id{0};
                
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for ( auto nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    double val2 = gref_.N( (*nit) )->Read (u_key);
                    double dc   = (val2 -val1);

                    // computing temporary variables for least squares calculation; needed for RHS of LGS
                    // NOTE: id's are the same as in CalculateLeastSquares since list is gone thorugh in same order
                    dcxyz.Component( 0, dcxyz[0]+dc*distance[ node_id ][id][0]);
                }

                if(det[ node_id]!=0.0){
                    grad(0)  = dcxyz[0]/det[ node_id ];
                    if ( fabs( grad(0) ) < tolerance ) grad(0) = 0.0;
                }else grad=0.0;

        }
        // gradient is zero if nodes are flagged DIRICH or NEUMANN
        else grad = 0.;
        
        // store the gradient
        // ------------------
        (*cvit)->Store( grad_key, grad );
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
    VectorVariable<DIM>  xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);

    for( vector<Node<DIM>*>::const_iterator cvit=gref_.NodesBegin();
         cvit!=gref_.NodesEnd(); cvit++){

        size_t node_id = (*cvit)->Idx();

        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );

        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );

        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        size_t id{0};
        for ( auto nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {

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
    VectorVariable<DIM>  grad(PLAIN,0.0);
    VectorVariable<DIM>  dcxyz(PLAIN,0.0);

    for( std::vector<Node<DIM>*>::const_iterator
           cvit=gref_.NodesBegin();
           cvit!=gref_.NodesEnd(); cvit++){

           VARIABLE_FLAG status  = (*(*cvit)).Status( u_key );
           size_t        node_id = (*(*cvit)).Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY) {

                // get neighbor ids at every cv and read concentration
                double val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.0;
                size_t id{0};
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for (auto nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    double val2 = gref_.N( (*nit) )->Read (u_key);
                    double dc   = (val2 -val1);

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

        (*(*cvit)).Store( grad_key, grad );
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
    VectorVariable<DIM>  xyz1(PLAIN,0.0), xyz2(PLAIN,0.0), dxyz(PLAIN,0.0);

    for( vector<Node<DIM>*>::const_iterator cvit=gref_.NodesBegin();
         cvit!=gref_.NodesEnd(); cvit++){

        const size_t node_id = (*cvit)->Idx();

        // get neighbor ids at every segment and fv barycenter
        GenericCenterOfMass( node_id, xyz1 );

        // stores the distance to the neighboring FVs; Used in CalculateNodalGradient
        // for the RHS

        distance[ node_id ].resize( neighbors_[ node_id ].size() );
        
        // loop over all neighbor ids, get coordinate values, sum up
        // -----------------------------------------------------------
        size_t id = 0;
        for ( auto nit = neighbors_[ node_id].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ) {

            GenericCenterOfMass( (*nit), xyz2 );

            // subtract x,y,z-coordinate_neighbour_center - x,y,z-coordinate_current_fv
            for(auto i{0};i<DIM;i++)
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
    VectorVariable<DIM>  grad(PLAIN,0.0);
    VectorVariable<DIM>  dcxyz(PLAIN,0.0);

    for( vector<Node<DIM>*>::const_iterator
         cvit=gref_.NodesBegin();
         cvit!=gref_.NodesEnd(); cvit++){

           VARIABLE_FLAG status = (*(*cvit)).Status( u_key );
           const size_t node_id = (*(*cvit)).Idx();

           // gradient is zero if nodes are flagged DIRICH or NEUMANN
           if ( status == PLAIN || status == ANY ) {

                // get neighbor ids at every cv and read concentration
                double val1 = (*(*cvit)).Read(u_key );

                // set temporary variables to zero
                dcxyz = 0.0;
                size_t id  = 0;
                // loop over all neighbor ids, get coordinate values and concentrations, sum up
                // -----------------------------------------------------------------------------
                for ( auto nit = neighbors_[ node_id ].begin(); nit != neighbors_[ node_id ].end(); nit++,id++ ){
                    // read at neighbour node!
                    double val2 = gref_.N(*nit)->Read(u_key);
                    double dc   = (val2 -val1);

                    // computing temporary variables for least squares calculation; needed for RHS of LGS
                    // NOTE: id's are the same as in CalculateLeastSquares since list is gone thorugh in same order
                    for(auto i{0};i<DIM;i++)
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
        (*cvit)->Store( grad_key, grad );
    }
} // end CalculateNodalGradient


template class GenericNodePropertyGradient<1U>;
template class GenericNodePropertyGradient<2U>;
template class GenericNodePropertyGradient<3U>;

} // end namespace













































