#include "LineElementMesher.h"
#include "Box.h"
#include "IsoparametricLinearLineElement.h"

#include "CSMP_mathUtilities.h"

#include "TextFileIO.h"

using namespace std;

namespace csmp {

template<size_t dim>
LineElementMesher<dim>::LineElementMesher()
  {
  }


/// create an uniform mesh
template<size_t dim>
void LineElementMesher<dim>
::BuildUniformMesh( VSet<dim>& vset,
                    double64 length,
                    size_t n_elements,
                    const std::vector<double64>& splitnode_coordinates,
                    const Point<dim>& origin,
                    const Point<dim>& destination )
 {
   // 1.0 Distribute points
   UniformMesh( vset,
                length, n_elements );

   // 2.0 Complete mesh
   CompleteMesh( vset,
                 splitnode_coordinates,
                 origin, destination );

 } // end Initialize







/// create a mesh size gradient
template<size_t dim>
void LineElementMesher<dim>
::BuildRefinedMesh( VSet<dim>& vset,
                    double64 length,
                    double64 dx_min,
                    double64 dx_max,
                    double64 width_of_transition_zone,
                    MeshDensity* density,
                    const std::vector<double64>& splitnode_coordinates,
                    const Point<dim>& origin,
                    const Point<dim>& destination )
{
   // 1.0 Distribute points
   RefinedMesh( vset,
                length, dx_min, dx_max, width_of_transition_zone, density );

   // 2.0 Complete mesh
   CompleteMesh( vset,
                 splitnode_coordinates,
                 origin, destination );

} // end Initialize (exponential refinement)

template<size_t dim>
void LineElementMesher<dim>
::BuildCustomMesh( VSet<dim>& vset,
                   const std::vector<double64>& node_coordinates,
                   const std::vector<double64>& splitnode_coordinates,
                   const Point<dim>& origin,
                   const Point<dim>& destination )
{
    // 1.0 Distribute points
    CustomMesh( vset,
                node_coordinates );

    // 2.0 Complete mesh
    CompleteMesh( vset,
                  splitnode_coordinates,
                  origin, destination );

} // end Initialize

template<size_t dim>
void LineElementMesher<dim>
::BuildCustomMesh( VSet<dim>& vset,
                   const std::vector<Point<dim> >& node_coordinates,
                   const std::vector<Point<dim> >& splitnode_coordinates,
                   const Point<dim>& origin,
                   const Point<dim>& destination )
{
    // 1.0 Distribute points
    CustomMesh( vset,
                node_coordinates );

    // 2.0 Complete mesh
    CompleteMesh( vset,
                  splitnode_coordinates,
                  origin, destination );

} // end Initialize









// WORK AROUND ON ALREADY EXISTING MESH


/** Modifies the mesh based on the new coordinates of corner points
    without changing the connectivity of previously created mesh

  @author R. Manasipov
*/

template<size_t dim>
void LineElementMesher<dim>::AssignCornerPoints( VSet<dim>& vset,
                                        const Point<dim>& origin,
                                        const Point<dim>& destination )
{
    if( vset.Vertices() > 0 )
    {
        const size_t  n_vertices( vset.Vertices() );

        // direction of new line model
        Point<dim> direction( destination - origin );
        double64   length( direction.Length() );
        direction.NormalizeLengthTo( 1.0 );

        // identify the direction of old line model
        std::set<Point<dim> > nodes;
        Point<dim> p( 0.0 );
        for( size_t i =0U; i< n_vertices; i++ ) {
            for( size_t dimension=0U; dimension<dim; dimension++ )
                p[ dimension ] = vset.P(dimension,i);
            nodes.insert( p );
        }
        Point<dim> old_origin( *nodes.begin()  );
        Point<dim> old_destination( *nodes.rbegin() );
        Point<dim> old_direction( old_destination - old_origin );
        double64   old_length( old_direction.Length() );
        old_direction.NormalizeLengthTo( 1.0 );

        // VSet - assigning new coordinates to all points
        double64   distance( 0.0 );
        double64   length_factor( length / old_length );
        for( size_t i=0U; i<n_vertices; i++ ) {
             for( size_t dimension=0U; dimension<dim; dimension++ )
                 p[ dimension ] = vset.P( dimension, i );
             distance  = (p - old_origin).Length();
             distance *= length_factor;
             p  = origin;
             p += distance*direction;

             // assign new coordinate
             for( size_t dimension=0U; dimension<dim; dimension++ )
                 vset.P( dimension, i, p[ dimension ] );
          }
    }

} // end Initialize

template<size_t dim>
void LineElementMesher<dim>::InsertSplitNodes( VSet<dim>& vset,
                                      const std::vector<double64>& splitnodes_x )
{
    if( !splitnodes_x.empty() && ( vset.Vertices() > 0 ) )
    {
        // convert splitnode coordinates along x axis to generic Point format
        std::set<Point<dim> > splitnodes;
        Point<dim> p( 0.0 );
        for( std::vector<double64>::const_iterator it = splitnodes_x.begin(); it != splitnodes_x.end(); it++ )
        {
            p[ 0U ] = (*it);
            for( size_t  dimension = 1U; dimension<dim; dimension++)
                p[ dimension ] = 0.0;
            splitnodes.insert( p );
        }

        InsertSplitNodes( vset, splitnodes );
    }

} // InsertSplitNodes


template<size_t dim>
void LineElementMesher<dim>::InsertSplitNodes( VSet<dim>& vset,
                                      std::set<Point<dim> >& splitnodes )
{
    // add splitnodes and reestablish connectivety
    EstablishConnectivity( vset, splitnodes );

} // InsertSplitNodes

template<size_t dim>
void LineElementMesher<dim>::CompleteMesh( VSet<dim>& vset,
                                  const std::vector<double64>& splitnode_coordinates,
                                  const Point<dim>& origin,
                                  const Point<dim>& destination )
{
    if( splitnode_coordinates.empty() )
        EstablishConnectivity( vset );
    else
        InsertSplitNodes( vset, splitnode_coordinates );

    if( origin != destination )
        AssignCornerPoints( vset, origin, destination );
}

template<size_t dim>
void LineElementMesher<dim>::CompleteMesh( VSet<dim>& vset,
                                  const std::vector<Point<dim> >& splitnode_coordinates,
                                  const Point<dim>& origin,
                                  const Point<dim>& destination )
{
    if( splitnode_coordinates.empty() )
        EstablishConnectivity( vset );
    else
    {
        std::set<Point<dim> > splitnodes( splitnode_coordinates.begin(), splitnode_coordinates.end() );
        InsertSplitNodes( vset, splitnodes );
    }

    if( origin != destination )
        AssignCornerPoints( vset, origin, destination );
}















// BUILDING BLOCKS

/// create an uniform mesh
template<size_t dim>
void LineElementMesher<dim>::UniformMesh( VSet<dim>& vset,
                                          double64 length,
                                          size_t n_elements )
 {
   const size_t  n_vertices( n_elements + 1U );
   const double  dx( length/static_cast<double>( n_elements ) );

   //1.0 VSet - refitting VSet
   vset.Resize( IsoparametricLinearLineElement().Nodes(),
                IsoparametricLinearLineElement().Neighbors(),
                IsoparametricLinearLineElement().ElementType(),
                n_vertices,
                n_elements );

   //1.1 VSet - assigning point x-coordinates
   for ( size_t i=0U; i<n_vertices; i++ ){
       vset.P( 0U, i, dx * i );
       for( size_t dimension=1U; dimension<dim; dimension++ )
           vset.P( dimension, i, 0.0 );
   }

 } // end UniformMesh


/// create a refined mesh
template<size_t dim>
void LineElementMesher<dim>
::RefinedMesh( VSet<dim>& vset,
               double64 length,
               double64 dx_min,
               double64 dx_max,
               double64 width_of_transition_zone,
               MeshDensity* density
              )
{
   // Initial checks
   if ( dx_min > dx_max ) {
       std::cerr << "\nLineElementMesher::RefinedMesh: ";
       std::cerr << "Min element size (arg1) must not exceed max-size (arg2)." << endl;
       return;
     }
   if ( dx_max > width_of_transition_zone ) {
       std::cerr << "\nLineElementMesher::RefinedMesh: ";
       std::cerr << "Max element size (arg2) should not exceed the width of variably ";
       std::cerr << "refined mesh region (arg3)." << endl;
       return;
     }
   if ( dx_max > length ) {
       std::cerr << "\nLineElementMesher::RefinedMesh: ";
       std::cerr << "Max element size (arg2) must not exceed model size (arg4)." << endl;
       return;
     }

   // 0. creating a list of nodes and element lenght
   // ------------------------------------------------------------------
   std::vector<double64>  element_length;
   double64 x0(0.0);
   double64 delta_x(dx_min);
   double64 x(dx_min);
   double64 scale_factor(dx_max);
   double64 x_scale_factor( 1.0/width_of_transition_zone );
   bool     first_call(true);

   while( x < length ) {

        // the last element
        if ( x + delta_x > length ) {
             delta_x = length - x;
             if ( delta_x < dx_min ) {
                  std::cerr <<"\nLineElementMesher::Initialize: last element is short compare to prescribed minimum: ";
                  std::cerr << delta_x <<" vs. "<< dx_min << endl;
                  std::cerr << "\nThus it will be equally distributed over the all elements, which might slighly enlarge the elements!";
               }
             const size_t elements( element_length.size() );
             delta_x /= static_cast<double64>( elements );
             for( size_t i = 0; i<elements; i++ )
                element_length[ i ] += delta_x;
             break;
          }

        // storing the coordinate
        element_length.push_back( delta_x );

        // calculating the new finite-element length
        // if this is the first call, an x-value must be found for which
        // the error function yields a delta_x that is >= min_size
        if ( first_call ) {
             delta_x  = scale_factor * (*density)( ( x0 + x ) * x_scale_factor );
             while ( delta_x < dx_min ) {
                   x0     += dx_min;
                  delta_x  = scale_factor * (*density)( ( x0 + x ) * x_scale_factor );
               }
             first_call = false;
          }
        else
            delta_x = scale_factor * (*density)( ( x0 + x ) * x_scale_factor );

        x += delta_x;
     }

    const size_t  n_vertices( element_length.size() + 1U );
    const size_t  n_elements( element_length.size()      );

    //1.0 VSet - refitting VSet
    vset.Resize( IsoparametricLinearLineElement().Nodes(),
                 IsoparametricLinearLineElement().Neighbors(),
                 IsoparametricLinearLineElement().ElementType(),
                 n_vertices,
                 n_elements );

    //1.1 VSet - assigning point x-coordinates
    size_t  n(0U);
    x = 0.;
    for ( std::vector<double64>::const_iterator
          it=element_length.begin(); it!=element_length.end(); it++ )
      {
         vset.P( 0U, n, x );
         for( size_t dimension=1U; dimension<dim; dimension++ )
             vset.P( dimension, n, 0.0 );
         n++;
         x += (*it);
      }
    // the last node
    vset.P( 0U, n, length );
    for( size_t dimension=1U; dimension<dim; dimension++ )
        vset.P( dimension, n, 0.0 );
}

template<size_t dim>
void LineElementMesher<dim>::CustomMesh( VSet<dim>& vset,
                                const std::vector<double64>& node_coordinates )
{
    const size_t  n_vertices( node_coordinates.size() );
    const size_t  n_elements( n_vertices - 1 );

    //1.0 VSet - refitting VSet
    vset.Resize( IsoparametricLinearLineElement().Nodes(),
                 IsoparametricLinearLineElement().Neighbors(),
                 IsoparametricLinearLineElement().ElementType(),
                 n_vertices,
                 n_elements );

    //1.1 VSet - assigning point x-coordinates
    size_t counter(0U);
    std::vector<double64>::const_iterator node_coordinate( node_coordinates.begin());
    std::vector<double64>::const_iterator last_coordinate( node_coordinates.end() );
    while( node_coordinate != last_coordinate ) {
         vset.P( 0U, counter, *node_coordinate );
         for( size_t dimension=1U; dimension<dim; dimension++ )
             vset.P( dimension, counter, 0.0 );
         counter++;
         node_coordinate++;
      }

} // end CustomMesh

template<size_t dim>
void LineElementMesher<dim>::CustomMesh( VSet<dim>& vset,
                                const std::vector<Point<dim> >& node_coordinates )
{
    const size_t  n_vertices( node_coordinates.size() );
    const size_t  n_elements( n_vertices - 1 );

    //1.0 VSet - refitting VSet
    vset.Resize( IsoparametricLinearLineElement().Nodes(),
                 IsoparametricLinearLineElement().Neighbors(),
                 IsoparametricLinearLineElement().ElementType(),
                 n_vertices,
                 n_elements );

    //1.1 VSet - assigning point x-coordinates
    size_t counter(0U);
    typename std::vector<Point<dim> >::const_iterator node_coordinate( node_coordinates.begin());
    typename std::vector<Point<dim> >::const_iterator last_coordinate( node_coordinates.end() );
    while( node_coordinate != last_coordinate ) {
         for( size_t dimension=0U; dimension<dim; dimension++ )
             vset.P( dimension, counter, (*node_coordinate)[ dimension ] );
         counter++;
         node_coordinate++;
      }

} // end CustomMesh



/**
    Creates line element mesh
    where the corners = endpoints are flagged:

    First node = CNR1
    Last  node = CNR2

    in correspondance with the labeling of box-shaped models.

    @attention for a 1D model which is oriented vertically the assignment of the
    second corner is not strictly correct.
*/


template<size_t dim>
void LineElementMesher<dim>::EstablishConnectivity( VSet<dim>& vset )
{
    if( vset.Vertices() > 0 )
    {
        const size_t n_vertices( vset.Vertices()      );
        const size_t n_elements( vset.Vertices() - 1U );

        //1.0 PList - establish
        deque<vector<size_t> >  plist( n_elements, vector<size_t>(2) );

        //1.1 PList - assigning content
        for ( size_t i=0; i<n_elements; ++i )
          for ( size_t j=0; j<2; ++j )
            plist[i][j] = i + j;

        //2.0 PFVerts - establish
        deque<vector<long64> > pfvert( n_elements, vector<long64>(2) );

        //2.1 PFVerts - assigning content
        for ( size_t i=0; i<n_elements; ++i ) {
             pfvert[i][0] = static_cast<long64>(i) - 1;
             pfvert[i][1] = static_cast<long64>(i) + 1;
          }

        //2.3 PFVerts - assigning boundary flags
        pfvert[0][0]             = CNR1;
        pfvert[n_elements-1U][1] = CNR2;

        //3.0 VSet - adding PList, PFVerts
        vset.AddPlist( plist.begin(), plist.end() );
        vset.AddPfverts( pfvert.begin(), pfvert.end() );

        //4.0 VSet - adding boundary flags
        vset.AddBFlag( 0,            CNR1 );
        vset.AddBFlag( n_vertices-1, CNR2 );
    }

} // EstablishConnectivity

/**
    Inserting SplitNodes to 1D mesh, reestablishing the connectivity

    @author R. Manasipov
*/

template<size_t dim>
void LineElementMesher<dim>::EstablishConnectivity( VSet<dim>& vset,
                                           std::set<Point<dim> >& splitnodes )
{
    if( vset.Vertices() > 0 )
    {
        const size_t  n_vertices( vset.Vertices() );

        // collect the coordinates of  nodes
        std::set<Point<dim> > nodes;
        Point<dim> p( 0.0 );
        for( size_t i =0U; i< n_vertices; i++ ) {
            for( size_t dimension=0U; dimension<dim; dimension++ )
                p[ dimension ] = vset.P(dimension,i);
            nodes.insert( p );
        }
        Point<dim> min_coord( *nodes.begin()  );
        Point<dim> max_coord( *nodes.rbegin() );

        // check whether all the splitnodes are lying inside the model
        const size_t n_splitnodes( splitnodes.size() );
        size_t splitnode_counter( 0U );
        typename std::set<Point<dim> >::iterator spit( splitnodes.begin() );
        while( splitnode_counter < n_splitnodes )
        {
            if( (*spit < min_coord) || (*spit > max_coord) || (*spit == min_coord) || (*spit == max_coord) )
                splitnodes.erase( spit );
            else
                spit++;
            splitnode_counter++;
        }

        if( splitnodes.size() > 0 )
        {
            // add splitnodes to the original nodes
            nodes.insert( splitnodes.begin(), splitnodes.end() );

            const size_t n_vertices( nodes.size() + splitnodes.size() );
            const size_t n_elements( nodes.size() - 1U );

            //1.0 VSet - refitting VSet
            vset.Resize( IsoparametricLinearLineElement().Nodes(),
                         IsoparametricLinearLineElement().Neighbors(),
                         IsoparametricLinearLineElement().ElementType(),
                         n_vertices,
                         n_elements );

            //2.0 PFVerts - establish
            deque<vector<long64> > pfvert( n_elements, vector<long64>(2) );

            //2.1 PFVerts - assigning content
            for ( size_t i=0U; i<n_elements; i++ ) {
                 pfvert[i][0] = static_cast<int32>(i) - 1U;
                 pfvert[i][1] = static_cast<int32>(i) + 1U;
              }

            //2.3 PFVerts - assigning boundary flags
            pfvert[0][0]             = CNR1;
            pfvert[n_elements-1U][1] = CNR2;

            //3.1 VSet  - assigning point x-coordinates
            //    PList - establish
            //    VSet  - adding boundary flags
            deque<vector<size_t> >  plist( n_elements, vector<size_t>(2) );

            size_t node_counter( 0U );
            size_t element_counter( 0U );
            typename std::set<Point<dim> >::const_iterator node_coordinate( nodes.begin() );
            typename std::set<Point<dim> >::const_iterator last_node_coordinate( nodes.end()   );
            while( node_coordinate != last_node_coordinate )
            {
                if( splitnodes.find( *node_coordinate ) == splitnodes.end() )
                {
                    for( size_t dimension=0U; dimension<dim; dimension++ )
                        vset.P( dimension, node_counter, (*node_coordinate)[dimension] );
                    node_counter++;
                }
                else
                {
                    for( size_t dimension=0U; dimension<dim; dimension++ )
                        vset.P( dimension, node_counter, (*node_coordinate)[dimension] );
                    node_counter++;

                    for( size_t dimension=0U; dimension<dim; dimension++ )
                        vset.P( dimension, node_counter, (*node_coordinate)[dimension] );
                    node_counter++;

                    // assign boundary flags to splitnode
                    vset.AddBFlag( node_counter - 2U, REGION_BOUNDARY );
                    vset.AddBFlag( node_counter - 1U, REGION_BOUNDARY );

                    // assign boundary flags to elements across the interface
                    pfvert[element_counter][1U]     = REGION_BOUNDARY;
                    pfvert[element_counter+1U][0U]  = REGION_BOUNDARY;
                }

                // identify local indexes of nodes inside the element
                if( element_counter != n_elements )
                {
                    plist[ element_counter ][ 0 ] = node_counter -1U;
                    plist[ element_counter ][ 1 ] = node_counter;
                }

                element_counter++;
                node_coordinate++;
            }

            //3.2 VSet - adding boundary flags
            vset.AddBFlag( 0,            CNR1 );
            vset.AddBFlag( n_vertices-1, CNR2 );

            //3.3 VSet - adding PList, PFVerts
            vset.AddPlist( plist.begin(), plist.end() );
            vset.AddPfverts( pfvert.begin(), pfvert.end() );
        }
    }

} // EstablishConnectivity

// MESH DENSITY FUNCTIONS

MeshDensity::MeshDensity(){}
MeshDensity::~MeshDensity(){}
double64 MeshDensity::operator()( double64 x ) { return x; }

LinDensity::LinDensity( double64 a, double64 b ):a_(a),b_(b){}
LinDensity::~LinDensity(){}
double64 LinDensity::operator()( double64 x ){ return a_ * x + b_; }
void LinDensity::SetA( double64 a ){ a_ = a; }
void LinDensity::SetB( double64 b ){ b_ = b; }

ExpDensity::ExpDensity( double64 a, double64 b ):a_(a),b_(b){}
ExpDensity::~ExpDensity(){}
double64 ExpDensity::operator()( double64 x ){ return b_ * std::exp( a_ * x ); }
void ExpDensity::SetA( double64 a ){ a_ = a; }
void ExpDensity::SetB( double64 b ){ b_ = b; }

ErfDensity::ErfDensity( double64 a, double64 b ):a_(a),b_(b){}
ErfDensity::~ErfDensity(){}
double64 ErfDensity::operator()( double64 x )
{
#ifdef _MSC_VER
    #if _MSC_VER < 1800  // detecting versions older than 2013
        return b_ * csmp::erf_Chebyshev( a_ * x );
    #else
        return b_ * std::erf( a_ * x );
    #endif
#else
    return b_ * std::erf( a_ * x );
#endif
}
void ErfDensity::SetA( double64 a ){ a_ = a; }
void ErfDensity::SetB( double64 b ){ b_ = b; }

template class LineElementMesher<1U>;
template class LineElementMesher<2U>;
template class LineElementMesher<3U>;

} // end namespace csmp
