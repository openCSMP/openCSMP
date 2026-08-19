#include "MeshDiagnostics.h"
#include "ANSYS_ElementSpecifications.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "PropertyConstraints.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "Standard_IO_Handler.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Internal helpers
// ----------------------------------------------------------------------------

namespace {

/**
    @brief Derives coordinate validity bounds from the model bounding box.

    The bounds are set to 100x the bounding box diagonal in each direction,
    centred on the model. This is large enough to flag genuinely erratic
    coordinates while remaining meaningful for models ranging from
    centimetre-scale to continental-scale.

    @param xyz_min  Minimum corner of the model bounding box.
    @param xyz_max  Maximum corner of the model bounding box.
    @param[out] valid_min  Lower validity bound per coordinate.
    @param[out] valid_max  Upper validity bound per coordinate.
*/
template<uint32_t dim>
void deriveValidityBounds( const Point<dim>& xyz_min,
                           const Point<dim>& xyz_max,
                           Point<dim>&       valid_min,
                           Point<dim>&       valid_max )
{
    for ( uint32_t d{ 0U }; d < dim; ++d ) {
        const double centre = ( xyz_min[d] + xyz_max[d] ) / 2.;
        const double half   = ( xyz_max[d] - xyz_min[d] ) / 2.;
        // 100x the half-extent, but at least 1 m to handle degenerate flat models
        const double margin = max( 100. * half, 1.0 );
        valid_min[d] = centre - margin;
        valid_max[d] = centre + margin;
    }
}


/**
    @brief Derives a minimum plausible volume/area/length threshold.

    Set to the cube/square/length of 1e-4 times the shortest bounding box
    dimension. This allows cm-scale objects to be resolved in km-scale models
    while still flagging genuinely degenerate elements.

    @param xyz_min  Minimum corner of the model bounding box.
    @param xyz_max  Maximum corner of the model bounding box.
    @return         Minimum plausible measure for an element of type @c dim.
*/
template<uint32_t dim>
double deriveMinVolume( const Point<dim>& xyz_min,
                        const Point<dim>& xyz_max )
{
    double shortest = numeric_limits<double>::max();
    for ( uint32_t d{ 0U }; d < dim; ++d ) {
        const double extent = xyz_max[d] - xyz_min[d];
        if ( extent > 0. && extent < shortest )
            shortest = extent;
    }
    // 1e-4 of the shortest extent, raised to the power of dim
    const double unit = 1.0e-4 * shortest;
    double result = unit;
    for ( uint32_t d{ 1U }; d < dim; ++d ) result *= unit;
    return result;
}


/**
    @brief Derives a maximum plausible volume/area/length threshold.

    Set to the cube/square/length of the longest bounding box dimension.

    @param xyz_min  Minimum corner of the model bounding box.
    @param xyz_max  Maximum corner of the model bounding box.
    @return         Maximum plausible measure for an element of type @c dim.
*/
template<uint32_t dim>
double deriveMaxVolume( const Point<dim>& xyz_min,
                        const Point<dim>& xyz_max )
{
    double longest = 0.;
    for ( uint32_t d{ 0U }; d < dim; ++d ) {
        const double extent = xyz_max[d] - xyz_min[d];
        if ( extent > longest ) longest = extent;
    }
    double result = longest;
    for ( uint32_t d{ 1U }; d < dim; ++d ) result *= longest;
    return result;
}


/**
    @brief Returns a human-readable label for the measure of an element.

    @param e  Pointer to the element.
    @return   "length", "area", or "volume" as appropriate.
*/
template<uint32_t dim, template<uint32_t> class CELL>
const char* measureLabel( const CELL<dim>* e )
{
    if ( e->IsLine()    ) return "length";
    if ( e->IsSurface() ) return "area";
    return "volume";
}

} // anonymous namespace


// ----------------------------------------------------------------------------
//  ScrutinizeMesh
// ----------------------------------------------------------------------------

template<uint32_t dim>
MeshDiagnosticsResult
MeshDiagnostics<dim>::ScrutinizeMesh( Model<dim>& model,
                                      double aspect_ratio_warn,
                                      double coincidence_tol ) const
{
    MeshDiagnosticsResult result;

    Region<dim>& sgroup( model.Region( "Model" ) );
    sgroup.UpdateMemberIndexes();

    // -------------------------------------------------------------------------
    // 0. Bounding box and derived thresholds
    // -------------------------------------------------------------------------
    Point<dim> xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );

    Point<dim> valid_min, valid_max;
    deriveValidityBounds<dim>( xyz_min, xyz_max, valid_min, valid_max );

    const double vol_min = deriveMinVolume<dim>( xyz_min, xyz_max );
    const double vol_max = deriveMaxVolume<dim>( xyz_min, xyz_max );

    cout << "\nMeshDiagnostics<" << dim << ">::ScrutinizeMesh: "
         << "Model bounding box:\n";
    for ( uint32_t d{ 0U }; d < dim; ++d )
        cout << "\t[" << d << "] " << xyz_min[d] << " to " << xyz_max[d] << "\n";
    cout << "\tDerived volume/area/length thresholds: ["
         << vol_min << ", " << vol_max << "]\n";


    // -------------------------------------------------------------------------
    // 1. Node coordinate check
    // -------------------------------------------------------------------------
    cout << "\nMeshDiagnostics<" << dim << ">::ScrutinizeMesh: "
         << "Verifying node coordinates...\n";

    for ( auto nit = sgroup.NodesBegin(); nit != sgroup.NodesEnd(); ++nit ) {
        const Node<dim>* n = *nit;
        for ( uint32_t d{ 0U }; d < dim; ++d ) {
            const double coord = n->Coordinate()[d];
            if ( coord < valid_min[d] || coord > valid_max[d] ) {
                cerr << "\tErratic coordinate[" << d << "] on node "
                     << n->Idx() << ": " << coord << "\n";
                ++result.erratic_coordinates;
                result.problems = true;
            }
        }
    }


    // -------------------------------------------------------------------------
    // 2. Per-element checks: coincident nodes, segment lengths, volume
    // -------------------------------------------------------------------------
    cout << "\nMeshDiagnostics<" << dim << ">::ScrutinizeMesh: "
         << "Checking elements...\n";

    vector<size_t> negative_volume_ids;
    vector<size_t> zero_volume_ids;
    vector<size_t> erratic_volume_ids;
    vector<size_t> large_aspect_ratio_ids;
    vector<size_t> coincident_node_ids;

    // pre-reserve to avoid repeated allocation
    negative_volume_ids.reserve( 1024U );
    zero_volume_ids.reserve( 1024U );
    erratic_volume_ids.reserve( 1024U );
    large_aspect_ratio_ids.reserve( 1024U );
    coincident_node_ids.reserve( 1024U );

    result.sl_min         =  numeric_limits<double>::max();
    result.sl_max         =  0.;
    result.aspect_ratio_min =  numeric_limits<double>::max();
    result.aspect_ratio_max =  0.;

    bool first_element = true;

    for ( auto eit = sgroup.CellsBegin(); eit != sgroup.CellsEnd(); ++eit ) {
        const auto* e = *eit;
        const size_t eidx = e->Idx();

        // --- 2a. Coincident node pairs ---
        for ( uint32_t i{ 0U }; i < e->Nodes(); ++i ) {
            const Point<dim> pi( e->N(i)->Coordinate() );
            for ( uint32_t j{ i + 1U }; j < e->Nodes(); ++j ) {
                const Point<dim> pj( e->N(j)->Coordinate() );
                if ( pi.CoincidesWithWithinTolerance( pj, coincidence_tol ) ) {
                    cerr << "\tElement " << eidx
                         << ": nodes " << e->N(i)->Idx()
                         << " (local " << i << ") and "
                         << e->N(j)->Idx()
                         << " (local " << j << ") are coincident\n";
                    ++result.coincident_node_pairs;
                    result.problems = true;
                    coincident_node_ids.push_back( eidx );
                }
            }
        }

        // --- 2b. Segment lengths and per-element aspect ratio ---
        vector<double> lengths;
        e->SegmentLengths( lengths );

        double e_sl_min =  numeric_limits<double>::max();
        double e_sl_max =  0.;

        for ( const double len : lengths ) {
            if ( len < result.sl_min ) result.sl_min = len;
            if ( len > result.sl_max ) result.sl_max = len;
            if ( len < e_sl_min      ) e_sl_min      = len;
            if ( len > e_sl_max      ) e_sl_max      = len;
        }

        // aspect ratio only meaningful when element has more than one segment
        if ( lengths.size() > 1U && e_sl_min > 0. ) {
            const double ratio = e_sl_max / e_sl_min;

            if ( first_element ) {
                result.aspect_ratio_min = result.aspect_ratio_max = ratio;
                first_element = false;
            }
            if ( ratio > result.aspect_ratio_max ) result.aspect_ratio_max = ratio;
            if ( ratio < result.aspect_ratio_min ) result.aspect_ratio_min = ratio;

            if ( ratio > aspect_ratio_warn ) {
                cerr << "\tLarge aspect ratio in element " << eidx
                     << ": " << ratio
                     << " (segments: " << e_sl_min << " to " << e_sl_max << ")\n";
                ++result.large_aspect_ratios;
                result.problems = true;
                large_aspect_ratio_ids.push_back( eidx );
            }
        }

        // --- 2c. Volume / area / length ---
        double volume = numeric_limits<double>::quiet_NaN();
        try {
            volume = e->Volume();
        }
        catch ( const csmp::Exception& ex ) {
            cerr << "\tVolume() threw for element " << eidx
                 << ": " << ex.what() << "\n";
            result.problems = true;
            negative_volume_ids.push_back( eidx );
            continue;
        }
        catch ( ... ) {
            cerr << "\tVolume() threw unknown exception for element "
                 << eidx << "\n";
            result.problems = true;
            negative_volume_ids.push_back( eidx );
            continue;
        }

        if ( isnan( volume ) ) continue;

        const char* label = measureLabel( e );

        if ( volume < 0. ) {
            cerr << "\tNegative " << label << " in element "
                 << eidx << ": " << volume << "\n";
            ++result.negative_volumes;
            result.problems = true;
            negative_volume_ids.push_back( eidx );
        }
        else if ( volume < vol_min ) {
            cerr << "\tNear-zero " << label << " in element "
                 << eidx << ": " << volume
                 << " (threshold: " << vol_min << ")\n";
            ++result.zero_volumes;
            result.problems = true;
            zero_volume_ids.push_back( eidx );
        }
        else if ( volume > vol_max ) {
            cerr << "\tImplausibly large " << label << " in element "
                 << eidx << ": " << volume
                 << " (threshold: " << vol_max << ")\n";
            ++result.erratic_volumes;
            result.problems = true;
            erratic_volume_ids.push_back( eidx );
        }
    }


    // -------------------------------------------------------------------------
    // 3. Form diagnostic regions and write VTK output
    // -------------------------------------------------------------------------
    auto formAndOutput = [&]( vector<size_t>&  ids,
                               const string&    region_name,
                               const string&    vtk_filename )
    {
        if ( ids.empty() ) return;

        // deduplicate — an element may appear in multiple sub-checks
        sort( ids.begin(), ids.end() );
        ids.erase( unique( ids.begin(), ids.end() ), ids.end() );

        // FormRegionFrom( const char* regionname, std::vector<size_t>& elmt_ids, bool is_unique=true );
        model.FormRegionFrom( region_name.c_str(), ids, false );

        VTK_Interface<dim> vtk_out;
        vtk_out.OutputDataToVTK( model, region_name, vtk_filename,
                                 "permeability", 0, true );

        cout << "\tRegion '" << region_name << "' formed with "
             << ids.size() << " element(s) — written to " << vtk_filename << "\n";
    };

    cout << "\nMeshDiagnostics<" << dim << ">::ScrutinizeMesh: "
         << "Forming diagnostic regions...\n";

    formAndOutput( negative_volume_ids,
                   "diag: negative volume",
                   "DIAG_NEGATIVE_VOLUME" );

    formAndOutput( zero_volume_ids,
                   "diag: zero volume",
                   "DIAG_ZERO_VOLUME" );

    formAndOutput( erratic_volume_ids,
                   "diag: erratic volume",
                   "DIAG_ERRATIC_VOLUME" );

    formAndOutput( large_aspect_ratio_ids,
                   "diag: large aspect ratio",
                   "DIAG_LARGE_ASPECT_RATIO" );

    formAndOutput( coincident_node_ids,
                   "diag: coincident nodes",
                   "DIAG_COINCIDENT_NODES" );


    // -------------------------------------------------------------------------
    // 4. Summary report
    // -------------------------------------------------------------------------
    cout << "\nMeshDiagnostics<" << dim << ">::ScrutinizeMesh: Summary\n";
    cout << "\tBounding box min:            " << xyz_min << "\n";
    cout << "\tBounding box max:            " << xyz_max << "\n";
    cout << "\tShortest segment:            " << result.sl_min << "\n";
    cout << "\tLongest segment:             " << result.sl_max << "\n";
    cout << "\tMin per-element aspect ratio:" << result.aspect_ratio_min << "\n";
    cout << "\tMax per-element aspect ratio:" << result.aspect_ratio_max << "\n";
    cout << "\tCoincident node pairs:       " << result.coincident_node_pairs << "\n";
    cout << "\tNegative volumes:            " << result.negative_volumes << "\n";
    cout << "\tNear-zero volumes:           " << result.zero_volumes << "\n";
    cout << "\tImplausibly large volumes:   " << result.erratic_volumes << "\n";
    cout << "\tErratic coordinates:         " << result.erratic_coordinates << "\n";
    cout << "\tLarge aspect ratio elements: " << result.large_aspect_ratios << "\n";
    cout << "\tOverall problems found:      "
         << ( result.problems ? "YES" : "no" ) << "\n";

    return result;

} // end ScrutinizeMesh




template<uint32_t dim>
void MeshDiagnostics<dim>::ElementVolumeRange( const Model<dim>& sg,
                                               double& vmin, double& vmax ) const
 {                                             
    double    vol;
    bool  first_element(true); 
    
    const Region<dim>& sgref(sg.Region("Model"));
    
    for ( auto eit=sgref.CellsBegin();
          eit!=sgref.CellsEnd(); eit++ ) 
      {
         if ( first_element ) {
              vmin = vmax = (*eit)->Volume();
              first_element = false;
           }
         vol = (*eit)->Volume();
         if ( vol < vmin ) vmin = vol;
         if ( vol > vmax ) vmax = vol;
  
      } // end for all elements
    
    cout <<"\nMeshDiagnostics::ElementVolumeRange: "<< vmin <<" to "<< vmax << endl;

 } // end ElementVolumeRange








template<uint32_t dim>
void MeshDiagnostics<dim>::FixFiniteElementNeighborOrientationOfSurfaceMeshes( Model<dim>& sg ) const
  {
    std::vector<Point<dim> > bc_vec(0);
    std::vector<size_t>      id_vec(0);
    double                   sign;
    size_t                   id, cntr(0);
    
    cout <<"\nMeshDiagnostics::FixFiniteElementNeighborOrientationOfSurfaceMeshes: Checking for any clockwise ordered element neighbor IDs..." << endl;
    
    Region<dim>& sgref(sg.Region("Model"));

    for ( auto eit=sgref.CellsBegin(); eit!=sgref.CellsEnd(); eit++ )
      { 
        if ( bc_vec.size() != (*eit)->Neighbors() ) bc_vec.resize((*eit)->Neighbors());
        sign = 0.;
        // construct a polygon of with vertices equal to the barycenters of the neighbor FEs
        // if there is no neighbor FE, use barycenter of current FE.
        for ( uint32_t i{0}; i<(*eit)->Neighbors(); i++ ) {
            if ( (*eit)->Neighbor(i) != NULL ) bc_vec[i] = (*eit)->Neighbor(i)->BaryCenter();
            else                               bc_vec[i] = (*eit)->BaryCenter(); 
          }
        for ( uint32_t i{0}; i<(bc_vec.size()-1U); i++ ) {
            // calculate sign of polygon determinant z = x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + xn * y1 + x1 * yn
            sign += bc_vec[i][0] * bc_vec[i+1U][1] - bc_vec[i+1][0] * bc_vec[i][1];
          } 
        // last entry  
        sign += bc_vec[bc_vec.size()-1U][0] * bc_vec[0][1] - bc_vec[0][0] * bc_vec[bc_vec.size()-1U][1];
        
        // negative sign means that neighbor FEs are ordered in a clock-wise fashion
        if ( sign < 0. ) {
            cntr++;
            if ( id_vec.size() != (*eit)->Neighbors() ) id_vec.resize((*eit)->Neighbors());
            id = (*eit)->Neighbors()-1;
            // reorder element ids
            for ( uint32_t i{0U}; i<(*eit)->Neighbors(); i++ ) {
                if ( (*eit)->Neighbor(i) != NULL ) id_vec[id-i] = (*eit)->Neighbor(i)->Idx();
                else                               id_vec[id-i] = 0;
              }
            for ( uint32_t i{0U}; i<id_vec.size(); i++ ) {
                if ( id_vec[i] > 0 ) (*eit)->Assign( i, sgref.E(id_vec[i]-1) );
                else                 (*eit)->Assign( i, static_cast<Element<dim>*>(nullptr) );
              }
          }
    
      }
    
    cout <<"\nMeshDiagnostics::FixFiniteElementNeighborOrientationOfSurfaceMeshes: ";
    cout <<"\nChanged order of element neighbor IDs to counter clockwise for " << cntr << " elements " << endl;
    
 } // end FixFiniteElementNeighborOrientationOfSurfaceMeshes





/**

Method tries to test the current finite-element mesh for internal
consistency. Thus it checks whether there are no gaps in the mesh object
numbering (e.g., node numbers which are larger than the actual number of
nodes etc.), no zero segment lengths of elements, no NAN coordinates etc.
The method will not report duplicate Node or IntegrationPoint coordinates. 

@section implementation Implementation

Thus far, scrutinizeMesh() assumes a fixed maximum dimension of 40'000km.

@section application Application

scrutinizeMesh() is recommended if meshes from external tools are
used in CSMP computations. Its performance is fast, and while it does not
guarantee for the detection of all possible errors many very obvious 
errors are reported. 

@section messages Messages 

At the end of the test the method will write a few mesh
diagnostics to stdout: 
 
"horizontal min, max coordinate
"vertical   min, max coordinate
"length of longest segment
"length of shortest segment
"maximum segment length ratio
"minimum segment length ratio

SKM fix 18/10/2014: additions to catch exceptions associated with negative Jacobians.
      
 */
/* buggy old version

template<uint32_t dim>
bool MeshDiagnostics<dim>::ScrutinizeMesh( Model<3U>& sg ) const
 {
    Region<3>&  sgroup(sg.Region("Model"));
    sgroup.UpdateMemberIndexes();
    double l_segm, segm_length_ratio(DBL_MAX),
           segm_length_ratio_min(DBL_MAX), segm_length_ratio_max(0.),
           sl_min, sl_max,
           sl_min_e(-4e7), sl_max_e(4e7), // the equator of the earth
           h_min(-4.0e7), h_max(4.0e7),
           v_min(-4.0e7), v_max(4.0e7),
           h_valmin, v_valmin,
           h_valmax, v_valmax;
    bool   problems = false;
           
    DenseMatrix<DM_MIN>    XY(3,3);
    vector<uint32_t>       nids(3);
    Standard_IO_Handler    stdio;

    // Node coordinates
    cout <<"\nMeshDiagnostics<dim>::ScrutinizeMesh: Verifying node coordinates..."<< endl;
    //cout <<"\nplease enter intended horizontal min, max coordinates: ";
    //cin >> h_min >> h_max;
    //cout <<"\nplease enter intended vertical min, max coordinates: ";
    //cin >> v_min >> v_max;
 
    h_valmin = h_valmax = sgroup.N(0)->x();
    v_valmin = v_valmax = sgroup.N(0)->y();    
 
    for ( auto eit=sgroup.CellsBegin(); eit!=sgroup.CellsEnd(); eit++ )
      {
         (*eit)->NodeCoordinateMatrix( XY );
         for ( uint32_t i{0}; i<(*eit)->Nodes(); i++ )
           {
              if ( XY(i,0) > h_valmax ) h_valmax = XY(i,0);
              if ( XY(i,1) > v_valmax ) v_valmax = XY(i,1);
              if ( XY(i,0) < h_valmin ) h_valmin = XY(i,0);
              if ( XY(i,1) < v_valmin ) v_valmin = XY(i,1);
              if ( XY(i,0) < h_min ||  XY(i,0) > h_max )  
                cout <<"\nErratic horizontal coordinate of node "<< nids[i] <<": "<< XY(i,0);
              if ( XY(i,1) < v_min ||  XY(i,1) > v_max )  
                cout <<"\nErratic vertical coordinate of node   "<< nids[i] <<": "<< XY(i,1);
           }
         if ( XY(0,0) == XY(1,0) && XY(0,1) == XY(1,1) )
           cout <<"\nFirst 2 coordinates of element "<< (*eit)->Idx() <<" are identical";
         if ( XY(1,0) == XY(2,0) && XY(1,1) == XY(2,1) )
           cout <<"\nSecond 2 coordinates of element "<< (*eit)->Idx() <<" are identical";
      }

    // Segment length
    vector<double> lengths;

    cout <<"\n\nMeshDiagnostics<dim>::ScrutinizeMesh: Testing length of element segments..."<< endl;
    
    (*sgroup.CellsBegin())->SegmentLengths( lengths );
    sl_min = sl_max = lengths[0];

    for ( auto eit=sgroup.CellsBegin(); eit!=sgroup.CellsEnd(); eit++ )
      {
         (*eit)->SegmentLengths( lengths );
         for ( uint32_t i{0}; i<lengths.size(); i++ )
           {
              l_segm = lengths[i];
              if ( l_segm > sl_max )   sl_max   = l_segm;
              if ( l_segm < sl_min )   sl_min   = l_segm;
              if ( l_segm > sl_max_e ) sl_max_e = l_segm;
              if ( l_segm < sl_min_e ) sl_min_e = l_segm;
           }
         segm_length_ratio = sl_max_e / sl_min_e;
         if ( eit == sgroup.CellsBegin() ) segm_length_ratio_max = segm_length_ratio_min = segm_length_ratio;
         if ( segm_length_ratio > 10. )
           cout <<"\nlarge segment length ratio in element "<< (*eit)->Idx() <<": "<< segm_length_ratio;
          
         if ( segm_length_ratio > segm_length_ratio_max ) segm_length_ratio_max = segm_length_ratio;
         if ( segm_length_ratio < segm_length_ratio_min ) segm_length_ratio_min = segm_length_ratio;
      }


    // Element volume
    double          volume;
    bool            repeat(true);
    vector<size_t>  element_numbers;
   
    cout <<"\n\nMeshDiagnostics<dim>::ScrutinizeMesh: Verifying element volumes/areas/lengths..."<< endl;
    for ( auto eit=sgroup.CellsBegin(); eit!=sgroup.CellsEnd(); eit++ )
      {
         // catching exceptions that might originate from negative Jacobian calculations
         try {
              volume=(*eit)->Volume();
           }
         catch(...) {
             // getting ready for trouble
             element_numbers.reserve(10000);
             // nothing is done here
           }
         if ( volume < 0. && repeat == true )
           {
              cerr <<"\nMeshDiagnostics<dim>::ScrutinizeMeshErratic area/volume of element ";
              cerr << (*eit)->Idx() <<": "<< volume;
              repeat = true; // stdio.YesNo("Would you like to form region of elements with negative area");
              if ( repeat ) element_numbers.push_back( (*eit)->Idx() );
           }
         if ( volume == 0. || fabs(volume) < 1.0e-6 || volume > 1.0e+10 )  
           cerr <<"\nErratic area/volume of element "<< (*eit)->Idx() <<": "<< volume;
      }
   
    // creation of new region "erratic elements" and output to VTK
    if ( repeat && element_numbers.size() > 0U ) {
         sg.FormRegionFrom( "erratic elements", element_numbers );
         VTK_Interface<3U>  vtk_out;
         vtk_out.OutputDataToVTK( sg, "erratic elements", "NEGATIVE_JACOBIAN_ELEMENTS",
                                      "permeability", 0, true );
      }

    // output all the diagnosed values
    cout <<"\n\nModel<3U>::ScrutinizeMesh: The mesh test found the characteristics:"<< endl;
    cout <<"\thorizontal min, max coordinate: "<< h_valmin <<","<< h_valmax << endl;
    cout <<"\tvertical   min, max coordinate: "<< v_valmin <<","<< v_valmax << endl;
    cout <<"\tlength of longest segment:      "<< sl_max << endl;
    cout <<"\tlength of shortest segment:     "<< sl_min << endl;
    cout <<"\tmaximum segment length ratio:   "<< segm_length_ratio_max << endl;
    cout <<"\tminimum segment length ratio:   "<< segm_length_ratio_min << endl;         
    cout <<"\nModel<dim>::ScrutinizeMesh: tests completed.";
    
    // determining return value
    if ( h_valmin < h_min || h_valmax > h_max || v_valmin < v_min || v_valmax > v_max ) problems = true;
    if ( segm_length_ratio_min < 0.1 || segm_length_ratio_max > 10.0 || sl_min == 0.0 ) problems = true;
    
    return problems;
 
 } // end scrutinize mesh
*/



/** Creates a range of quality variables, computes them and outputs them from the target region to VTU; creates diagnostic regions for this purpose.
 
    Diagnostics obtained: cells:
         
         1.  "cell volume" - cells that are smaller than 1e-5 of average cell size are reported in region "small cells"
           (a special treatment is applied to lower dimensional elements)
         2.  "Jacobian determinant" - a negative value indicates a problem with the node-numbering of the cell or twistedness
         3.  "aspect ratio" - longest over shortest segment length
         4. "minimum inter-segment angle" - indicating splinter elements
         5. "width over height ratio" - suitable for the detection of non-manifold vertices (where all nodes lie in the same plane)

    Diagnostics obtained: nodes:
    
         1. "node connections" - number of nodes the node is connected to. A highly variable number indicates a poor quality mesh
         2. "shortest distance to node" - any of the other nodes surrounding the node

         
    Cells that do not suffer from any of the checked criteria are stored in a new region that fullfil the criterion "fit for computation" = 1.
*/
/*
template<uint32_t dim>
bool MeshDiagnostics<dim>::ComputeQualityMetricsAndOutputToVTU( Model<dim>& model, const std::string& region )
 {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    
    if ( !model.ContainsRegion(region) ) {
         csmp_error.Note( ERROR, "MeshDiagnostics<dim>::ComputeQualityMetricsAndOutputToVTU",
                         "Region to diagnose not found; no diagnostics obtained");
         return false;
      }
    Region<dim>& subdomain = model.Region( region );
    bool         mesh_suitable_for_computation{ true };
    
    // 0. creating the properties needed to store the diagnostic values
    const csmp::Index  evol_key = model.CreateProperty( "cell volume", "m3", SCALAR, ELEMENT );
    const csmp::Index  mJac_key = model.CreateProperty( "Jacobian determinant", "none", SCALAR, ELEMENT );
    const csmp::Index  los_key = model.CreateProperty( "aspect ratio", "X", SCALAR, ELEMENT );
    
    // 1. obtaining and storing diagnostic values
    for ( auto& it : subdomain.CellVector() )
      {
         // cell volume
         const double cell_volume = it->Volume();
         it->Store( evol_key, makeScalar(ANY,cell_volume) );
         // determinant of Jacobian matrix (gets minimum value)
         double min_determinant = it->det_J_AtIntegrationPoint(0U);
         for ( auto i{1u}; i<it->IntegrationPoints(); i++ )
           min_determinant = min( min_determinant, it->det_J_AtIntegrationPoint(i) );
         it->Store( mJac_key, makeScalar(ANY,min_determinant) );
         // aspect ratio = longest over shortest segment length
         it->Store( los_key, makeScalar(ANY, it->AspectRatio() ) );
           

      }
      
    // 2. creating regions from value ranges of particular diagnostic values
     
    // 3. diagnostics applied to nodes
    const csmp::Index  conn_key  = model.CreateProperty( "node connections", "none", SCALAR, NODE );
    const csmp::Index  ndist_key = model.CreateProperty( "shortest distance to node", "m", SCALAR, NODE );
        
    // obtaining and storing diagnostic values
    for ( auto& nit : subdomain.NodeVector() )
      {
         // node connections
         const auto n_nbors{ nit->Neighbors() };
         nit->Store( conn_key, makeScalar(ANY,static_cast<double>(n_nbors)) );
         // shortest distance to node
         const Point<dim> ncoord{ nit->Coordinate() };
         double           min_distance{ 1.0e30 };
         for ( auto i{0u}; i< n_nbors; i++ )
           min_distance = min( min_distance, ncoord.DistanceTo( nit->Neighbor(i)->Coordinate() ) );
         nit->Store( ndist_key, makeScalar(ANY,min_distance) );
      }


    return mesh_suitable_for_computation;
     
 } // end ComputeQualityMetricsAndOutputToVTU
*/



/** Creates a range of quality variables, computes them and outputs them from the target region to VTU; creates diagnostic regions for this purpose.

    Diagnostics obtained: cells:

         1.  "cell volume" - cells that are smaller than 1e-5 of average cell size are reported in region "small cells"
           (a special treatment is applied to lower dimensional elements)

         2.  "Jacobian determinant" - a negative value indicates a problem with the node-numbering of the cell or twistedness

         3.  "aspect ratio" - longest over shortest segment length

         4. "minimum inter-segment angle" - indicating splinter elements

         5. "width over height ratio" - suitable for the detection of non-manifold vertices (where all nodes lie in the same plane)

    Diagnostics obtained: nodes:

         1. "node connections" - number of nodes the node is connected to. A highly variable number indicates a poor quality mesh

         2. "shortest distance to node" - any of the other nodes surrounding the node

    Cells that do not suffer from any of the checked criteria are stored in a new region that fullfil the criterion "fit for computation" = 1.

*/

template<uint32_t dim>
bool MeshDiagnostics<dim>::ComputeQualityMetricsAndOutputToVTU( Model<dim>& model, const std::string& region )
  {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( !model.ContainsRegion(region) ) {
      csmp_error.Note( ERROR, "MeshDiagnostics<dim>::ComputeQualityMetricsAndOutputToVTU",
                       "Region to diagnose not found; no diagnostics obtained");
      return false;
    }

    Region<dim>& subdomain = model.Region( region );
    bool         mesh_suitable_for_computation{ true };

    // 0. creating the properties needed to store the diagnostic values
    const csmp::Index  evol_key = model.CreateProperty( "cell volume", "Ve", "m3", SCALAR, ELEMENT );
    const csmp::Index  mJac_key = model.CreateProperty( "Jacobian determinant", "Jinv", "none", SCALAR, ELEMENT );
    const csmp::Index  los_key = model.CreateProperty( "aspect ratio", "AR", "X", SCALAR, ELEMENT );

    // 1. obtaining and storing diagnostic values
    for ( auto& it : subdomain.CellVector() )
    {
      // cell volume
      const double cell_volume = it->Volume();
      it->Store( evol_key, makeScalar(ANY,cell_volume) );

      // determinant of Jacobian matrix (gets minimum value)
      double min_determinant = it->det_J_AtIntegrationPoint(0U);
      for ( auto i{1u}; i<it->IntegrationPoints(); i++ )
        min_determinant = min( min_determinant, it->det_J_AtIntegrationPoint(i) );
      it->Store( mJac_key, makeScalar(ANY,min_determinant) );

      // aspect ratio = longest over shortest segment length
      it->Store( los_key, makeScalar(ANY, it->AspectRatio() ) );
    }


    // 2. diagnostics applied to nodes
    const csmp::Index  conn_key  = model.CreateProperty( "node connections", "NC", "none", SCALAR, NODE );
    const csmp::Index  ndist_key = model.CreateProperty( "shortest distance to node", "dist-nd", "m", SCALAR, NODE );

    // obtaining and storing diagnostic values
    for ( auto& nit : subdomain.NodeVector() )
    {
      // node connections
      const auto n_nbors{ nit->Neighbors() };
      nit->Store( conn_key, makeScalar(ANY,static_cast<double>(n_nbors)) );

      // shortest distance to node
      const Point<dim> ncoord{ nit->Coordinate() };
      double           min_distance{ 1.0e30 };
      for ( auto i{0u}; i< n_nbors; i++ )
        min_distance = min( min_distance, ncoord.DistanceTo( nit->Neighbor(i)->Coordinate() ) );
      nit->Store( ndist_key, makeScalar(ANY,min_distance) );

    }

    // 3. creating regions from value ranges of particular diagnostic values
    VTK_Interface<dim> vtk_output;
    //cell volume
    double cell_volume_threshold;
    cout << "\nPlease type in a threshold value of cell volume (m3) below which the cells will be output to VTU"<<endl;
    cin >> cell_volume_threshold;
    model.FormRegionFrom("cells_with_small_volumes", "cell volume", 0., cell_volume_threshold, false);
    vtk_output.OutputDataToVTK(model, "cells_with_small_volumes", "", "cell volume", 0., true);
    //Jacobian_determinant
//    const double jacobian_determinant_threshold(0.);
    model.FormRegionFrom("cells_with_negative_Jacobian_determinant", "Jacobian determinant", -1.0e30, 0., false);
    vtk_output.OutputDataToVTK(model, "cells_with_negative_Jacobian_determinant", "", "Jacobian determinant", 0., true);
    //aspect ratio
    double aspect_ratio_threshold;
    cout << "\nPlease type in a threshold value of aspect ratio above which the cells will be output to VTU"<<endl;
    cin >> aspect_ratio_threshold;
    model.FormRegionFrom("cells_with_large_aspect_ratios", "aspect ratio", aspect_ratio_threshold, 1.0e30, false);
    vtk_output.OutputDataToVTK(model, "cells_with_large_aspect_ratios", "", "aspect ratio", 0., true);
    
    //shortest distance to node
    double shortest_distance_to_node_threshold;
    cout << "\nPlease type in a threshold value of the shortest distance between nodes, below which the cells will be output to VTU"<<endl;
    cin >> shortest_distance_to_node_threshold;
    model.FormRegionFrom("cells_with_close_neighboring_nodes", "shortest distance to node", 0., shortest_distance_to_node_threshold, false);
    vtk_output.OutputDataToVTK(model, "cells_with_close_neighboring_nodes", "", "shortest distance to node", 0., true);

    
    //create property constraints for "fit for computation" cells
    PropertyConstraints prop_constraints(model.Database(), "cell volume", cell_volume_threshold, 1.0e30);
    prop_constraints.AddConstraint("Jacobian determinant", 0., 1.0e30);
    prop_constraints.AddConstraint("aspect ratio", 0., aspect_ratio_threshold);
    prop_constraints.AddConstraint("shortest distance to node", shortest_distance_to_node_threshold, 1.0e30);
    //form a region with defined property constraints
    model.FormRegionFrom("region_with_high_quality_cells", prop_constraints);
    // RegionInterface::FormRegionFrom PRopertyConstraints erases region if it
    // is empty
    if( model.ContainsRegion( "region_with_high_quality_cells" ) ) {
        VTU_Interface<dim>  vtu_output(model);
        list<string> output_properties;
        output_properties.emplace_back( "cell volume" );
        output_properties.emplace_back( "Jacobian determinant" );
        output_properties.emplace_back( "aspect ratio" );
        output_properties.emplace_back( "shortest distance to node" );
        vtu_output.OutputDataToVTU( model.Name(), output_properties, "region_with_high_quality_cells", 0 );
    }

    return mesh_suitable_for_computation;
    
  } // end ComputeQualityMetricsAndOutputToVTU







/**
    tests for negative element volume / area / length
*/
template<uint32_t dim>
bool MeshDiagnostics<dim>::DetectPotentiallyMisnumberedElements( const Model<dim>& model ) const
 {
    const Region<dim> model_domain(model.Region("Model"));
    bool discovered_negative_element_volume(false);
   
   
    for ( auto eit=model_domain.CellsBegin(); eit!=model_domain.CellsEnd(); ++eit )
      if ( (*eit)->Volume() < 0. ) {
           discovered_negative_element_volume = true;
           cerr <<"\nMeshDiagnostics: Deteceted element with negative volume: ";
           (*eit)->Out();
        }
    
    return false;
 
 } // end




/**
    detects whether different Dirichlet conditions have been assigned to a single finite element.
    
    @attention method has been implemented only for scalar and vector variables.
*/
template<uint32_t dim>
bool MeshDiagnostics<dim>::DetectConflictingDirichletConditions( const Model<dim>& model, const char* variable_of_interest, VARIABLE_FLAG status ) const
 {
    const csmp::Index key = model.Database().StorageKey( variable_of_interest );
    const Region<dim> model_domain(model.Region("Model"));
    assert( key.place == NODE );
    assert( key.type == SCALAR || key.type == VECTOR );
   
    bool duplicate_constraits(false);
   
    // scalar node variables
    if ( key.type == SCALAR )
      for ( auto eit=model_domain.CellsBegin(); eit!=model_domain.CellsEnd(); ++eit ) {
           double value(numeric_limits<double>::quiet_NaN());
           bool     detected_status(false);
           for ( uint32_t i{0}; i<(*eit)->Nodes(); ++i ) {
                // finding status-flagged nodes and reading their stored values
                if ( !detected_status && (*eit)->N(i)->Status(key) == status ) {
                     value           = (*eit)->N(i)->Read(key);
                     detected_status = true;
                  }
                // if such nodes were already discovered, a comparison with previous values is made
                if ( detected_status && (*eit)->N(i)->Status(key) == status ) {
                     double next_value = (*eit)->N(i)->Read(key);
                     if ( next_value != value ) {
                          cerr <<"\nMeshDiagnostics<dim>::DetectConflictingDirichletConditions: detected conflicting constrained ScalarVariable values ";
                          cerr << value <<" vs. "<< next_value <<" ";
                          cerr <<"for variable '"<< variable_of_interest <<"' in Element: "<< (*eit)->Idx() <<"\n";
                          (*eit)->Out();
                          duplicate_constraits = true;
                       }
                  }
             }
        }
   
    if ( key.type == VECTOR ) {
          for ( auto eit=model_domain.CellsBegin(); eit!=model_domain.CellsEnd(); ++eit ) {
               VectorVariable<dim> vc;
               bool     detected_status(false);
               for ( uint32_t i{0}; i<(*eit)->Nodes(); ++i ) {
                    // recovering the variable
                    (*eit)->N(i)->Read( key, vc );
                 
                    // for each variable component
                    for ( auto j{0U}; j<dim; j++ )
                      {
                         double value(numeric_limits<double>::quiet_NaN());
                         // finding status-flagged nodes and reading their stored values
                         if ( !detected_status && vc.Flag(j) == status ) {
                              value = vc[j];
                              detected_status = true;
                           }
                         // if such nodes were already discovered, a comparison with previous values is made
                         if ( detected_status && vc.Flag(j) == status ) {
                              double next_value = vc[j];
                              if ( next_value != value ) {
                                   cerr <<"\nMeshDiagnostics<dim>::DetectConflictingDirichletConditions: detected conflicting constrained VectorVariable component values, component ";
                                   cerr << j <<": "<< value <<" vs. "<< next_value <<" ";
                                   cerr <<"for variable '"<< variable_of_interest <<"' in Element: "<< (*eit)->Idx() <<"\n";
                                   (*eit)->Out();
                                   duplicate_constraits = true;
                                }
                            }
                      }
                }
            }
     } // end for
 
    return duplicate_constraits;
   
 } // end DetectConflictingDirichletConditions
 
 
 
  
/**
    checks whether the any value of a computed node variable (P,T,C) lies outside of the range of the values in its neighborhood
*/
template<uint32_t dim>
bool MeshDiagnostics<dim>::DetectNonMonotonicity( const Model<dim>& model, const char* variable_of_interest ) const
 {
    const csmp::Index key = model.Database().StorageKey( variable_of_interest );

    const Region<dim> model_domain(model.Region("Model"));
    assert( key.place == NODE );
    assert( key.type == SCALAR );
   
    bool local_peak_values(false);
   
    // scalar node variables
    if ( key.type == SCALAR )
      for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
           double n_value = (*nit)->Read(key);
           double val_min(1e30), val_max(-1e30);
           for ( uint32_t i{0u}; i<(*nit)->Neighbors(); ++i ) {
                val_min = std::min( val_min, (*nit)->Neighbor(i)->Read(key) );
                val_max = std::max( val_max, (*nit)->Neighbor(i)->Read(key) );
             }
           // checking that the value lies within range of node neighbors
           if ( n_value < val_min || n_value > val_max ) {
                cerr <<"\nvalue: "<< n_value <<" vs. "<< val_min <<"-"<< val_max;
                (*nit)->Out();
                local_peak_values = true;
             }
        }
 
    return local_peak_values;
   
 } // end CheckMonotonicity





/// elements where all nodes have a status constraint so that they do not participate in the computation
template<uint32_t dim>
bool MeshDiagnostics<dim>::DetectOverConstrainedElements( const Model<dim>& model, const char* variable_of_interest, VARIABLE_FLAG status ) const
{
    const csmp::Index key = model.Database().StorageKey( variable_of_interest );

    const Region<dim> model_domain(model.Region("Model"));
    assert( key.place == NODE );
  
    bool over_constrained_elmts(false);
   
    for ( auto eit=model_domain.CellsBegin(); eit!=model_domain.CellsEnd(); ++eit ) {
           const size_t nodes = (*eit)->Nodes();
           size_t status_constraints(0U);
           if ( key.type == SCALAR ) {
                for ( uint32_t i{0u}; i<nodes; ++i )
                  if ( (*eit)->N(i)->Status(key) == status )
                    status_constraints++;
              }
           if ( status_constraints == nodes ) {
                over_constrained_elmts = true;
                (*eit)->Out();
             }
      }
 
   return over_constrained_elmts;
  
} // end CheckMonotonicity
  



template class MeshDiagnostics<2>;
template class MeshDiagnostics<3>;


// ==================================================================================================
//
//   NON-MEMBER FUNCTIONS
//
// ==================================================================================================



/**
    Tries to identify duplicate elements by collocated barycentres true if there are duplicate elements in the model
*/
template<uint32_t dim>
bool detectCollocatedElements( const Model<dim>& m )
 {
    set<csmp::Point<dim> >  element_barycenters;
    bool                    coincident_barycenters(false);
    
    const Region<dim>&  gref(m.Region("Model"));
    
    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ ) {
          csmp::Point<dim> bc = (*it)->BaryCenter();
          pair<typename set<csmp::Point<dim> >::iterator,bool>
          bc_it = element_barycenters.insert( bc );
          // if we were not able to insert this barycenter point
          // this means it is already in the set and we are dealing
          // with a duplicate element
          if ( !bc_it.second ) {
               cerr <<"\ndetectDuplicateElements: detected element whose barycentre coincides with other element.\n";
               (*it)->Out();
               coincident_barycenters = true;
            }
      
      }
 
    return coincident_barycenters;
 
 } // detectCollocatedElements
 
 template bool detectCollocatedElements( const Model<1U>& );
 template bool detectCollocatedElements( const Model<2U>& );
 template bool detectCollocatedElements( const Model<3U>& );








} // end namespace csmp
