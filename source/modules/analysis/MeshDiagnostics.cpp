#include "MeshDiagnostics.h"
#include "ANSYS_ElementSpecifications.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "Standard_IO_Handler.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {

template<size_t dim>
void MeshDiagnostics<dim>::ElementVolumeRange( const Model<dim>& sg,
                                               double64& vmin, double64& vmax ) const
 {                                             
    double64    vol;
    bool  first_element(true); 
    
    const Region<dim>& sgref(sg.Region("Model"));
    
    for ( typename vector<Element<dim>*>::const_iterator
          eit=sgref.ElementsBegin(); 
          eit!=sgref.ElementsEnd(); eit++ ) 
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








template<size_t dim>
void MeshDiagnostics<dim>::FixFiniteElementNeighborOrientationOfSurfaceMeshes( Model<dim>& sg ) const
  {
    std::vector<Point<dim> > bc_vec(0);
    std::vector<size_t >     id_vec(0);
    double64                 sign;
    size_t                   id, cntr(0);
    
    cout <<"\nMeshDiagnostics::FixFiniteElementNeighborOrientationOfSurfaceMeshes: Checking for any clockwise ordered element neighbor IDs..." << endl;
    
    const Region<dim>& sgref(sg.Region("Model"));

    for ( typename vector<Element<dim>*>::const_iterator
          eit=sgref.ElementsBegin(); 
          eit!=sgref.ElementsEnd(); eit++ )
      { 
        if ( bc_vec.size() != (*eit)->Neighbors() ) bc_vec.resize((*eit)->Neighbors());
        sign = 0.;
        // construct a polygon of with vertices equal to the barycenters of the neighbor FEs
        // if there is no neighbor FE, use barycenter of current FE.
        for ( size_t i=0U; i<(*eit)->Neighbors(); i++ ) {
            if ( (*eit)->Neighbor(i) != NULL ) bc_vec[i] = (*eit)->Neighbor(i)->BaryCenter();
            else                               bc_vec[i] = (*eit)->BaryCenter(); 
          }
        for ( size_t i=0U; i<(bc_vec.size()-1U); i++ ) {
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
            for ( size_t i=0; i<(*eit)->Neighbors(); i++ ) {
                if ( (*eit)->Neighbor(i) != NULL ) id_vec[id-i] = (*eit)->Neighbor(i)->Idx();
                else                               id_vec[id-i] = 0;
              }
            for ( size_t i=0; i<id_vec.size(); i++ ) {
                if ( id_vec[i] > 0 ) (*eit)->Assign( i, sgref.E(id_vec[i]-1) );
                else                 (*eit)->Assign( i, static_cast<Element<dim>*>(NULL) );
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
template<size_t dim>
bool MeshDiagnostics<dim>::ScrutinizeMesh( Model<3U>& sg ) const
 {
    Region<3>&  sgroup(sg.Region("Model"));
    sgroup.UpdateMemberIndexes();
    double64 l_segm, segm_length_ratio(DBL_MAX),
           segm_length_ratio_min(DBL_MAX), segm_length_ratio_max(0.),
           sl_min, sl_max,
           sl_min_e(-4e7), sl_max_e(4e7), // the equator of the earth
           h_min(-4.0e7), h_max(4.0e7),
           v_min(-4.0e7), v_max(4.0e7),
           h_valmin, v_valmin,
           h_valmax, v_valmax;
    bool   problems = false;
           
    DenseMatrix<DM_MIN>    XY(3,3);
    vector<size_t>         nids(3);
    Standard_IO_Handler    stdio;

    // Node coordinates
    cout <<"\nMeshDiagnostics<dim>::ScrutinizeMesh: Verifying node coordinates..."<< endl;
    //cout <<"\nplease enter intended horizontal min, max coordinates: ";
    //cin >> h_min >> h_max;
    //cout <<"\nplease enter intended vertical min, max coordinates: ";
    //cin >> v_min >> v_max;
 
    h_valmin = h_valmax = sgroup.N(0)->x();
    v_valmin = v_valmax = sgroup.N(0)->y();    
 
    for ( vector<Element<3U>*>::const_iterator
          eit=sgroup.ElementsBegin(); eit!=sgroup.ElementsEnd(); eit++ )
      {
         (*eit)->NodeCoordinateMatrix( XY );
         for ( size_t i=0U; i<(*eit)->Nodes(); i++ )
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
    vector<double64> lengths;

    cout <<"\n\nMeshDiagnostics<dim>::ScrutinizeMesh: Testing length of element segments..."<< endl;
    
    (*sgroup.ElementsBegin())->SegmentLengths( lengths );
    sl_min = sl_max = lengths[0];

    for ( vector<Element<3U>*>::const_iterator
          eit=sgroup.ElementsBegin(); eit!=sgroup.ElementsEnd(); eit++ )
      {
         (*eit)->SegmentLengths( lengths );
         for ( size_t i=0U; i<lengths.size(); i++ )
           {
              l_segm = lengths[i];
              if ( l_segm > sl_max )   sl_max   = l_segm;
              if ( l_segm < sl_min )   sl_min   = l_segm;
              if ( l_segm > sl_max_e ) sl_max_e = l_segm;
              if ( l_segm < sl_min_e ) sl_min_e = l_segm;
           }
         segm_length_ratio = sl_max_e / sl_min_e;
         if ( eit == sgroup.ElementsBegin() ) segm_length_ratio_max = segm_length_ratio_min = segm_length_ratio;
         if ( segm_length_ratio > 10. )
           cout <<"\nlarge segment length ratio in element "<< (*eit)->Idx() <<": "<< segm_length_ratio;
          
         if ( segm_length_ratio > segm_length_ratio_max ) segm_length_ratio_max = segm_length_ratio;
         if ( segm_length_ratio < segm_length_ratio_min ) segm_length_ratio_min = segm_length_ratio;
      }


    // Element volume
    double64        volume;
    bool            repeat(true);
    vector<size_t>  element_numbers;
   
    cout <<"\n\nMeshDiagnostics<dim>::ScrutinizeMesh: Verifying element volumes/areas/lengths..."<< endl;
    for ( vector<Element<3U>*>::const_iterator
           eit=sgroup.ElementsBegin(); eit!=sgroup.ElementsEnd(); eit++ )
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
              cerr <<"\nMeshDiagnostics<dim>::ScrutinizeMeshErratic area/volume of element ",
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


template class MeshDiagnostics<2>;
template class MeshDiagnostics<3>;




///@return true if there are duplicate elements in the model
template<size_t dim>
bool detectDuplicateElement( const Model<dim>& m )
 {
    set<csmp::Point<dim> >  element_barycenters;
    
    const Region<dim>&  gref(m.Region("Model"));
    
    for ( typename vector<Element<dim>*>::const_iterator
        it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ ) {
          csmp::Point<dim> bc = (*it)->BaryCenter();
          pair<typename set<csmp::Point<dim> >::iterator,bool>
          bc_it = element_barycenters.insert( bc );
          // if we were not able to insert this barycenter point
          // this means it is already in the set and we are dealing
          // with a duplicate element
          if ( !bc_it.second ) return true;
      }
 
    return false;
 
 } // detectDuplicateElement
 
 template bool detectDuplicateElement( const Model<1U>& );
 template bool detectDuplicateElement( const Model<2U>& );
 template bool detectDuplicateElement( const Model<3U>& );


} // end namespace csmp
