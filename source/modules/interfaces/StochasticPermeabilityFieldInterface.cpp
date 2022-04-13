#include <sstream>
#include "StochasticPermeabilityFieldInterface.h"
#include "MJL_Point.h"
#include "InterFace.h"
#include "Region.h"
#include "Model.h"
#include "FemFromGridVisitor.h"
#include "PropertyHandle.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
StochasticPermeabilityFieldInterface<dim>::StochasticPermeabilityFieldInterface( Model<dim>& sg ) 
: nodal_perm( sg, "nodal permeability", SCALAR, NODE ),
  log_perm( sg, "log nodal permeability", SCALAR, NODE )
{}

template<uint32_t dim>
StochasticPermeabilityFieldInterface<dim>::~StochasticPermeabilityFieldInterface() {}


/**
 
Reads in a HYDROGEN (Bellin and Rubin, 1996; Fortran source code in folder 
/csp/varia/hydrogen) generated stochastic permeability field for a regular 
2D finite element grid where the node-spacing in x and y direction is constant
and the elements have the same volume. HYDROGEN input file must have the 
ending '*.perm' and be in x-y-z column format. Permeability values are
assigned to the nodes of the FE mesh and the nodal permeability values are
interpolated to the element center. If node coordinates and stochastic 
permeability coordinates do not coincide, the method terminates without
mapping the stochastic permeability field to the Model. A variable 'nodal 
logarithmic permeability' is created such that the log values of the 
permeability at the nodes can be visualized from the main() file.  

@section arguments Input Arguments 

A reference to the Model and the name of the HYDROGEN file that stores 
the permeability field.

@section application Application

Call the function from your main() file to assing the stochastic permeability
values to the Model.  

 
@test tested: O.K. */
template<uint32_t dim>
bool StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid( Model<dim>& sg, const char* fname ) 
 {
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid: Mapping 2D Permeability Field to Nodes" << endl;
    
    if ( dim != 2 ) {
        throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                     "Method works only for 2D meshes, nothing is done" );
        return false;
      }             
    
    typename std::map<mjl::Point, size_t>  coor_map, node_map;

    mjl::Point  xy;
    double      xval, yval, perm;
    char        file[200];
    std::string text_line;
    typename std::vector<std::string>  tokens;
    strcpy( file, fname );
    strcat( file, ".perm");
    ifstream ifs( file );
    ScalarVariable k, log_k;
    double vol;

    Region<dim>&  model_domain(sg.Region("Model"));
    
    // loop over all elements and check if all ahve the same volume. if yes, it is most likely that grid is regular
    auto eit = model_domain.CellsBegin();
    vol = (*eit)->Volume();
    for ( eit = model_domain.CellsBegin(); eit != model_domain.CellsEnd(); eit++ ) {
        if ( (*eit)->Volume() != vol ) {
            throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                     "Element volume varies, mesh appears to be not regular, exiting function, nothing is done" );
            cout << endl;
            return false;
          }
      }
    
    // store xy-coordinates and node ID in a map to find node id using xy-coordinates from input file
    for ( auto it = model_domain.NodesBegin(); it !=  model_domain.NodesEnd(); it++ ) {
        xy(0) = (*it)->x();
        xy(1) = (*it)->y();
        node_map.insert( make_pair( xy, (*it)->Idx() ) );
      }

    // read in file header and echo the information
    if ( !ifs.is_open() ) {
        throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                     "\nFile could not be opened... nothing is done" );
        cout << endl;
        return false;
      }
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid: \nEchoing File Header: " << endl << endl;
    for ( int i=0; i<6; i++ ) {
        getline( ifs, text_line );
        cout << text_line << endl;
        if ( i == 2 ) {
            if ( strcmp( text_line.c_str(), " the data are stored in column format (x y z)" ) != 0 ) {
                throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                         "\nHYDROGEN output file appears not to be in x y z column format, nothing is done" );
                cout << endl;
                return false;
              }
          }
      }
    
    // read in xy coordinates and log k from file, find node with same xy coordinates, store log k and k
    for ( size_t i{0U}; i<model_domain.Nodes(); i++ ) {
        ifs >> xval >> yval >> perm;
        xy(0) = xval;
        xy(1) = yval;
        auto mit = node_map.find( xy );
        if ( xy != (*mit).first ) {
            cout << "\nCoordinates x: " << xval << ", y: " << yval << endl;
            throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                     "\nCould not find the node coordinates, exiting function, nothing is done" );
            cout << endl;
            return false;
          }
        coor_map.insert( make_pair( xy, (*mit).second ) );
        k = std::pow( 10., perm );
        log_k = perm;
        model_domain.N( (*mit).second )->Store( nodal_perm.Key(), k );
        model_domain.N( (*mit).second )->Store( log_perm.Key(), log_k );
        //cout << "\nNode: " << (*mit).second << ", x: " << (*mit).first[0] << ", y: " << (*mit).first[1] << ", k: " << k();
        //cout << "\nNode: " << (*mit).second << ", x: " << (*mit).first[0] << ", y: " << (*mit).first[1] << ", log k: " << log k();
        
      }  
    
    // now check that each node as a permeability assigned to it, if this is not the case, assign NAN to each node for nk, logk  
    for ( auto it = model_domain.NodesBegin(); it !=  model_domain.NodesEnd(); it++ ) {
        xy(0) = (*it)->x();
        xy(1) = (*it)->y();
        auto mit = coor_map.find( xy );
        if ( xy != (*mit).first ) {
            cout << "\nCoordinates x: " << xy(0) << ", y: " << xy(1) << endl;
            throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid", 
                     "\nNodes where found that do not have a stochastic permeability assigned to them, nothing is done" );
            cout << endl;
            ScalarVariable nan(PLAIN,std::strtod("NAN",NULL));         
            for ( it =model_domain.NodesBegin(); it !=  model_domain.NodesEnd(); it++ ) {
                (*it)->Store( nodal_perm.Key(), nan );
                (*it)->Store( log_perm.Key(), nan );
              }
            return false;
          }
      }

                                                    
    model_domain.InterpolateNodeToCellProperty( "nodal permeability", "permeability" );
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityFieldOnRegularGrid: Successfully mapped permeability field to nodes and elements" << endl;
    
    return true;
 
 }


/**
 
Same as above bu maps the stochastic permeability field to an irregular
finite element mesh. A FiniteDifferenceGrid that holds the stochastic
permeability values is created and the FemFromGridVisitor is employed to
transfer the grid data to the finite element mesh. If the stochastic 
permeability field extend and the dimensions of the Model do not coincide,
the field is rescaled, in which case some of the accuracy of the stochastic field
is lost. Further, if finite elements are smaller than the grid resolution of
the FiniteDifferenceGrid such that no new permeability values can be assigned
to these elements, the values of their connected neighbor elements are averaged
to compute the new element permeability. In this case, accuracy of the stochastic
permeability field is lost too. If more than 1/3 of the finite elements have no
new permeability assigned, the method terminates and the original permeability
field is retained.   

 
@test tested: O.K. */
template<uint32_t dim>
bool StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField( Model<dim>& sg, const char* fname ) 
 {
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: Mapping 2D Permeability Field to Nodes" << endl;
    
    if ( dim != 2 ) {
        throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField", 
                     "\nMethod works only for 2D meshes, nothing is done" );
        return false;
      }             
      
    Region<dim>&  model_domain(sg.Region("Model"));
    
    typename std::vector<double>  xvec, yvec, kvec;
    mjl::Point xy;                                    
    double dx=std::numeric_limits<double>::quiet_NaN(), dy=std::numeric_limits<double>::quiet_NaN(), x, y, k, xmax(0.), ymax(0.);
    size_t index(0);
    char    file[200];
    std::string   text_line;
    typename std::vector<std::string>  tokens;
    strcpy( file, fname );
    strcat( file, ".perm");
    ifstream ifs( file );
    ScalarVariable nk, log_k, perm;
    PropertyHandle<dim> k_backup( sg, "backup permeability", SCALAR, ELEMENT );
    

    // read in file header and echo the information
    if ( !ifs.is_open() ) {
            throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField", 
                     "\nFile could not be opened... nothing is done" );
        return false;
      }
    
    // backup k if anything fails later down
    sg.CopyReplace( "permeability", "backup permeability" );

    // read in the file, echo header, extract dx and dy at 4th line
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: \nEchoing File Header: " << endl << endl;
    while ( getline( ifs, text_line ) ) {
        if ( index < 6 ) cout << text_line << endl;
        if ( index == 2 ) {
            if ( strcmp( text_line.c_str(), " the data are stored in column format (x y z)" ) != 0 ) {
                throw csmp::Exception( ERROR, "StochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField", 
                     "\nHYDROGEN output file appears not to be in x y z column format, nothing is done" );
                cout << endl;
                return false;
              }
          }
        if ( index == 3 ) {
            Tokenize( text_line, tokens );
            dx = atof(tokens[2].c_str());
            dy = atof(tokens[3].c_str());
          }
        // now read the x and y coordinates and k value  
        if ( index >= 6 ) {  
            Tokenize( text_line, tokens );
            x = atof(tokens[0].c_str());  
            y = atof(tokens[1].c_str());
            k = atof(tokens[2].c_str());
            if ( x > xmax ) xmax = x;
            if ( y > ymax ) ymax = y;
            xvec.push_back(x);
            yvec.push_back(y);
            kvec.push_back(k);
            
          }
        index++; 
      }
    
    double min, max;
    sg.MinMaxOf( "inner radius", min, max );
    if ( isnan(min) || isnan(max) ) {
        sg.AssignCellCharacteristicsTo( "inner radius", "inner radius" ); 
        sg.MinMaxOf( "inner radius", min, max );
      }
    if ( min < dx || min < dy ) {
        cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField:";
        cout << "\nInner radius of finite elements smaller than dx or dy of k-field, some elements may have no new permeability assigned to them " << endl;
      }
      
    // compare dimensions of k-field with Model dimensions and scale x-y resolution if necessary
    double xres, yres;
    Point<dim>  xyz_min, xyz_max;
    sg.MinMaxCoordinates( xyz_min, xyz_max );
    if ( xmax != xyz_max[0] || ymax != xyz_max[1] ) {
        cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: ";
        cout << "\nDimensions of Model and stochastic k-field are not the same: ";
        cout << "\nHorizontal extension: "<< xyz_min[0] << " to " << xmax;
        cout << "\nVertical extension: "<< xyz_min[1] << " to " << ymax;
        cout << "\nResizing k-field, correlation lengths may change now!" << endl;
        xres = xyz_max[1] / ( xmax / dx );
        yres = xyz_max[3] / ( ymax / dy );
      }
    else {
        xres = dx;
        yres = dy;
      }  
    
    // construct a finite difference grid with the Model's horizontal and vertical extension
    // and scaled (if necessary) resolution. Map all uniformly gridded k-values to the FD Grid
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: \nTransferring k-data to FiniteDifferenceGrid" << endl;
    FiniteDifferenceGrid k_field( xyz_max[0], xyz_max[1],  xres, yres );
    for ( size_t i{0U}; i<kvec.size(); i++ ) {
        k_field(static_cast<int>(xvec[i]/dx),static_cast<int>(yvec[i]/dy) ) = pow( 10., kvec[i] );
      }
    
    // now use the FemFromGridVisitor to map grid k-data to finite elements
    cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: \nTransferring FiniteDifferenceGrid to Finite Elements" << endl;
    FemFromGridVisitor<dim> fem_visitor( sg.Database(), k_field, "permeability", model_domain.Cells() );
    sg.Accept( fem_visitor ); 
    
    // find the min and max k and count the elements that have no new k assigned and are nan
    csmp::Index perm_key(sg.Database().StorageKey("permeability"));  
    double mink(1.0e+20), maxk(0.0);
    size_t nan_count(0);
    for ( auto eit = model_domain.CellsBegin(); eit !=  model_domain.CellsEnd(); eit++ ) {
        (*eit)->Read( perm_key, perm );
        if ( !isnan(perm()) ) {
            if ( perm() > maxk ) maxk = perm();
            if ( perm() < mink ) mink = perm();
          }
        else nan_count++;  
      }
    // if many elements have no random k assigned, it makes no sense defining a stochastic k-field,
    // use the originial k field instead and exit function
    if ( nan_count > model_domain.Cells() / 3 ) {
        cout << "\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField:";
        cout << "\n " << nan_count << " Elements of " << model_domain.Cells() << " have no new stochastic k value assigned";
        cout << "\nNo sense to employ a stochastic k-field any longer, re-using original k-field instead and exiting function, nothing is done" << endl;
        sg.CopyReplace( "backup permeability", "permeability" );
        return false;
      }  
    
   
    // check if any elements were too small such that no FD grid points are assigned and their k's are nan
    // if k is nan, average the k's of the neighbor elements. If no neighbor elements are given or k's are all nan
    // as well, assign random k
    double k_avg, counter;
    for ( auto eit = model_domain.CellsBegin(); eit !=  model_domain.CellsEnd(); eit++ ) {
        (*eit)->Read( perm_key, perm );
        if ( isnan( perm() ) ) {
            k_avg = counter = 0.;
            for ( auto i{0U}; i<(*eit)->Neighbors(); i++ ) {
                if ( (*eit)->Neighbor(i) != nullptr ) {
                    (*eit)->Neighbor(i)->Read( perm_key, perm );
                    if ( !isnan( perm() ) ) {
                        k_avg += perm();
                        counter += 1.;
                      }  
                  }
              }
            if ( k_avg > 0.0 ) perm = k_avg / counter;
            else               perm = ( maxk + mink ) / 0.5;
            (*eit)->Store( perm_key, perm );
          }
      }

    // extrapolate element permeability to nodes, compute log k's for visualization
    model_domain.ExtrapolateCellToNodeProperty( "permeability", "nodal permeability" );

    // store log k values for visualzation
    for ( auto it = model_domain.NodesBegin(); it !=  model_domain.NodesEnd(); it++ ) {
        (*it)->Read( nodal_perm.Key(), nk );
        log_k = std::log10( nk() );
        (*it)->Store( log_perm.Key(), log_k );
        
      }                                                
    
    cout <<"\nStochasticPermeabilityFieldInterface<dim>::Read2DStochasticPermeabilityField: ";
    cout <<"\nSuccessfully mapped permeability field to nodes and elements" << endl;
    
    return true;
 
 }

template<uint32_t dim>
void StochasticPermeabilityFieldInterface<dim>::Tokenize( const std::string& str, std::vector<std::string>& tkns ) 
 {
   std::string buf;
   tkns.erase( tkns.begin(), tkns.end() );
   std::stringstream ss(str);
   while ( ss >> buf ) tkns.push_back(buf);
 }

template class StochasticPermeabilityFieldInterface<2U>;
 
} // end namespace csmp 
 
 
 
 
 
 
 
