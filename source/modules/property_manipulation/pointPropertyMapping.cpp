//
//  pointPropertyMapping.cpp
//
//  Created by Stephan Matthai on 1/6/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//
//  Edited by Lukas Mosser on 24.01.2013
/* 
Edit 24.01.2013:
  Fixed a bug where the SAMG did not converge do to overconstraining.
  Fixed by including a check if the number of dirichlet nodes is equal to 
  the number of nodes in the region.

ToDo: List:
    Remove unnecessary maps and vector structures
    Trim code for higher speed.
    Combine Volumetric and Surface Mapping
    Port Volumetric Data Model to the surface model

    Lukas Mosser 2013 */

#include "pointPropertyMapping.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "Standard_IO_Handler.h"
#include "VTU_Interface.h"

#include "PDE_Integrator.h"
#include "LinearSolver.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#include "GaussJordan_Solver.h"
#include "LUdcmp_Solver.h"

#include "DenseMatrix.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "PDE_Integrator.h"

using namespace std;
using           std::ofstream;

namespace csmp {

const double ZERO_DATA_VALUE(1e30);

void propertiesToRegions( Model<3U>&, const char* model_name, const char* point_data_file, bool debug=false );

/**
   Mapping property values specified at points stored in textfile to named surfaces and volumes of input model
   
   Input file format (no blank lines nor whitepace allowed)
   
   NAME,X,Y,Z,porosity
   REGION_NAME,699874.918945,-4332.125000,-1042574.062500,0.018301
   ...
   
   @attention some properties have already been assigned and must not be overwritten).

*/   
void pointPropertyMapping( const string& model_name )
 {
    // 1. Reading the CSMP model onto which the properties shall be mapped
    // -------------------------------------------------------------------
    cout <<"\n\npointPropertyMapping: Enter name of textfile with point property values that shall be mapped onto model '"<< model_name <<"' (ommit '.txt'): ";
    string  property_file;
    cin >> property_file;
    property_file += ".txt";
    
    ifstream  my_file(property_file.c_str());
    while ( !my_file.good() ) {
           cerr <<"\n\tfile does not exist; enter name again: ";
           cin >> property_file;
           my_file.open(property_file.c_str());
       }
    my_file.close();

    csmp::Model<3U>  model(model_name);

    // 1. mapping property point data to unstructured grid
    // ---------------------------------------------------
    VTU_Interface<3U>  vtu_output( model, "pointPropertyMapping_test");
    const bool debug(false);
    propertiesToRegions( model, model_name.c_str(), property_file.c_str(), debug );
   
    const double phi_min = printRangeOfVariable( model, "porosity", false );
    if ( phi_min <= 0. ) {
         vtu_output.OutputDataToVTU( (model_name + "-erratic_porosity"), string("porosity"), string("Model"), 0 );
         throw csmp::Exception( ERROR, "pointPropertyMapping:", "negative/zero porosity values detected; terminating mapping process.");
      }
   
   
    // 2. checking whether permeability data may have been entered in mD rather than m2 and offering to convert it
    // -----------------------------------------------------------------------------------------------------------
    Standard_IO_Handler  stdio( (model_name + "pointPropertyMapping").c_str() );
    const double k_max = printRangeOfVariable( model, "permeability", true );
    if ( k_max > 1.0e-5 )
      if ( stdio.YesNo("pointPropertyMapping: Were the 'permeability' data entered in [mD]? - CSMP needs [m2]; convert to m2") ) {
           csmp::Index k_key(model.Database().StorageKey("permeability"));
           Region<3U>&  mref(model.Region("Model"));
           for ( auto it=mref.CellsBegin(); it!=mref.CellsEnd(); it++ )
              //                                                                 md -> m2
              (*it)->Store( k_key, makeScalar( (*it)->Status(k_key), (*it)->Read(k_key) * 1.0e-15 ) );
         }
   
   
    // 3. allowing the user to check the assigned property distribution
    // ----------------------------------------------------------------
    if ( model.Database().IsDefined("rocktype") ) {
         // converting rock-type values to integers
         const csmp::Index rkey(model.Database().StorageKey("rocktype"));
         Region<3U>& model_domain(model.Region("Model"));
         for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it )
           (*it)->Store( rkey, makeScalar( (*it)->Status(rkey), round((*it)->Read(rkey))) );
         vtu_output.OutputDataToVTU( (model_name + "-mapped_rocktype"), string("rocktype"), string("Model"), 0 );
      }
    if ( model.Database().IsDefined("porosity") )
      vtu_output.OutputDataToVTU( (model_name + "-mapped_porosity"), string("porosity"), string("Model"), 0 );
    if ( model.Database().IsDefined("permeability") )
      vtu_output.OutputDataToVTU( (model_name + "-mapped_permeability"), string("permeability"), string("Model"), 0 );
   
   
    // 4. saving model to disk if requested
    // ------------------------------------
    if ( stdio.YesNo("point_property_mapping: Do you want to save the model with the newly mapped values") )
       model.OutputToBinaryFile( model_name.c_str() );
    
    cout <<"\n\npointPropertyMapping: That's it!"<< endl;

} // end pointPropertyMapping


// SPECIAL DESIGN: ALL THE FOLLOWING UTILITY FUNCTION ARE KNOWN ONLY TO THIS CPP FILE


/// assign scalar variable value to all nodes of element
static void AssignToNodes( const csmp::Index& prop_key, Element<3U>* elmt, ScalarVariable sp )
{
  const size_t nodes(elmt->Nodes());
  for (uint32_t i = 0; i < nodes; ++i )
    elmt->N(i)->Store(prop_key,sp);
}


/// set variable at all nodes of the element to Dirichlet
static void MakeNodesDirichlet( const csmp::Index& prop_key, Element<3U>* elmt )
{
  const size_t nodes(elmt->Nodes());
  for(uint32_t i = 0;i < nodes;++i)
    elmt->N(i)->Status(prop_key,DIRICH);
}


//Added Check if the zero data value comes up then dont use for averaging
static pair<double,size_t> VectorSum(const vector<double>& vector){
  double sum = 0;
  size_t length = 0;
  const size_t vec_size(vector.size());
  for(size_t i = 0;i<vec_size;++i){
    if(vector[i]!=ZERO_DATA_VALUE){
      sum += vector[i];
      length++;
    }
  }
  return make_pair(sum,length);
}

static double VectorAverage(const vector<double>& vector){
  double avg;
  pair<double,size_t> sum_length = VectorSum(vector);
  avg = sum_length.first/sum_length.second;
  return avg;
}

static vector<double> AverageVectorOfVectors(const vector<vector<double > >& vecvecs){
  vector<double> averages;
  const size_t vec_size(vecvecs.size());
  for(size_t i = 0;i < vec_size;++i){
    averages.push_back(VectorAverage(vecvecs[i]));
  }
  return averages;
}



/** 
  
  point data input from text file
  
  Features Property Descriptors from first line of config File
  Return Property Model Map
 
  @attentionThis function supercedes all other File Input Functions except for the Region Indicators Function.

  LM 2013
  */
static void readPointDataFromFile( const char* filename,                      ///< text file containing ','-delimited values with a headerline with property names
                            vector<string>& property_descriptors, ///< names of the properties read from file
                            map<string,map<csmp::Point<3U>,vector<double> > >&  PropertyModel )
                            ///< target region name        point from file  array of property values
{
  cout <<"\n\nreadPointDataFromFile: reading point property data from file '"<< filename <<"'...\n";
  if ( !PropertyModel.empty() ) PropertyModel.clear();
  
  ifstream myfile(filename); //Ifstream definition
  string line; //Container for ifstream
  string item; //Container for item in line in IFStream

  size_t line_index = 0; //Line Number 
  size_t val_index = 0; //Value Index

  //Temporary Storage Variables for Points and Region Name
  string region_name;
  double xcord(numeric_limits<double>::quiet_NaN());
  double ycord(numeric_limits<double>::quiet_NaN());
  double zcord(numeric_limits<double>::quiet_NaN());
  double proptemp;
  vector<double> props; //Temporary Storage Vector for Point Properties
  map<csmp::Point<3U>,vector<double> > emptyMap; //An Empty Map to insert as a placeholder if a Region is not yet in the list
  //Iterator for the PropertyModel
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      //Clear the line variables
      props.clear();
      emptyMap.clear();
      val_index = 0;
      //Get the next line and copy it's values to line
      getline( myfile, line);
      if ( line.empty() ) break; //If a line is empty => Break
      //Create a Stringstream object that reads the items saved in "line"
      //Access "line" elements with getline(ss,item,',')
      stringstream ss(line);
      //Check if the line index is zero to get PropertyDescriptors
      if(line_index == 0){
        //Get the First Line Descriptors
        val_index = 0;
        while(getline(ss,item, ',')){
          //Descriptors Start at position 4 REGION, X, Y, Z, SCALAR1, etc.
          if(val_index > 3){
            //Iterate over the elements in this line
            //Check for numeric 0
            property_descriptors.push_back(item);
          }
          val_index++;
        }
        //increment to access the next element in the current line
      }
      else{
        //do This for each line except the first: Line Format: Fault Name, X Coord, Y Coord, Z, Coord, Scalar 1, Scalar 2 etc.
        //Iterate over the current lines elements
        while(getline(ss,item,',')){
          //If it's the first element in the current line, this is the region
          if(val_index==0){
            region_name = item; //Region Name
          }
          else if(val_index==1) xcord = atof(item.c_str()); // string_to_double(item); //XCoord
          else if(val_index==2) ycord = atof(item.c_str()); // string_to_double(item);//YCoord
          else if(val_index==3) zcord = atof(item.c_str()); // string_to_double(item);//ZCoord
          else{
            proptemp = atof(item.c_str()); // string_to_double(item);
            props.push_back(proptemp);//Scalars pushed into props vector for storage
          }
          val_index++;
        }
        //All the elements in the current line have been stored
        //Now check the map if the Region name has already been found.
        map<string,map<Point<3U>,vector<double> > >::iterator iter = PropertyModel.find(region_name);
        
        if(iter==PropertyModel.end()){
          //iterator didn't find region name inside map
          //Create a new pair
          emptyMap.insert(make_pair(csmp::Point<3U>(xcord,ycord,zcord),props));
          //Insert the newly created map into the property model
          PropertyModel.insert(make_pair(region_name,emptyMap));
        }
        else{
          //Found the fault name, add the pair to the map of the fault
          iter->second.insert(make_pair(csmp::Point<3U>(xcord,ycord,zcord),props));
        }
      }
      //Evaluation of Current Line in DataFile is Done
      line_index++; //Increment the Line_Index
    }
    myfile.close();
  }
  
  cout <<"\nreadPointDataFromFile: read point data from input file successfully.\n";
}




/**
    any of the input property data are checked to ascertain that the values to be assigned 
    to the model are indeed within the range specified in the property database.
    
    @author SKM
*/
static bool checkInputPropertyRanges( const PropertyDatabase<3U>& database,
                               const vector<string>& property_descriptors, ///< names
                               const map<string,map<csmp::Point<3U>,vector<double> > >&  pdata ) ///< region-by-region data
 {
    // property counter
    size_t property(0U);
    bool   failed(false);
   
    for ( auto dit=property_descriptors.begin(); dit!=property_descriptors.end(); ++dit ) {
         // property as identified by key
         // has the property been discretised on the model?
         if ( !database.IsDefined((*dit).c_str()) ) {
              cerr <<"\ncheckInputPropertyRanges: point property '"<< (*dit) <<"' is not defined in the property database.\n";
              throw csmp::Exception( ERROR, "checkInputPropertyRanges:", (*dit).c_str(), "not defined in database.");
           }
      
         cout <<"\ncheckInputPropertyRanges: checking point property: '"<< (*dit) <<"'\n";
         // are its values within range
         double pmin, pmax;
         database.RangeOf( (*dit).c_str(), pmin, pmax );
         // looping over the regions
         for ( map<string,map<csmp::Point<3U>,vector<double> > >::const_iterator
               rit=pdata.begin(); rit!=pdata.end(); ++rit ) {
              // looping over the regional data values of the point properties
              bool range_check_failed(false);
              for ( map<csmp::Point<3U>,vector<double> >::const_iterator
                    it=(*rit).second.begin(); it!=(*rit).second.end(); ++it ) {
                   // range check
                   if ( (*it).second[property] < pmin or
                        (*it).second[property] > pmax ) {
                        range_check_failed=true;
                        break;
                     }
                }
              // reporting failure
              if ( range_check_failed ) {
                   cerr <<"\n\tcheckInputPropertyRanges: value of '"<< (*dit).c_str();
                   cerr <<"' in region '"<< (*rit).first <<"' out of range specified in database";
                   cerr <<" ("<< pmin <<","<< pmax <<").\n";
                   failed = true;
                }
            }
         property++;
         
      }
   
    return failed;
   
 } // end checkInputPropertyRanges





/// changed return value to avoid copying
static void  AveragePropertiesOnElements( size_t, map<Element<3U>*,vector<vector<double> > >& elementPropertyMap,
                                          bool, map<Element<3U>*,vector<double> >&  elementPropertyAveragedMap ){
  
  if ( !elementPropertyAveragedMap.empty() ) elementPropertyAveragedMap.clear();
  
  for( map<Element<3U>*,vector<vector<double> > >::const_iterator it=elementPropertyMap.begin();it!=elementPropertyMap.end();++it){
    
    elementPropertyAveragedMap.insert(make_pair(it->first,AverageVectorOfVectors(it->second)));
  
  }
} // end  CollocateProperties





/// assigns element variable values to nodes, fixing their flags to DIRICH
static void SpreadElementPropertyToNodeProperty( map<Element<3U>*,vector<double> >& elementPropertyMap,
                                                  const csmp::Index elmtProperty,
                                                  const csmp::Index& nodeProperty )
 {
  for( map<Element<3U>*,vector<double> >::iterator
       it=elementPropertyMap.begin(); it!=elementPropertyMap.end(); ++it )
    {
        ScalarVariable  sp(DIRICH,it->first->Read(elmtProperty));
        for( auto i{0U}; i<it->first->Nodes(); ++i )
          it->first->N(i)->Store(nodeProperty,sp);
    }
}






static void AssignElementProperties( map<Element<3U>*,vector<double> >& elementPropertyMap,
                                      const list<string>&,
                                      const vector<csmp::Index>& property_keys )
{
  for( map<Element<3U>*,vector<double> >::iterator it=elementPropertyMap.begin();it!=elementPropertyMap.end(); ++it )
    for( auto i{0U}; i < property_keys.size(); ++i )
      if ( it->second[i] != ZERO_DATA_VALUE )
        it->first->Store( property_keys[i], makeScalar(PLAIN,it->second[i]) ); //Flag Plain if element has value, Will have ANY Flag if it does not.
}






/// returns an entire copy of the vector!
static void  InitializeEmptyElementVector( const csmp::Region<3U>& gref, const vector<csmp::Index>& property_keys,
                                    vector<vector<Element<3U>*> >& noElementPropertyVectors )
 {
  if ( !noElementPropertyVectors.empty() ) noElementPropertyVectors.clear();
  noElementPropertyVectors.resize(property_keys.size());
  vector<Element<3U>*>  noPropElements;
  
  for ( auto i{0U}; i<property_keys.size(); ++i )
   {
      noPropElements.clear();
      noPropElements.reserve(gref.Cells());
      for ( size_t j{0U}; j<gref.Cells(); ++j ){
        if( gref.E(j)->Status(property_keys[i]) == ANY){
          noPropElements.push_back(gref.E(j));
        }
      }
      noElementPropertyVectors[i].assign( noPropElements.begin(), noPropElements.end() );
   }
   
} // end InitializeEmptyElementVector





static void AssignDefaultValuesToElements( const vector<vector<Element<3U>*> >& elementVectors,
                                    const vector<csmp::Index>& property_keys,
                                    const vector<double>& defaultValues )
{
  Element<3U>* elmt(NULL);
  ScalarVariable sp;
  vector<Element<3U>*> vector;
  for(size_t i = 0;i<property_keys.size();++i){
    vector = elementVectors[i];
    sp() = defaultValues[i];
    for(size_t j = 0;j<vector.size();++j){
      elmt = vector[j];
      elmt->Store(property_keys[i],sp);
    }
  }
}





static void AssignElementFromNodalProperty( const vector<vector<Element<3U>*> >&  elementVectors,
                                            const vector<csmp::Index> nodal_keys,vector<csmp::Index>&  property_keys )
{
  Element<3U>* elmt(NULL);
  vector<Element<3U>*> vector;
  csmp::ScalarVariable propertyVal;
  double value;
  for(size_t i = 0;i<property_keys.size();++i){
    vector = elementVectors[i];
    for(size_t j = 0;j<vector.size();++j){
      value = 0;
      elmt = vector[j];
      for(uint32_t k = 0;k<elmt->Nodes();++k){
        value += elmt->N(k)->Read(nodal_keys[i]);
      }
      propertyVal() = (value/(elmt->Nodes()));
      elmt->Store(property_keys[i],propertyVal);
    }
  }
}




static void AssignElementPropertyFromElementList( const vector<vector<Element<3U>*> >& elementVectors,
                                                  const vector<csmp::Index>& temp_keys, const vector<csmp::Index>& property_keys  )
{
  Element<3U>* elmt(NULL);
  vector<Element<3U>*> vector;
  csmp::ScalarVariable propertyVal;
  for(size_t i = 0;i<property_keys.size();++i){
    vector = elementVectors[i];
    for(size_t j = 0;j<vector.size();++j){
      elmt = vector[j];
      propertyVal() = elmt->Read(temp_keys[i]);
      elmt->Store(property_keys[i],propertyVal);
    }
  }
}





static void AssignElementPropertyFromElementList( const vector<Element<3U>*>& elementVector,
                                           const csmp::Index& temp_key, const csmp::Index& property_key ,const csmp::Index& empty_elmt_key)
{
  Element<3U>* elmt(NULL);
  csmp::ScalarVariable propertyVal;
  for(size_t j = 0;j<elementVector.size();++j){
    elmt = elementVector[j];
    propertyVal() = elmt->Read(temp_key);
    elmt->Store(property_key,propertyVal);
    elmt->Store(empty_elmt_key, ScalarVariable(PLAIN,2.0));
  }
}




// replace with TriangularFacet - unit normal operation
static csmp::Point<3U> normalVectorElementUnitized( const vector<csmp::Point<3U> >& points )
{
  assert( points.size() >= 3 );
  csmp::Point<3U> p_1 = points[0];
  csmp::Point<3U> p_2 = points[1];
  csmp::Point<3U> p_3 = points[2];
  
  csmp::Point<3U> normalvector;
  vector<double> p1c = p_1.Coordinates();
  vector<double> p2c = p_2.Coordinates();
  vector<double> p3c = p_3.Coordinates();
  
  csmp::Point<3U> p21 = p_2 - p_1;
  csmp::Point<3U> p31 = p_3 - p_1;
  
  csmp::Point<3U> normal = crossProduct(p21,p31);
  normal.NormalizeLengthTo();
  return normal;
}



/**
     this implementation is deadly expensive - use full storage matrix and fill directly!
*/
static vector<double> computeSurfaceBarycentricCoordinates( const Element<3U>* elmt, csmp::Point<3U> point )
{
  assert( elmt != nullptr );
  vector<csmp::Point<3U> > elmtNs;
  csmp::Point<3U> coords;
  for( uint32_t i{0U};i<elmt->Nodes();++i){
    coords = elmt->N(i)->Coordinate();
    elmtNs.push_back(coords);
  }
  csmp::Point<3U> normalVector = normalVectorElementUnitized(elmtNs);
  //SAMG_Solver solver;
  GaussJordan_Solver solver;
  size_t matrixdim = 4;
  SparseMatrix LHS(matrixdim);
  uint32_t dimension = 3;
  vector<double>  RHS(dimension+1,0.0),RESULT(dimension+1,0.0);
  vector<double> n;
  vector<double> V0;
  vector<double> V1;
  vector<double> V2;
  
  for( uint32_t i = 0;i<dimension+1;++i){
    if(i!=dimension){
      n.push_back(normalVector[i]);
      V0.push_back(elmtNs[0][i]);
      V1.push_back(elmtNs[1][i]);
      V2.push_back(elmtNs[2][i]);
    }
    else{
    n.push_back(0);
    V0.push_back(1);
    V1.push_back(1);
    V2.push_back(1);
    }
  }
  for( uint32_t i = 0;i<dimension+1;i++){
    for(size_t j = 0;j<dimension+1;j++){
      if(j==0){
        LHS.Add(i,j,-1*n[i]);
      }
      else if(j==1){
        LHS.Add(i,j,V0[i]);
      }
      else if(j==2){
      
        LHS.Add(i,j,V1[i]);
      }
      else if(j==3){
        LHS.Add(i,j,V2[i]);
      }
      else{
        cout<<"Error in Assembly"<<endl;
      }
    }
  }//End Left Hand Side Assembly

  //Right hand side Assembly
  for( uint32_t i = 0;i<dimension+1;++i){
    if(i!=dimension){
      RHS[i] = point[i];
    }
    else{
      RHS[i] = 1;
    }
  }//End RHS Assembly
  
  solver.Solve(LHS,RHS,RESULT);
  
  //cout<<RESULT[0]<<endl;
//  double t = RESULT[0];
  double u = RESULT[1];
  double v = RESULT[2];
  double w = RESULT[3];
  
  //double d = u*u+v*v+w*w;
  csmp::Point<3U> barycentricCoords(u,v,w);
  //vector<double>
  //barycentricCoords.Out();
  return RESULT;
}

static Element<3U>* GetRandomElement(const csmp::Region<3U>& gref){
  size_t emax = gref.Cells();
  auto rand_value{ rand() % emax };
  Element<3U>* elmt = gref.E(rand_value);
  return elmt;
}


static void OutputVector(vector<double>& vec){
  const auto vec_size{vec.size()};
  for( auto i{0U}; i < vec_size;++i ){
    cout << "Position: " << i << " Value: "<<vec[i]<<" ";
  }
  cout << endl;
}


static double SumVec(const vector<double>& vec){
  double sum = 0;
  const auto vec_size(vec.size());
  for(auto  i{0U}; i < vec_size; ++i){
    sum = sum + vec[i];
  }
  return sum;
}



/** This function evaluates the barycentric coordinates to find the direction of the next element.
The next element as well as a boolean is returned as a pair.
The boolean indicates if the point is inside the current element.

LM 2013 */
static pair<Element<3U>*,bool>  EvaluateBarycentricCoordinates( Element<3U>* element, vector<double>& RESULT, size_t nodes )
{
  double e = 1E-4;
  vector<double> pos_position;
  pair<Element<3U>*,bool> return_vals;
  Element<3U>* next_element(NULL);
  vector<double> coordinates;
  
  for(auto i{1U}; i <=nodes;++i){
    coordinates.push_back(RESULT[i]);
  }
  for(auto i{0U}; i < coordinates.size();++i){
    if(coordinates[i] >= (0.0-e) && coordinates[i] <= (1.0+e) ){
      pos_position.push_back(i);
    }
  }
  if(pos_position.size() == nodes){
    return_vals = make_pair(element,true);
  }else{
    next_element = element->Neighbor( static_cast<uint32_t>(distance(coordinates.begin(),min_element(coordinates.begin(),coordinates.end()))));
    return_vals = make_pair(next_element, false);
  }
  return return_vals;
}






/// comment on type vs. operator
static Element<3U>* FindPointIn3DVolumetricRegion(csmp::Point<3U> pXYZ, const Region<3U>&, Element<3U>* elmt )
{
  const uint32_t dimension = 3;
  
  vector<double>  RHS(dimension+1,0.0),RESULT(dimension+1,0.0);
  //SAMG_Solver solver;
  LUdcmp_Solver solver;
  Element<3U>* element(elmt);
  Element<3u>* next_element(NULL);
  size_t maxiter = 2000;
  bool stop = false;
  Node<3U> *node;
  pair<Element<3U>*,Node<3U>*> return_vals;
  size_t num_nodes;
  map<double,csmp::Element<3U>*> face_list;
  size_t iteration = 0;
  
  vector<double> pos_position;
  ScalarVariable sp;
  while(!stop){
    if(!element) {
          cerr << "\nNULL\n";
          break;
      }
    num_nodes = element->Nodes();
    SparseMatrix LHS(num_nodes);
    vector<double> coords = pXYZ.Coordinates();
    for(size_t i{0U};i <= dimension;i++){
      if(i<3){
        RHS[i] = coords[i];
      }
      else{
        RHS[i] = 1.0;
      }      
    }
    for( size_t i{0U}; i<=dimension; i++ ) RESULT[i] = 0.;
    
    //create LHS vector from tetrahedron nodes coordinates in node_coordinates as  well as add the barycenter constraint Epsilon(Lambda_i) = 1
    for(size_t i=0U;i < dimension+1;i++){
      for(auto j{0U}; j < num_nodes;j++){//finish setting up the matrix
        if(i==0U){
        // TODO:  this access of the element's node coordinates with a global node number looks like an ERROR ?
          LHS.Add(i,j,element->N(j)->x());
        } 
        else if(i==1U){
          
          LHS.Add(i,j,element->N(j)->y()); 
          
        } 
        else if(i==2U){
          
          LHS.Add(i,j,element->N(j)->z());  
          
        }
        else if(i==3U){
          
          LHS.Add(i,j,1.0);
        
        }
        else{
          cout<<"Error in LHS construction"<<endl;
        }
      
      }
    }
    //solve system of linear equations resulting the barycentric coordinates of the current point 
    //relative to the current element
    solver.Solve(LHS,RHS,RESULT);
    //check if all three barycentric coordinates are positive => the point is inside the element    
    pos_position.clear();
    for(size_t i{0U}; i<RESULT.size();i++){
      if(RESULT[i]>=0.0){
        pos_position.push_back(i);
      }
      else{
      }
    }
    for( uint32_t i{0U}; i<RESULT.size(); ++i ){
 // TODO: fix error: this uses i to access neighbors which will definitely create overflows
      if ( i >= element->Neighbors() )
        throw csmp::Exception( ERROR, "FindPointIn3DVolumetricRegion","out of range access of element neighbor.");
      face_list.insert(make_pair(RESULT[i],element->Neighbor(i)));
    
    }
    
    for(map<double,csmp::Element<3U>*>::iterator p = face_list.begin();p!=face_list.end();++p){
      //cout<<"Barycentric Coordinate: "<<p->first<<endl;
    }
    if (pos_position.size() == 4){
      auto n_assign = distance(RESULT.begin(), max_element (RESULT.begin(),RESULT.end()));
      node = element->N( static_cast<uint32_t>(n_assign) );
      return_vals = make_pair(element,node);
      //cout << "point found!"<<endl;
      stop = true;
    }
    else{
      next_element = element->Neighbor( static_cast<uint32_t>(distance(RESULT.begin(), min_element (RESULT.begin(),RESULT.end()))) );
      
      if(!next_element){
        return_vals = make_pair(element,node);
       
        return element;
       
      }else{
        element = next_element;
      }
    //throw csmp::Exception( ERROR, "FindNearestNode", "Null neighbor result" ); //Null ¸berpr¸fung ausserhalb der schleife und dort dann continue bei while setzen
    }
    iteration++;
    if(iteration==maxiter){
      element = NULL;
      return_vals = make_pair(element,node);
      stop=true;
    }
  } // while
  return element;
}








static Element<3U>*  FindPointIn3DSurfaceTetraMesh( const csmp::Point<3U>& point,
                                                    csmp::Region<3U>& surface,
                                                    Element<3U>* elmt )
{
  const bool debug = false;
  //surface.InputPropertyValue(iter_prop,ScalarVariable(PLAIN,0.0),COMPLETE);
  Element<3U>* nullElmt(nullptr);
  vector<Element<3U>* > path;
  Element<3U>* lastElement{nullptr};
  vector<double> maxiter_last;
  vector<double> maxiter_current;
  Element<3U>* containingElement{nullptr}; //Will contain the Element that contains the point being searched for
  pair<Element<3U>*,bool> traversingValues;
  csmp::Point<3U> barycentricCoordinates;
  vector<double> RESULT;
  vector<uint32_t> iterations;
  const size_t maxiter = 250; //Do a maximum of 1000 Steps inside the mesh
  size_t i = 0;
  
  if(!elmt){//If the specified elmt vector was NULL then choose a random element in the region to start with.
   elmt = GetRandomElement(surface);
  }

  traversingValues = make_pair(elmt,false);
  while(traversingValues.second == false){//Calls the following functions as long as the Evaluation of B.Coords returns false, meaning point outside current elmt.
    RESULT = computeSurfaceBarycentricCoordinates(traversingValues.first,point);
    traversingValues = EvaluateBarycentricCoordinates(traversingValues.first, RESULT,traversingValues.first->Nodes());
    //DEBUG BLOCK
    containingElement = traversingValues.first;
    if(traversingValues.first){
      
      //containingElement->Store(iter_key,ScalarVariable(PLAIN,i));
    }
    else{}
    //DEBUG BLOCK END
    if(!(traversingValues.first)){//If a Null Point is returned, return a Null Element, and the Point could not be found.
      //Here probably the problem arises  with points not being found. Check this, maybe restart improves the quality.
      if(debug){
        cout << "Null Element Returned " << "Iterations to here: " << i << endl;
        //vtu.OutputDataToVTU("Debug_Null",iter_prop,surface);
      }
      return nullElmt;
    }
    i++;
    
    if(i==maxiter){
      if(debug){
        maxiter_last = computeSurfaceBarycentricCoordinates(lastElement,point);
        maxiter_current = computeSurfaceBarycentricCoordinates(containingElement,point);
        cout << "Max Iter Reached" << endl;
        cout << "Outputing Barycentric Coordinates" << endl;
        cout << "Barycentric Coordinates of Last Element" <<endl;
        OutputVector(maxiter_last);
        cout << "Barycentric Coordinates of Current Element" <<endl;
        OutputVector(maxiter_current);
        cout << endl;
        //vtu.OutputDataToVTU("Debug_Maxiter",iter_prop,surface);
      }
      return nullElmt;
    }
    lastElement = containingElement;
  }
  
  //The search was succesfull and the element is returned.
  containingElement = traversingValues.first;
  //cout << "Iterations needed: " << i <<endl;
  //vtu.OutputDataToVTU("Debug_Normal",iter_prop,surface);
  return containingElement;
}







/** Output a list of points to a csv file
Lukas Mosser 2013*/
static void OutputPointsToCSV( vector<csmp::Point<3U> >& point_vec )
{
  if ( point_vec.empty() ) {
       cerr <<"\nOutputPointsToCSV: dataset does not contain any points; nothing was done.\n";
       return;
    }

  ofstream outdata("points_not_found.csv");
  
  outdata<<"x coord, y coord, z coord"<<endl;
  //for (it=point_vec.begin();it!=point_vec.end();it++){
  for(size_t i = 0;i<point_vec.size();++i){
    auto x = point_vec[i].Coordinates()[0];
    auto y = point_vec[i].Coordinates()[1];
    auto z = point_vec[i].Coordinates()[2];
    outdata << x;
    outdata <<',';
    outdata << y;
    outdata << ',';
    outdata << z <<endl;
  }
  outdata.close();
}



static void OutputVectorToCSV(vector<size_t>& vec,string file_name){
  ofstream outdata; // outdata is like cin
  outdata.open(file_name.c_str());// opens the file
  for( size_t i = 0;i<vec.size();++i){
    outdata << vec[i]<<endl;
  }
  outdata.close();
}





static void OutputMapToCSV(map<double,size_t>& dtmap,string file_name){
  ofstream outdata; // outdata is like cin
  outdata.open(file_name.c_str());
  for(auto it = dtmap.begin();it!=dtmap.end();++it){
    outdata << it->first << "," << it->second << endl;
  }
  outdata.close();
}






/** This function Remaps a vector<double> Structure to a vector<vector<double>> structure
  LM 2013
*/
static vector<vector<double> > RemapToVectorOfProperties(const vector<double>& vec){
  vector<vector<double> > vecvecs;
  vector<double> temp;
  for(size_t i = 0;i < vec.size();++i){
    temp.clear();
    temp.push_back(vec[i]);
    vecvecs.push_back(temp);
  }
  return vecvecs;
}






/** 

This function will find the points in points variable inside the mesh and add the element pointer
to the element that includes that point to its argument map. The properties assigned to this point will be 
assigned to this element in the map. Multiple points inside one element each contribute their vector.
processing this vector of property data vectors will be handled by another function.

Lukas Mosser 2013
*/
static void  MapPointsTo3DVolumetricRegion( const map<csmp::Point<3U>,vector<double> >& points,
                                     csmp::Region<3U>&  gref, map<Element<3U>*,vector<vector<double> > >& elementMap )
{
  
  Element<3U>* elmt(0);
  vector<Point<3U> > points_not_found;
  
  int not_found = 0;
  
  Point<3U> currentP;
  bool restart = false;
  
  auto points_size = points.size();
  cerr<<"total points to map = "<<points_size<<endl;
  size_t count = 0, i=0;
  
  for(auto pts = points.begin(); pts!=points.end() ; ++pts){
    //monitoring the progress
    if (count>=points_size/10*i){
        cerr<<"processed "<<round(double(count)/double(points_size)*100.0)<<" %"<<endl;
        i++;
    }
        
    if(pts==points.begin()||restart==true){
      elmt = GetRandomElement(gref);
      restart = false;
    }
    currentP = pts->first;
    
    //Find the point
    elmt = FindPointIn3DVolumetricRegion(currentP,gref,elmt);
    
    //If It returned a null vector start again from someplace random and go to the next point. 
    if(!elmt){
      not_found++;
      points_not_found.push_back(currentP);
      restart = true;
      count++;
      continue;
    }//End NULL If
    else{
      //If the point is found do this
      //Check if the Element is already in the elementPropertyMap
      map<Element<3U>*,vector<vector<double> > >::iterator it;
      it = elementMap.find(elmt);
      
      //The element was found in the elementPropertyMap
      if(it!=elementMap.end()){
        
        //Element has been found
        //Each point has a vector with its properties This needs to be split into the structure PropertyVector=>Vector of Properties for this element.
        for( i = 0; i < pts->second.size();++i){
          it->second[i].push_back(pts->second[i]);
        }
      }else{
        //Element hasn't been found, Insert the element as a new pair
        elementMap.insert(make_pair(elmt,RemapToVectorOfProperties(pts->second)));
      }
    }
    count++;
  }
  
  if ( not_found > 0 )
    cout << "\n\nMapPointsToElements3DSurfaceTetraMesh: "<< not_found <<" points of "<< points.size()<<" could not be found."<<endl;
  else 
    cout << "\n\nMapPointsToElements3DSurfaceTetraMesh: "<< points.size()<<" mapped successfully."<< endl;
  
  OutputPointsToCSV(points_not_found);
} // end FindPointsInRegion




/**
    This function will find the points in points variable inside the mesh and add the element pointer
  to the element that includes that point to a map. The properties assigned to this point will be placed
  assigned to this element in the map. Multiple points inside one element each contribute their vector.
  processing this vector of property data vectors will be handled by another function.
  Lukas Mosser 2012

*/
static void  MapPointsToElements3DSurfaceTetraMesh( const map<csmp::Point<3U>, vector<double> >& points,
                                             Region<3U>& matrix, map<Element<3U>*,vector<vector<double> > >& elementMap)
{
  
  if ( !elementMap.empty() ) elementMap.clear();
  
  //Container Vector for an Elements Properties
  vector<vector<double> > elmtPropVec;
  
  //Container for Points that aren't found
  vector<csmp::Point<3U> > points_not_found;
  Element<3U>* elmt(NULL);
  csmp::Point<3U> currentP;
  
  int not_found = 0;
  
  bool restart = false;
  
  //Iterator iterates over all points in in the point map
  for( map<csmp::Point<3U>,vector<double> >::const_iterator pts = points.begin(); pts!=points.end() ; ++pts){
    
    if(pts==points.begin()||restart==true){
      elmt = GetRandomElement(matrix);
      //restart = false;
    }
    
    //End first iteration if
    currentP = pts->first;
    elmt = FindPointIn3DSurfaceTetraMesh(currentP, matrix, elmt);
    
    //If the point is not found do this !elmt checks for Null elmt, which means the point lies outside of the current region.
    if(!elmt){
      not_found++;
      points_not_found.push_back(currentP);
      //restart=true;
      continue;
    }
    else{
      
      //If the point is found do this
      //Check if the Element is already in the elementPropertyMap
      map<Element<3U>*,vector<vector<double> > >::iterator it;
      it = elementMap.find(elmt);
      
      //The element was found in the elementPropertyMap
      if(it!=elementMap.end()){
        
        //Element has been found
        //Each point has a vector with its properties This needs to be split into the structure PropertyVector=>Vector of Properties for this element.
        
        for(size_t i = 0; i < pts->second.size();++i){
          it->second[i].push_back(pts->second[i]);
        }
      }else{
        //Element hasn't been found, Insert the element as a new pair
        elementMap.insert(make_pair(elmt,RemapToVectorOfProperties(pts->second)));
      }
    }
  }
  
  if ( not_found > 0 )
    cout << "\n\nMapPointsToElements3DSurfaceTetraMesh: "<< not_found <<" points of "<< points.size()<<" could not be found."<<endl;
  else 
    cout << "\n\nMapPointsToElements3DSurfaceTetraMesh: "<< points.size()<<" mapped successfully."<< endl;
  
  OutputPointsToCSV(points_not_found);
  
} // MapPointsToElements3DSurfaceTetraMesh






/// establishes the dimensionality of the regions that are to be mapped
static void ReadRegionIndicatorsFromFile( const char* filename, map<string,bool>& indicatorMap )
{
  if ( !indicatorMap.empty() ) indicatorMap.clear();
  
  
  ifstream myfile(filename);
  string line;
  string item;
  
  size_t line_index = 0;
  size_t value_index = 0;
  
  bool indicator_bool;
  
  string name;
  
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      //do This for each line, Line Format: Fault Name, X Coord, Y Coord, Z, Coord, Scalar 1, Scalar 2 etc.
      getline( myfile, line);
      if ( line.empty() ) break;
      stringstream ss(line);
      value_index = 0;
      //Iterate over the current lines values
      while(getline(ss,item,',')){
        
        if(value_index==0){
          //First Value is the Region Name
          name = item; 
        }
        else if(value_index==1){
          //Second Value is the Indicator: 0 For Volumetric and 1 For Surface in 3D Mesh
          indicator_bool = (atoi(item.c_str()) == 1) ? true : false;
        }
        value_index++;
      }
      indicatorMap.insert(make_pair(name,indicator_bool));
      line_index++;
    }
    myfile.close();
  }
}





static void OutPutPropertyModel( const map<string,map<csmp::Point<3U>,vector<double> > >& PropertyModel )
{
  for( map<string,map<csmp::Point<3U>,vector<double> > >::const_iterator
       it = PropertyModel.begin();it!=PropertyModel.end();++it ){
    cout<<"\n\nOutPutPropertyModel: Outputing Data for Fault "<<it->first << endl;
    cout<<"\n\tNumber of points for this surface: "<<it->second.size()<<endl;
  }
}




/**

  Checks if the number of Dirichlet Nodes is equal to the number of nodes in the region.
  If so, the region is fully constrained and no property interpolation is necessary.
 
*/
static size_t degreesOfFreedom( const Region<3U>& gref, csmp::Index prop_key )
 {
   size_t num = 0;
   for ( size_t j{0U}; j<gref.Nodes(); ++j )
    if(gref.N(j)->Status(prop_key) == DIRICH) num++;

   cout << "\ndegreesOfFreedom: N nodes: "<< gref.Nodes() <<" vs. Dirichlet nodes: "<< num << endl;
   return gref.Nodes() - num;
}





/**
    Adaption of Lukas 2D mapping main file

  Tested with:
   
       modelname = "doimoi_reverse";
       variablesfile = "RSP-variables.txt";
       property_file = "doimoi_distance.csv";
  Edited by Lukas Mosser 2013
*/
/*
list<string> ConvertToList(const vector<string>& property_descriptors){
  list<string> name_list;
  for(auto i = 0; i < property_descriptors.size();++i){
    name_list.push_back(property_descriptors.at(i));
  }
  return name_list;
}
*/

static vector<csmp::Index> DatabaseKeysFromNameList( Model<3U>& model, const vector<string>& property_descriptors ){
  vector<csmp::Index> keys;
  csmp::Index property_key;
  for(size_t i = 0; i < property_descriptors.size();++i){
    property_key = model.Database().StorageKey(property_descriptors[i].c_str());
    keys.push_back(property_key);
  }
  return keys;
}


/**
    Zeroing out property values in the target model subregions.
*/
static void InitializePropertiesInDatabase( Model<3U>& model, const char* region_name, const vector<string>& property_descriptors )
 {
    for ( size_t i{0U}; i<property_descriptors.size(); ++i ) {
         Region<3U>& ref = model.Region(region_name);
         auto property_key = model.Database().StorageKey(property_descriptors.at(i).c_str());
         switch (property_key.type) {
             case SCALAR:
             {
               ref.InputPropertyValue( property_descriptors.at(i).c_str(), makeScalar( ANY, 0. ) );
               break;
             }

             case TENSOR:
             {
               TensorVariable<3u> var(ANY, 0.);
               ref.InputPropertyValue( property_descriptors.at(i).c_str(), var );
               break;
             }

             default:
             {
                throw csmp::Exception( ERROR, "InitializePropertiesInDatabase", "Can't initialise property", property_descriptors.at(i) );
             }
         }
      }
 }



static void GetElementsAndCurrentPropertyMap( const size_t& index,map<Element<3U>*,vector<double > >& elementAveragedPropertyMap,
                                              map<Element<3U>*,double>& elementAndProperty)
{
  for(map<Element<3U>*,vector<double > >::iterator it = elementAveragedPropertyMap.begin(); it != elementAveragedPropertyMap.end();++it){
    elementAndProperty.insert(make_pair(it->first,it->second[index]));
  }
}



/**
    Overwrites the element property values in the specified target regions with the data
    supplied through the point-data file and the variables specified therein.
    
    ToDo: propertiesToRegions: break this spaghetti up into reasonable parts
*/
void propertiesToRegions( Model<3U>& model, const char*, const char* point_data_file, bool debug )
{
  // ---------------------------------------------------------
  // 1. reading point property data, and setting up storage
  //    for temporary variables
  // ---------------------------------------------------------
  const char* scalar_node = "scalar nodal container";
  const char* scalar_elmt = "scalar element temp container";
  const char* empty_elmt  = "empty element";
  
  // creating temporary property values (they have function scope so that variables will be deleted when function exits)
  PropertyHandle<3U>  scn( model, scalar_node,  SCALAR, NODE );
  PropertyHandle<3U>  sec( model, scalar_elmt, SCALAR, ELEMENT );
  PropertyHandle<3U>  sem(model,empty_elmt,SCALAR, ELEMENT);
  csmp::Index scalar_node_key = model.Database().StorageKey(scalar_node);
  csmp::Index scalar_elmt_key = model.Database().StorageKey(scalar_elmt);
  csmp::Index empty_elmt_key = model.Database().StorageKey(empty_elmt);

  model.InputPropertyValue( scalar_node, makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( scalar_elmt, makeScalar( PLAIN, 0.0 ) );
  model.InputPropertyValue( empty_elmt, makeScalar(PLAIN, 0.0) );
  
  // read file, initializing propertyPointData descriptors and keys 
  vector<string>                                property_descriptors;
  map<string,map<csmp::Point<3U>,vector<double> > >  propertyPointData;
  
  readPointDataFromFile( point_data_file, property_descriptors, propertyPointData );
  
  checkInputPropertyRanges( model.Database(), property_descriptors, propertyPointData );
  
  vector<csmp::Index>  property_keys = DatabaseKeysFromNameList(model, property_descriptors);
  list<string>              property_names(property_descriptors.begin(),property_descriptors.end());
  //set default values for given properties to zero
  vector<double>          defaultValues( property_descriptors.size(), 0. );
  
  if(debug){
    //Outputs the Property descriptors as well as the Property Model
    for(size_t i = 0; i < property_descriptors.size(); ++i) cout << property_descriptors.at(i) << endl;
    OutPutPropertyModel(propertyPointData);
    for( map<string,map<csmp::Point<3U>,vector<double> > >::iterator
         it = propertyPointData.begin(); it!=propertyPointData.end(); ++it ) {
      cout << "\npropertiesToRegions: OutPuting Property Model for Region: "<<it->first<<endl;
      for(map<Point<3U>,vector<double> >::iterator iter = it->second.begin();iter!=it->second.end();++iter){
        cout << "Current Point Info: " << endl;
        iter->first.Out();
        for(size_t i = 0;i < iter->second.size();++i)
          cout << "Property Value At Current Point: " << iter->second[i] <<endl;
      }
    }
  }
  
  
  // ---------------------------------------------------------
  // 2. property mapping
  // ---------------------------------------------------------
  // LHS matrix for Laplacian interpolation
  NumIntegral_dNT_dN_dV<3U>  stiffness_matrix(model.Database(),scalar_node,scalar_node);
  // RHS vector integral for RHS = 0 (homogeneous boundary conditions)
  NumIntegral_SetRHS_to_Zero<3U>  source(model.Database(),scalar_node);

  const size_t MIN_DOF_TO_INVOKE_SOLVER(100);
  
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings settings;
  settings.SetSolverInstance(1);
  settings.Set_eps(0.);
  settings.Set_iswit(5);
  
  SAMG_Solver solver( &settings );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  solver;
#endif

  // for all regions for which mapping needs to be carried out
  for ( map<string,map<csmp::Point<3U>,vector<double> > >::iterator
        it=propertyPointData.begin(); it!=propertyPointData.end(); ++it )
    {
       cout << "\nProcessing region: "<< (*it).first.c_str() << endl;
       Region<3U>&  gref(model.Region((*it).first.c_str()));
       map<csmp::Point<3U>,vector<double> >&  pointsAndProperties(it->second);
       // zeroing property values in the target region
       InitializePropertiesInDatabase( model, (*it).first.c_str(), property_descriptors );    
    
      // ---------------------------------
      // PROCESSING THE VOLUMETRIC REGIONS
      // ---------------------------------
      // Identifying whether we are dealing with a volumetric region
      // (assuming that regions only consist of a single dimensionality of elements)
      if ( (*gref.CellsBegin())->FE()->IsVolume() )
        {
           map<Element<3U>*,vector<vector<double> > >  mapped;
           MapPointsTo3DVolumetricRegion( pointsAndProperties, gref, mapped ); // ToDo: Bug: properties are not put on diferent vectors in map??

           if (debug){
             cout << "Outputing Non-Averaged Property Map"<<endl;
             map<Element<3U>*,vector<vector<double> > >::iterator map_it;
             for(map_it = mapped.begin();map_it !=mapped.end();++map_it){
               cout << "Outputing Values for this Element" << endl;
            for(size_t i = 0;i<map_it->second.size();++i){
              cout << "Ouputing Values for this Property"<<endl;
              for(size_t j = 0;j<map_it->second[i].size();++j)
                cout << "Value: "<<map_it->second[i][j]<<endl;
            }
          }
        }
     
      // averaging property values on the elements
      map<csmp::Element<3U>*,vector<double > > elementAveragedPropertyMap;
      AveragePropertiesOnElements(property_descriptors.size(),mapped,true,elementAveragedPropertyMap);
      if (debug){
            cout << "Output of averaged property map: " <<endl;
            for( map<Element<3U>*,vector<double> >::iterator
                 avg_it = elementAveragedPropertyMap.begin(); avg_it !=elementAveragedPropertyMap.end();++avg_it )
            for(size_t i = 0;i<avg_it->second.size();++i)
              cout << "value: "<<avg_it->second[i]<<endl;
        }
      
      AssignElementProperties(elementAveragedPropertyMap,property_names,property_keys);
      
      vector<vector<Element<3U>*> > noPropertyElementVectors;
      InitializeEmptyElementVector( gref, property_keys, noPropertyElementVectors );
      AssignDefaultValuesToElements(noPropertyElementVectors, property_keys, defaultValues);
      
      cout << "\n\tStarting the interpolation of "<< property_keys.size() <<" properties onto the volumetric regions...\n";
      for( size_t i{0U}; i<property_keys.size(); ++i )
        {
          cout <<"Interpolating property: "<<property_descriptors[i]<<endl;
          gref.InputPropertyValue(scalar_node, makeScalar(PLAIN,0.0));
          gref.InputPropertyValue(scalar_elmt, makeScalar(PLAIN,0.0));
            
          cout << "Number of elements in surface mesh: "<< gref.Cells()<<endl;
          cout << "Number of elements without current property assigned."<<noPropertyElementVectors[i].size()<< endl;
          
          // overwrites property values at element nodes with the value of the corresponding property stored on the element 
          SpreadElementPropertyToNodeProperty( elementAveragedPropertyMap, property_keys[i], scalar_node_key );

          const size_t dof = degreesOfFreedom(gref,scalar_node_key);
          if ( dof >= MIN_DOF_TO_INVOKE_SOLVER )
            {
               cout << "\n\tUsing SAMG solver for Laplace interpolation in volumetric domain..." << endl;
               PDE_Integrator<3U,Element> property_int_samg( solver );
               property_int_samg.Add(&stiffness_matrix);
               property_int_samg.Add(&source);            
               // Laplace interpolatation of unknown values
               property_int_samg.IntegrateOver( model, gref );
            }
          //Interpolate nodal values to elements using the element temp variable container
          gref.InterpolateNodeToCellProperty(scalar_node,scalar_elmt);
            
          /*Assign the known and interpolated values*/
          AssignElementPropertyFromElementList(noPropertyElementVectors[i],scalar_elmt_key,property_keys[i],empty_elmt_key);
          
          // doing a nearest neighbour fill to eliminate potential left-over no-data values
          if ( dof >= 1 ) nearestNeighborFill( model, (*it).first.c_str(), property_descriptors[i].c_str(), ZERO_DATA_VALUE );
       }
    }


      // ------------------------------------
      // PROCESSING REGIONS THAT ARE SURFACES
      // ------------------------------------
      else if ( (*gref.CellsBegin())->FE()->IsSurface() )
        {
        cout << "\nMapping point data to surface mesh." << endl;
        map<Element<3U>*,vector<vector<double> > >  mapped;
        MapPointsToElements3DSurfaceTetraMesh( pointsAndProperties, gref, mapped ); // ToDo: Bug: properties are not put on diferent vectors in map
        
        if (debug) {
             cout << "\nOutput of non-averaged property map."<< endl;
             for ( map<Element<3U>*,vector<vector<double> > >::iterator
                   map_it=mapped.begin(); map_it !=mapped.end(); ++map_it ) {
                  cout << "Outputing Values for this Element" << endl;
                  for ( size_t i{0U}; i<map_it->second.size(); ++i ) {
                       cout << "Ouputing Values for this Property"<<endl;
                       for ( size_t j{0U}; j<map_it->second[i].size(); ++j )                
                         cout << "Value: "<<map_it->second[i][j]<<endl;              
                    }
              }
          }
       
        // value averaging inside of the elements
        map<csmp::Element<3U>*,vector<double > > elementAveragedPropertyMap;
        AveragePropertiesOnElements(property_descriptors.size(),mapped,true,elementAveragedPropertyMap);
        if(debug){
          cout << "Outputing Averaged Property Map" <<endl;
          map<Element<3U>*,vector<double> >::iterator avg_it;
          for(avg_it = elementAveragedPropertyMap.begin();avg_it !=elementAveragedPropertyMap.end();++avg_it)
            for(size_t i = 0;i<avg_it->second.size();++i)
              cout << "Value: "<<avg_it->second[i]<<endl;
        }

        AssignElementProperties(elementAveragedPropertyMap,property_names,property_keys);
        
        vector<vector<Element<3U>*> > noPropertyElementVectors;
        InitializeEmptyElementVector( gref, property_keys, noPropertyElementVectors );
        AssignDefaultValuesToElements(noPropertyElementVectors, property_keys, defaultValues);
        
        cout << "\n\tStarting the interpolation of "<< property_keys.size() <<" properties onto the surface regions...\n";
        for( size_t i{0U}; i<property_keys.size(); ++i )
          {
            cout <<"Interpolating Property: "<<property_descriptors[i]<<endl;
            gref.InputPropertyValue(scalar_node, makeScalar(PLAIN,0.0));
            gref.InputPropertyValue(scalar_elmt, makeScalar(PLAIN,0.0));
            
            SpreadElementPropertyToNodeProperty( elementAveragedPropertyMap, property_keys[i], scalar_node_key );
            
            const size_t dof = degreesOfFreedom(gref,scalar_node_key);
            if ( dof >= MIN_DOF_TO_INVOKE_SOLVER )
              {
                cout << "\n\tUsing SAMG solver for Laplace interpolation of property values in surface domain..." << endl;
                PDE_Integrator<3U,Element> property_int_samg( solver );
                property_int_samg.Add(&stiffness_matrix);
                property_int_samg.Add(&source);
                property_int_samg.IntegrateOver(gref);
              }
            gref.InterpolateNodeToCellProperty(scalar_node,scalar_elmt);
            AssignElementPropertyFromElementList(noPropertyElementVectors[i],scalar_elmt_key,property_keys[i],empty_elmt_key);
            if ( dof >= 1 ) nearestNeighborFill( model, (*it).first.c_str(), property_descriptors[i].c_str(), ZERO_DATA_VALUE );
          }
      }

  }
  
} // end propertiesToRegions





/** 
  Finds Empty Elements based on the Flag of the property in that Element.
  If Value in Field then set PLAIN
  Otherwise initialise Values with PLAIN ????
  LM 2013
*/
static void FindEmptyElements(Region<3U>& gref, csmp::Index property_key, vector<Element<3U>*> emptyElements){
  if ( !emptyElements.empty() ) emptyElements.clear();
  vector<Element<3U>*>  noPropElements;
  emptyElements.clear();
  emptyElements.reserve(gref.Cells());
  for ( size_t j{0U}; j<gref.Cells(); ++j ){
    if( gref.E(j)->Status(property_key)==ANY) emptyElements.push_back(gref.E(j));
  }
}


} // end csmp
