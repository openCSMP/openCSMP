//
//  SKUA_Interface.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/20/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#include "SKUA_Interface.h"
#include "Model.h"
#include "RegionInterface.h"
#include "Point.h"

using namespace std;

namespace csmp {


/**  VariableToPointCloud

   Outputs any element variable from CSMP to SKUA point cloud (ASCII table) format.
   Variables discretized on volume, surface or line elements will be output 
   as point data placed on the barycenter of respective elements.  
   
   @attention methods outputs coordinates in SKUA format where y=-z  and  z=-y
   
*/   
void SKUA_Interface::VariableToPointCloud( const Model<3U>& model,
                                           const set<string>& regions_of_interest,
                                           const char* filename, const char* var ) const
 {
    const csmp::Index  var_key(model.Database().StorageKey(var));
    if ( var_key.place != ELEMENT )
      throw csmp::Exception( CSMP_ERROR, "SKUA_Interface::VariableToPointCloud", "This method works only for 'Element' variables" );
   
    string var_name(var);
    replaceWhiteSpaceBy( var_name, '-' );
   
    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::VariableToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t";
    if ( var_key.type == SCALAR ) ofs << var_name;
    else if ( var_key.type == VECTOR ) ofs << var_name <<"[0]\t"<< var_name <<"[1]\t"<< var_name <<"[2]\n";
    else if ( var_key.type == TENSOR ) {
         int counter(0);
         for ( size_t i=0; i<3U; i++ )
           for ( size_t j=0; j<3U; j++ ) ofs << var_name <<"["<< counter++ <<"]\t";
      }
    else { // ARRAY variable
         for ( size_t i=0; i<var_key.index; i++ ) ofs << var_name <<"["<< i <<"]\t";
      }
    ofs <<"\n";
   
    // for all model regions
    for (  set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( vector<Element<3U>*>::const_iterator
               it=gref.ElementsBegin(); it!=gref.ElementsEnd(); ++it )
           {
              ofs << (*rt) <<"\t";
              Point<3U> xyz((*it)->BaryCenter());
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";
              if ( var_key.type == SCALAR ) ofs << (*it)->Read( var_key ) <<"\n";
              else if ( var_key.type == VECTOR ) {
                   VectorVariable<3U> vc;
                   (*it)->Read( var_key, vc );
                   ofs << vc[0] <<"\t"<< vc[1] <<"\t"<< vc[2] <<"\n";
                }
              else if ( var_key.type == TENSOR ) {
                   TensorVariable<3U> ts;
                   (*it)->Read( var_key, ts );
                   for ( size_t i=0; i<3U; i++ )
                     for ( size_t j=0; j<3U; j++ ) ofs << ts(i,j) <<"\t";
                   ofs <<"\n";
                }
              else { // ARRAY variable
                   ArrayVariable  ary;
                   (*it)->Read( var_key, ary );
                   for ( size_t i=0; i<ary.Size(); i++ ) ofs << ary[i] <<"\t";
                   ofs <<"\n";
                }
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::VariableToPointCloud: variable '";
    cout << var <<"' written successfully to file "<< filename << endl;
   
 } // end VariableToPointCloud







/// like previous method, but for multiple output variables discretized on the elements
void SKUA_Interface::VariablesToPointCloud( const Model<3U>& model,
                                            const set<string>& regions_of_interest,
                                            const char* filename, const set<string>& element_vars ) const
 {
    // processing variable names
    map<string,csmp::Index>  var_keys;
    for ( set<string>::const_iterator i=element_vars.begin(); i!=element_vars.end(); i++ )
      {
         const csmp::Index  var_key(model.Database().StorageKey((*i).c_str()));
         if ( var_key.place != ELEMENT )
            throw csmp::Exception( CSMP_ERROR, "SKUA_Interface::VariablesToPointCloud", "This method works only for 'Element' variables" );
         if ( var_key.type != SCALAR )
            throw csmp::Exception( CSMP_ERROR, "SKUA_Interface::VariablesToPointCloud", "This method works only for SCALAR variables" );
         string var_name((*i).c_str());
         replaceWhiteSpaceBy( var_name, '_' );
         var_keys.insert( make_pair(var_name,var_key) );
      }
   
    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::VariablesToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t";
    for ( map<string,csmp::Index>::const_iterator vit=var_keys.begin(); vit!=var_keys.end(); vit++ )
      ofs << (*vit).first <<"\t";
    ofs <<"\n";
   
    // for all model regions
    for ( set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( vector<Element<3U>*>::const_iterator
               it=gref.ElementsBegin(); it!=gref.ElementsEnd(); ++it )
           {
              ofs << (*rt) <<"\t";
              Point<3U> xyz((*it)->BaryCenter());
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";
              // writing the values of all the enlisted scalar element variables
              for ( map<string,csmp::Index>::const_iterator vit=var_keys.begin(); vit!=var_keys.end(); vit++ )
                ofs << (*it)->Read( (*vit).second ) <<"\t";
              ofs <<"\n";
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::VariablesToPointCloud: variables: ";
    for ( set<string>::const_iterator i=element_vars.begin(); i!=element_vars.end(); i++ ) cout <<"'"<< (*i) <<"'  ";
    cout <<"written successfully to file "<< filename << endl;
   
 } // end VariablesToPointCloud








/// creates point clouds around surface of interest using the 'thickness' attribute
// TODO: use UnitNormal() of surface elements to construct normals
void SKUA_Interface::SurfaceArrayVariableToPointCloud( const Model<3U>& model,
                                                       const set<string>& regions_of_interest,
                                                       const char* filename, const char* var ) const
 {
    const csmp::Index  var_key(model.Database().StorageKey(var));
    if ( var_key.place != ELEMENT )
      throw csmp::Exception( CSMP_ERROR, "SKUA_Interface::SurfaceArrayVariableToPointCloud", "This method works only for 'Element' variables" );
    assert( var_key.type == ARRAY );
    string var_name(var);
    replaceWhiteSpaceBy( var_name, '-' );

    const csmp::Index  fth_key(model.Database().StorageKey("thickness"));
    assert( fth_key.type == SCALAR );
    Point<3U>      nrml;
    ArrayVariable  ary(6);

    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::SurfaceArrayVariableToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t"<< var_name <<"\n";

    // for all model regions
    for (  set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( vector<Element<3U>*>::const_iterator
               it=gref.ElementsBegin(); it!=gref.ElementsEnd(); ++it )
           {
              // each array variable entry is output as a singe line 
              assert( (*it)->FE()->IsSurfaceElement() );
             
              // getting the normal along which the output points will be created
              if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
                nrml = normalAtFacetCenter( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(),
                                            (*it)->N(2)->Coordinate(), (*it)->N(3)->Coordinate() );
              else
                if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
                  nrml = normalOfTriangle( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(), (*it)->N(2)->Coordinate() );
              else
              throw Exception( CSMP_ERROR, "SKUA_Interface::SurfaceArrayVariableToPointCloud",
                              "attempt to compute normal on volume element rather than fault surface element");

              // reading the ARRAY variable and the thickness attribute
              (*it)->Read( var_key, ary );
              // writing out the first value representing center of surface
              Point<3U> xyz = (*it)->BaryCenter();
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << (*rt) <<"\t" << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";

              // doing all subsequent points, assuming that values are symmetrically distributed around surface
              // and using the unit normal
              const double64 dx = (*it)->Read( fth_key ) / static_cast<double64>(ary.Size()*2);
              for ( size_t i=1U; i<ary.Size(); i++ )
                {
                   Point<3U> out_pt =  (dx * i) * nrml;
                   Point<3U> ins_pt = (-dx * i) * nrml;
                   Point<3U> out = xyz + out_pt;
                   Point<3U> inp = xyz + ins_pt;
                   ofs << (*rt) <<"\t" << out[0] <<"\t"<< out[1] <<"\t"<< out[2] <<"\t"<< ary[i] <<"\n";
                   ofs << (*rt) <<"\t" << inp[0] <<"\t"<< inp[1] <<"\t"<< inp[2] <<"\t"<< ary[i] <<"\n";
                }
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::SurfaceArrayVariableToPointCloud: variable '";
    cout << var <<"' written successfully to file "<< filename << endl;

 } // end SurfaceArrayVariableToPointCloud





} // end csmp
