#include "MatlabInterface.h"
#include "Region.h"

using namespace std;


namespace csmp {

MatlabInterface::MatlabInterface() {}
MatlabInterface::~MatlabInterface() {};


void MatlabInterface::ExtractAndWrite1DVariableProfileAlongX( Model<2U>& mdl, double distance, const char* name, const char* variable, long step )
 {
    const Index        key(mdl.Database().StorageKey(variable));
    size_t             size(0);
    const Region<2>&   mref = mdl.Region("Model");

    for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
        if ( (*it)->y() == distance ) size++;
      }
    vector<pair<double,double> >                  vec1(size);
    vector<pair<double,pair<double,double> > >  vec2(size);
    if ( key.type == SCALAR ) {
            double  val;
            size_t    i(0);
            for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
                if ( (*it)->y() == distance ) {
                    val = (*it)->Read( key );
                    vec1[i].first  = (*it)->x();
                    vec1[i].second = val;
	            i++;
	          }
	      }
        sort( vec1.begin(), vec1.end() );
      }
    else {
       VectorVariable<2U>  vals;
       size_t              i(0);
       for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
            if ( (*it)->y() == distance ) {
                (*it)->Read( key, vals );
                vec2[i].first  = (*it)->x();
                vec2[i].second.first  = vals(0);
                vec2[i].second.second = vals(1);
                i++;
              }
          }
        sort( vec2.begin(), vec2.end() );
      }
    
    string  outfile(name), num1, num2;
    stringstream sstream1, sstream2;
    sstream1 << step;
    sstream1 >> num1;
    sstream2 << static_cast<long>(distance);
    sstream2 >> num2;
    outfile += "_1d_x_";
    outfile += num2;
    outfile += "_t_";
    outfile += num1;
    outfile += ".txt";

    ofstream ofs;
    ofs.open( outfile.c_str(), ios::out|ios::trunc );
    
    if ( key.type == SCALAR )
      for ( auto it2=vec1.begin(); it2!=vec1.end(); it2++ ) ofs << it2->first << "\t" << it2->second << endl;
    else
      for ( auto it2=vec2.begin(); it2!=vec2.end(); it2++ )
        ofs << it2->first << "\t" << it2->second.first <<  "\t" << it2->second.second << endl;
 }


void MatlabInterface::ExtractAndWrite1DVariableProfileAlongY( Model<2U>& mdl, double distance, const char* name, const char* variable, long step )
 {
    const Index       key(mdl.Database().StorageKey(variable));
    size_t            size(0);
    const Region<2>&  mref = mdl.Region("Model");
    
    for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
        if ( (*it)->x() == distance ) size++;
      }
        vector<pair<double,double> >                  vec1(size);
        vector<pair<double,pair<double,double> > >  vec2(size);
    if ( key.type == SCALAR ) {
            double  val;
            size_t    i(0);
            for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
                if ( (*it)->x() == distance ) {
                    val = (*it)->Read( key );
                    vec1[i].first  = (*it)->y();
                    vec1[i].second = val;
	            i++;
	          }
	      }
        sort( vec1.begin(), vec1.end() );
      }
    else {
        VectorVariable<2U>  vals;
        size_t              i(0);
        for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
            if ( (*it)->x() == distance ) {
                (*it)->Read( key, vals );
                vec2[i].first  = (*it)->y();
                vec2[i].second.first  = vals(0);
                vec2[i].second.second = vals(1);
                i++;
              }
          }
        sort( vec2.begin(), vec2.end() );
      }
    
    string  outfile(name), num1, num2;
    stringstream sstream1, sstream2;
    sstream1 << step;
    sstream1 >> num1;
    sstream2 << static_cast<long>(distance);
    sstream2 >> num2;
    outfile += "_1d_y_";
    outfile += num2;
    outfile += "_t_";
    outfile += num1;
    outfile += ".txt";

    ofstream ofs;
    ofs.open( outfile.c_str(), ios::out|ios::trunc );
    
    if ( key.type == SCALAR )
      for ( vector<pair<double,double> >::iterator it2=vec1.begin(); it2!=vec1.end(); it2++ ) ofs << it2->first << "\t" << it2->second << endl;
    else
      for ( vector<pair<double,pair<double,double> > >::iterator it2=vec2.begin(); it2!=vec2.end(); it2++ ) ofs << it2->first << "\t" << it2->second.first <<  "\t" << it2->second.second << endl;
    

    
 } 


void MatlabInterface::WriteSavedTimeStepsFile( long step, const char* name )
 {
    static bool first_pass(true);
    static ofstream ofs;
    if ( first_pass ) { ofs.open( name, ios::out|ios::trunc ); first_pass = false; }
    ofs << step << endl;
 }

void MatlabInterface::Write2DMatlabFile( Model<2U>& mdl, const char* file_name, const char* variable_name, long step, const char* region_name )
 {
    const Index       key(mdl.Database().StorageKey(variable_name));
    static bool       written(false);
    const Region<2>&   mref = mdl.Region(region_name);
    
    // generate name of variable file with/without sub region
    string outputfile(file_name), regionfile(region_name), number;
    stringstream sstream;
    if ( regionfile != "Model" ) {
        outputfile += "_";
        outputfile += regionfile;
    }
    sstream << step;
    sstream >> number;
    outputfile += number;
    outputfile += ".txt";

    // only write node coordinates and node IDs for each element once
    if ( !written ) {
        ofstream nodes, elements;
        nodes.open( "nodes.txt", ios::out|ios::trunc );
        elements.open( "elements.txt", ios::out|ios::trunc );

        // node coordinates
        for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
            nodes << (*it)->x() << "\t" << (*it)->y() << endl;
          }

        // get node-element connectivity
        for ( auto it2=mref.ElementsBegin(); it2!=mref.ElementsEnd(); it2++ ) {
            // regular triangles
            if ( (*it2)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE or
                 (*it2)->FE_Type() == LINEAR_TRIANGLE ) {
                    for ( auto i=0; i<(*it2)->Nodes(); i++ ) {
                        if ( i<((*it2)->Nodes()-1) )
                            elements << (*it2)->N(i)->Idx()+1 << "\t";
	                else
                            elements << (*it2)->N(i)->Idx()+1 << endl;
	              }
	          }
	        // quadrilaterals are split into 2 triangles  
                else if ( (*it2)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) {
	            // upper triangle
                    for ( auto i=1; i<(*it2)->Nodes(); i++ ) {
                        if ( i<((*it2)->Nodes()-1) )
                            elements << (*it2)->N(i)->Idx()+1 << "\t";
	                else
                            elements << (*it2)->N(i)->Idx()+1 << endl;
	              }
	            // lower triangle
                    for ( auto i=0; i<(*it2)->Nodes(); i++ ) {
                        if ( i<((*it2)->Nodes()-2) )
                            elements << (*it2)->N(i)->Idx()+1 << "\t";
                        if ( i==((*it2)->Nodes()-1) )
                            elements << (*it2)->N(i)->Idx()+1 << endl;
	              }
	          }
	        /*else {
	            cout << "\nwrite2DMatlabFile: This function works only for element types ";
	            cout << "\nISOPARAMETRIC_LINEAR_TRIANGLE, LINEAR_TRIANGLE and ISOPARAMETRIC_LINEAR_QUADRILATERAL  ";
	          }*/
            
          }	   
       nodes.close();
       elements.close();
       written = true;    	    
     }
    // save nodal values to file
    ofstream ofs;
    ofs.open( outputfile.c_str(), ios::out|ios::trunc );
    VectorVariable<2U> vec;
    TensorVariable<2U> ts;
    for ( auto it=mref.NodesBegin(); it!=mref.NodesEnd(); it++ ) {
        if ( key.type == SCALAR ) {
            ofs << (*it)->Read( key ) << endl;
          }
        else if ( key.type == VECTOR ) {
            (*it)->Read( key, vec );
            ofs << vec(0) << "\t" << vec(1) << endl;
          }
        else if ( key.type == TENSOR ) {
            (*it)->Read( key, ts );
            ofs << ts(0,0) << "\t" << ts(0,1) <<"\t" << ts(1,0) << "\t" << ts(1,1) << endl;
          }
      }
    ofs.close();    
    
    cout << "\nMatlabInterface::Write2DMatlabFile: File '" << outputfile << "' successfully written " << endl;
           
 } 
 
}

 
