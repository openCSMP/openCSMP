#include "MapleInterface.h"
#include "Region.h"
#include "Model.h"

using namespace std;

namespace csmp {

/**

writes Maple plot description which can be pasted directly into Maple  
worksheet.
*/
void writeVariableToMapleTextFile( const Model<1U>& sg, 
                                   const char* variable, uint32_t timestep, double time )
 {
    const Region<1>&  model_domain(sg.Region("Model"));
    char   num[30];
    snprintf( num, sizeof(num), "%u", timestep );
    string fname(variable);
    fname += "-maple-dataset";
    fname += num;
    fname += ".mpl";
    for ( string::iterator it=fname.begin(); it!=fname.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';

    ofstream  ofs( fname.c_str() );

    ofs <<"writeVariableToMapleTextFile: "<< fname <<", variable '"<< variable <<"'"<< endl;
    ofs << endl << endl;
    
    csmp::Index  prop_key = sg.Database().StorageKey(variable);
    assert( prop_key.type == SCALAR );
    // removing potential hyphens from the dataset name
    string dataset(variable);
    dataset += num;
    for ( string::iterator it=dataset.begin(); it!=dataset.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';
    
    // writing a weighted list for maple listname and assignment
    size_t  counter(1U);
    ofs << dataset <<" := [ ";

    if ( prop_key.place == NODE ) {
        const auto  nit_last(--model_domain.NodesEnd());
	      for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
	           ofs <<"["<< (*nit)->x();
	           ofs <<","<< (*nit)->Read( prop_key );
	           if ( nit != nit_last ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    else if ( prop_key.place == ELEMENT ) {
        const auto  eit_last(--model_domain.CellsEnd());
	      for ( auto eit=model_domain.CellsBegin(); eit!=model_domain.CellsEnd(); eit++ ) {
	           Point<1U>  x((*eit)->BaryCenter());
	           ofs <<"["<< x[0];
	           ofs <<","<< (*eit)->Read( prop_key );
	           if ( eit != eit_last ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    // writing the plot function (looks like this for a single variable)
    Point<1U>  xyz_min, xyz_max;
    sg.MinMaxCoordinates( xyz_min, xyz_max );
    double var_min, var_max; 
    sg.MinMaxOf( variable, var_min, var_max );
    ofs << endl;
    ofs <<"plot( "<< dataset <<", a="<< xyz_min[0] <<".."<< xyz_max[0];
    ofs <<", y="<< var_min <<".."<< var_max <<","<< endl;
  	ofs <<"labels=[\"distance (m)\",\""<< variable <<"\"],"<< endl;
  	ofs <<"font=[HELVETICA,18],"<< endl;
  	ofs <<"labelfont=[HELVETICA,20],"<< endl;
  	ofs <<"titlefont=[HELVETICA,20],"<< endl;
  	ofs <<"axes=BOXED,"<< endl;
  	ofs <<"style=[POINT],"<< endl;
  	ofs <<"symbol=[BOX],"<< endl;
  	ofs <<"symbolsize=12,"<< endl;
  	ofs <<"linestyle=[SOLID],"<< endl; 
  	ofs <<"legend=[\"'"<< variable <<"' profile, t="<< time <<" secs.\"],"<< endl;
  	ofs <<"color=[black],"<< endl;
  //	ofs <<"xtickmarks=5,"<< endl;
  	ofs <<"resolution=2000,"<< endl;
  	ofs <<"thickness=[2] );"<< endl;    
      
    ofs.close();
    cout <<"\nwriteVariableToMapleTextFile: file '"<< fname <<"' written successfully." << endl;
    
 } // end writeVariableToMapleTextFile (model)








void writeVariablesToMapleTextFile( const Model<1U>& sg, 
                                    const char* variable1, const char* variable2, 
                                    uint32_t timestep, double time )
 {
    const Region<1>&  model_domain(sg.Region("Model"));
    char   num[30];
    snprintf( num, sizeof(num), "%u", timestep );
    string fname(variable1);
    fname += "-";
    fname += variable2;
    fname += "-maple-dataset";
    fname += num;
    fname += ".mpl";
    for ( string::iterator it=fname.begin(); it!=fname.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';

    ofstream  ofs( fname.c_str() );

    ofs <<"writeVariableToMapleTextFile: "<< fname <<", variables '";
    ofs << variable1 <<"', '"<< variable2 <<"'"<< endl;
    ofs << endl << endl;
    
    // 1. writing the first dataset
    // ----------------------------
    csmp::Index  prop_key(sg.Database().StorageKey(variable1));
    assert( prop_key.type == SCALAR );
    // removing potential hyphens from the dataset name
    string dataset1(variable1);
    dataset1 += num;
    for ( string::iterator it=dataset1.begin(); it!=dataset1.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';
    
    // writing a weighted list for maple listname and assignment
    size_t  counter(1U);
    ofs << dataset1 <<" := [ ";
    
    // NB: nodes will not necessarily be in the order of their coordinates
    //     so a line graph will not work
    if ( prop_key.place == NODE ) {
	      for ( auto i{0}; i<model_domain.Nodes(); i++ ) {
	           ofs <<"["<< model_domain.N(i)->x();
	           ofs <<","<< model_domain.N(i)->Read( prop_key );
	           if ( i < model_domain.Nodes()-1U ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    else if ( prop_key.place == ELEMENT ) {
	      for ( size_t i{0U}; i<model_domain.Cells(); i++ ) {
	           Point<1U> x(model_domain.E(i)->BaryCenter());
	           ofs <<"["<< x[0];
	           ofs <<","<< model_domain.E(i)->Read( prop_key );
	           if ( i < model_domain.Cells()-1U ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    // 2. writing the second dataset
    // -----------------------------
    csmp::Index prop_key2(sg.Database().StorageKey(variable2));
    assert( prop_key2.type == SCALAR );
    // removing potential hyphens from the dataset name
    string dataset2(variable2);
    dataset2 += num;
    for ( string::iterator it=dataset2.begin(); it!=dataset2.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';
    
    // writing a weighted list for maple listname and assignment
    counter = 1U;
    ofs << dataset2 <<" := [ ";
    
    if ( prop_key2.place == NODE ) {
	      for ( size_t i{0U}; i<model_domain.Nodes(); i++ ) {
	           ofs <<"["<< model_domain.N(i)->x();
	           ofs <<","<< model_domain.N(i)->Read( prop_key2 );
	           if ( i < model_domain.Nodes()-1U ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    else if ( prop_key2.place == ELEMENT ) {
	      for ( size_t i{0U}; i<model_domain.Cells(); i++ ) {
	           Point<1U>  x(model_domain.E(i)->BaryCenter());
	           ofs <<"["<< x[0];
	           ofs <<","<< model_domain.E(i)->Read( prop_key2 );
	           if ( i < model_domain.Cells()-1U ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }

    // writing the plot function (looks like this for a single variable)
    Point<1U>  xyz_min, xyz_max;
    sg.MinMaxCoordinates( xyz_min, xyz_max );
    double var_min, var_max; 
    sg.MinMaxOf( variable1, var_min, var_max );
    ofs << endl;
    ofs <<"plot( ["<< dataset1 <<","<< dataset2 <<"],"; 
    ofs <<" a="<< xyz_min[0] <<".."<< xyz_max[0] <<", y="<< var_min <<".."<< var_max <<","<< endl;
  	ofs <<"labels=[\"distance (m)\",\"var1, var2\"],"<< endl;
  	ofs <<"font=[HELVETICA,18],"<< endl;
  	ofs <<"labelfont=[HELVETICA,20],"<< endl;
  	ofs <<"titlefont=[HELVETICA,20],"<< endl;
  	ofs <<"axes=BOXED,"<< endl;
  	ofs <<"style=[LINE,LINE],"<< endl;
  	ofs <<"symbol=[BOX, BOX],"<< endl;
  	ofs <<"symbolsize=12,"<< endl;
  	ofs <<"linestyle=[SOLID,SOLID],"<< endl; 
  	ofs <<"legend=[\"'"<< variable1 <<"'\",\"'"<< variable2 <<"' profile, t="<< time <<" secs.\"],"<< endl;
  	ofs <<"color=[red,black],"<< endl;
  //	ofs <<"xtickmarks=5,"<< endl;
  	ofs <<"resolution=2000,"<< endl;
  	ofs <<"thickness=[2,2] );"<< endl;    
      
    ofs.close();
    cout <<"\nwriteVariableToMapleTextFile: file '"<< fname <<"' written successfully." << endl;
    
 } // end writeVariablesToMapleTextFile (model)





void writeVariableToMapleTextFile( const Model<1U>& sg, const char* group, 
                                   const char* variable, uint32_t timestep, double time )
 {
    const Region<1>& gref = sg.Region(group);
 
    char   num[30];  
    snprintf( num, sizeof(num), "%u", timestep );
    string fname(group); fname+="-"; fname+=variable; fname+="-maple-dataset"; fname+=num; fname+=".mpl";

    for ( string::iterator it=fname.begin(); it!=fname.end(); it++ )
      if ( *it == '-' || *it == ' ' || *it == '\t' || *it == '\n' ) *it = '_';

    ofstream  ofs( fname.c_str() );

    ofs <<"writeVariableToMapleTextFile: "<< fname <<", variable '"<< variable <<"'"<< endl;
    ofs << endl << endl;
    
    csmp::Index  prop_key = sg.Database().StorageKey(variable);
    assert( prop_key.type == SCALAR );
    // removing potential hyphens from the dataset name
    string dataset(group); dataset+="-"; dataset+=variable; dataset+=num;

    // writing a weighted list for maple listname and assignment
    size_t  counter(1U);
    ofs << dataset <<" := [ ";
    
    if ( prop_key.place == NODE ) {
          auto end_it=gref.NodesEnd(); end_it--;
          for ( auto nit=gref.NodesBegin(); nit!=gref.NodesEnd(); nit++ ) {
	           ofs <<"["<< (*nit)->x();
	           ofs <<","<< (*nit)->Read( prop_key );
	           if ( nit != end_it ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    else if ( prop_key.place == ELEMENT ) {
        vector<double>  x;
        auto end_it=gref.CellsEnd();
        end_it--;
	      for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ ) {
	           x = (*eit)->BaryCenter().Coordinates();
	           ofs <<"["<< x[0];
	           ofs <<","<< (*eit)->Read( prop_key );
	           if ( eit != end_it ) ofs <<"],";
	           else ofs <<"]";
	           if ( counter++ == 5U ) { ofs << endl; counter=1; }
	        }
	      ofs <<" ]:"<< endl;
      }
      
    // writing the plot function (looks like this for a single variable)
    double var_min, var_max, xmin, xmax; 
    auto nit=gref.NodesBegin();
    xmin = xmax = (*nit)->x();
    while (  nit!=gref.NodesEnd() ) {
         xmin = std::min ( xmin, (*nit)->x() ); 
         xmax = std::max ( xmax, (*nit)->x() ); 
         nit++;
      }
    gref.MinMaxOf( variable, var_min, var_max );
    ofs << endl;
    ofs <<"plot( "<< dataset <<", a="<< xmin <<".."<< xmax <<", y="<< var_min <<".."<< var_max <<","<< endl;
	  ofs <<"labels=[\"distance (m)\",\""<< variable <<"\"],"<< endl;
	  ofs <<"font=[HELVETICA,18],"<< endl;
	  ofs <<"labelfont=[HELVETICA,20],"<< endl;
	  ofs <<"titlefont=[HELVETICA,20],"<< endl;
	  ofs <<"axes=BOXED,"<< endl;
	  ofs <<"style=[POINT],"<< endl;
	  ofs <<"symbol=[BOX],"<< endl;
	  ofs <<"symbolsize=12,"<< endl;
	  ofs <<"linestyle=[SOLID],"<< endl; 
	  ofs <<"legend=[\"'"<< variable <<"' profile, t="<< time <<" secs.\"],"<< endl;
	  ofs <<"color=[black],"<< endl;
//	ofs <<"xtickmarks=5,"<< endl;
	  ofs <<"resolution=2000,"<< endl;
	  ofs <<"thickness=[2] );"<< endl;    
      
    ofs.close();
    cout <<"\nwriteVariableToMapleTextFile: file '"<< fname <<"' written successfully." << endl;
    
 } // end writeVariableToMapleTextFile (region)



} // end csmp
