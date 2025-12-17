#include "FluidPropertyLookupTable.h"
#include "Exception.h"

using namespace std;

namespace csmp {

  /**
  Fluid Property Lookup Table Container
  Created by LJM 03.2013

  Imports data from a BlackOil File Format 
  and creates a Lookup Table that can be accessed via indices.

  First Create a Lookup Table Object via:
  FluidPropertyLookupTable table(const char* file_name);

  then use: table.BOLookup[i][j] to access the Lookup Table values

  To get a property at a certain pressure use:
  IntFluidProp(double& p, size_t& prop_index)
  
  Indices:
  Pressure -> 0
  GOR -> 1
  Oil FVF -> 2
  Oil Viscosity -> 3
  Oil Density -> 4
  Oil Compressibility -> 5
  Gas FVF -> 6
  Gas Viscosity -> 7
  Water Viscosity -> 8
  Water Compressibility -> 9
  Z Factor -> 10
  Gas Density -> 11
  CGR -> 12
  Vapour CGR -> 13 

  @todo: Create a handler for when the pressure exceeds the tabulated data
         Make Multi Table Compatible
  */
FluidPropertyLookupTable::FluidPropertyLookupTable(const char* file_name)
  {
    LoadTableFromFile(file_name);
    FlushToScreen();
  }



void FluidPropertyLookupTable::FlushToScreen(){
    cout << "PRES\tGOR\tOFVF\tOVIS\tODEN\tOCOM\tGFVF\tGVIS\tWVIS\tWCOM\tZFAC\tGDEN\tCGR\tVCGR" << endl;
    for(size_t i  = 0;i < BOLookup[0].size();++i){
      for(size_t j = 0; j < BOLookup.size();++j){
        cout << BOLookup[j][i] << "\t";
      }
      cout << endl;
    }
    cout << "------------------------------End of Fluid Data---------------------------------------------------------------------" << endl;
    cout << "\n------------------------------Start of Extrapolation Matrix---------------------------------------------------------------------" << endl;
    cout << "PRES\tGOR\tOFVF\tOVIS\tODEN\tOCOM\tGFVF\tGVIS\tWVIS\tWCOM\tZFAC\tGDEN\tCGR\tVCGR" << endl;
    for(size_t i  = 0;i < 2;++i){
      for(size_t j = 0; j < ExtraPolMat.size();++j){
        cout << ExtraPolMat[j][i] << "\t";
      }
      cout << endl;
    }
    cout << "------------------------------End of Extrapolation Matrix---------------------------------------------------------------------" << endl;
  }



  
void FluidPropertyLookupTable::LoadTableFromFile(const char* fluid_property_file)
  {
    ifstream myfile(fluid_property_file);
    string line;
    string item;
  
    long line_index = 1;
    long value_index = 0;
  
  
    string name;
  
    if (myfile.is_open())
    {
      while ( myfile.good() )
      {
        //do This for each line, Line Format: Fault Name, X Coord, Y Coord, Z, Coord, Scalar 1, Scalar 2 etc.
        getline( myfile, line);
        //cout << "Line Number: "<< line_index << "Line: "<<line<<endl;
        
        stringstream ss(line);
        value_index = 0;
        //Iterate over the current lines values
        while(getline(ss,item,' ')){
          if(line_index == 55){ //Number of Vals, Tres, Psat reading
            if(value_index == 1){
              vals = atoi(item.c_str());
            }
            else if(value_index == 2){
              tres = string_to_double(item);
            }
            else if(value_index == 3){
              psat = string_to_double(item);
            }
            value_index++;
          }
          if( line_index >=56 && line_index <= (56+vals-1) ){ //Read all the values in the lookup table file
            if(value_index == 0) pres.push_back(string_to_double(item));
            else if(value_index == 1) gor.push_back(string_to_double(item));
            else if(value_index == 2) ofvf.push_back(string_to_double(item));
            else if(value_index == 3) ovis.push_back(string_to_double(item));
            else if(value_index == 4) oden.push_back(string_to_double(item));
            else if(value_index == 5) ocom.push_back(string_to_double(item));
            else if(value_index == 6) gfvf.push_back(string_to_double(item));
            else if(value_index == 7) gvis.push_back(string_to_double(item));
            else if(value_index == 8) zfac.push_back(string_to_double(item));
            else if(value_index == 9) gden.push_back(string_to_double(item));
            else if(value_index == 10) cgr.push_back(string_to_double(item));
            else if(value_index == 11) vcgr.push_back(string_to_double(item));
            value_index++;
          }
        }
        line_index++;
      }
      myfile.close();
    }
    //Build LookupTable
    BOLookup.push_back(pres);
    AddToLookup(gor);
    AddToLookup(ofvf);
    AddToLookup(ovis);
    AddToLookup(oden);
    AddToLookup(ocom);
    AddToLookup(gfvf);
    AddToLookup(gvis);
    AddToLookup(wvis);
    AddToLookup(wcom);
    AddToLookup(zfac);
    AddToLookup(gden);
    AddToLookup(cgr);
    AddToLookup(vcgr);

    //Set the number of values in pressure vector to size_t vals for later lookup
    vals = static_cast<long>(pres.size());
    ConvertFromFieldToSIUnits();
    CalculateExtrapolationMatrix();
  }



void FluidPropertyLookupTable::ConvertFromFieldToSIUnits()
  {
    for(size_t i = 0;i < BOLookup.size();++i){
      for(size_t j = 0; j < BOLookup[0].size();++j){
            if(i == 0) BOLookup[i][j] = BOLookup[i][j]*6894.757;
            else if(i == 1) /*GOR Do nothing for now*/;
            else if(i == 2) /*OFVF Do nothing for now*/;
            else if(i == 3) BOLookup[i][j] = BOLookup[i][j]*0.001; //cP to Pa*s
            else if(i == 4) BOLookup[i][j] = BOLookup[i][j]*16.0185;//lb/ft3 to kg/m3 is this true? or is it in API?
            else if(i == 5) BOLookup[i][j] = BOLookup[i][j]*145.03774; //1/psi to 1/Pa ? is this correct?
            else if(i == 6) /*GFVF Do nothing for now*/;
            else if(i == 7) BOLookup[i][j] = BOLookup[i][j]*0.001; //cP to Pa*s
            else if(i == 8) /*ZFAC Do nothing for now*/;
            else if(i == 9) /*GDEN Do nothing for now*/;
            else if(i == 10) /*CGR Do nothing for now*/;
            else if(i == 11) /*VCGR Do nothing for now*/;
            else
            throw csmp::Exception( ERROR, "FluidPropertyLookupTable::ConvertFromFieldToSIUnits:",
                                  "look-up table row index case not handled yet.");
      }
    }
  }



void FluidPropertyLookupTable::CalculateExtrapolationMatrix(){
    //Precomputes the Extrapolation Matrix ExtraPolMat that contains derivatives
    //to calculate extrapolated values. Precomputing these derivatives should increase 
    //speed dramatically.
    //Needs Testing
    ExtraPolMat.resize(BOLookup.size());
    for( size_t i = 0; i < BOLookup.size(); ++i ){
      ExtraPolMat[i].resize(2);
      for( int j=0; j <2; ++j ){
        if(j == 0){
          //Precompute Lower Extrapolation Derivative
          ExtraPolMat[i][j] = (BOLookup[i][1]-BOLookup[i][0])/(BOLookup[0][1]-BOLookup[0][0]);
        }
        else if(j == 1){
          //Precompute Upper Extrapolation Derivative
          ExtraPolMat[i][j] = (BOLookup[i][vals-1]-BOLookup[i][vals-2])/(BOLookup[0][vals-1]-BOLookup[0][vals-2]);
        }
      }
    }
  }



void FluidPropertyLookupTable::AddToLookup(std::vector<double>& prop){
    if(prop.empty()){
      prop.resize(pres.size());
      BOLookup.push_back(prop);
    }
    else{
      BOLookup.push_back(prop);
    }
  }



double FluidPropertyLookupTable::IntFluidProp(double p, size_t prop_index)
  {
      long close_p = FindClosestPressure(p);
      if(close_p > -1){
        return BOLookup[prop_index][close_p-1]+(BOLookup[prop_index][close_p]-BOLookup[prop_index][close_p-1])/(p-BOLookup[0][close_p-1]);
      }
      else if(close_p == -1){
        //Extrapolate to the pressure above the maximum value in the lookup table
        //Using the gradient of the last two values for the property in the table
        //Needs to be tested
        return BOLookup[prop_index][vals-1]+ExtraPolMat[prop_index][1]*(p-BOLookup[0][vals-1]);
      }
      else if(close_p == -2){
        //Extrapolate to the pressure below the minimum value in the lookup table
        //Using the gradient of the first two values for the property in the table
        //Needs to be tested
        return BOLookup[prop_index][0]-ExtraPolMat[prop_index][0]*(BOLookup[0][0]-p);
      }
    
    return numeric_limits<double>::quiet_NaN();
  }


/// returns line index of closest pressure for interpolation
long FluidPropertyLookupTable::FindClosestPressure(double& p){
    auto it = std::lower_bound(pres.begin(),pres.end(),p);
    long index = std::distance(pres.begin(), it);
    if(it == pres.end()) return -1;
    else if(it == pres.begin()) return -2;
    else return index;
  }
  
  

double FluidPropertyLookupTable::string_to_double( const std::string& s )
  {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
  }
  
  
}//csmp
