//
//  FluidPropertyLookupTable.h
//  CSMP_ReservoirSimulator
//
//  Created by Lukas Mosser on 3/25/13.
//  Provides a way to retrieve fluid properties
//  from a BlackOil Format Lookup Table
//
//  Copyright (c) 2013 Mosser . All rights reserved.
//

#ifndef CSMP_FLUID_PROPERTY_LOOKUP_TABLE_H
#define CSMP_FLUID_PROPERTY_LOOKUP_TABLE_H

#include "CSMP_definitions.h"

namespace csmp {

/** FluidPropertyLookupTable
   
    Loads fluid properties from a blackoil format lookup table
    and provides functions to get values for the black oil reservoir simulator.
    
    Interpolation between two data points is done linear.
    
    @author Lukas Mosser
    @date 25/3/2013
    
    @todo Add a Function to load the fluid properties into csmp
*/
class FluidPropertyLookupTable
  {
  public:
    FluidPropertyLookupTable( const char* file_name );

    void FlushToScreen();

    double IntFluidProp(double p,size_t prop_index);

    std::vector<std::vector<double > > BOLookup;

  private:
    void LoadTableFromFile( const char* fluid_property_file );
    void AddToLookup(std::vector<double>& prop);
    void CalculateExtrapolationMatrix();
    void ConvertFromFieldToSIUnits();
    
    long FindClosestPressure(double& p);

    std::vector<std::vector<double > > ExtraPolMat;

    double string_to_double(const std::string& s);

    double tres;
    double psat;
    long   vals;
    
    std::vector<double> pres;
    std::vector<double> gor;
    std::vector<double> ofvf;
    std::vector<double> ovis;
    std::vector<double> ocom;
    std::vector<double> oden;
    std::vector<double> gfvf;
    std::vector<double> gvis;
    std::vector<double> wvis;
    std::vector<double> wcom;
    std::vector<double> zfac;
    std::vector<double> gden;
    std::vector<double> cgr;
    std::vector<double> vcgr;
  };

} // end csmp

#endif // CSMP_FLUID_PROPERTY_LOOKUP_TABLE_H
