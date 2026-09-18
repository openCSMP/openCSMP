// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  LayeredCompositeProcessor1.cpp
//  CSMP_FECFVM_Simulator
//
//  Created by Stephan Matthai on 26/6/19.
//

#include "LayeredCompositeProcessor1.h"
#include "TextFileIO.h"
#include "ErrorHandler.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

/**
    TODO: cleanup - remove variables that need to come from input file
    use rocktype properly.
*/
LayeredCompositeProcessor1::LayeredCompositeProcessor1( const char* RRT_data )
 :  // material properties
     k_low_(3.4759e-14), k_high_(3.6075e-13),  // layer permeabilities
     LY_low_(0.5), LY_high_(0.5),        // cumulative layer thickness in the vertical direction (Y)
     Swi_low_(0.18), Swi_high_(0.159),     // irreducible saturations of the 2 different layers
     m_low_(0.5), m_high_(0.6),            // van Genuchten exponents for the 2 different layers
     pd_low_(20684280.), pd_high_(6894760.),        // capillary (drainage) entry pressure of low and high k layers - high = far too high!
     bcp_low_(1.), bcp_high_(1.5)         // Brooks-Corey 64' exponents for low and high k layers
 {
    ReadRockTypeData( RRT_data );
 }
 
 
/**

Reading the datafile:

@code
ARP_FSst_Slt_new_names.txt - Input data table to define Heterogeneity Aware Relative Permeability as a function of sw and Nc = k grad P / sigma, sigma = interfacial tension H2O-CO2
# comment line
1 # total number of rocktypes in file
1 # rocktype identifier 1: FSst-Slt (Fine Sandstone - Silt) Planar Bedding
30 6 # table size(rows vs. columns)
# F-flow    Sw_VL (cell average)  Sw_CL(cell average)  Sw_low_CL  Sw_high_CL    Sfactor
0.990000000000000  0.958327612404074  0.978851642918739  1  0.957703285837478  137051613.363363
0.956206896551724  0.935003844772692  0.964807268982292  1  0.929614537964584  28739845.1229135
@endcode

*/
void LayeredCompositeProcessor1::ReadRockTypeData( const char* datafile )
 {
    if ( !doesFileExist( datafile ) )
      throw csmp::Exception( ERROR, "LayeredCompositeProcessor1::ReadRockTypeData",
                             datafile, "does not exist; nothing was done." );
    ifstream ifs(datafile);
   
    // reading and echoing file header to console
    std::string header;
    readFileHeader( ifs, header );
    cout <<"\nLayeredCompositeProcessor1::ReadRockTypeData: reading '"<< datafile <<"' header:\n";
    cout << header << endl;
   
    // reading how many rock type records are contained in file
    size_t line_length = INFO_STRING;
    char   text_line[INFO_STRING];
    int    rocktypes(0);
    while ( !ifs.eof() )
      if ( !readLineTellIfBlank( ifs, text_line, line_length ) )
        if ( !isCommentLine(text_line) ) {
             // parse number of records
             rocktypes = atoi(text_line);
             assert ( rocktypes > 0 );
             break;
          }

    // reading the rocktype records
    int rocktype(0), rocktype_counter(0);
    RRT_.reserve(rocktypes);
    vector<int> RRT_identifiers;
    RRT_identifiers.reserve(rocktypes);
    double value;
   
    while ( rocktype_counter < rocktypes and !ifs.eof() )
      {
         // read and parse rocktype identifier
         if ( !readLineTellIfBlank( ifs, text_line, line_length ) ) rocktype = atoi(text_line);
         else cerr <<"\nLayeredCompositeProcessor1::ReadRockTypeData: failed to read rocktype.\n";
         assert( rocktype <= rocktypes );
         RRT_identifiers.push_back( rocktype );
        
         // reading the table dimensions and setting up the storage vector for the table
         int rows(0), columns(0);
         if ( !readLineTellIfBlank( ifs, text_line, line_length ) ) {
               stringstream table_dimension(text_line);
               table_dimension >> rows >> columns;
           }
         else cerr <<"\nLayeredCompositeProcessor1::ReadRockTypeData: failed to read rocktype table dimensions.\n";
         assert( rows > 0 );
         assert( columns >= 6 );

         const int extra_columns(3);
         vector<vector<double> > table(rows,vector<double>(columns + extra_columns));
      
        // reading the table, swallowing the heading
        advancePastCommentLine( ifs );
         for ( int i=0; i<rows; ++i ) {
             if ( !readLineTellIfBlank( ifs, text_line, line_length ) ) {
                  stringstream row(text_line);
                  for ( int j=0; j<columns; ++j ) {
                       row >> value;
                       assert( value > -1.0e30 && value < 1.0e30 );
                       table[i][j] = value;
                    }
               }
             else cerr <<"\nLayeredCompositeProcessor1::ReadRockTypeData: failed to read table row: "<< i <<"\n";
           }
         // storing the recorded table
         RRT_.emplace_back( table );
        
         rocktype_counter++;
      }
   
    cout <<"\nLayeredCompositeProcessor1::ReadRockTypeData: "<< rocktypes <<" rocktype tables read from file '";
    cout << datafile <<"' successfully.\n";

 } // end ReadRockTypeData
  
  
  
  
// KEY METHODS FOR THE USER


/**
    Costly prototype which computes values for each row in the rocktype specific rocktype table
 
    The following table entries are expected: f(sw), Sw_VL, Sw_CL, Sw_CL_low_k, Sw_CL_high_k, Sfactor
 
*/
double LayeredCompositeProcessor1::krw( double sw, double Nc, int rocktype )
{
   assert( rocktype >= 0 );
   assert( rocktype < RRT_.size() );
  
   // computing Nc-specific average saturation values for each rock in the rocktype table
   const size_t table_columns(6U); // including fractional flow
   const size_t table_rows(RRT_[rocktype].size());
   for ( int i=0U; i<table_rows; ++i )
     {
         const double sw_at_Nc        = SwAtNc( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_low_k  = SwAtNc_Low_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_low_k], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_high_k = SwAtNc_High_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_high_k], RRT_[rocktype][i][Sfactor] );
         const double krw_comp        = KrwComposite( sw_at_Nc_low_k, sw_at_Nc_high_k );
         //storing the result values in the 3 last columns of the table
         RRT_[rocktype][i][table_columns]    = sw_at_Nc;
         RRT_[rocktype][i][table_columns+1U] = krw_comp;
     }
  
    // interpolating between the rows of the table
    // dealing with the extrema - if sw is out-of-range, the nearest value in the table will be chosen
    // NOTE: this is a table for drainage: sw decreases with increasing row number
    if      ( sw >= RRT_[rocktype][0U][table_columns] ) return RRT_[rocktype][0U][table_columns+1U];
    else if ( sw <= RRT_[rocktype][table_rows-1U][table_columns] ) return RRT_[rocktype][table_rows-1U][table_columns+1U];
    // searching through the rows in reverse order until an sw_at_Nc > sw is found
    else {
        int lower_row(0U);
        for ( int i=static_cast<int>(table_rows-1U); i>0; --i )
          if ( sw < RRT_[rocktype][i][table_columns] ) { lower_row = i; break; }
        // interpolating the kri values
        const double delta_sw = RRT_[rocktype][lower_row][table_columns] - RRT_[rocktype][lower_row+1][table_columns];
        const double ds = (sw - RRT_[rocktype][lower_row+1][table_columns]) / delta_sw;
        // weighting the kri's to find the interpolated value
        return ds * RRT_[rocktype][lower_row+1U][table_columns+1U] + (1. - ds) * RRT_[rocktype][lower_row][table_columns+1U];
     }
  
   // in case an error occurred
   return std::numeric_limits<double>::signaling_NaN();
  
} // end krw



double LayeredCompositeProcessor1::krn( double sw, double Nc, int rocktype )
 {
   assert( rocktype >= 0 );
   assert( rocktype < RRT_.size() );
  
   // computing Nc-specific average saturation values for each rock in the rocktype table
   const size_t table_columns(6U); // including fractional flow
   const size_t table_rows(RRT_[rocktype].size());
   for ( int i=0U; i<table_rows; ++i )
     {
         const double sw_at_Nc        = SwAtNc( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_low_k  = SwAtNc_Low_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_low_k], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_high_k = SwAtNc_High_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_high_k], RRT_[rocktype][i][Sfactor] );
         const double krn_comp        = KrnComposite( sw_at_Nc_low_k, sw_at_Nc_high_k );
         //storing the result values in the 3 last columns of the table
         RRT_[rocktype][i][table_columns]    = sw_at_Nc;
         RRT_[rocktype][i][table_columns+2U] = krn_comp;
     }
  
    // interpolating between the rows of the table
    // dealing with the extrema - if sw is out-of-range, the nearest value in the table will be chosen
    // NOTE: this is a table for drainage: sw decreases with increasing row number
    if      ( sw >= RRT_[rocktype][0U][table_columns] ) return RRT_[rocktype][0U][table_columns+2U];
    else if ( sw <= RRT_[rocktype][table_rows-1U][table_columns] ) return RRT_[rocktype][table_rows-1U][table_columns+2U];
    // searching through the rows in reverse order until an sw_at_Nc > sw is found
    else {
        int lower_row(0U);
        for ( int i=static_cast<int>(table_rows-1U); i>0; --i )
          if ( sw < RRT_[rocktype][i][table_columns] ) { lower_row = i; break; }
        // interpolating the kri values
        const double delta_sw = RRT_[rocktype][lower_row][table_columns] - RRT_[rocktype][lower_row+1][table_columns];
        const double ds = (sw - RRT_[rocktype][lower_row+1][table_columns]) / delta_sw;
        // weighting the kri's to find the interpolated value
        return ds * RRT_[rocktype][lower_row+1U][table_columns+2U] + (1. - ds) * RRT_[rocktype][lower_row][table_columns+2U];
     }
  
   // in case an error occurred
   return std::numeric_limits<double>::signaling_NaN();
 }
 
 
 
 
 /**
     More efficient version that computes and returns both relperms in one go.
  
     @return pair of krw,krn
 */
 pair<double,double>  LayeredCompositeProcessor1::RelativePermeability( double sw, double Nc, int rocktype )
{
   assert( rocktype >= 0 );
   assert( rocktype < RRT_.size() );
  
   // computing Nc-specific average saturation values for each rock in the rocktype table
   const size_t table_columns(6U); // including fractional flow
   const size_t table_rows(RRT_[rocktype].size());
   for ( int i=0U; i<table_rows; ++i )
     {
         const double sw_at_Nc        = SwAtNc( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_low_k  = SwAtNc_Low_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_low_k], RRT_[rocktype][i][Sfactor] );
         const double sw_at_Nc_high_k = SwAtNc_High_k_Layer( Nc, RRT_[rocktype][i][Sw_VL], RRT_[rocktype][i][Sw_CL_high_k], RRT_[rocktype][i][Sfactor] );
         const double krw_comp        = KrwComposite( sw_at_Nc_low_k, sw_at_Nc_high_k );
         const double krn_comp        = KrnComposite( sw_at_Nc_low_k, sw_at_Nc_high_k );
         //storing the result values in the 3 last columns of the table
         RRT_[rocktype][i][table_columns]    = sw_at_Nc;
         RRT_[rocktype][i][table_columns+1U] = krw_comp;
         RRT_[rocktype][i][table_columns+2U] = krn_comp;
     }
  
    // interpolating between the rows of the table
    // dealing with the extrema - if sw is out-of-range, the nearest value in the table will be chosen
    // NOTE: this is a table for drainage: sw decreases with increasing row number
    if      ( sw >= RRT_[rocktype][0U][table_columns] )
      return make_pair( RRT_[rocktype][0U][table_columns+1U], RRT_[rocktype][0U][table_columns+2U] );
    else if ( sw <= RRT_[rocktype][table_rows-1U][table_columns] )
      return make_pair( RRT_[rocktype][table_rows-1U][table_columns+1U], RRT_[rocktype][table_rows-1U][table_columns+2U] );
    // searching through the rows in reverse order until an sw_at_Nc > sw is found
    else {
        int lower_row(0U);
        for ( int i=static_cast<int>(table_rows-1U); i>0; --i )
          if ( sw < RRT_[rocktype][i][table_columns] ) { lower_row = i; break; }
        // interpolating the kri values
        const double delta_sw = RRT_[rocktype][lower_row][table_columns] - RRT_[rocktype][lower_row+1][table_columns];
        const double ds = (sw - RRT_[rocktype][lower_row+1][table_columns]) / delta_sw;
        // weighting the kri's to find the interpolated value
        return make_pair( ds * RRT_[rocktype][lower_row+1U][table_columns+1U] + (1. - ds) * RRT_[rocktype][lower_row][table_columns+1U],
                          ds * RRT_[rocktype][lower_row+1U][table_columns+2U] + (1. - ds) * RRT_[rocktype][lower_row][table_columns+2U] );
     }
  
   // in case an error occurred
   return make_pair( std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN() );
  
} // end RelativePermeability



/**
    writes textfile with sw, krw(sw,Nc), krn(sw,Nc), and pc(sw) values computed for (composite) rocktype in 0.05 saturation increments
*/
void LayeredCompositeProcessor1::WriteRelativePermeabilityTable( const char* filename, int rocktype, double Nc )
 {
    ofstream  ofs( string(filename) + ".txt" );
    assert( rocktype < RRT_.size() );
    assert( Nc > 0. );
   
    ofs <<"sw(Nc="<< Nc <<")\t krw(sw,Nc,rocktype="<< rocktype <<")\t krnw(sw,Nc)\n";
    for ( double sw(0.); sw<=1.0; sw+=0.05 ) {
         ofs << sw << "\t"<< krw( sw, Nc, rocktype );
         ofs <<"\t"<< krn( sw, Nc, rocktype );
//         ofs <<"\t"<< pc( sw, rocktype );
         ofs << endl;
      }
   
 } // end WriteRelativePermeabilityTable
  





void LayeredCompositeProcessor1::Out() const
 {
     cout <<"\nLayeredCompositeProcessor1::Out: ";
     cout <<"\n\tstored number of rocktypes: "<< RRT_.size();
     cout <<"\n\trock-type data tables:";
  
     for ( size_t rrt=0U; rrt<RRT_.size(); rrt++ ) {
         cout <<"\n\t\trocktype "<< rrt <<", data table: "<< RRT_[rrt].size() <<" x "<< RRT_[rrt][0U].size();
         cout <<"\n\t\t"<<"f(sw) Sw_VL  Sw_CL SW_CL_low_k_layer  Sw_CL_high_k_layer  Sfactor  Sw(Nc)\n\t\t";
         const size_t rows    = RRT_[rrt].size();
         const size_t columns = RRT_[rrt][0U].size();
         for ( size_t i{0U}; i<rows; ++i ) {
              for ( size_t j{0U}; j<columns; ++j )
                cout << RRT_[rrt][i][j] <<"  ";
              cout <<"\n\t\t";
           }
       }
     cout << endl;
     
 } // end Out


 
/**
    Best approach still subject of discussion. Capillary limit!
 
    Horizontal case.
    - entry pressure is that of the high-k layer
    - capillary pressure curve is stepped
      - first the high-k layer gets saturated
      - now the low key layer gets saturated once its entry pressure is overcome
*/
double LayeredCompositeProcessor1::pc( double sw, int rocktype ) const
 {
    double PV_low  = LY_low_ * phi_low_;
    double PV_high = LY_high_ * phi_high_;
    double sw_threshold = PV_high / (PV_high + PV_low); // saturation at which low-k-layer will start saturating
   
    // Brooks-Corey relations for the layers, pc = pd * sw_eff^(-1/lambda)
    if ( sw < sw_threshold ){
        return pd_high_ + pow( (sw - Swi_high_)/( 1. - Swi_high_), -1./bcp_high_ );
      }
   
    return pd_low_ + pow( ((sw-sw_threshold) - Swi_low_)/( 1. - Swi_low_), -1./bcp_low_ );
 }
 
 

// permeability averages

/// vertical thickness weighted (harmonic) mean of the vertical layer permeabilities
double LayeredCompositeProcessor1::k_AverageY() const
 {
    return (LY_low_ + LY_high_) / ( LY_low_/k_low_  + LY_high_/k_high_ );
 }
 
 
 
/// horizontal thickness-weighted average of the horizontal layer permeabilities
double LayeredCompositeProcessor1::k_AverageX() const
 {
    return k_low_ * LY_low_ + k_high_ * LY_high_;
 }
 
 
 
/// water saturation in cell at given capillary number; calculated from entries in table
double LayeredCompositeProcessor1::SwAtNc(  double Nc, double Sw_VL, double Sw_CL, double Sfactor ) const
 {
    return (Nc * Sfactor * Sw_VL + Sw_CL) / (Nc * Sfactor + 1.);
 }
 
 
 
/// corresponding sw in low-k laminations
double LayeredCompositeProcessor1::SwAtNc_Low_k_Layer( double Nc, double Sw_VL, double Sw_low_CL, double Sfactor ) const
 {
    return (Nc * Sfactor * Sw_VL + Sw_low_CL) / (Nc * Sfactor + 1.);
 }
 
 
 
/// corresponding sw in high-k laminations
double LayeredCompositeProcessor1::SwAtNc_High_k_Layer( double Nc, double Sw_VL, double Sw_high_CL, double Sfactor ) const
 {
    return (Nc * Sfactor * Sw_VL + Sw_high_CL) / (Nc * Sfactor + 1.);
 }
 
 
 
/// effective sw in low-k laminations
double LayeredCompositeProcessor1::SwStarLow( double Sw_low_at_Nc ) const
 {
    return (Sw_low_at_Nc - Swi_low_) / (1. - Swi_low_);
 }
 
 
 
/// effective sw in low-k laminations
double LayeredCompositeProcessor1::SwStarHigh( double Sw_high_at_Nc ) const
 {
    return (Sw_high_at_Nc - Swi_high_) / (1. - Swi_high_);
 }
 
 

/**
    krw relperm of low-k layer at given saturation and capillary number
    (uses van Genuchten model - with m parameter)
*/
double LayeredCompositeProcessor1::KrwLow( double SwStar_low ) const
 {
    double t1 = sqrt(SwStar_low);
    double t4 = pow(SwStar_low, 1. / m_low_);
    double t6 = pow(1. - t4, m_low_);
    double t8 = pow(1. - t6, m_low_);
    return(t8 * t1);
 }
 
 


/**
    krw relperm of high-k layer at given saturation and capillary number
    (uses van Genuchten model - with m parameter)
*/
double LayeredCompositeProcessor1::KrwHigh( double SwStar_high ) const
 {
    double t1 = sqrt(SwStar_high);
    double t4 = pow(SwStar_high, 0.1e1 / m_high_);
    double t6 = pow(0.1e1 - t4, m_high_);
    double t8 = pow(0.1e1 - t6, m_high_ );
    return(t8 * t1);
 }
 
 


/**
    krn relperm of low-k layer at given saturation and capillary number
    (uses Brooks-Corey model - with bcp parameter)
*/
double LayeredCompositeProcessor1::KrnLow( double SwStar_low ) const
 {
    double t2 = pow(1. - SwStar_low, bcp_low_);
    double t3 = SwStar_low * SwStar_low;
    return((0.1e1 - t3) * t2);
 }
 
 
 
/**
    krn relperm of high-k layer at given saturation and capillary number
    (uses Brooks-Corey model - with bcp parameter)
*/
double LayeredCompositeProcessor1::KrnHigh( double SwStar_high ) const
 {
    double t2 = pow(1. - SwStar_high, bcp_high_);
    double t3 = SwStar_high * SwStar_high;
    return((0.1e1 - t3) * t2);
 }
 
 
 
/**
    relative permeability at the given water saturation and capillary number
    (uses van Genuchten model - with m parameter)
*/
double LayeredCompositeProcessor1::KrwComposite( double swAtNc_Low_k_Layer, double swAtNc_High_k_Layer ) const
 {
    double t2 = SwStarLow(swAtNc_Low_k_Layer);
    double t3 = KrwLow(t2);
    double t6 = SwStarHigh(swAtNc_High_k_Layer);
    double t7 = KrwHigh(t6);
    return(1. / k_AverageX() / (LY_low_ + LY_high_) * (LY_high_ * k_high_ * t7 + LY_low_ * k_low_ * t3));

 }
 
 
 
double LayeredCompositeProcessor1::KrnComposite( double swAtNc_Low_k_Layer, double swAtNc_High_k_Layer ) const
 {
    double t2 = SwStarLow(swAtNc_Low_k_Layer);
    double t3 = KrnLow(t2);
    double t6 = SwStarHigh(swAtNc_High_k_Layer);
    double t7 = KrnHigh(t6);
    return(1. / k_AverageX() / (LY_low_ + LY_high_) * (LY_high_ * k_high_ * t7 + LY_low_ * k_low_ * t3));
 }
 
 
 




} // end csmp
