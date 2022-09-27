//
//  SandPropertiesFor_VE_Model.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 26/8/2022.
//

#include "SandPropertiesFor_VE_Model.h"
#include "CSMP_mathUtilities.h"

#include "ErrorHandler.h"
#include "Model.h"
#include "Region.h"
#include "Node.h"
#include "compareFloats.h"

// RULE: Prefer std::string_view over std::string when you need a read-only string, especially for function parameters.

using namespace std;

namespace csmp {

SandPropertiesFor_VE_Model::SandPropertiesFor_VE_Model( const Model<2U>& model,
                                                        const string& dimM1_region_name,
                                                        double bcp, double swr, double snr )
 : target_region_(dimM1_region_name), bcp_(bcp), swr_(swr), snr_(snr),
   thi_key_(model.Database().StorageKey("thickness")),
   vt_key_(model.Database().StorageKey("total velocity")),
   dip_key_(model.Database().StorageKey("dip vector")),
   sw_key_(model.Database().StorageKey("saturation aqueous phase")),
   bcp_key_(model.Database().StorageKey("brooks corey parameter")),
   swr_key_(model.Database().StorageKey("residual saturation aqueous phase")),
   snr_key_(model.Database().StorageKey("residual saturation carbonic phase")),
   mobw_key_(model.Database().StorageKey("mobility aqueous phase")),
   mobn_key_(model.Database().StorageKey("mobility carbonic phase"))
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    if ( !model.ContainsRegion( dimM1_region_name ) )
      csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                      "target region does not exist");
                      
    constexpr uint32_t dim{ 2U };
    const Region<dim>& subdomain = model.Region( dimM1_region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 ) {
         csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                                 "target region is not lower dimensional (dim region != dim-1)");
         target_region_consists_of_line_elmts_ = false;
      }
                      
    // checking that the parameter in the target region which will be modified are indeed
    // single valued
    double vmin, vmax;
    subdomain.MinMaxOf( "brooks corey parameter", vmin, vmax );
    if ( !approximatelyEqual(vmin,vmax) ) {
         cout <<"\nSandPropertiesFor_VE_Model: original range of 'brooks corey parameter'"<< endl;
         printRangeOfVariable( model, target_region_.c_str(), "brooks corey parameter" );
         cout <<"new value that will be assigned: "<< bcp << endl;
      }
    subdomain.MinMaxOf( "residual saturation aqueous phase", vmin, vmax );
    if ( !approximatelyEqual(vmin,vmax) ) {
         cout <<"\nSandPropertiesFor_VE_Model: original range of 'residual saturation aqueous phase'"<< endl;
         printRangeOfVariable( model, target_region_.c_str(), "residual saturation aqueous phase" );
         cout <<"new value that will be assigned: "<< swr << endl;
      }
    subdomain.MinMaxOf( "residual saturation carbonic phase", vmin, vmax );
    if ( !approximatelyEqual(vmin,vmax) ) {
         cout <<"\nSandPropertiesFor_VE_Model: original range of 'residual saturation carbonic phase'"<< endl;
         printRangeOfVariable( model, target_region_.c_str(), "residual saturation carbonic phase" );
         cout <<"new value that will be assigned: "<< snr << endl;
      }
    // fluid viscosities
    subdomain.MinMaxOf( "viscosity aqueous phase", vmin, vmax );
    if ( !approximatelyEqual(vmin,vmax) ) {
         cout <<"\nSandPropertiesFor_VE_Model: original range of 'viscosity aqueous phase'"<< endl;
         muw_ = printRangeOfVariable( model, target_region_.c_str(), "viscosity aqueous phase" );
         cout <<"this value will be assigned: "<< muw_ << endl;
      }
    subdomain.MinMaxOf( "viscosity carbonic phase", vmin, vmax );
    if ( !approximatelyEqual(vmin,vmax) ) {
         cout <<"\nSandPropertiesFor_VE_Model: original range of 'viscosity carbonic phase'"<< endl;
         mun_ = printRangeOfVariable( model, target_region_.c_str(), "viscosity carbonic phase" );
         cout <<"this value will be assigned: "<< mun_ << endl;
      }
    
 } // end constructor





/**
    Computes multiphase flow parameters for the finite volumes associated with the line elements of the lower-dimensional
    region associated with the sand properties object.
    
    @attention method must be called after the standard multiphase properties were calculated because it
    uses these property values.
*/
void SandPropertiesFor_VE_Model::Compute2PhaseFlowPropertiesForSandLayer( Model<2U>& model )
 {
    constexpr uint32_t dim{ 2U };
    Region<dim>& subdomain = model.Region( target_region_ );

    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    // looping over the elements computing the input parameters for the relperm adjustment
    VectorVariable<dim> vt;
    double              sw;
    
    for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); ++it )
      {
         // 1. computation of FV diameter, L, and vertical extent, H (direction given by node numbers)
         auto ncoord1 = (*it)->N(0U)->Coordinate();
         auto ncoord2 = (*it)->N(1U)->Coordinate();
         const double vE = fabs( ncoord1[1] - ncoord2[1] ); // vertical extent = height
         // only a linear relperm model needs to be used if the FV is horizontal
         if ( approximatelyEqual(vE,0) ) {
              (*it)->Store( bcp_key_, makeScalar(ANY,0.) );
              // calculation of linear relative permeability-model parameters
              (*it)->Read( vt_key_, vt ); // total velocity of fluid mixture needed to find upstream sw
              const bool right_pointing = ( (*it)->N(0)->x() < (*it)->N(0)->x() ) ? true : false;
              if ( right_pointing )
                   sw = ( vt[0] >= 0. ) ? (*it)->N(0)->Read(sw_key_) : (*it)->N(1)->Read(sw_key_);
              else sw = ( vt[0] >= 0. ) ? (*it)->N(1)->Read(sw_key_) : (*it)->N(0)->Read(sw_key_);
              const double sw_eff  = min( max( (sw - swr_) / (1. - swr_ - snr_), 0. ), 1. );
              const double krw = sw_eff;
              const double krn = 1. - sw_eff;
              (*it)->Store( mobw_key_, makeScalar( (*it)->Status(mobw_key_), krw / muw_ ) );
              //(*it)->Store( krw_key_, makeScalar( (*it)->Status(krw_key_), krw ) );
              (*it)->Store( mobn_key_, makeScalar( (*it)->Status(mobn_key_), krn / mun_ ) );
              //(*it)->Store( krn_key_, makeScalar( (*it)->Status(krn_key_), krn ) );
              continue;
           }
         // else, the facet is tilted and the relative permeabilities need to be scaled
         const double H    = (*it)->Read( thi_key_ );
         const double hCO2 = CO2LevelAtFacet( (*it) ); // returns 'h'
         // TODO: mob_co2 must be the mobility at the end-point saturation
         const double mob_co2 = (hCO2 / H) * (*it)->Read( mobn_key_ );
         // saving the new mobility of the carbonic phase & relperm
         (*it)->Store( mobn_key_, makeScalar( (*it)->Status(mobn_key_), mob_co2 ) );
         //(*it)->Store( krn_key_, makeScalar( (*it)->Status(krn_key_), mob_co2 * mun_ ) );
         // water mobility & relperm
         // TODO: mob_h2o must be the mobility at the end-point saturation
         const double mob_h2o = ((H - hCO2) / H) * (*it)->Read( mobw_key_ );
         (*it)->Store( mobw_key_, makeScalar( (*it)->Status(mobw_key_), mob_h2o ) );
         //(*it)->Store( krw_key_, makeScalar( (*it)->Status(krw_key_), mob_h2o * muw_ ) );
         
      } // end for lower-dim elements in target region
      
 } // end Compute2PhaseFlowPropertiesForSandLayer
 
 
 
 
 
/**
    Finds h/H at upstream or downstream facet of inclided FV from the updip or down-dip water saturation and the flow direction.
    In the case of updip flow, the saturation of the lower FV is used and the opposite for down-dip flow
    
*/
double SandPropertiesFor_VE_Model::CO2LevelAtFacet( const Element<2U>* const eptr ) const
  {
     assert( eptr != nullptr );
     csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );

     // 1. computation of FV diameter, L, and vertical extent, H (direction given by node numbers)
     // ------------------------------------------------------------------------------------------
     auto ncoord1 = eptr->N(0U)->Coordinate();
     auto ncoord2 = eptr->N(1U)->Coordinate();
     const double V = fabs( ncoord1[1] - ncoord2[1] ); // vertical extent of upper triangle of FV
     const double L = ncoord1.DistanceTo( ncoord2 );   // length = diameter of FV
     const double thickness = eptr->Read( thi_key_ );  // H = height of facet = thickness attribute of lower-dim element
     if ( thickness >= L ) {
         cout <<"\n"<<"cell length vs thickness: "<< L <<" "<< thickness;
         csmp_error.Note( WARNING, "SandPropertiesFor_VE_Model::CO2WaterColumnHeightRatioAtFacet",
                         "thickness of cell greater than diameter; not likely to yield correct results");
       }
     // for essentially horizontal FVs, V = 0., and the input saturation is equivalent to the CO2 / water column ratio
     if ( V <= numeric_limits<double>::epsilon() * 100. ) {
          const double sw = (eptr->N(0U)->Read( sw_key_ ) + eptr->N(1U)->Read( sw_key_ )) / 2.;
          return (1. - sw) / sw;
       }
     // if there is a dip, the dip angle of element (in radians) in calculated
     ncoord1[1] = ncoord2[1] = 0.;
     const double L_horizontal = ncoord1.DistanceTo( ncoord2 );
     const double dip_angle = atan( V / L_horizontal );
     const double dip = radiansToDegrees( dip_angle ); // angle in degrees
     assert( dip <= 90. );
       
       
     // 2. flow direction
     // ------------------------------------------------------------------------------------------
     VectorVariable<2U> vt;
     eptr->Read( vt_key_, vt ); // total velocity of fluid mixture in FV from previous time step
     // since the flow is aligned with the dim-m1 element it can be used to diagnose flow mode
     enum { DOWN_DIP_FLOW, UP_DIP_FLOW } flow = (vt[1] >= 0.) ? UP_DIP_FLOW : DOWN_DIP_FLOW;
     
     
     // 3. Computing CO2 column height vs facet height from the water saturation
     // ------------------------------------------------------------------------------------------
     // the facet is on the upper end of the finite volume
     // and the facet is within this current element with its thickness attribute
     const double fv_vol = thickness * L;
     
     if ( flow == UP_DIP_FLOW ) {
          // in up-dip flow the upstream saturation comes from the FV below the facet and vice versa
          const double sw_FV  = ( vt[0] >= 0. ) ? eptr->N(0U)->Read(sw_key_) : eptr->N(1U)->Read(sw_key_);
          // if dip angle 45o, only the upper triangle of the FV needs to be filled for
          // the whole facet to be covered with CO2
          // (the following applies to all dip angles up to 90o)
          // the ratio hCO2 / H <= 1 of the upstream facet is given by
          return min( 1., sqrt( 2. * fv_vol * (1. - sw_FV) / tan(dip_angle) ) );
       }
     else { // DOWN_DIP_FLOW
          // the down-dip face will only start to see CO2 once most of the FV is filled
          // in up-dip flow the upstream saturation comes from the FV below the facet and vice versa
          const double sw_FV  = ( vt[0] < 0. ) ? eptr->N(0U)->Read(sw_key_) : eptr->N(1U)->Read(sw_key_);          
          // if the water saturation is greater than what it takes to fill the lower triangular area
          // the ratio hCO2 / H will be zero, H = thickness
          const double V_crit = 0.5 * thickness * thickness / tan(dip_angle);
          if ( fv_vol * sw_FV >= V_crit ) return 0.;
          else {
               // the ratio hCO2 / H scales positively with the water saturation between sw=0 and sw_crit
               return max( 0., (1. - sqrt( 2. * fv_vol * sw_FV / tan( dip_angle) )) );
            }
       }

 } // end CO2LevelAtFacet

 
 
 //template SandPropertiesFor_VE_Model<2U>; // only 2D so far (special treatment)


} // end csmp

