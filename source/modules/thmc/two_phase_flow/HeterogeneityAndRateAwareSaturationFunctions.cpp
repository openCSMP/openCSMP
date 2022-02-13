#include "CSMP_mathUtilities.h"
#include "HeterogeneityAndRateAwareSaturationFunctions.h"
#include "FlowFunctionsModule.h"
#include "Element.h"
#include "ErrorHandler.h"

//#define TURN_OFF_RATE_AWARE

using namespace std;

namespace csmp {
  
/**
   The default constructor of Brooks Corey Saturation Functions class
*/
template<uint32_t dim, template<uint32_t> class USER> 
HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::HeterogeneityAndRateAwareSaturationFunctions()
{
}


template<uint32_t dim, template<uint32_t> class USER> 
HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::~HeterogeneityAndRateAwareSaturationFunctions()
 {
 }



template<uint32_t dim, template<uint32_t> class USER>
void HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::InitilizeRockProperties( Element<dim>* const e ) const
{
    if(e->Status( User()->key_srH2O ) != DIRICH)
      e->Store(User()->key_srH2O, makeScalar( e->Status( User()->key_srH2O ), GetSwr(e) ));
      
    if(e->Status( User()->key_srCO2 ) != DIRICH)
      e->Store(User()->key_srCO2, makeScalar( e->Status( User()->key_srCO2 ), GetSnr(e) ));
      
    if(e->Status( User()->key_k ) != DIRICH)
      e->Store(User()->key_k, makeScalar( e->Status( User()->key_k ), GetK(e) ));
    if(e->Status( User()->key_kV ) != DIRICH)
      e->Store(User()->key_kV, makeScalar( e->Status( User()->key_kV ), GetKV(e) ));
    if(e->Status( User()->key_phi ) != DIRICH)
      e->Store(User()->key_phi, makeScalar( e->Status( User()->key_phi ), GetPhi(e) ));
    
    if(e->Status( User()->key_pd ) != DIRICH)
      e->Store(User()->key_pd, makeScalar( e->Status( User()->key_pd ), GetPd(e) ));
    if(e->Status( User()->key_bcp ) != DIRICH)
      e->Store(User()->key_bcp, makeScalar( e->Status( User()->key_bcp ), GetBcp(e) ));
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSwr( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double swr(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: swr = get<0>(Otway_.rocktype_).Swi_; break;
        case 1: swr = get<1>(Otway_.rocktype_).Swi_; break;
        case 2: swr = get<2>(Otway_.rocktype_).Swi_pc_; break;
        case 3: swr = get<3>(Otway_.rocktype_).Swi_; break;
        case 4: swr = get<4>(Otway_.rocktype_).Swi_; break;
        case 5: swr = get<5>(Otway_.rocktype_).Swi_pc_; break;
        case 6: swr = get<6>(Otway_.rocktype_).Swi_; break;
        case 7: swr = get<7>(Otway_.rocktype_).Swi_pc_; break;
        case 8: swr = get<8>(Otway_.rocktype_).Swi_; break;
        case 9: swr = get<9>(Otway_.rocktype_).Swi_pc_; break;
        case 10: swr = get<10>(Otway_.rocktype_).Swi_; break;
        case 11: swr = get<11>(Otway_.rocktype_).Swi_; break;
        case 12: swr = get<12>(Otway_.rocktype_).Swi_; break;
        case 13: swr = get<13>(Otway_.rocktype_).Swi_; break;
        case 14: swr = get<14>(Otway_.rocktype_).Swi_; break;
        case 15: swr = get<15>(Otway_.rocktype_).Swi_; break;
        case 16: swr = get<16>(Otway_.rocktype_).Swi_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSwr", "rocktype not recognized");
    }

    return swr;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSnr( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double snr(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: snr = get<0>(Otway_.rocktype_).Sgr_; break;
        case 1: snr = get<1>(Otway_.rocktype_).Sgr_; break;
        case 2: snr = get<2>(Otway_.rocktype_).Sgr_; break;
        case 3: snr = get<3>(Otway_.rocktype_).Sgr_; break;
        case 4: snr = get<4>(Otway_.rocktype_).Sgr_; break;
        case 5: snr = get<5>(Otway_.rocktype_).Sgr_; break;
        case 6: snr = get<6>(Otway_.rocktype_).Sgr_; break;
        case 7: snr = get<7>(Otway_.rocktype_).Sgr_; break;
        case 8: snr = get<8>(Otway_.rocktype_).Sgr_; break;
        case 9: snr = get<9>(Otway_.rocktype_).Sgr_; break;
        case 10: snr = get<10>(Otway_.rocktype_).Sgr_; break;
        case 11: snr = get<11>(Otway_.rocktype_).Sgr_; break;
        case 12: snr = 0.; break;
        case 13: snr = get<13>(Otway_.rocktype_).Sgr_; break;
        case 14: snr = get<14>(Otway_.rocktype_).Sgr_; break;
        case 15: snr = get<15>(Otway_.rocktype_).Sgr_; break;
        case 16: snr = get<16>(Otway_.rocktype_).Sgr_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSnr", "rocktype not recognized");
    }

    return snr;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPhi( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double phi(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: phi = get<0>(Otway_.rocktype_).phi_; break;
        case 1: phi = get<1>(Otway_.rocktype_).phi_; break;
        case 2: phi = get<2>(Otway_.rocktype_).phi_; break;
        case 3: phi = get<3>(Otway_.rocktype_).phi_; break;
        case 4: phi = get<4>(Otway_.rocktype_).phi_; break;
        case 5: phi = get<5>(Otway_.rocktype_).phi_; break;
        case 6: phi = get<6>(Otway_.rocktype_).phi_; break;
        case 7: phi = get<7>(Otway_.rocktype_).phi_; break;
        case 8: phi = get<8>(Otway_.rocktype_).phi_; break;
        case 9: phi = get<9>(Otway_.rocktype_).phi_; break;
        case 10: phi = get<10>(Otway_.rocktype_).phi_; break;
        case 11: phi = get<11>(Otway_.rocktype_).phi_; break;
        case 12: phi = get<12>(Otway_.rocktype_).phi_; break;
        case 13: phi = get<13>(Otway_.rocktype_).phi_; break;
        case 14: phi = get<14>(Otway_.rocktype_).phi_; break;
        case 15: phi = get<15>(Otway_.rocktype_).phi_; break;
        case 16: phi = get<16>(Otway_.rocktype_).phi_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPhi", "rocktype not recognized");
    }

    return phi;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetK( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double k(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: k = get<0>(Otway_.rocktype_).k_; break;
        case 1: k = get<1>(Otway_.rocktype_).k_; break;
        case 2: k = PermeabilityParallelToLaminations(e); break;
        case 3: k = get<3>(Otway_.rocktype_).k_; break;
        case 4: k = PermeabilityParallelToLaminations(e); break;
        case 5: k = PermeabilityParallelToLaminations(e); break;
        case 6: k = PermeabilityParallelToLaminations(e); break;
        case 7: k = PermeabilityParallelToLaminations(e); break;
        case 8: k = PermeabilityParallelToLaminations(e); break;
        case 9: k = PermeabilityParallelToLaminations(e); break;
        case 10: k = get<10>(Otway_.rocktype_).k_; break;
        case 11: k = PermeabilityParallelToLaminations(e); break;
        case 12: k = PermeabilityParallelToLaminations(e); break;
        case 13: k = get<13>(Otway_.rocktype_).k_; break;
        case 14: k = PermeabilityParallelToLaminations(e); break;
        case 15: k = get<15>(Otway_.rocktype_).k_; break;
        case 16: k = get<16>(Otway_.rocktype_).k_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetK", "rocktype not recognized");
    }

    return k;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKV( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double k(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: k = get<0>(Otway_.rocktype_).k_; break;
        case 1: k = get<1>(Otway_.rocktype_).k_; break;
        case 2: k = PermeabilityPerpendicularToLaminations(e); break;
        case 3: k = get<3>(Otway_.rocktype_).k_; break;
        case 4: k = PermeabilityPerpendicularToLaminations(e); break;
        case 5: k = PermeabilityPerpendicularToLaminations(e); break;
        case 6: k = PermeabilityPerpendicularToLaminations(e); break;
        case 7: k = PermeabilityPerpendicularToLaminations(e); break;
        case 8: k = PermeabilityPerpendicularToLaminations(e); break;
        case 9: k = PermeabilityPerpendicularToLaminations(e); break;
        case 10: k = get<10>(Otway_.rocktype_).k_; break;
        case 11: k = PermeabilityPerpendicularToLaminations(e); break;
        case 12: k = PermeabilityPerpendicularToLaminations(e); break;
        case 13: k = get<13>(Otway_.rocktype_).k_; break;
        case 14: k = PermeabilityPerpendicularToLaminations(e); break;
        case 15: k = get<15>(Otway_.rocktype_).k_; break;
        case 16: k = get<16>(Otway_.rocktype_).k_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKV", "rocktype not recognized");
    }

    return k;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSwiPc( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double swi_pc(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: swi_pc = get<0>(Otway_.rocktype_).Swi_pc_; break;
        case 1: swi_pc = get<1>(Otway_.rocktype_).Swi_pc_; break;
        case 2: swi_pc = get<2>(Otway_.rocktype_).Swi_pc_; break;
        case 3: swi_pc = get<3>(Otway_.rocktype_).Swi_pc_; break;
        case 4: swi_pc = get<4>(Otway_.rocktype_).Swi_pc_; break;
        case 5: swi_pc = get<5>(Otway_.rocktype_).Swi_pc_; break;
        case 6: swi_pc = get<6>(Otway_.rocktype_).Swi_pc_; break;
        case 7: swi_pc = get<7>(Otway_.rocktype_).Swi_pc_; break;
        case 8: swi_pc = get<8>(Otway_.rocktype_).Swi_pc_; break;
        case 9: swi_pc = get<9>(Otway_.rocktype_).Swi_pc_; break;
        case 10: swi_pc = get<10>(Otway_.rocktype_).Swi_pc_; break;
        case 11: swi_pc = get<11>(Otway_.rocktype_).Swi_pc_; break;
        case 12: swi_pc = get<12>(Otway_.rocktype_).Swi_pc_; break;
        case 13: swi_pc = get<13>(Otway_.rocktype_).Swi_pc_; break;
        case 14: swi_pc = get<14>(Otway_.rocktype_).Swi_pc_; break;
        case 15: swi_pc = get<15>(Otway_.rocktype_).Swi_pc_; break;
        case 16: swi_pc = get<16>(Otway_.rocktype_).Swi_pc_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetSwiPc", "rocktype not recognized");
    }

    return swi_pc;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetmVG( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double m_VG(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: m_VG = get<0>(Otway_.rocktype_).m_; break;
        case 1: m_VG = get<1>(Otway_.rocktype_).m_; break;
        case 2: m_VG = get<2>(Otway_.rocktype_).m_; break;
        case 3: m_VG = get<3>(Otway_.rocktype_).m_; break;
        case 4: m_VG = get<4>(Otway_.rocktype_).m_; break;
        case 5: m_VG = get<5>(Otway_.rocktype_).m_; break;
        case 6: m_VG = get<6>(Otway_.rocktype_).m_; break;
        case 7: m_VG = get<7>(Otway_.rocktype_).m_; break;
        case 8: m_VG = get<8>(Otway_.rocktype_).m_; break;
        case 9: m_VG = get<9>(Otway_.rocktype_).m_; break;
        case 10: m_VG = get<10>(Otway_.rocktype_).m_; break;
        case 11: m_VG = get<11>(Otway_.rocktype_).m_; break;
        case 12: m_VG = get<12>(Otway_.rocktype_).m_; break;
        case 13: m_VG = get<13>(Otway_.rocktype_).m_; break;
        case 14: m_VG = get<14>(Otway_.rocktype_).m_; break;
        case 15: m_VG = get<15>(Otway_.rocktype_).m_; break;
        case 16: m_VG = get<16>(Otway_.rocktype_).m_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetmVG", "rocktype not recognized");
    }

    return m_VG;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetmLow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double m_low(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: m_low = get<2>(Otway_.rocktype_).m_low_; break;
        case 4: m_low = get<4>(Otway_.rocktype_).m_low_; break;
        case 5: m_low = get<5>(Otway_.rocktype_).m_low_; break;
        case 6: m_low = get<6>(Otway_.rocktype_).m_low_; break;
        case 7: m_low = get<7>(Otway_.rocktype_).m_low_; break;
        case 8: m_low = get<8>(Otway_.rocktype_).m_low_; break;
        case 9: m_low = get<9>(Otway_.rocktype_).m_low_; break;
        case 11: m_low = get<11>(Otway_.rocktype_).m_low_; break;
        case 12: m_low = get<12>(Otway_.rocktype_).m_low_; break;
        case 14: m_low = get<14>(Otway_.rocktype_).m_low_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetmLow", "rocktype not supported");
    }

    return m_low;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPd( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double pd(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: pd = get<0>(Otway_.rocktype_).pd_; break;
        case 1: pd = get<1>(Otway_.rocktype_).pd_; break;
        case 2: pd = get<2>(Otway_.rocktype_).pd_low_; break;
        case 3: pd = get<3>(Otway_.rocktype_).pd_; break;
        case 4: pd = get<4>(Otway_.rocktype_).pd_low_; break;
        case 5: pd = get<5>(Otway_.rocktype_).pd_low_; break;
        case 6: pd = get<6>(Otway_.rocktype_).pd_low_; break;
        case 7: pd = get<7>(Otway_.rocktype_).pd_low_; break;
        case 8: pd = get<8>(Otway_.rocktype_).pd_low_; break;
        case 9: pd = get<9>(Otway_.rocktype_).pd_low_; break;
        case 10: pd = get<10>(Otway_.rocktype_).pd_; break;
        case 11: pd = get<11>(Otway_.rocktype_).pd_low_; break;
        case 12: pd = get<12>(Otway_.rocktype_).pd_low_; break;
        case 13: pd = get<13>(Otway_.rocktype_).pd_; break;
        case 14: pd = get<14>(Otway_.rocktype_).pd_low_; break;
        case 15: pd = get<15>(Otway_.rocktype_).pd_; break;
        case 16: pd = get<16>(Otway_.rocktype_).pd_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPd", "rocktype not recognized");
    }

    return pd;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPd_VG( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double pd(numeric_limits<double>::quiet_NaN());
    VectorVariable<dim> vt; 
    e->Read( User()->key_vt, vt); 
    double vt_magnitude_x = fabs(vt[0]);
    if(dim==3U) vt_magnitude_x = pow( (pow(vt[0], 2.) + pow(vt[2], 2.)) , 0.5);

   
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude_x = 1.0e-12; //set to below capillary limit
#endif   
    switch( rock_type) {
        case 0: pd = get<0>(Otway_.rocktype_).pd_; break;
        case 1: pd = get<1>(Otway_.rocktype_).pd_; break;
        case 2: pd = get<2>(Otway_.rocktype_).pd_; break;
        case 3: pd = get<3>(Otway_.rocktype_).pd_; break;
        case 4: pd = get<4>(Otway_.rocktype_).Pd( vt_magnitude_x ); break; 
        case 5: pd = get<5>(Otway_.rocktype_).pd_; break;
        case 6: pd = get<6>(Otway_.rocktype_).Pd( vt_magnitude_x ); break; 
        case 7: pd = get<7>(Otway_.rocktype_).pd_; break;
        case 8: pd = get<8>(Otway_.rocktype_).Pd( vt_magnitude_x ); break; 
        case 9: pd = get<9>(Otway_.rocktype_).pd_; break;
        case 10: pd = get<10>(Otway_.rocktype_).pd_; break;
        case 11: pd = get<11>(Otway_.rocktype_).Pd( vt_magnitude_x ); break; 
        case 12: pd = get<12>(Otway_.rocktype_).Pd( vt_magnitude_x ); break; 
        case 13: pd = get<13>(Otway_.rocktype_).pd_; break;
        case 14: pd = get<14>(Otway_.rocktype_).pd_; break;
        case 15: pd = get<15>(Otway_.rocktype_).pd_; break;
        case 16: pd = get<16>(Otway_.rocktype_).pd_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPd", "rocktype not recognized");
    }

    return pd;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPdLow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double pd_low(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: pd_low = get<2>(Otway_.rocktype_).pd_low_; break;
        case 4: pd_low = get<4>(Otway_.rocktype_).pd_low_; break;
        case 5: pd_low = get<5>(Otway_.rocktype_).pd_low_; break;
        case 6: pd_low = get<6>(Otway_.rocktype_).pd_low_; break;
        case 7: pd_low = get<7>(Otway_.rocktype_).pd_low_; break;
        case 8: pd_low = get<8>(Otway_.rocktype_).pd_low_; break;
        case 9: pd_low = get<9>(Otway_.rocktype_).pd_low_; break;
        case 11: pd_low = get<11>(Otway_.rocktype_).pd_low_; break;
        case 12: pd_low = get<12>(Otway_.rocktype_).pd_low_; break;
        case 14: pd_low = get<14>(Otway_.rocktype_).pd_low_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetPdLow", "rocktype not supported");
    }

    return pd_low;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetBcp( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double bcp(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: bcp = Bcp(e); break;
        case 1: bcp = Bcp(e); break;
        case 2: bcp = Bcp_low(e); break;
        case 3: bcp = Bcp(e); break;
        case 4: bcp = Bcp_low(e); break;
        case 5: bcp = Bcp_low(e); break;
        case 6: bcp = Bcp_low(e); break;
        case 7: bcp = Bcp_low(e); break;
        case 8: bcp = Bcp_low(e); break;
        case 9: bcp = Bcp_low(e); break;
        case 10: bcp = Bcp(e); break;
        case 11: bcp = Bcp_low(e); break;
        case 12: bcp = Bcp_low(e); break;
        case 13: bcp = Bcp(e); break;
        case 14: bcp = Bcp_low(e); break;
        case 15: bcp = Bcp(e); break;
        case 16: bcp = Bcp(e); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetBcp", "rocktype not recognized");
    }

    return bcp;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrw( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double krw(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: krw = get<0>(Otway_.rocktype_).Krw(Sw); break;
        case 1: krw = get<1>(Otway_.rocktype_).Krw(Sw); break;
        case 3: krw = get<3>(Otway_.rocktype_).Krw(Sw); break;
        case 10: krw = get<10>(Otway_.rocktype_).Krw(Sw); break;
        case 13: krw = get<13>(Otway_.rocktype_).Krw(Sw); break;
        case 15: krw = get<15>(Otway_.rocktype_).Krw(Sw); break;
        case 16: krw= get<16>(Otway_.rocktype_).Krw(Sw); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrw", "rocktype not supported");
    }

    return krw;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrn( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double krn(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 0: krn = get<0>(Otway_.rocktype_).Krn(Sw); break;
        case 1: krn = get<1>(Otway_.rocktype_).Krn(Sw); break;
        case 3: krn = get<3>(Otway_.rocktype_).Krn(Sw); break;
        case 10: krn = get<10>(Otway_.rocktype_).Krn(Sw); break;
        case 13: krn = get<13>(Otway_.rocktype_).Krn(Sw); break;
        case 15: krn = get<15>(Otway_.rocktype_).Krn(Sw); break;
        case 16: krn= get<16>(Otway_.rocktype_).Krn(Sw); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrn", "rocktype not supported");
    }

    return krn;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwParallelDrainage( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double krw_parallel_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: krw_parallel_drainage = get<2>(Otway_.rocktype_).Krw_ParallelDrainage( Sw ); break;
        case 5: krw_parallel_drainage = get<5>(Otway_.rocktype_).Krw_ParallelDrainage( Sw ); break;
        case 7: krw_parallel_drainage = get<7>(Otway_.rocktype_).Krw_ParallelDrainage( Sw ); break;
        case 9: krw_parallel_drainage = get<9>(Otway_.rocktype_).Krw_ParallelDrainage( Sw ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwParallelDrainage", "rocktype not supported");
    }

    return krw_parallel_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwParallelDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude = 1.0e-12; //set to below capillary limit
#endif 
    const int32_t rock_type = RockType(e);
    double krw_parallel_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 4: krw_parallel_drainage = get<4>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        case 6: krw_parallel_drainage = get<6>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        case 8: krw_parallel_drainage = get<8>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        case 11: krw_parallel_drainage = get<11>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        case 12: krw_parallel_drainage = get<12>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        case 14: krw_parallel_drainage = get<14>(Otway_.rocktype_).Krw_ParallelDrainage( Sw, vt_magnitude ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwParallelDrainage", "rocktype not supported");
    }

    return krw_parallel_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwCrossDrainage( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    const int32_t rock_type = RockType(e);
    double krw_cross_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: krw_cross_drainage = get<2>(Otway_.rocktype_).Krw_CrossDrainage( Sw ); break;
        case 5: krw_cross_drainage = get<5>(Otway_.rocktype_).Krw_CrossDrainage( Sw ); break;
        case 7: krw_cross_drainage = get<7>(Otway_.rocktype_).Krw_CrossDrainage( Sw ); break;
        case 9: krw_cross_drainage = get<9>(Otway_.rocktype_).Krw_CrossDrainage( Sw ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwCrossDrainage", "rocktype not supported");
    }

    return krw_cross_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwCrossDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude = 1.0e-12; //set to below capillary limit
#endif 
    const int32_t rock_type = RockType(e);
    double krw_cross_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 4: krw_cross_drainage = get<4>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        case 6: krw_cross_drainage = get<6>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        case 8: krw_cross_drainage = get<8>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        case 11: krw_cross_drainage = get<11>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        case 12: krw_cross_drainage = get<12>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        case 14: krw_cross_drainage = get<14>(Otway_.rocktype_).Krw_CrossDrainage( Sw, vt_magnitude ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrwCrossDrainage", "rocktype not supported");
    }

    return krw_cross_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnParallelDrainage( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double Krn_parallel_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: Krn_parallel_drainage = get<2>(Otway_.rocktype_).Krn_ParallelDrainage( Sw ); break;
        case 5: Krn_parallel_drainage = get<5>(Otway_.rocktype_).Krn_ParallelDrainage( Sw ); break;
        case 7: Krn_parallel_drainage = get<7>(Otway_.rocktype_).Krn_ParallelDrainage( Sw ); break;
        case 9: Krn_parallel_drainage = get<9>(Otway_.rocktype_).Krn_ParallelDrainage( Sw ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnParallelDrainage", "rocktype not supported");
    }

    return Krn_parallel_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnParallelDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude = 1.0e-12; //set to below capillary limit
#endif 
    const int32_t rock_type = RockType(e);
    double Krn_parallel_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 4: Krn_parallel_drainage = get<4>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        case 6: Krn_parallel_drainage = get<6>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        case 8: Krn_parallel_drainage = get<8>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        case 11: Krn_parallel_drainage = get<11>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        case 12: Krn_parallel_drainage = get<12>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        case 14: Krn_parallel_drainage = get<14>(Otway_.rocktype_).Krn_ParallelDrainage( Sw, vt_magnitude ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnParallelDrainage", "rocktype not supported");
    }

    return Krn_parallel_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnCrossDrainage( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double Krn_cross_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: Krn_cross_drainage = get<2>(Otway_.rocktype_).Krn_CrossDrainage( Sw ); break;
        case 5: Krn_cross_drainage = get<5>(Otway_.rocktype_).Krn_CrossDrainage( Sw ); break;
        case 7: Krn_cross_drainage = get<7>(Otway_.rocktype_).Krn_CrossDrainage( Sw ); break;
        case 9: Krn_cross_drainage = get<9>(Otway_.rocktype_).Krn_CrossDrainage( Sw ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnCrossDrainage", "rocktype not supported");
    }

    return Krn_cross_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnCrossDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude = 1.0e-12; //set to below capillary limit
#endif 
    const int32_t rock_type = RockType(e);
    double Krn_cross_drainage(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 4: Krn_cross_drainage = get<4>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        case 6: Krn_cross_drainage = get<6>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        case 8: Krn_cross_drainage = get<8>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        case 11: Krn_cross_drainage = get<11>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        case 12: Krn_cross_drainage = get<12>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        case 14: Krn_cross_drainage = get<14>(Otway_.rocktype_).Krn_CrossDrainage( Sw, vt_magnitude ); break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKrnCrossDrainage", "rocktype not supported");
    }

    return Krn_cross_drainage;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetLYLow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double LY_low(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: LY_low = get<2>(Otway_.rocktype_).LY_low_; break;
        case 4: LY_low = get<4>(Otway_.rocktype_).LY_low_; break;
        case 5: LY_low = get<5>(Otway_.rocktype_).LY_low_; break;
        case 6: LY_low = get<6>(Otway_.rocktype_).LY_low_; break;
        case 7: LY_low = get<7>(Otway_.rocktype_).LY_low_; break;
        case 8: LY_low = get<8>(Otway_.rocktype_).LY_low_; break;
        case 9: LY_low = get<9>(Otway_.rocktype_).LY_low_; break;
        case 11: LY_low = get<11>(Otway_.rocktype_).LY_low_; break;
        case 12: LY_low = get<12>(Otway_.rocktype_).LY_low_; break;
        case 14: LY_low = get<14>(Otway_.rocktype_).LY_low_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetLYLow", "rocktype not supported");
    }

    return LY_low;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetLYHigh( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double LY_high(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: LY_high = get<2>(Otway_.rocktype_).LY_high_; break;
        case 4: LY_high = get<4>(Otway_.rocktype_).LY_high_; break;
        case 5: LY_high = get<5>(Otway_.rocktype_).LY_high_; break;
        case 6: LY_high = get<6>(Otway_.rocktype_).LY_high_; break;
        case 7: LY_high = get<7>(Otway_.rocktype_).LY_high_; break;
        case 8: LY_high = get<8>(Otway_.rocktype_).LY_high_; break;
        case 9: LY_high = get<9>(Otway_.rocktype_).LY_high_; break;
        case 11: LY_high = get<11>(Otway_.rocktype_).LY_high_; break;
        case 12: LY_high = get<12>(Otway_.rocktype_).LY_high_; break;
        case 14: LY_high = get<14>(Otway_.rocktype_).LY_high_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetLYHigh", "rocktype not supported");
    }

    return LY_high;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKLow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double k_low(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: k_low = get<2>(Otway_.rocktype_).k_low_; break;
        case 4: k_low = get<4>(Otway_.rocktype_).k_low_; break;
        case 5: k_low = get<5>(Otway_.rocktype_).k_low_; break;
        case 6: k_low = get<6>(Otway_.rocktype_).k_low_; break;
        case 7: k_low = get<7>(Otway_.rocktype_).k_low_; break;
        case 8: k_low = get<8>(Otway_.rocktype_).k_low_; break;
        case 9: k_low = get<9>(Otway_.rocktype_).k_low_; break;
        case 11: k_low = get<11>(Otway_.rocktype_).k_low_; break;
        case 12: k_low = get<12>(Otway_.rocktype_).k_low_; break;
        case 14: k_low = get<14>(Otway_.rocktype_).k_low_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetkLow", "rocktype not supported");
    }

    return k_low;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKHigh( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const int32_t rock_type = RockType(e);
    double k_high(numeric_limits<double>::quiet_NaN());
    switch( rock_type) {
        case 2: k_high = get<2>(Otway_.rocktype_).k_high_; break;
        case 4: k_high = get<4>(Otway_.rocktype_).k_high_; break;
        case 5: k_high = get<5>(Otway_.rocktype_).k_high_; break;
        case 6: k_high = get<6>(Otway_.rocktype_).k_high_; break;
        case 7: k_high = get<7>(Otway_.rocktype_).k_high_; break;
        case 8: k_high = get<8>(Otway_.rocktype_).k_high_; break;
        case 9: k_high = get<9>(Otway_.rocktype_).k_high_; break;
        case 11: k_high = get<11>(Otway_.rocktype_).k_high_; break;
        case 12: k_high = get<12>(Otway_.rocktype_).k_high_; break;
        case 14: k_high = get<14>(Otway_.rocktype_).k_high_; break;
        default:
            cerr <<"\n\t rocktype: "<< rock_type;
            csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::GetKHigh", "rocktype not supported");
    }

    return k_high;
}


/// get seff at the element barycentre
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::EffectiveSaturation( Element<dim>* const e ) const
{
    double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    double Swi_pc = GetSwiPc(e);
    return seffL(Sw,Swi_pc);
}



/// get seff from the supplied saturation value
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::EffectiveSaturation_at( Element<dim>* const e, double Sw ) const
{
    assert( Sw >= 0. );
    assert( Sw <= 1. );

    double Swi_pc = GetSwiPc(e);
    return seffL(Sw,Swi_pc);
}



template<uint32_t dim, template<uint32_t> class USER>
int32_t HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::RockType( Element<dim>* const e ) const
{
     const double rock_type = e->Read(User()->key_RRT);
     assert( !isnan(rock_type) );
     return static_cast<int32_t>( rock_type );    
}



template<uint32_t dim, template<uint32_t> class USER>
bool HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::IsComposite( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    bool is_composite = false;
    const int  rock_type = RockType(e);
    if(rock_type==2 || rock_type==4 || rock_type==5 || rock_type==6 || rock_type==7 || rock_type==8 || \
       rock_type==9 || rock_type==11 || rock_type==12 || rock_type==14)
       is_composite = true;
    if(rock_type<0 || rock_type>16) {
       cerr <<"\n\t rocktype: "<< rock_type;
       csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::IsComposite", "rocktype not recognized");    
    }
    
    return is_composite;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::Bcp( Element<dim>* const e ) const
{
    const int  rock_type = RockType(e);
    if(rock_type==0) return 0.; //WELL
    
    // VG m parameter
    double m_VG= GetmVG(e);
    // Brooks-Corey lambda parameter from VG m parameter, Lenhard et al. (1989)
    double bcp = 0.5e-1 * std::exp(5.8 * m_VG) + 1.;

    return bcp;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::Bcp_low( Element<dim>* const e ) const
{    
    // VG m parameter
    double m_low = GetmLow(e);
    // Brooks-Corey lambda parameter from VG m parameter, Lenhard et al. (1989)
    double bcp_low= 0.5e-1 * std::exp(5.8 * m_low) + 1.;
        
    return bcp_low;
}



template<uint32_t dim, template<uint32_t> class USER>
VectorVariable<dim> HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::InitializeVelocity( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    VectorVariable<dim> vt, vt_normalised;
    e->Read( User()->key_vt, vt ); 
    
    bool has_nan(false);
    for(int i=0U;i<dim;i++)
        if(isnan(vt[i])) {has_nan = true; break;}
    
    if ( has_nan )
      csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::InitializeVelocity ", "the velocity variable has not been initialised.");
    
    double vt_magnitude = vt.Length();
    // in zero velocity case, the horizontal relative permeability is set to dominate
    if ( vt_magnitude <= numeric_limits<double>::epsilon() * 100. ) {
         if(dim == 2U) { 
             vt_normalised(0) = 1.; 
             vt_normalised(1) = 0.;
         }
         else if(dim == 3U) { //x^2 + y^2 + z^2 = 1; x = z; y = 0
             vt_normalised(0) = pow(0.5, 0.5); 
             vt_normalised(1) = 0.;
             vt_normalised(2) = pow(0.5, 0.5);
         }
    }
    else vt_normalised = vt / vt_magnitude;
    
    return vt_normalised;    
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc( Element<dim>* const e ) const
{
    const int  rock_type = RockType(e); 
    
    //only rocktype 16 uses Brooks Corey
    if(rock_type == 16) return pc_BrookCorey(e);
     
    return pc_VanGenuchten(e);
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc_at( Element<dim>* const e, double Sw ) const
{
    const int  rock_type = RockType(e); 
    
    //only rocktype 16 uses Brooks Corey
    if(rock_type == 16) return pc_at_BrookCorey(e, Sw);

    return pc_at_VanGenuchten(e, Sw);
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc_BrookCorey( Element<dim>* const e ) const
{
    bool is_composite = IsComposite(e);
 
    if ( !is_composite ) {
         // if BC-lambda = 0, pc is assumed to be pd, ie. constant
         double bcp = Bcp(e);
         assert(!isnan(bcp));
         if ( bcp <= numeric_limits<double>::epsilon() ) return GetPd(e);
         double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
         double Swi_pc = GetSwiPc(e);
         double pd = GetPd(e);
         return min( pc_BC( Sw, Swi_pc, pd, bcp ), MaxCapillaryPressure() );
      }
        
    // for horizontal flow, weighted average is used
    double bcp_low = Bcp_low(e);
    double pd_low = GetPdLow(e);
    double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    double Swi_pc = GetSwiPc(e);
    
    const double pc_low = (bcp_low  == 0.) ? pd_low : pc_BC( Sw, Swi_pc, pd_low, bcp_low );
    
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    
    double pc = pc_low * (vt_normalised[0]*vt_normalised[0]) + pc_low * (vt_normalised[1] * vt_normalised[1]);
    
    if(dim==3U) pc += pc_low * (vt_normalised[2]*vt_normalised[2]); 
                        
    return min( pc, MaxCapillaryPressure()); 
    
} // end pc_BrookCorey
  
 

template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc_at_BrookCorey( Element<dim>* const e, double Sw ) const
{
    assert( Sw >= 0. );
    assert( Sw <= 1. );
    
    bool is_composite = IsComposite(e);
 
    if ( !is_composite ) {
         // if BC-lambda = 0, pc is assumed to be pd, ie. constant
         double bcp = Bcp(e);
         assert(!isnan(bcp));
         if ( bcp <= numeric_limits<double>::epsilon() ) return GetPd(e);
         double Swi_pc = GetSwiPc(e);
         double pd = GetPd(e);
         return min( pc_BC( Sw, Swi_pc, pd, bcp ), MaxCapillaryPressure() );
      }
        
    // for horizontal flow, weighted average is used
    double bcp_low = Bcp_low(e);
    double pd_low = GetPdLow(e);
    double Swi_pc = GetSwiPc(e);
    
    const double pc_low = (bcp_low  == 0.) ? pd_low : pc_BC( Sw, Swi_pc, pd_low, bcp_low );

    if(dim==1U) return min( pc_low, MaxCapillaryPressure()); 

    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    
    double pc = pc_low * (vt_normalised[0]*vt_normalised[0]) + pc_low * (vt_normalised[1] * vt_normalised[1]);

    if(dim==3U) pc += pc_low * (vt_normalised[2]*vt_normalised[2]); 
                        
    return min( pc, MaxCapillaryPressure());     
} // end pc_at_BrookCorey




template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc_VanGenuchten( Element<dim>* const e ) const
{ 

    const double Pc_MAX(1.0e+7);
    double Swi_pc = GetSwiPc(e);
    double pd = GetPd_VG(e);
    double m_VG = GetmVG(e);    
    double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    const double PC_HIGH_SW_LIMIT = 1.0 - 0.01;

    // linear regularization for lower part of sw range
    if ( Sw < PC_LOW_SW_LIMIT ){
        const double pc_lim = pc_VG( PC_LOW_SW_LIMIT, pd, m_VG, Swi_pc );
        const double seff_mult( 1.0/ (1.0 - Swi_pc ) );
        const double dpcds_lim = dpcds_at_VG( e, PC_LOW_SW_LIMIT )/seff_mult;
        const double pc = pc_lim + dpcds_lim*( Sw - PC_LOW_SW_LIMIT );
        return std::min(pc, Pc_MAX);    
    
    // linear regularization for higher part of sw range
    }else if ( Sw > PC_HIGH_SW_LIMIT ){
        const double pc_lim = pc_VG( PC_HIGH_SW_LIMIT, pd, m_VG, Swi_pc );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT ) );
        const double pc = dpcds_lim*( Sw - 1.0 );
        return std::max(pc, 0.);
    }

    return std::min( pc_VG( Sw, pd, m_VG, Swi_pc ), Pc_MAX);    
    
} // end pc_VanGenuchten



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::pc_at_VanGenuchten( Element<dim>* const e, double Sw ) const
{

    const double Pc_MAX(1.0e+7);
    double Swi_pc = GetSwiPc(e);
    double pd = GetPd_VG(e);
    double m_VG = GetmVG(e);    
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    const double PC_HIGH_SW_LIMIT = 1.0 - 0.01;

    // linear regularization for lower part of sw range
    if ( Sw < PC_LOW_SW_LIMIT ){
        const double pc_lim = pc_VG( PC_LOW_SW_LIMIT, pd, m_VG, Swi_pc );
        const double seff_mult( 1.0/ (1.0 - Swi_pc ) );
        const double dpcds_lim = dpcds_at_VG( e, PC_LOW_SW_LIMIT )/seff_mult;
        const double pc = pc_lim + dpcds_lim*( Sw - PC_LOW_SW_LIMIT );
        return std::min(pc, Pc_MAX);    
    
    // linear regularization for higher part of sw range
    }else if ( Sw > PC_HIGH_SW_LIMIT ){
        const double pc_lim = pc_VG( PC_HIGH_SW_LIMIT, pd, m_VG, Swi_pc );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT ) );
        const double pc = dpcds_lim*( Sw - 1.0 );
        return std::max(pc, 0.);
    }

    return std::min( pc_VG( Sw, pd, m_VG, Swi_pc ), Pc_MAX);    
} //end pc_at_VanGenuchten



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds( Element<dim>* const e ) const
{
    const int  rock_type = RockType(e); 
    
    //use constant derivative at PC_LOW_SW_LIMIT for lower part of sw range
    const double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    const double Swi_pc = GetSwiPc(e);
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    if(Sw < PC_LOW_SW_LIMIT) {
        //only rocktype 16 uses Brooks Corey
        if(rock_type == 16) return dpcds_at_BC(e, PC_LOW_SW_LIMIT);
        return dpcds_at_VG(e, PC_LOW_SW_LIMIT);
    }
    
    //only rocktype 16 uses Brooks Corey
    if(rock_type == 16) return dpcds_BC(e);
    
    return dpcds_VG(e); 
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_at( Element<dim>* const e, double sw ) const
{
    const int  rock_type = RockType(e); 
    
    //use constant derivative at PC_LOW_SW_LIMIT for lower part of sw range
    const double Swi_pc = GetSwiPc(e);
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    if(sw < PC_LOW_SW_LIMIT) {
        //only rocktype 16 uses Brooks Corey
        if(rock_type == 16) return dpcds_at_BC(e, PC_LOW_SW_LIMIT);
        return dpcds_at_VG(e, PC_LOW_SW_LIMIT);
    }
    
    //only rocktype 16 uses Brooks Corey
    if(rock_type == 16) return dpcds_at_BC(e, sw);
    
    return dpcds_at_VG(e, sw); 
}

  
   
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_BC( Element<dim>* const e ) const
{
    double h(0.001);
    return dpcds_Numerical_BC(e, h);  
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_at_BC( Element<dim>* const e, double sw ) const
{
    double h(0.001);
    return dpcds_at_Numerical_BC(e, sw, h);  
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_VG( Element<dim>* const e ) const
{
    double h(0.001);
    
    //use constant derivative at PC_LOW_SW_LIMIT for lower part of sw range
    double Swi_pc = GetSwiPc(e);
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    if(sw < PC_LOW_SW_LIMIT) return dpcds_at_Numerical_VG(e, PC_LOW_SW_LIMIT, h );
        
    return dpcds_Numerical_VG(e, h);  
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_at_VG( Element<dim>* const e, double sw ) const
{
    double h(0.001);
    
    //use constant derivative at PC_LOW_SW_LIMIT for lower part of sw range
    double Swi_pc = GetSwiPc(e);
    const double PC_LOW_SW_LIMIT = Swi_pc + 0.01;
    if(sw < PC_LOW_SW_LIMIT) return dpcds_at_Numerical_VG(e, PC_LOW_SW_LIMIT, h );

    return dpcds_at_Numerical_VG(e, sw, h);  
}

 
 
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel( Element<dim>* const e, size_t direction ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if(direction != 0 && direction != 2) 
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel", "wrong input flow direction, mush be 0 or 2");
    
    double parallel_krw(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        parallel_krw = GetKrwParallelDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ //composite and rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[direction]);
        parallel_krw = GetKrwParallelDrainage( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel", "this function only supports composite rock type");
    }
    
    return parallel_krw;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_crossflow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    double crossflow_krw(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        crossflow_krw = GetKrwCrossDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ //composite and rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[1]);
        crossflow_krw = GetKrwCrossDrainage( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_crossflow_x", "this function only supports composite rock type");
    }
    
    return crossflow_krw;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel_at( Element<dim>* const e, double Sw, size_t direction ) const 
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if(direction != 0 && direction != 2) 
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel_at", "wrong input flow direction, mush be 0 or 2");    
    
    double krw_parallel(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        krw_parallel = GetKrwParallelDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[direction]);
        krw_parallel = GetKrwParallelDrainage( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_parallel_at", "this function only supports composite rock type");
    }
    
    return krw_parallel;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_crossflow_at( Element<dim>* const e, double Sw  ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    double krw_crossflow(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        krw_crossflow = GetKrwCrossDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[1]);
        krw_crossflow = GetKrwCrossDrainage ( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_crossflow_at", "this function only supports composite rock type");
    }
    
    return krw_crossflow;
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::PermeabilityPerpendicularToLaminations( Element<dim>* const e ) const
{
    bool is_composite = IsComposite(e);
    assert( is_composite );
    
    double L_low = GetLYLow(e);
    double L_high = GetLYHigh(e);
    const double sum_of_weights = L_low + L_high;
    
    double k_high = GetKHigh(e);
    double k_low = GetKLow(e);
    
    return sum_of_weights / (L_high / k_high + L_low / k_low);    
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::PermeabilityParallelToLaminations( Element<dim>* const e ) const
{
    bool is_composite = IsComposite(e);
    assert( is_composite );
    
    double L_low = GetLYLow(e);
    double L_high = GetLYHigh(e);
    
    double k_high = GetKHigh(e);
    double k_low = GetKLow(e);
   
    return (L_high * k_high + L_low * k_low) / (L_high + L_low);
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::PermeabilityInFlowDirection( Element<dim>* const e,
                                                                                              const TensorVariable<dim>& KK) const
{
    bool is_composite = IsComposite(e);
    assert( is_composite );
    assert( KK(0,0) > 0. );
    assert( KK(1,1) > 0. );
    if(dim==3U) assert( KK(2,2) > 0. );
    
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    
    if ( fabs(vt_normalised.Length() - 1.) <= numeric_limits<double>::epsilon() )
      // horizontal permeability
      return KK(0,0);

    // finding the permeability in the direction of the velocity vector
    VectorVariable<dim>  vc = KK * vt_normalised;

    // the length of vc_ is equal to the magnitude of permeability in the target direction
    return vc.Length();
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::K_reduction_in_flow_direction( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    bool is_composite = IsComposite(e);
    double K_reduction(numeric_limits<double>::quiet_NaN());
    
    if(is_composite) {
        const double k_crossflow = PermeabilityPerpendicularToLaminations(e);
        const double k_parallel  = PermeabilityParallelToLaminations(e);
        // assuming that layers are horizontal and that the stored K is the horizontal one
        TensorVariable<dim> KK;
        KK = 0.; // all entries = zero
        KK(0,0) = k_parallel;
        KK(1,1) = k_crossflow;
        if(dim==3U) KK(2,2) = k_parallel; 
        // uses KK(0,0) internally
        double K_flow_direction = PermeabilityInFlowDirection( e, KK );
        K_reduction = K_flow_direction / k_parallel;         
        
    } else {
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::K_reduction_in_flow_direction", "this function only supports composite rock type");
    }
    
    return K_reduction;
}
  

  
  
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw( Element<dim>* const e ) const
{
    bool is_composite = IsComposite(e);
    if ( !is_composite ) {
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        return GetKrw(e, Sw); // krw_VG( Sw_, m_VG_ );
    }
    
    // to get the ensemble krw for the composite, the parallel and perpendicular values are blended
    // taking into account the flow direction 
    // --------------------------------------
    double parallel_krw_x = krw_parallel(e,0); //x direction 
    double crossflow_krw = krw_crossflow(e); //y direcction
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    
    double krw = parallel_krw_x * (vt_normalised[0]*vt_normalised[0]) + crossflow_krw * (vt_normalised[1]*vt_normalised[1]); 
    if(dim==3U) {
        double parallel_krw_z = krw_parallel(e,2); //z direction 
        krw += parallel_krw_z * (vt_normalised[2]*vt_normalised[2]); 
    }
    // scaling the relative permeability by the vertical permeability
    double K_reduction = K_reduction_in_flow_direction(e);
    
    return max( krw * K_reduction, 0. );
}
  
  
   
 
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krw_at( Element<dim>* const e, double Sw ) const
{
    bool is_composite = IsComposite(e);
    if ( !is_composite ) {
        return GetKrw(e, Sw); // krw_VG( Sw_, m_VG_ );
    }
    
    // to get the ensemble krw for the composite, the parallel and perpendicular values are blended
    // taking into account the flow direction 
    // --------------------------------------
    double parallel_krw_x = krw_parallel_at(e, Sw, 0); //x direction 
    double crossflow_krw = krw_crossflow_at(e, Sw); //y direction
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    double krw = parallel_krw_x * (vt_normalised[0]*vt_normalised[0]) + crossflow_krw * (vt_normalised[1]*vt_normalised[1]); 

    if(dim==3U) { 
        double parallel_krw_z = krw_parallel_at(e, Sw, 2); //z direction 
        krw += parallel_krw_z * (vt_normalised[2]*vt_normalised[2]); 
    } 

    // scaling the relative permeability by the vertical permeability
    double K_reduction = K_reduction_in_flow_direction(e);
    
    return max( krw * K_reduction, 0. );
}




template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel( Element<dim>* const e, size_t direction ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if(direction != 0 && direction != 2) 
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel", "wrong input flow direction, mush be 0 or 2");     
    
    double parallel_krn(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        parallel_krn = GetKrnParallelDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[direction]);
        parallel_krn = GetKrnParallelDrainage( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel", "this function only supports composite rock type");
    }
    
    return parallel_krn;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_crossflow( Element<dim>* const e ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    double crossflow_krn(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        crossflow_krn = GetKrnCrossDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[1]);
        crossflow_krn = GetKrnCrossDrainage ( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_crossflow", "this function only supports composite rock type");
    }
    
    return crossflow_krn;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel_at( Element<dim>* const e, double Sw, size_t direction ) const 
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if(direction != 0 && direction != 2) 
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel_at", "wrong input flow direction, mush be 0 or 2");       
    
    double krn_parallel(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        krn_parallel = GetKrnParallelDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[direction]);
        krn_parallel = GetKrnParallelDrainage( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_parallel_at", "this function only supports composite rock type");
    }
    
    return krn_parallel;
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_crossflow_at( Element<dim>* const e, double Sw ) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    double krn_crossflow(numeric_limits<double>::quiet_NaN());
    const int  rock_type = RockType(e);

    if(rock_type==2 || rock_type==5 || rock_type==7 || rock_type==9) { //composite but not rate-dependent
        krn_crossflow = GetKrnCrossDrainage( e, Sw );
    } else if(rock_type==4 || rock_type==6 || rock_type==8 || rock_type==11 || rock_type==12 || rock_type==14){ ////composite and rate-dependent
        VectorVariable<dim> vt;
        e->Read( User()->key_vt, vt);
        double vt_magnitude = fabs(vt[1]); 
        krn_crossflow = GetKrnCrossDrainage ( e, Sw, vt_magnitude );
    } else { //not composite
        csmp_error.notice( ERROR, "HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_crossflow_at", "this function only supports composite rock type");
    }
    
    return krn_crossflow;
}


  
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn( Element<dim>* const e ) const
{
    bool is_composite = IsComposite(e);
    if ( !is_composite ) {
        double Sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
        return GetKrn(e, Sw); // krn_VG( Sw_, m_VG_ );
    }
    
    double parallel_krn_x = krn_parallel(e, 0); //x direction
    double crossflow_krn = krn_crossflow(e); //y direction
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    double krn = parallel_krn_x * (vt_normalised[0]*vt_normalised[0]) + crossflow_krn * (vt_normalised[1]*vt_normalised[1]);
    
    if(dim==3U) { 
        double parallel_krn_z = krn_parallel(e, 2); //z direction 
        krn += parallel_krn_z * (vt_normalised[2]*vt_normalised[2]);
    }
    
    // scaling the relative permeability by the vertical permeability
    double K_reduction = K_reduction_in_flow_direction(e);
    
    return max( krn * K_reduction, 0. );
}




template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::krn_at( Element<dim>* const e, double Sw ) const
{
    assert( Sw >= 0. );
    assert( Sw <= 1. );
    
    bool is_composite = IsComposite(e);
    if ( !is_composite ) {
        return GetKrn(e, Sw); // krn_VG( Sw_, m_VG_ );
    }
    
    // to get the ensemble krw for the composite, the parallel and perpendicular values are blended
    // taking into account the flow direction 
    // --------------------------------------
    double parallel_krn_x = krn_parallel_at(e, Sw, 0); //x direction
    double crossflow_krn = krn_crossflow_at(e, Sw); //y direction
    VectorVariable<dim> vt_normalised = InitializeVelocity(e);
    double krn = parallel_krn_x * (vt_normalised[0]*vt_normalised[0]) + crossflow_krn * (vt_normalised[1]*vt_normalised[1]);
    if(dim==3U) {
        double parallel_krn_z = krn_parallel_at(e, Sw, 2); //z direction
        krn += parallel_krn_z * (vt_normalised[2]*vt_normalised[2]);
    }

    // scaling the relative permeability by the vertical permeability
    double K_reduction = K_reduction_in_flow_direction(e);
    
    return max( krn * K_reduction, 0. );
}


  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrwds( Element<dim>* const e ) const
{
    double h(0.001);
    return dkrwds_Numerical(e, h);  
}

  

/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrwds_at( Element<dim>* const e, double sw ) const
{
    double h(0.001);
    return dkrwds_at_Numerical(e, sw, h);  
}

  
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrnds( Element<dim>* const e ) const
{
    double h(0.001);
    return dkrnds_Numerical(e, h);   
}
   
  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrnds_at( Element<dim>* const e, double sw ) const
{
    double h(0.001);
    return dkrnds_at_Numerical(e, sw, h);    
}

  

// Numerical derivative
  
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrwds_Numerical( Element<dim>* const e, double h ) const
{
    double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double krw1 = krw_at(e, 1.);
         double krw2 = krw_at(e, 1. - h);
         return (krw1 - krw2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double krw1 = krw_at(e, h);
         double krw2 = krw_at(e, 0.);
         return (krw1 - krw2) / h;
    }

    // water saturation between the endpoints
    double krw1 = krw_at(e, sw + h);
    double krw2 = krw_at(e, sw - h);   
    return (krw1 - krw2) / (2. * h);
}
  
  
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrwds_at_Numerical( Element<dim>* const e, double sw, double h ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double krw1 = krw_at(e, 1.);
         double krw2 = krw_at(e, 1. - h);
         return (krw1 - krw2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double krw1 = krw_at(e, h);
         double krw2 = krw_at(e, 0.);
         return (krw1 - krw2) / h;
    }

    // water saturation between the endpoints
    double krw1 = krw_at(e, sw + h);
    double krw2 = krw_at(e, sw - h);   
    return (krw1 - krw2) / (2. * h);
}


template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrnds_Numerical( Element<dim>* const e, double h ) const
{
    double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double krn1 = krn_at(e, 1.);
         double krn2 = krn_at(e, 1. - h);
         return (krn1 - krn2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double krn1 = krn_at(e, h);
         double krn2 = krn_at(e, 0.);
         return (krn1 - krn2) / h;
    }

    // water saturation between the endpoints
    double krn1 = krn_at(e, sw + h);
    double krn2 = krn_at(e, sw - h);   
    return (krn1 - krn2) / (2. * h);
}
 
  
template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dkrnds_at_Numerical( Element<dim>* const e, double sw, double h ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double krn1 = krn_at(e, 1.);
         double krn2 = krn_at(e, 1. - h);
         return (krn1 - krn2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double krn1 = krn_at(e, h);
         double krn2 = krn_at(e, 0.);
         return (krn1 - krn2) / h;
    }

    // water saturation between the endpoints
    double krn1 = krn_at(e, sw + h);
    double krn2 = krn_at(e, sw - h);   
    return (krn1 - krn2) / (2. * h);
}
 
  

template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_Numerical_BC( Element<dim>* const e, double h ) const
{
    double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double pc1 = pc_at_BrookCorey(e, 1.);
         double pc2 = pc_at_BrookCorey(e, 1. - h);
         return (pc1 - pc2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double pc1 = pc_at_BrookCorey(e, h);
         double pc2 = pc_at_BrookCorey(e, 0.);
         return (pc1 - pc2) / h;
    }

    // water saturation between the endpoints
    double pc1 = pc_at_BrookCorey(e, sw + h);
    double pc2 = pc_at_BrookCorey(e, sw - h);   
    return (pc1 - pc2) / (2. * h);
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_at_Numerical_BC( Element<dim>* const e, double sw, double h ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double pc1 = pc_at_BrookCorey(e, 1.);
         double pc2 = pc_at_BrookCorey(e, 1. - h);
         return (pc1 - pc2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double pc1 = pc_at_BrookCorey(e, h);
         double pc2 = pc_at_BrookCorey(e, 0.);
         return (pc1 - pc2) / h;
    }

    // water saturation between the endpoints
    double pc1 = pc_at_BrookCorey(e, sw + h);
    double pc2 = pc_at_BrookCorey(e, sw - h);   
    return (pc1 - pc2) / (2. * h);
}




template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_Numerical_VG( Element<dim>* const e, double h ) const
{
    double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
    assert( sw >= 0. );
    assert( sw <= 1. );
    
    h = 0.001;
    
    const double pd = GetPd_VG(e);
    const double m_VG = GetmVG(e);
    const double Swi_pc = GetSwiPc(e);     
    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double pc1 = pc_VG( 1., pd, m_VG, Swi_pc );
         double pc2 = pc_VG( 1.-h, pd, m_VG, Swi_pc );
         return (pc1 - pc2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double pc1 = pc_VG( h, pd, m_VG, Swi_pc );
         double pc2 = pc_VG( 0., pd, m_VG, Swi_pc );
         return (pc1 - pc2) / h;
    }

    // water saturation between the endpoints
    double pc1 = pc_VG( sw + h, pd, m_VG, Swi_pc );
    double pc2 = pc_VG( sw - h, pd, m_VG, Swi_pc );
    return (pc1 - pc2) / (2. * h);
}



template<uint32_t dim, template<uint32_t> class USER>
double HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::dpcds_at_Numerical_VG( Element<dim>* const e, double sw, double h ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );

    h = 0.001;
    
    const double pd = GetPd_VG(e);
    const double m_VG = GetmVG(e);
    const double Swi_pc = GetSwiPc(e); 

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double pc1 = pc_VG( 1., pd, m_VG, Swi_pc );
         double pc2 = pc_VG( 1.-h, pd, m_VG, Swi_pc );
         return (pc1 - pc2) / h;
    }
  
    // low water saturation
    if ( sw <= h ) {
         double pc1 = pc_VG( h, pd, m_VG, Swi_pc );
         double pc2 = pc_VG( 0., pd, m_VG, Swi_pc );
         return (pc1 - pc2) / h;
    }

    // water saturation between the endpoints
    double pc1 = pc_VG( sw + h, pd, m_VG, Swi_pc );
    double pc2 = pc_VG( sw - h, pd, m_VG, Swi_pc );
    return (pc1 - pc2) / (2. * h);
}



/**
    writes textfile with sw, krw(sw,Nc), krn(sw,Nc), and pc(sw) values computed for (composite) rocktype in 0.005 saturation increments
*/
template<uint32_t dim, template<uint32_t> class USER>
void HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::WriteRelativePermeabilityTable (const char* filename, long RT, Element<dim>* const e)
{

    std::string text_file = "RT-";
    text_file += to_string(RT);
    std::string file_name(filename);
    text_file += ("-" + file_name);

    ofstream  ofs( string(text_file) + ".txt" );
   
    const string   rocktype(parseRockType(RT));
    
    ofs <<"rocktype="<< rocktype <<"="<< RT <<endl;
    ofs <<"sw\t krw \t krn \t pc \t dpcds"<<endl;
    for ( double sw(0.); sw<1.005; sw+=0.005 )
      {
         // dynamic parameters
         ofs << sw << "\t"<< krw_at(e, sw);
         ofs <<"\t"<< krn_at(e, sw);
         ofs <<"\t"<< pc_at(e, sw);
         ofs <<"\t"<< dpcds_at(e, sw);
         ofs << endl;
      }
} // end WriteRelativePermeabilityTable
 


template<uint32_t dim, template<uint32_t> class USER>
void HeterogeneityAndRateAwareSaturationFunctions<dim,USER>::OutputTestingResults(Element<dim>* const e)
{

    VectorVariable<dim> vt; 
    vt=0.;
    double vt_magnitude(0.);

    //1. horizontal flows
    //1.1 vt(0) < CL
    vt_magnitude = 4.0e-8;
    vt(0) = vt_magnitude; vt(1) = 0.;
    if(dim==3U) vt(2) = vt_magnitude;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Horizontal_CL_4e-8", i, e );
    }
    
    //1.2 CL < vt(0) < VL: 3 velocities
    vt_magnitude = 5.0e-5;
    vt(0) = vt_magnitude ; vt(1) = 0.;
    if(dim==3U) vt(2) = vt_magnitude;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Horizontal_CL_VL_5.0e-5", i, e );
    }
    
    vt_magnitude = 5.0e-4;
    vt(0) = vt_magnitude; vt(1) = 0.;
    if(dim==3U) vt(2) = vt_magnitude;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Horizontal_CL_VL_5.0e-4", i, e );
    }
    
    vt_magnitude = 5.0e-3;
    vt(0) = vt_magnitude; vt(1) = 0.;
    if(dim==3U) vt(2) = vt_magnitude;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Horizontal_CL_VL_5.0e-3", i, e );
    }
    
    //1.3 vt(0) > VL
    vt_magnitude = 6.0e-1;
    vt(0) = vt_magnitude; vt(1) = 0.;
    if(dim==3U) vt(2) = vt_magnitude;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Horizontal_VL_6.0e-1", i, e );
    }   

    //2. vertical flows
    //2.1 vt(1) < CL
    vt_magnitude = 4.0e-8;
    vt(0) = 0; vt(1) = vt_magnitude;
    if(dim==3U) vt(2) = 0.;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Vertical_CL_4.0e-8", i, e );
    } 
    
    //2.2 CL < vt(1) < VL: 3 velocities
    vt_magnitude = 5.0e-6;
    vt(0) = 0.; vt(1) = vt_magnitude;
    if(dim==3U) vt(2) = 0.;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Vertical_CL_VL_5.0e-6", i, e );
    } 
    
    vt_magnitude = 5.0e-5;
    vt(0) = 0.; vt(1) = vt_magnitude;
    if(dim==3U) vt(2) = 0.;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Vertical_CL_VL_5.0e-5", i, e );
    } 
    
    vt_magnitude = 5.0e-4;
    vt(0) = 0.; vt(1) = vt_magnitude;
    if(dim==3U) vt(2) = 0.;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Vertical_CL_VL_5.0e-4", i, e );
    } 
    
    //2.3 vt(1) > VL
    vt_magnitude = 6.0e-2;
    vt(0) = 0.; vt(1) = vt_magnitude;
    if(dim==3U) vt(2) = 0.;
    e->Store(User()->key_vt, vt);
    for ( long i=0; i<=16; i++ ) { // for all rocktypes
        e->Store(User()->key_RRT, makeScalar( e->Status( User()->key_RRT ), i));   
        WriteRelativePermeabilityTable( "Vertical_VL_6.0e-2", i, e );
    }    
}
  

template class HeterogeneityAndRateAwareSaturationFunctions<1U,FlowFunctionsModule7>;
template class HeterogeneityAndRateAwareSaturationFunctions<2U,FlowFunctionsModule7>;
template class HeterogeneityAndRateAwareSaturationFunctions<3U,FlowFunctionsModule7>;


} // csmp
