#ifndef EXPERIMENTAL_RT_H
#define EXPERIMENTAL_RT_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"
#include <vector>

namespace csmp {


template<size_t dim>
class ExperimentalRT : public TwoPhaseModel<dim> {

public:

    ExperimentalRT( const PropertyDatabase<dim>& database,
                    const char* rt_file_name,
                    const char* rt_number,
                    const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false
                           
    ExperimentalRT( const PropertyDatabase<dim>& database,
                    const char* rt_file_name,
                    const char* rt_number,
                    const char* permeability,
                    const char* viscosity_nw, const char* viscosity_w,
                    const char* density_nw, const char* density_w,
                    const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                    const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    virtual ~ExperimentalRT();
    
    void ConstructRTs(const char* rt_file_name);
    void Initialize(const Element<dim>& e);

    // relative permeabilities
    virtual double64 krn_Phase() const;
    virtual double64 krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double64 dkrnds_Phase() const;
    virtual double64 dkrwds_Phase() const;

    // capillary pressure
    virtual double64 pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double64 dpcds_Phase( ) const;

private:
    std::vector<csmp::CubicSpline> krn_, krw_, pc_;
    csmp::Index rt_key_;
    unsigned int rt_number_;


};


// Input file structure for rock types
// To use experimental data for rock types you need to define  
// RELPERM_MODEL in Numerical.txt file as EXPERIMENTAL
// Input file for rock types should be put in the execution folder
// and the file name that contains relative permeability data and
// capillary values should be 'rock type.txt'
// File structure is as follow:
//2 ->number of tables
//21 ->number of entries for table 0 (first table) excluding the first row
//-1.99513	-0.00488	0.000125	3.709875	-10557.3	-519.567 ->derivatives at start and end for kro, krw and pc
//0	1	0	5000 -> data (sw  kro krw pc)
//0.05	0.90024375	0.00000625	4472.135955
//0.1	0.8019	0.0001	3162.27766
//0.15	0.70624375	0.00050625	2581.988897
//0.2	0.6144	0.0016	2236.067977
//0.25	0.52734375	0.00390625	2000
//0.3	0.4459	0.0081	1825.741858
//0.35	0.37074375	0.01500625	1690.308509
//0.4	0.3024	0.0256	1581.13883
//0.45	0.24124375	0.04100625	1490.711985
//0.5	0.1875	0.0625	1414.213562
//0.55	0.14124375	0.09150625	1348.399725
//0.6	0.1024	0.1296	1290.994449
//0.65	0.07074375	0.17850625	1240.347346
//0.7	0.0459	0.2401	1195.228609
//0.75	0.02734375	0.31640625	1154.700538
//0.8	0.0144	0.4096	1118.033989
//0.85	0.00624375	0.52200625	1084.652289
//0.9	0.0019	0.6561	1054.092553
//0.95	0.00024375	0.81450625	1025.978352
//1	0	1	1000
//6 ->number of entries for table 1 (second table) excluding the first row
//-1	-1	1	1	0	0 ->derivatives
//0	1	0	0 ->data
//0.2	0.8	0.2	0
//0.4	0.6	0.4	0
//0.6	0.4	0.6	0
//0.8	0.2	0.8	0
//1	0	1	0


} // end namespace csmp

#endif
