// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ExperimentalRT.h"


namespace csmp {


template<uint32_t dim>
ExperimentalRT<dim>::ExperimentalRT( const PropertyDatabase<dim>& database,
                                     const char* rt_file_name,
                                     const char* rt_number,
                                     const bool sw_ro_mu_placement )

    : TwoPhaseModel<dim>(database, "permeability",
                         "viscosity oil", "viscosity water",
                         "density oil", "density water", "saturation water",
                         "residual saturation non-wetting phase",
                         "residual saturation wetting phase",
                         sw_ro_mu_placement ),
      rt_key_(database.StorageKey(rt_number))
{

    ConstructRTs(rt_file_name);

}

template<uint32_t dim>
ExperimentalRT<dim>::ExperimentalRT( const PropertyDatabase<dim>& database,
                                     const char* rt_file_name,
                                     const char* rt_number,
                                     const char* permeability,
                                     const char* viscosity_nw, const char* viscosity_w,
                                     const char* density_nw, const char* density_w,
                                     const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                                     const bool sw_ro_mu_placement )

    :    TwoPhaseModel<dim>(database, permeability,
                            viscosity_nw, viscosity_w,
                            density_nw, density_w,
                            sat_w, res_sat_nw, res_sat_w,
                            sw_ro_mu_placement ),
      rt_key_(database.StorageKey(rt_number))
{

    ConstructRTs(rt_file_name);

}


template<uint32_t dim>
ExperimentalRT<dim>::~ExperimentalRT()
{
}


template<uint32_t dim>
void ExperimentalRT<dim>::ConstructRTs(const char* rt_file_name)
{
    std::cout << "\nExperimentalRT::ConstructRTs()\n";
    
    std::ifstream rt_file;
    rt_file.open(rt_file_name);

    if (!rt_file)
    {
        std::cout << "\nExperimentalRT::ConstructRTs(): " << rt_file_name << " file doesn't exist.\n";
        std::exit(1);
    }
    else
    {
        unsigned int number_of_tables;
        rt_file >> number_of_tables;
        
        std::cout << "\nExperimentalRT::ConstructRTs(): Number of tables: " << number_of_tables << "\n";
        
        krn_.resize(number_of_tables);
        krw_.resize(number_of_tables);
        pc_.resize(number_of_tables);

        for (unsigned int i = 0; i < number_of_tables; i++)
        {
            std::cout << "\nTable " << i << "\n";
            
            double kro_start_derivative, kro_end_derivative,
                    krw_start_derivative, krw_end_derivative,
                    pc_start_derivative, pc_end_derivative;
            std::vector<double> sw, kro, krw, pc;
            double sw_value, kro_value, krw_value, pc_value;
            unsigned int number_of_entries;
            
            rt_file >> number_of_entries;
            
            rt_file >> kro_start_derivative >> kro_end_derivative
                    >> krw_start_derivative >> krw_end_derivative
                    >> pc_start_derivative >> pc_end_derivative;
            
            std::cout << "Derivatives:\t" << kro_start_derivative << "\t" << kro_end_derivative
                      << "\t" << krw_start_derivative << "\t" << krw_end_derivative
                      << "\t" << pc_start_derivative << "\t" << pc_end_derivative << "\n";
            
            std::cout << "sw\tkro\tkrw\tpc\n";
            
            for (unsigned int n = 0; n < number_of_entries; n++)
            {
                rt_file >> sw_value >> kro_value >> krw_value >> pc_value;
                
                std::cout << sw_value << "\t" << kro_value << "\t" << krw_value << "\t" << pc_value <<"\n";
                
                sw.push_back(sw_value);
                kro.push_back(kro_value);
                krw.push_back(krw_value);
                pc.push_back(pc_value);

            }

            std::cout << "\n";
            
            krn_[i].Initialize(sw, kro, kro_start_derivative, kro_end_derivative);
            krw_[i].Initialize(sw, krw, krw_start_derivative, krw_end_derivative);
            pc_[i].Initialize(sw, pc, pc_start_derivative, pc_end_derivative);
            
        }
        
    }
    
}






/*
template<uint32_t dim>
void ExperimentalRT<dim>::ConstructRTs(const char* rt_file_name)
{
    std::cout << "\nExperimentalRT::ConstructRTs()\n";
    
    std::ifstream rt_file;
    rt_file.open(rt_file_name);

    if (!rt_file)
    {
        std::cout << "\nExperimentalRT::ConstructRTs(): " << rt_file_name << " file doesn't exist.\n";
        std::exit(1);
    }
    else
    {
        std::string temp;
        int position;
        int tables_counter(0);
        
        while (!rt_file.eof())
        {
            position = rt_file.tellg();
            getline(rt_file, temp);

            if (int(temp[0]) < 48 || int(temp[0]) > 57)
            {
                continue;
            }
            else
            {
                break;
            }

        }

        rt_file.seekg(position);
        
        unsigned int number_of_tables;
        
        rt_file >> number_of_tables;
        
        std::cout << "\nExperimentalRT::ConstructRTs(): Number of tables: " << number_of_tables << "\n";
        
        kro_.resize(number_of_tables);
        krw_.resize(number_of_tables);
        pc_.resize(number_of_tables);
        
        // to read the whole line and go to the next
        rt_file.seekg(position);
        getline(rt_file, temp);
        
        for (unsigned int i = 0; i < number_of_tables; i++)
        {
            std::cout << "Table " << i << "\n";
            std::cout << "sw\tkro\tkrw\tpc\n";
            
            position = rt_file.tellg();
            getline(rt_file, temp);
            while (int(temp[0]) < 48 || int(temp[0]) > 57)
            {
                position = rt_file.tellg();
                getline(rt_file, temp);
            }
            
            rt_file.seekg(position);

            double kro_start_derivative, kro_end_derivative,
                     krw_start_derivative, krw_end_derivative,
                     pc_start_derivative, pc_end_derivative;
            std::vector<double> sw, kro, krw, pc;
            double sw_value, kro_value, krw_value, pc_value;
            
            position = rt_file.tellg();
            
            rt_file >> kro_start_derivative >> kro_end_derivative
                    >> krw_start_derivative >> krw_end_derivative
                    >> pc_start_derivative >> pc_end_derivative;
                    
            std::cout << "\n" << kro_start_derivative << "\t" << kro_end_derivative
                      << "\t" << krw_start_derivative << "\t" << krw_end_derivative
                      << "\t" << pc_start_derivative << "\t" << pc_end_derivative << "\n";

            // to read the whole line and go to the next
            rt_file.seekg(position);
            getline(rt_file, temp);
            
            position = rt_file.tellg();
            getline(rt_file, temp);
            while (!(int(temp[0]) < 48 || int(temp[0]) > 57))
            {
                rt_file.seekg(position);
                rt_file >> sw_value >> kro_value >> krw_value >> pc_value;
                
                std::cout << sw_value << "\t" << kro_value << "\t" << krw_value << "\t" << pc_value <<"\n";
                
                sw.push_back(sw_value);
                kro.push_back(kro_value);
                krw.push_back(krw_value);
                pc.push_back(pc_value);

                // to read the whole line and go to the next
                rt_file.seekg(position);
                getline(rt_file, temp);

                position = rt_file.tellg();
                getline(rt_file, temp);

            }
            
            std::cout << "\n";
            
            kro_[i].Initialize(sw, kro, kro_start_derivative, kro_end_derivative);
            krw_[i].Initialize(sw, krw, krw_start_derivative, krw_end_derivative);
            pc_[i].Initialize(sw, pc, pc_start_derivative, pc_end_derivative);

        }

    }
    
    
}

*/


template<uint32_t dim>
void ExperimentalRT<dim>::Initialize(const Element<dim>& e)
{
    TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
    TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );
    rt_number_ = e.Read(rt_key_);

    if( TwoPhaseModel<dim>::tensor_permeability_){
        e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );
    }

} // end Initialize



template<uint32_t dim>
double ExperimentalRT<dim>::krn_Phase() const
{
    if (TwoPhaseModel<dim>::sat_ <= TwoPhaseModel<dim>::swr_) return krn_[rt_number_].Value(TwoPhaseModel<dim>::swr_);
    if (TwoPhaseModel<dim>::sat_ >= 1. - TwoPhaseModel<dim>::snr_) return krn_[rt_number_].Value(1. - TwoPhaseModel<dim>::snr_);
    return krn_[rt_number_].Value(TwoPhaseModel<dim>::sat_);
}


template<uint32_t dim>
double ExperimentalRT<dim>::krw_Phase() const
{
    if (TwoPhaseModel<dim>::sat_ <= TwoPhaseModel<dim>::swr_) return krw_[rt_number_].Value(TwoPhaseModel<dim>::swr_);
    if (TwoPhaseModel<dim>::sat_ >= 1. - TwoPhaseModel<dim>::snr_) return krw_[rt_number_].Value(1. - TwoPhaseModel<dim>::snr_);
    return krw_[rt_number_].Value(TwoPhaseModel<dim>::sat_);
}


template<uint32_t dim>
double ExperimentalRT<dim>::pc_Phase( ) const
{
    if (TwoPhaseModel<dim>::sat_ <= TwoPhaseModel<dim>::swr_) return pc_[rt_number_].Value(TwoPhaseModel<dim>::swr_);
    if (TwoPhaseModel<dim>::sat_ >= 1. - TwoPhaseModel<dim>::snr_) return pc_[rt_number_].Value(1. - TwoPhaseModel<dim>::snr_);
    return pc_[rt_number_].Value(TwoPhaseModel<dim>::sat_);
}


template<uint32_t dim>
double ExperimentalRT<dim>::dkrnds_Phase() const
{
    if (TwoPhaseModel<dim>::sat_ <= TwoPhaseModel<dim>::swr_) return krn_[rt_number_].Derivative(TwoPhaseModel<dim>::swr_);
    if (TwoPhaseModel<dim>::sat_ >= 1. - TwoPhaseModel<dim>::snr_) return krn_[rt_number_].Derivative(1. - TwoPhaseModel<dim>::snr_);
    return krn_[rt_number_].Derivative(TwoPhaseModel<dim>::sat_);
}


template<uint32_t dim>
double ExperimentalRT<dim>::dkrwds_Phase() const
{
    if (TwoPhaseModel<dim>::sat_ <= TwoPhaseModel<dim>::swr_) return krw_[rt_number_].Derivative(TwoPhaseModel<dim>::swr_);
    if (TwoPhaseModel<dim>::sat_ >= 1. - TwoPhaseModel<dim>::snr_) return krw_[rt_number_].Derivative(1. - TwoPhaseModel<dim>::snr_);
    return krw_[rt_number_].Derivative(TwoPhaseModel<dim>::sat_);
}

// capillary pressure derivatives
template<uint32_t dim>
double ExperimentalRT<dim>::dpcds_Phase() const
{
    return pc_[rt_number_].Derivative( TwoPhaseModel<dim>::sat_);
}


template class ExperimentalRT<1U>;
template class ExperimentalRT<2U>;
template class ExperimentalRT<3U>;


}
