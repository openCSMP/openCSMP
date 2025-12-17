#include "IMPES_Setup.h"
#include <iostream>
#include <fstream>
#include <map>
#include "PropertyDatabase.h"
#include <limits>

namespace csmp{

template<uint32_t dim>
IMPES_Setup<dim>::IMPES_Setup()
 : with_capillary_forces_(false),
   with_gravitational_forces_(false)
 { DefaultSetting(); }


template<uint32_t dim>
 IMPES_Setup<dim>::IMPES_Setup(bool with_gravity, bool with_capillary)

 : with_capillary_forces_(with_capillary),
   with_gravitational_forces_(with_gravity)
 { DefaultSetting(); }


 template<uint32_t dim>
 IMPES_Setup<dim>::~IMPES_Setup()
{ }



 template<uint32_t dim>
 void IMPES_Setup<dim>::WriteVariableFile( const char* file_name )
 {
     std::ifstream variable_file_check;

     variable_file_check.open(file_name);

     if (!variable_file_check)
     {
         variable_file_check.close();

         std::ofstream variable_file(file_name);

         std::cout << "\nIMPES_Setup::WriteVariableFile() " << file_name << " file doesn't exist. It will be created automatically.\n";
         
         variable_file << "name\tMParameter\tunit\tindex\tmin.\tmax.\tplace\n";
         variable_file << "\n";
         variable_file << "fluid pressure\t\t\t\t\tPF\t\tPa\t\t\t1\t0.E+10\t1.E+10\tnode\n";
         variable_file << "permeability\t\t\t\t\tK\t\tm2\t\t\t1\t1.E-25\t1.E+02\telement\n";
         variable_file << "oil volume source\t\t\t\t\tQOV\t\tm3 m-3 s-1\t\t1\t-1.E+02\t1.E+02\telement\n";
         variable_file << "fluid volume source\t\t\t\tQFV\t\tm3 m-3 s-1\t\t1\t-1.E+02\t1.E+02\telement\n";
         variable_file << "water volume source\t\t\t\tQWV\t\tm3 m-3 s-1\t\t1\t-1.E+02\t1.E+02\telement\n";
         variable_file << "nodal fluid volume source\t\t\tNQV\t\tm3 s-1\t\t1\t-1.E+02\t1.E+02\tnode\n";
         variable_file << "nodal oil volume source\t\t\t\tNOV\t\tm3 s-1\t\t1\t-1.E+02\t1.E+02\tnode\n";
         variable_file << "nodal water volume source\t\t\tNWV\t\tm3 s-1\t\t1\t-1.E+02\t1.E+02\tnode\n";
         variable_file << "total mobility\t\t\t\t\tTM\t\tm2\t\t\t1\t-1.E+02\t1.E+02\telement\n";
         variable_file << "density oil\t\t\t\t\t\tRHOO\t\tkg m-3\t\t1\t1.E-01\t2.E+03\telement\n";
         variable_file << "density water\t\t\t\t\tRHOW\t\tkg m-3\t\t1\t1.E-01\t2.E+03\telement\n";
         variable_file << "viscosity oil\t\t\t\t\tMUO\t\tPa s\t\t\t1\t1.E-05\t1.00E-01\telement\n";
         variable_file << "viscosity water\t\t\t\t\tMUW\t\tPa s\t\t\t1\t1.E-05\t1.60E-01\telement\n";
         variable_file << "saturation oil\t\t\t\t\tSO\t\tX\t\t\t1\t-5.00E-02\t1.05E+00\telement\n";
         variable_file << "saturation water\t\t\t\t\tSW\t\tX\t\t\t1\t-5.00E-02\t1.05E+00\telement\n";
         variable_file << "brooks corey parameter\t\t\t\tBC\t\tNA\t\t\t1\t0.0\t\t5.0\t\telement\n";
         variable_file << "entry pressure\t\t\t\t\tPD\t\tPa\t\t\t1\t0.0E+00\t5.0E+07\telement\n";
         variable_file << "residual saturation wetting phase\t\tSWR\t\tX\t\t\t1\t0.00E+00\t5.00E-01\telement\n";
         variable_file << "residual saturation non-wetting phase\tSNR\t\tX\t\t\t1\t0.00E+00\t5.00E-01\telement\n";
         variable_file << "gravity term\t\t\t\t\tGT\t\tkg Pa-1 s-3\t\t2\t-1.E+05\t1.E+05\telement\n";
         variable_file << "capillary term\t\t\t\t\tCT\t\tm s-1\t\t\t2\t-1.E+05\t1.E+05\telement\n";
         variable_file << "porosity\t\t\t\t\t\tPHI\t\tX\t\t\t1\t1.0E-05\t1.00E+00\telement\n";
         variable_file << "velocity\t\t\t\t\t\tVT\t\tm s-1\t\t\t2\t-1.0E+02\t1.00E+02\telement\n";
         variable_file << "velocity oil\t\t\t\t\tVO\t\tm s-1\t\t\t2\t-1.0E+02\t1.00E+02\telement\n";
         variable_file << "saturation\t\t\t\t\t\tSO\t\tX\t\t\t1\t-5.00E-02\t1.05E+00\tnode\n";
         variable_file << "boundary source\t\t\t\t\tNS\t\tm3 m-2 s-1\t\t1\t-1.E+02\t1.E+02\tnode\n";
         variable_file << "divergence\t\t\t\t\t\tDG\t\tX\t\t\t1\t-1.0E+02\t1.00E+02\telement\n";
         variable_file << "volume modifier\t\t\t\t\tVM\t\tX\t\t\t1\t0.\t\t1.00E+02\telement\n";
         variable_file << "rock type\t\t\t\t\tRT\t\tX\t\t\t1\t0.\t\t1.00E+02\telement\n";

         variable_file.close();
     }
     else
     {
         variable_file_check.close();
     }

 }



 template<uint32_t dim>
 void IMPES_Setup<dim>::DefaultSetting()
 {
    if (with_gravitational_forces_)
    {
        normal_saturation_changes_ = 0.01;
        maximum_saturation_changes_ = 0.025;
        CFL_multiplier_ = 0.95;
        maximum_number_transport_timesteps_ = 30U;
        maximum_transport_timestep_size_ = std::numeric_limits<double>::max();
        with_divergence_correction_ = true;
        formulation_ = GLOBAL;
        relperm_model_ = BROOKSCOREY;
        solution_ = IMPES;
        activation_criteria_ = 1.e-10;
        limiter_ = 1.e-3;
        deactivate_ = 0U;
        samg_reconstruction_ = 10U;
    }
    else if (with_capillary_forces_)
    {
        normal_saturation_changes_ = 0.025;
        maximum_saturation_changes_ = 0.05;
        CFL_multiplier_ = 0.95;
        maximum_number_transport_timesteps_ = 50U;
        maximum_transport_timestep_size_ = std::numeric_limits<double>::max();
        with_divergence_correction_ = true;
        formulation_ = GLOBAL;
        relperm_model_ = BROOKSCOREY;
        solution_ = IMPES;
        activation_criteria_ = 1.e-10;
        limiter_ = 1.e-3;
        deactivate_ = 0U;
        samg_reconstruction_ = 10U;
    }
    else
    {
        normal_saturation_changes_ = 0.05;
        maximum_saturation_changes_ = 0.1;
        CFL_multiplier_ = 0.95;
        maximum_number_transport_timesteps_ = 75U;
        maximum_transport_timestep_size_ = std::numeric_limits<double>::max();
        with_divergence_correction_ = true;
        formulation_ = GLOBAL;
        relperm_model_ = BROOKSCOREY;
        solution_ = IMPES;
        activation_criteria_ = 1.e-10;
        limiter_ = 1.e-3;
        deactivate_ = 0U;
        samg_reconstruction_ = 10U;
    }
    
    output_once_.insert("permeability");
    
    output_always_.insert("saturation oil");
    output_once_.insert("fluid pressure");
    output_always_.insert("velocity");
    
 
 }

 template<uint32_t dim>
 void IMPES_Setup<dim>::ReadSettingsFromFile(const PropertyDatabase<dim>& database)
 {
    std::cout << "\nIMPES_Setup::ReadSettingsFromFile()\n";
    
    // numerical settings always should come in a file with following name
    const char* file_name("Numerical.txt");
    
    std::ifstream numerical_file;
    
    numerical_file.open(file_name);

    if (!numerical_file)
    {
        std::cout << "\nIMPES_Setup::ReadSettingsFromFile() 'Numerical.txt' file doesn't exist. Default settings will be applied.\n";
    }
    else
    {
        std::map<std::string, double> number_variables;
        std::map<std::string, std::string> text_variables;
        std::set<std::string> output_once;
        std::set<std::string> output_always;
        bool output_once_keyword(false);
        bool output_always_keyword(false);
        
        std::set<std::string> all_keywords;
        
        all_keywords.insert("DIVERGENCE_CORRECTION");
        all_keywords.insert("NORM_SAT_CHANGE");
        all_keywords.insert("MAX_SAT_CHANGE");
        all_keywords.insert("CFL_MULTIPLIER");
        all_keywords.insert("MAX_TRANSPORT_STEPS");
        all_keywords.insert("MAX_TRANSPORT_TIME_SIZE");
        all_keywords.insert("RELPERM_MODEL");
        all_keywords.insert("FORMULATION");
        all_keywords.insert("OUTPUT_ONCE");
        all_keywords.insert("OUTPUT_ALWAYS");
        all_keywords.insert("SOLUTION");
        all_keywords.insert("ACTIVATION_CRITERIA");
        all_keywords.insert("LIMITER");
        all_keywords.insert("DEACTIVATION");
#ifdef CSMP_WITH_SAMG_SOLVER
        all_keywords.insert("SAMG_RECONSTRUCTION");
#endif
        // creating the maps
        number_variables["DIVERGENCE_CORRECTION"] = (int) with_divergence_correction_;
        number_variables["NORM_SAT_CHANGE"] = normal_saturation_changes_;
        number_variables["MAX_SAT_CHANGE"] = maximum_saturation_changes_;
        number_variables["CFL_MULTIPLIER"] = CFL_multiplier_;
        number_variables["MAX_TRANSPORT_STEPS"] = maximum_number_transport_timesteps_;
        number_variables["MAX_TRANSPORT_TIME_SIZE"] = maximum_transport_timestep_size_;
        number_variables["ACTIVATION_CRITERIA"] = activation_criteria_;
        number_variables["LIMITER"] = limiter_;
        number_variables["DEACTIVATION"] = deactivate_;
#ifdef CSMP_WITH_SAMG_SOLVER
        number_variables["SAMG_RECONSTRUCTION"] = samg_reconstruction_;
#endif
        switch (relperm_model_)
        {
            case LINEARMODEL:
                text_variables["RELPERM_MODEL"] = "LINEARMODEL";
                break;
            
            case BROOKSCOREY:
                text_variables["RELPERM_MODEL"] = "BROOKSCOREY";
                break;

            case VANGENUCHTEN:
                text_variables["RELPERM_MODEL"] = "VANGENUCHTEN";
                break;

            case EXPERIMENTAL:
                text_variables["RELPERM_MODEL"] = "EXPERIMENTAL";
                break;
                
            default:
                text_variables["RELPERM_MODEL"] = "Unknown!";
        }

        switch (formulation_)
        {
            case GLOBAL:
                text_variables["FORMULATION"] = "GLOBAL";
                break;
            
            case PHASE:
                text_variables["FORMULATION"] = "PHASE";
                break;
                
            default:
                text_variables["FORMULATION"] = "Unknown!";
        }

        switch (solution_)
        {
        case IMPES:
            text_variables["SOLUTION"] = "IMPES";
            break;

        case QEIMPES:
            text_variables["SOLUTION"] = "QEIMPES";
            break;

        case QFIMPES:
            text_variables["SOLUTION"] = "QFIMPES";
            break;

        default:
            text_variables["SOLUTION"] = "Unknown!";
        }

        std::string temp;
        long position;

        std::cout << "\n";
        
        // reading file and putting values in the maps
        while (!numerical_file.eof())
        {
            position = numerical_file.tellg();
            
            getline(numerical_file, temp);
            if (temp[0] == '#')
            {
                continue;
            }
            else
            {
                numerical_file.seekg(position);
                
                numerical_file >> temp;
                
                if (number_variables.find(temp) != number_variables.end())
                {
                    std::cout << temp << "\t";
                    double value;
                    
                    numerical_file >> value;
                    
                    number_variables[temp] = value;
                    
                    std::cout << value << "\n";
                    
                    // to read the whole line and go to the next
                    numerical_file.seekg(position);
                    getline(numerical_file, temp);
                    
                    
                }
                else if (text_variables.find(temp) != text_variables.end())
                {
                    std::cout << temp << "\t";
                    std::string value;
                    
                    numerical_file >> value;
                    
                    text_variables[temp] = value;
                    
                    std::cout << value << "\n";

                    // to read the whole line and go to the next
                    numerical_file.seekg(position);
                    getline(numerical_file, temp);
                    
                
                }
                else if (temp == "OUTPUT_ONCE")
                {
                    output_once_keyword = true;
                    
                    numerical_file.seekg(position);
                    getline(numerical_file, temp);
                    
                    std::cout << temp << "\t";
                    
                    do
                    {
                        position = numerical_file.tellg();
                        getline(numerical_file, temp);
                        
                        if (temp[0] == '#')
                        {
                            continue;
                        }
                        else
                        {
                            if (all_keywords.find(temp) == all_keywords.end())
                            {
                                output_once.insert(temp);
                            }
                            
                        }
                       
                        
                    } while (!numerical_file.eof() && (all_keywords.find(temp) == all_keywords.end()));
                    
                    for (std::set<std::string>::iterator sit = output_once.begin(); sit != output_once.end(); sit++)
                    {
                        std::cout << (*sit) << "\t";
                    }
                    std::cout << "\n";
                    
                    if (all_keywords.find(temp) != all_keywords.end())
                    {
                        numerical_file.seekg(position);
                        continue;
                    }                    
                
                }
                else if (temp == "OUTPUT_ALWAYS")
                {
                    output_always_keyword = true;
                    
                    numerical_file.seekg(position);
                    getline(numerical_file, temp);
                    
                    std::cout << temp << "\t";
                    
                    do
                    {
                        position = numerical_file.tellg();
                        getline(numerical_file, temp);
                        
                        if ( temp.empty() || temp[0] == '#')
                        {
                            continue;
                        }
                        else
                        {
                            if (all_keywords.find(temp) == all_keywords.end())
                            {
                                output_always.insert(temp);
                            }
                            
                        }
                       
                        
                    } while (!numerical_file.eof() && (all_keywords.find(temp) == all_keywords.end()));
                    
                    for (std::set<std::string>::iterator sit = output_always.begin(); sit != output_always.end(); sit++)
                    {
                        std::cout << (*sit) << "\t";
                    }
                    std::cout << "\n";
                    
                    if (all_keywords.find(temp) != all_keywords.end())
                    {
                        numerical_file.seekg(position);
                        continue;
                    }

                }
                else
                {
                    numerical_file.seekg(position);
                    getline(numerical_file, temp);
                    std::cout << "Unrecognized keyword! Skip the line...  " << temp << "\n";
                }
                
            }
          
        }
        
        // assigning the values to the variables
        with_divergence_correction_ = (bool) number_variables["DIVERGENCE_CORRECTION"];
        normal_saturation_changes_ = number_variables["NORM_SAT_CHANGE"];
        maximum_saturation_changes_ = number_variables["MAX_SAT_CHANGE"];
        CFL_multiplier_ = number_variables["CFL_MULTIPLIER"];
        maximum_number_transport_timesteps_ = (size_t) number_variables["MAX_TRANSPORT_STEPS"];
        maximum_transport_timestep_size_ = number_variables["MAX_TRANSPORT_TIME_SIZE"];
        activation_criteria_ = number_variables["ACTIVATION_CRITERIA"];
        limiter_ = number_variables["LIMITER"];
        deactivate_ = (size_t) number_variables["DEACTIVATION"];
#ifdef CSMP_WITH_SAMG_SOLVER
        samg_reconstruction_ = (size_t) number_variables["SAMG_RECONSTRUCTION"];
#endif
        if (text_variables["RELPERM_MODEL"] == "LINEARMODEL")
        {
            relperm_model_ = LINEARMODEL;
        }
        else if (text_variables["RELPERM_MODEL"] == "BROOKSCOREY")
        {
            relperm_model_ = BROOKSCOREY;
        }
        else if (text_variables["RELPERM_MODEL"] ==  "VANGENUCHTEN")
        {
            relperm_model_ = VANGENUCHTEN;
        }
        else if (text_variables["RELPERM_MODEL"] ==  "EXPERIMENTAL")
        {
            relperm_model_ = EXPERIMENTAL;
        }
        // SKM_FIX: new rate dependent heterogeneity aware relative permeability model
        else if (text_variables["RELPERM_MODEL"] ==  "HETEROGENEITY_AWARE")
        {
            relperm_model_ = HETEROGENEITY_AWARE;
        }

        else
        {
            std::cout << "Unknown RELPERM_MODEL. It will be set to BROOKSCOREY\n";
            relperm_model_ = BROOKSCOREY;
        }

        if (text_variables["FORMULATION"] == "GLOBAL")
        {
            formulation_ = GLOBAL;
        }
        else if (text_variables["FORMULATION"] == "PHASE")
        {
            formulation_ = PHASE;
        }
        else
        {
            std::cout << "Unknown FORMULATION. It will be set to GLOBAL\n";
            formulation_ = GLOBAL;
        }

        if (text_variables["SOLUTION"] == "IMPES")
        {
            solution_ = IMPES;
        }
        else if (text_variables["SOLUTION"] == "QEIMPES")
        {
            solution_ = QEIMPES;
        }
        else if (text_variables["SOLUTION"] ==  "QFIMPES")
        {
            solution_ = QFIMPES;
        }
        else
        {
            std::cout << "Unknown SOLUTION. It will be set to IMPES\n";
            solution_ = IMPES;
        }

        if (output_once_keyword)
        {
            output_once_.clear();
            
            for (std::set<std::string>::iterator sit = output_once.begin(); sit != output_once.end(); sit++)
            {
                if (database.IsDefined((*sit).c_str()))
                {
                    output_once_.insert((*sit));
                }
                else
                {
                    std::cout << "\nProperty '" << (*sit) << "' is not defined in the database. It will be discarded.\n";
                }
                
            }
            
            
        }

        if (output_always_keyword)
        {
            output_always_.clear();
            
            for (std::set<std::string>::iterator sit = output_always.begin(); sit != output_always.end(); sit++)
            {
                if (database.IsDefined((*sit).c_str()))
                {
                    output_always_.insert((*sit));
                }
                else
                {
                    std::cout << "\nProperty '" << (*sit) << "' is not defined in the database. It will be discarded.\n";
                }
                
            }
            
            
        }

        
    
    }
    
    numerical_file.close();
    
 }



 template<uint32_t dim>
 void IMPES_Setup<dim>::ReportSetting()
 {
    std::cout << "\nIMPES_Setup::ReportSetting()\n";
    
    std::cout << "\n*******************************************\n";
    std::cout << "           IMPES Settings Report\n";
    std::cout << "*******************************************\n\n";
    
    std::cout << "Gravitational forces:         ";   
    if (with_gravitational_forces_)
        std::cout << "ON\n";
    else
        std::cout << "OFF\n";
        
    std::cout << "Capillary forces:             ";
    if (with_capillary_forces_)
        std::cout << "ON\n";
    else
        std::cout << "OFF\n";
        
    std::cout << "DIVERGENCE_CORRECTION:        ";
    if (with_divergence_correction_)
        std::cout << "YES\n";
    else
        std::cout << "NO\n";
        
    std::cout << "FORMULATION                   ";
    switch (formulation_)
    {
        case GLOBAL:
            std::cout << "GLOBAL\n";
            break;
        
        case PHASE:
            std::cout << "PHASE\n";
            break;
            
        default:
            std::cout << "Unknown!\n";
    }

    std::cout << "RELPERM_MODEL:                ";
    switch (relperm_model_)
    {
        case LINEARMODEL:
            std::cout << "LINEARMODEL\n";
            break;
        
        case BROOKSCOREY:
            std::cout << "BROOKSCOREY\n";
            break;

        case VANGENUCHTEN:
            std::cout << "VANGENUCHTEN\n";
            break;

        case EXPERIMENTAL:
            std::cout << "EXPERIMENTAL\n";
            break;
            
        case HETEROGENEITY_AWARE:
            std::cout << "HETEROGENEITY_AWARE\n";
            break;
        
        default:
            std::cout << "Unknown!\n";
    }

    std::cout << "SOLUTION:                     ";
    switch (solution_)
    {
    case IMPES:
        std::cout << "IMPES\n";
        break;

    case QEIMPES:
        std::cout << "QEIMPES\n";
        break;

    case QFIMPES:
        std::cout << "QFIMPES\n";
        break;

    default:
        std::cout << "Unknown!\n";
    }

    std::cout << "NORM_SAT_CHANGE:              " << normal_saturation_changes_ << "\n";    
    std::cout << "MAX_SAT_CHANGE:               " << maximum_saturation_changes_ << "\n";
    std::cout << "CFL_MULTIPLIER:               " << CFL_multiplier_ << "\n";
    std::cout << "MAX_TRANSPORT_STEPS:          " << maximum_number_transport_timesteps_ << "\n";
    std::cout << "MAX_TRANSPORT_TIME_SIZE:      " << maximum_transport_timestep_size_ << "\n";
#ifdef CSMP_WITH_SAMG_SOLVER
    std::cout << "SAMG_RECONSTRUCTION:          " << samg_reconstruction_ << "\n";
#endif
    std::cout << "ACTIVATION_CRITERIA:          " << activation_criteria_ << "\n";
    std::cout << "LIMITER:                      " << limiter_ << "\n";
    std::cout << "DEACTIVATE:                   " << deactivate_ << "\n";
    
    std::cout << "OUTPUT_ONCE:                  \n";
    for (std::set<std::string>::const_iterator sit = output_once_.begin(); sit != output_once_.end(); sit++)
    {
        std::cout << "                              " << (*sit) << "\n";
    }
    std::cout << "\n";

    std::cout << "OUTPUT_ALWAYS:                \n";
    for (std::set<std::string>::const_iterator sit = output_always_.begin(); sit != output_always_.end(); sit++)
    {
        std::cout << "                              " << (*sit) << "\n";
    }
    std::cout << "\n";
    
 }
 
 
 template class IMPES_Setup<2>;
} // end csmp
