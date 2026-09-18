// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IMPES_Setup.h"
#include "Generic2P2D_IMPES_Simulator.h"
#include "ErrorHandler.h"
#include "ModelTime.h"
#include "GlobalVerbose.h"
#include "SignalHandler.h"
// 2019 stuff
#include "LayeredCompositeProcessor1.h"

using namespace std;
using namespace csmp;

/**
     SKM Trial bed for new functionality 27/9/2019.
 
     1. Reading the rocktypes from file, including the table for the upscaling of variables.
 
     2. Computation rate-dependent relperm
 
     3. Population of and use of rate-dependent relperms in look-up table
 
     X. Use of relperms in Shao's simular
 
     X. Speed-up of Shaho's simulator
 
*/
int main( int argc, char* argv[] )
{
    bool cmdLineModelName(argc != 1);
    bool gravity(false);
    bool capillary(false);
    bool restart(false);
    int choice(9999);

    ModelTime::Instance().modelTime = 0.;
    
    signal(SIGINT, &signalHandler);
    
    cout <<"\n";
    cout <<"*****************************************************************************\n";
    cout <<"***   Element Based IMPES Reservoir Simulator 2D (version 2011.10 Beta)   ***\n";
    cout <<"*****************************************************************************\n";
    string model_name("undefined");
    if( !cmdLineModelName )
    {
      cout << "\nEnter name of simulation model: ";
      cin >> model_name;
    }
    else
    {
      const char* cache(argv[1]);
      model_name = cache;
      choice = atoi(argv[2]);
      
      cout << "\nName of simulation model: " << model_name << "\n";      
    }
    
    
    // checking if restart data exist
    string restart_geo_file_name(model_name);
    restart_geo_file_name += ".rsg";
    string restart_dat_file_name(model_name);
    restart_dat_file_name += ".rsd";
    
    ifstream restart_g_file;
    ifstream restart_d_file;
    
    restart_g_file.open(restart_geo_file_name.c_str());
    restart_d_file.open(restart_dat_file_name.c_str());
    
    if (restart_g_file && restart_d_file)
    {
        double restart_time;
        
        restart_d_file >> restart_time;
        
        cout << "\nRestart data found at time " << restart_time << 
                 " sec. \nDo you want to continue simulation from restart data (y/n)? ";
        
        char choice;
        cin.clear();
        cin >> choice;
        switch (choice)
        {
            case 'y':
            case 'Y':
                restart = true;
                cout << "\nSimulation runs from time = " << restart_time << " sec\n";
                break;
            default:
                cout << "\nSimulation runs from time = 0.\n";
        }
        
    }
    
    restart_g_file.close();
    restart_d_file.close();


    // check whether one wants to use csmp-binary file
    bool ansys_model_true_or_csmp_binary_false( true );
    char csmp_binary('n');
    cout << "\nDo you want to load ANSYS model from files ( "<< model_name<< ".asc, "<< model_name<< ".dat ) (y/n)? ";
    cout << "\nIn case if you choose 'no' you should have mesh files in csmp native binary format:";
    cin.clear();
    cin >> csmp_binary;
    switch (csmp_binary)
    {
        case 'n':
        case 'N':
            ansys_model_true_or_csmp_binary_false = false;
          break;
    }

    cout << "\n";
    cout << "Choose type of your simulation:\n\n";
    cout << "  (0) Effective permeability calculation\n";
    cout << "  (1) Viscous flow\n";
    cout << "  (2) Viscous flow with gravitational force\n";
    cout << "  (3) Viscous flow with capillary force\n";
    cout << "  (4) Viscous flow with gravitational and capillary forces\n\n";
    cout << "Your choice: ";
    
    if (choice == 9999)
    {
        cin >> choice;
    }
    else
    {
        cout << choice << "\n";
    }
    
    switch (choice)
    {
        case 0:
            break;

        case 1:
            break;
        
        case 2:
            gravity = true;
            break;
        
        case 3:
            capillary = true;
            break;
        
        case 4:
            gravity = true;
            capillary = true;
            break;
            
        default:
            cout << "\nUnrecognized option! Program terminates.";
            return 0;
        
    }

    try{

    IMPES_Setup<2> setup(gravity, capillary);
    setup.WriteVariableFile("IMPES-variables.txt");
    Generic2P2D_IMPES_Simulator simulator( model_name.c_str(),
                                           ansys_model_true_or_csmp_binary_false,
                                           setup );
    if (choice == 0U)
    {
        simulator.CalculateEffectivePermeability();
    } 
    else
    {
        simulator.Run(restart);
    }

    } 
    catch( csmp::Exception& exc )
    {
      cout << exc.What();
    }    
    catch( std::exception& exc )
    {
      cout << exc.what();
    }

    cout << "\nSimulation is done ...\n";
    return 0;

}

