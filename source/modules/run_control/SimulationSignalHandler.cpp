// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <iostream>
#include <string>
#include "SimulationSignalHandler.h"

using namespace std;

bool SimulationSignalHandler::outputVTUSignal_(false);
bool SimulationSignalHandler::quitSignal_(false);
bool SimulationSignalHandler::interactiveSignalRaised_(false);
bool SimulationSignalHandler::outputMonitorDataSignal_(false);
bool SimulationSignalHandler::outputRestartSignal_(false);

SimulationSignalHandler::SimulationSignalHandler()
{
}


SimulationSignalHandler::~SimulationSignalHandler()
{
}


void HandleSignal(int sig_code)
{
    SimulationSignalHandler sig;
    //    if (!sig.SignalRaised())
    //    {
    //        sig.SignalRaised(true);
    //    }
    //    else
    //    {
    //  sig.SignalRaised(false);

    if (sig_code == SIGINT)
    {
        int choice(0);
        while (choice < 1 || choice > 10){
            std::cout << "\n\nSimulationSignalHandler::CTRL+C signal raised. Signal handling options:\n\n";
            std::cout << "1. Continue\n";
            std::cout << "2. Write all files, and continue\n";
            std::cout << "3. Output Monitoring Data and continue\n";
            std::cout << "4. Output restart file and continue\n";
            std::cout << "5. Output to VTU and continue\n";

            std::cout << "6. Output Monitoring Data and quit\n";
            std::cout << "7. Output restart file and quit\n";
            std::cout << "8. Output to VTU and quit\n";

            std::cout << "9. Write all files, and quit\n";
            std::cout << "10. Quit\n\n";
            std::cout << "Your choice: ";
            sig.OutputMonitorDataSignal(false);
            sig.OutputRestartFileSignal(false);
            sig.OutputVTUSignal(false);
            sig.QuitSignal(false);
            std::string choice1;

            std::cin.clear();
            fflush(stdin);
            getline(cin,choice1);
            choice=std::atoi(choice1.c_str());
            switch (choice){
            case 1: {
                break;
            }case 2: {  // all + cont
                sig.OutputMonitorDataSignal(true);
                sig.OutputRestartFileSignal(true);
                sig.OutputVTUSignal(true);
                break;
            }case 3: {  // mon + cont
                sig.OutputMonitorDataSignal(true);
                break;
            }case 4: {  // restart + cont
                sig.OutputRestartFileSignal(true);
                break;
            }case 5: {  // vtu + cont
                sig.OutputVTUSignal(true);
                break;
            }case 6: {  // mon + quit
                sig.OutputMonitorDataSignal(true);
                sig.QuitSignal(true);
                break;
            }case 7: {  // restart + quit
                sig.OutputRestartFileSignal(true);
                sig.QuitSignal(true);
                break;
            }case 8: {  // vtu + quit
                sig.OutputVTUSignal(true);
                sig.QuitSignal(true);
                break;
            }case 9: {  // all + quit
                sig.OutputMonitorDataSignal(true);
                sig.OutputRestartFileSignal(true);
                sig.OutputVTUSignal(true);
                sig.QuitSignal(true);
                break;
            }case 10:{ // quit
                sig.QuitSignal(true);
                break;
            }default:{
                std::cout << "Unrecognized option! Please try again.\n";
                break;
            }

            }
        }
        cout<<" Output signal State:"<<endl;
        cout<<"    outputVTUSignal_         "<<sig.OutputVTUSignal()<<endl;
        cout<<"    quitSignal_              "<<sig.QuitSignal()<<endl;
        cout<<"    outputMonitorDataSignal_ "<<sig.OutputMonitorDataSignal()<<endl;
        cout<<"    outputRestartSignal_     "<<sig.OutputRestartFileSignal()<<endl;
    }
}

//} //namespace csmp{
