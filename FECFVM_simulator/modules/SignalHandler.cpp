// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SignalHandler.h"
#include <string>
#include "stdio.h"

bool SignalHandler::outputSignal_(false);

bool SignalHandler::restartSignal_(false);

bool SignalHandler::quitSignal_(false);

bool SignalHandler::signalRaised_(false);


SignalHandler::SignalHandler()
{
}


SignalHandler::~SignalHandler()
{
}


void signalHandler(int sig_code)
{
    SignalHandler sig;
    
    if (!sig.SignalRaised())
    {
        sig.SignalRaised(true);
    }
    else
    {
        sig.SignalRaised(false);
        
        if (sig_code == SIGINT)
        {
            std::cout << "\n\nsignalHandler::CTRL+C signal raised. Signal handling options:\n\n";
            std::cout << "1. Continue\n";
            std::cout << "2. Output to vtk and continue\n";
            std::cout << "3. Output to vtk and quit\n";
            std::cout << "4. Write restart file and quit\n";
            std::cout << "5. Write restart file, output to vtk and quit\n";
            std::cout << "6. Quit\n\n";
            std::cout << "Your choice: ";
            
            char choice;
            
            std::cin.clear();
            fflush(stdin);
            
            std::cin >> choice;
            
            switch (choice)
            {
                case '1':
                {
                    break;
                }
                
                case '2':
                {
                    sig.OutputSignal(true);
                    break;
                }

                case '3':
                {
                    sig.OutputSignal(true);
                    sig.QuitSignal(true);
                    break;
                }


                case '4':
                {
                    sig.RestartSignal(true);
                    sig.QuitSignal(true);
                    break;
                }

                case '5':
                {
                    sig.RestartSignal(true);
                    sig.OutputSignal(true);
                    sig.QuitSignal(true);
                    break;
                }

                case '6':
                {
                    sig.QuitSignal(true);
                    break;
                }
                
                default:
                {
                  std::cout << "Unrecognized option! Program continues.\n";
                  std::string all_string;
                  getline(std::cin, all_string);
                  std::cin.clear();
                  fflush(stdin);
                  break;
                }
                
            }

    
        }
        
    }
    
    std::signal(sig_code, &signalHandler);

}

