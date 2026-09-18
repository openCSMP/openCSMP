// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SIMULATION_SIGNAL_HANDLER_H
#define SIMULATION_SIGNAL_HANDLER_H

#include <csignal>
#include <iostream>

void HandleSignal(int sig_code);

/**
    Allows user to interrupt simulation by keystrokes.
    
    Used by FEMCFVM simulator.
    
    @author J. Mindel & S. Bazrafkan
    @date 2012
*/
class SimulationSignalHandler
{
public:
    SimulationSignalHandler();
    ~SimulationSignalHandler();


    void OutputVTUSignal(bool signal){outputVTUSignal_ = signal;}
    bool OutputVTUSignal(){return outputVTUSignal_;}

    // used for monitoring data output
    void OutputMonitorDataSignal(bool signal){outputMonitorDataSignal_= signal;}
    bool OutputMonitorDataSignal(){return outputMonitorDataSignal_;}

    // output restart file.
    void OutputRestartFileSignal(bool signal){outputRestartSignal_ = signal;}
    bool OutputRestartFileSignal(){return outputRestartSignal_;}

    void QuitSignal(bool signal){this->quitSignal_ = signal;}
    bool QuitSignal(){return this->quitSignal_;}

    void SignalRaised(bool signal){interactiveSignalRaised_ = signal;}
    bool SignalRaised(){return interactiveSignalRaised_;}
    //void HandleSignal(int sig_code);


private:

    static bool outputVTUSignal_;
    static bool outputRestartSignal_;
    static bool quitSignal_;
    static bool interactiveSignalRaised_;
    static bool outputMonitorDataSignal_;
};

#endif // SIMULATION_SIGNAL_HANDLER_H
