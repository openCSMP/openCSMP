// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SIGNAL_HANDLER_H
#define SIGNAL_HANDLER_H

#include <csignal>
#include <iostream>

void signalHandler(int sig_code);

class SignalHandler
{
public:
    SignalHandler();
    ~SignalHandler();
    
    void OutputSignal(bool signal);
    bool OutputSignal();
    void RestartSignal(bool signal);
    bool RestartSignal();
    void QuitSignal(bool signal);
    bool QuitSignal();
    void SignalRaised(bool signal);
    bool SignalRaised();
    
    
private:

    static bool outputSignal_;
    static bool restartSignal_;
    static bool quitSignal_;
    static bool signalRaised_;


};


inline
void SignalHandler::OutputSignal(bool signal)
{
    outputSignal_ = signal;
}


inline
bool SignalHandler::OutputSignal()
{
    return outputSignal_;
}


inline
void SignalHandler::RestartSignal(bool signal)
{
    restartSignal_ = signal;
}


inline
bool SignalHandler::RestartSignal()
{
    return restartSignal_;
}


inline
void SignalHandler::QuitSignal(bool signal)
{
    quitSignal_ = signal;
}


inline
bool SignalHandler::QuitSignal()
{
    return quitSignal_;
}


inline
void SignalHandler::SignalRaised(bool signal)
{
    signalRaised_ = signal;
}


inline
bool SignalHandler::SignalRaised()
{
    return signalRaised_;
}



#endif