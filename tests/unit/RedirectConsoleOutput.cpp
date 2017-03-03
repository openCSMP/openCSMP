//
//  RedirectConsoleOutput.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 8/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "RedirectConsoleOutput.hpp"


// actually define vars.
int StdCapture::m_pipe[2];
int StdCapture::m_oldStdOut;
int StdCapture::m_oldStdErr;
bool StdCapture::m_capturing;
std::mutex StdCapture::m_mutex;
std::string StdCapture::m_captured;
