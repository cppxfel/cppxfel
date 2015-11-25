//
//  PythonExt.h
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__PythonExt__
#define __cppxfel__PythonExt__

#include <string>
#include <vector>
#include "parameters.h"

class MtzRefiner;

void runScriptFromPython(std::string scriptName);
void runCommandLineArgs(int argc, vector<std::string> argv);
void runCommandLine(std::string fullArgs);

#endif /* defined(__cppxfel__PythonExt__) */
