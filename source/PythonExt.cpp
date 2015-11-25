//
//  PythonExt.cpp
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "PythonExt.h"
#include "InputFileParser.h"
#include "Logger.h"
#include "MtzRefiner.h"
#include "FileReader.h"

int new_main(int argc, char *argv[]);

void setupCppxfel()
{
    time_t startcputime;
    time(&startcputime);
    
    Logger::mainLogger = LoggerPtr(new Logger());
    boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);
   
    std::cout << "Welcome to Helen's XFEL tasks" << std::endl;
}

void runScriptFromPython(std::string scriptName)
{
    setupCppxfel();
    
    InputFileParser *parser = new InputFileParser(scriptName);
    
    parser->parse(true);
    
    delete parser;
}

void runCommandLine(std::string fullArgs)
{
    std::vector<std::string> strings = FileReader::split(fullArgs, ' ');
		
		runCommandLineArgs(strings.size(), strings);
}

void runCommandLineArgs(int argc, std::vector<std::string> stringArgv)
{
    std::cout << "Running cppxfel..." << std::endl;

    vector<char *> charVector;
    
    for (int i = 0; i < stringArgv.size(); i++)
    {
        charVector.push_back((char *)stringArgv[i].c_str());
    }
    
    new_main(argc, &(*(charVector.begin())));
}
