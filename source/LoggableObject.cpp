//
//  LoggableObject.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "LoggableObject.h"

void LoggableObject::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}