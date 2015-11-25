//
//  LoggableObject.h
//  cppxfel
//
//  Created by Helen Ginn on 12/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__LoggableObject__
#define __cppxfel__LoggableObject__

#include <stdio.h>
#include "Logger.h"

class LoggableObject
{
protected:
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);
};

#endif /* defined(__cppxfel__LoggableObject__) */
