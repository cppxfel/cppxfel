//
//  XManager.h
//  cppxfel
//
//  Created by Helen Ginn on 03/02/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#ifndef __cppxfel__XManager__
#define __cppxfel__XManager__

#include "MtzManager.h"
#include <vector>
#include <string>

#include <stdio.h>

class XManager : public MtzManager
{
private:
    vector<std::string> filenames;
    vector<int> lineSplitters;
public:
    XManager();
    ~XManager();
    
    void setFilenames(vector<std::string> newFiles);
    
    virtual void loadReflections(PartialityModel model);
};

#endif /* defined(__cppxfel__XManager__) */
