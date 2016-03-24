//
//  Beam.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Beam__
#define __cppxfel__Beam__

#include <stdio.h>
#include "parameters.h"

class Beam
{
private:
    
public:
    Beam();
    
    virtual double integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength) { return 0;};
    virtual double getNominalWavelength() { return 0;};
    
    virtual void addParameters(GetterSetterMapPtr map) {};
    virtual bool nonZeroPartialityExpected(double lowWavelength, double highWavelength) { return false;};
};

#endif /* defined(__cppxfel__Beam__) */
