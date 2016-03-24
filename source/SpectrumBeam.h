//
//  SpectrumBeam.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SpectrumBeam__
#define __cppxfel__SpectrumBeam__

#include <stdio.h>
#include "GaussianBeam.h"

class SpectrumBeam : public Beam
{
private:
    std::map<double, double> spectrum;
    double wavelength;
    double interpolateBetweenIterators(std::map<double, double>::iterator itLow, std::map<double, double>::iterator itHigh, double wavelength);
    
public:
    SpectrumBeam(GaussianBeamPtr gaussianStart, MtzPtr mtzStart);
    bool nonZeroPartialityExpected(double lowWavelength, double highWavelength);
    
    double integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength);
    double getNominalWavelength()
    {
        return wavelength;
    }
    
    void addParameters(GetterSetterMapPtr map) {};
};

#endif /* defined(__cppxfel__SpectrumBeam__) */
