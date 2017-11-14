//
//  GaussianBeam.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__GaussianBeam__
#define __cppxfel__GaussianBeam__
#include "Beam.h"

#include <stdio.h>

class GaussianBeam : public Beam
{
private:
    double optimisingWavelength;
    double optimisingBandwidth;
    double optimisingExponent;
    double stepSizeWavelength;
    double stepSizeBandwidth;
    double stepSizeExponent;
    double toleranceWavelength;
    double toleranceBandwidth;
    double toleranceExponent;

    std::vector<double> meanWavelengths;
    std::vector<double> bandwidths;
    std::vector<double> exponents;
    std::vector<double> heights;

public:
    GaussianBeam(double meanWavelength, double bandwidth, double exponent);
    double integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength);
    void addParameters(GetterSetterMapPtr map);
    double getNominalWavelength();
    double valueAtWavelength(double aWavelength);
    bool nonZeroPartialityExpected(double lowWavelength, double highWavelength);

    static void setWavelength(void *object, double aWavelength, int which)
    {
        static_cast<GaussianBeam *>(object)->meanWavelengths[which] = aWavelength;
    }

    static double getWavelength(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->meanWavelengths[which];
    }

    static void setBandwidth(void *object, double aBandwidth, int which)
    {
        static_cast<GaussianBeam *>(object)->bandwidths[which] = aBandwidth;
    }

    static double getBandwidth(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->bandwidths[which];
    }

    static void setExponent(void *object, double anExponent, int which)
    {
        static_cast<GaussianBeam *>(object)->exponents[which] = anExponent;
    }

    static double getExponent(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->exponents[which];
    }

    static void setHeight(void *object, double aHeight, int which)
    {
        static_cast<GaussianBeam *>(object)->heights[which] = aHeight;
    }

    static double getHeight(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->heights[which];
    }

    // individual parameter getters

    static double getWavelengthA(void *object)
    {
        return getWavelength(object, 0);
    }

    static double getWavelengthB(void *object)
    {
        return getWavelength(object, 1);
    }

    static double getBandwidthA(void *object)
    {
        return getBandwidth(object, 0);
    }

    static double getBandwidthB(void *object)
    {
        return getBandwidth(object, 1);
    }

    static double getExponentA(void *object)
    {
        return getExponent(object, 0);
    }

    static double getExponentB(void *object)
    {
        return getExponent(object, 1);
    }

    static double getHeightA(void *object)
    {
        return getHeight(object, 0);
    }

    static double getHeightB(void *object)
    {
        return getHeight(object, 1);
    }

    // individual parameter setters

    static void setWavelengthA(void *object, double aWavelength)
    {
        return setWavelength(object, aWavelength, 0);
    }

    static void setWavelengthB(void *object, double bWavelength)
    {
        return setWavelength(object, bWavelength, 1);
    }

    static void setBandwidthA(void *object, double aBandwidth)
    {
        return setBandwidth(object, aBandwidth, 0);
    }

    static void setBandwidthB(void *object, double aBandwidth)
    {
        return setBandwidth(object, aBandwidth, 1);
    }

    static void setExponentA(void *object, double aExponent)
    {
        return setExponent(object, aExponent, 0);
    }

    static void setExponentB(void *object, double aExponent)
    {
        return setExponent(object, aExponent, 1);
    }

    static void setHeightA(void *object, double aHeight)
    {
        return setHeight(object, aHeight, 0);
    }

    static void setHeightB(void *object, double aHeight)
    {
        return setHeight(object, aHeight, 1);
    }

};

#endif /* defined(__cppxfel__GaussianBeam__) */
