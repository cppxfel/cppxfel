//
//  GaussianBeam.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "GaussianBeam.h"
#include "Vector.h"
#include <math.h>
#include "GetterSetterMap.h"
#include "FileParser.h"

GaussianBeam::GaussianBeam(double aMeanWavelength, double aBandwidth, double anExponent)
{
    meanWavelengths.push_back(aMeanWavelength);
    bandwidths.push_back(aBandwidth);
    exponents.push_back(anExponent);
    heights.push_back(1);

    optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH",
                                              !OPTIMISED_WAVELENGTH);
    optimisingBandwidth = FileParser::getKey("OPTIMISING_BANDWIDTH",
                                             !OPTIMISED_BANDWIDTH);

    optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT",
                                            !OPTIMISED_EXPONENT);

    stepSizeWavelength = FileParser::getKey("STEP_SIZE_WAVELENGTH",
                                            MEAN_STEP);
    stepSizeBandwidth = FileParser::getKey("STEP_SIZE_BANDWIDTH",
                                           BANDWIDTH_STEP);

    stepSizeExponent = FileParser::getKey("STEP_SIZE_EXPONENT",
                                          EXPONENT_STEP);

    toleranceWavelength = FileParser::getKey("TOLERANCE_WAVELENGTH",
                                             MEAN_TOLERANCE);
    toleranceBandwidth = FileParser::getKey("TOLERANCE_BANDWIDTH",
                                            BANDWIDTH_TOLERANCE);

    toleranceExponent = FileParser::getKey("TOLERANCE_EXPONENT",
                                           EXPO_TOLERANCE);
}

double GaussianBeam::integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength)
{
    double pValue = valueAtWavelength(lowWavelength);
    double qValue = valueAtWavelength(highWavelength);

/*    for (int i = 0; i < meanWavelengths.size(); i++)
    {
        pValue += super_gaussian(lowWavelength, meanWavelengths[i], bandwidths[i], exponents[i]);
        qValue += super_gaussian(highWavelength, meanWavelengths[i], bandwidths[i], exponents[i]);
    }
    */
 //   double pX = (lowWavelength - meanWavelength) / bandwidth;
 //   double qX = (highWavelength - meanWavelength) / bandwidth;

    double width = fabs(highWavelength - lowWavelength);

    double area = (pValue + qValue) / 2 * width;

    return area;
}

double GaussianBeam::valueAtWavelength(double aWavelength)
{
    double value = 0;

    for (int i = 0; i < meanWavelengths.size(); i++)
    {
        value += super_gaussian(aWavelength, meanWavelengths[i], bandwidths[i], exponents[i]);
    }

    return value;
}

double GaussianBeam::getNominalWavelength()
{
    double sum = 0;
    double weights = 0;

    for (int i = 0; i < meanWavelengths.size(); i++)
    {
        sum += meanWavelengths[i] * heights[i];
        weights += heights[i];
    }

    sum /= weights;

    return sum;
}

void GaussianBeam::addParameters(GetterSetterMapPtr map)
{
    if (optimisingWavelength)
    {
        map->addParameter(this, getWavelengthA, setWavelengthA, stepSizeWavelength, toleranceWavelength);

        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getWavelengthB, setWavelengthB, stepSizeWavelength, toleranceWavelength);
        }
    }

    if (optimisingBandwidth)
    {
        map->addParameter(this, getBandwidthA, setBandwidthA, stepSizeBandwidth, toleranceBandwidth);

        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getBandwidthB, setBandwidthB, stepSizeBandwidth, toleranceBandwidth);
        }
    }

    if (optimisingExponent)
    {
        map->addParameter(this, getExponentB, setExponentB, stepSizeExponent, toleranceExponent);

        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getExponentB, setExponentB, stepSizeExponent, toleranceExponent);
        }
    }

    if (meanWavelengths.size() >= 2)
    {
        map->addParameter(this, getHeightA, setHeightA, stepSizeExponent, toleranceExponent);
        map->addParameter(this, getHeightB, setHeightB, stepSizeExponent, toleranceExponent);
    }
}

bool GaussianBeam::nonZeroPartialityExpected(double lowWavelength, double highWavelength)
{
    if (highWavelength < getNominalWavelength() - bandwidths[0] * 3)
        return true;
    if (lowWavelength > getNominalWavelength() + bandwidths[0] * 3)
        return true;

    return false;
}
