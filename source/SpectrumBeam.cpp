//
//  SpectrumBeam.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpectrumBeam.h"
#include "MtzManager.h"
#include "Holder.h"
#include "Miller.h"

SpectrumBeam::SpectrumBeam(GaussianBeamPtr gaussianStart, MtzPtr mtzStart)
{
   // double originalBandwidth = GaussianBeam::getBandwidthA(&*gaussianStart);
    double originalWavelength = GaussianBeam::getWavelengthA(&*gaussianStart);

    double lowRes = 1.9;
    double highRes = 1.6;

    for (int i = 0; i < mtzStart->reflectionCount(); i++)
    {
        Reflection *single = mtzStart->reflection(i);

        if (!single->acceptedCount())
            continue;

        if (!single->betweenResolutions(lowRes, highRes))
            continue;

        MtzManager *reference = MtzManager::getReferenceManager();
        Reflection *holder;

        reference->findReflectionWithId(single->getReflId(), &holder);

        if (!holder)
            continue;

        double merged = holder->meanIntensity();
        double unmerged = single->acceptedMiller(0)->getRawIntensity();

        double ratio = unmerged / merged;
        double wavelength = single->acceptedMiller(0)->getWavelength();

        if (!(wavelength == wavelength && ratio == ratio))
            continue;

        double value = gaussianStart->valueAtWavelength(wavelength);

        if (ratio < 0)
            ratio = 0;

        spectrum[wavelength] = value;// * ratio;
    //    std::cout << wavelength << "," << spectrum[wavelength] << std::endl;
    }
    /*
    for (double i = originalWavelength - 3 * originalBandwidth; i <= originalWavelength + 3 * originalBandwidth; i += originalBandwidth / 6)
    {
        double value = gaussianStart->valueAtWavelength(i);

      //  std::cout << value << ", ";

        spectrum[i] = value;
    }
    */
 //   std::cout << std::endl;

    wavelength = originalWavelength;
}

double SpectrumBeam::interpolateBetweenIterators(std::map<double, double>::iterator itLow, std::map<double, double>::iterator itHigh, double wavelength)
{
    double lowWav = itLow->first;
    double lowValue = itLow->second;

    double highWav = itHigh->first;
    double highValue = itHigh->second;

    double proportion = (wavelength - lowWav) / (highWav - lowWav);
    double value = proportion * (highValue - lowValue) + lowValue;

    return value;
}

bool SpectrumBeam::nonZeroPartialityExpected(double lowWavelength, double highWavelength)
{
    std::map<double, double>::iterator it = spectrum.begin();

    if (highWavelength < it->first)
    {
        return true;
    }

    it = spectrum.end();
    it--;

    if (lowWavelength > it->first)
    {
        return true;
    }

    return false;
}

double SpectrumBeam::integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength)
{
    bool foundLowWavelength = false;
    bool foundHighWavelength = false;
    double summation = 0;

    if (lowWavelength == highWavelength)
        return 0;

    std::map<double, double>::iterator it = spectrum.end();
    it--;

    if (highWavelength > it->first)
        return 0;

    for (std::map<double, double>::iterator it = spectrum.begin(); it != spectrum.end(); it++)
    {
        std::map<double, double>::iterator it2 = it;
        it2++;

        if (it2 == spectrum.end())
            break;

        double bottomWavelength = 0;
        double topWavelength = 0;
        double bottomValue = 0;
        double topValue = 0;

        if (lowWavelength > it->first && lowWavelength < it2->first)
        {
            bottomWavelength = lowWavelength;
        }

        if (lowWavelength < it->first)
        {
            bottomWavelength = it->first;
        }

        if (highWavelength > it->first && highWavelength < it2->first)
        {
            topWavelength = highWavelength;
        }

        if (highWavelength > it2->first)
        {
            topWavelength = it2->first;
        }

        if (bottomWavelength != 0 && topWavelength != 0)
        {
            bottomValue = interpolateBetweenIterators(it, it2, bottomWavelength);
            topValue = interpolateBetweenIterators(it, it2, topWavelength);

            assert(topWavelength > bottomWavelength);

            double averageValue = (topValue + bottomValue) / 2;
            double addition = averageValue * (topWavelength - bottomWavelength);

//            std::cout << averageValue << ", " << bottomWavelength << ", " << topWavelength << std::endl;

    //        std::cout << "Sum addition " << addition << " = " << summation << std::endl;

            summation += addition;
        }

        if (foundHighWavelength)
            return summation;
    }

    return summation;
}
