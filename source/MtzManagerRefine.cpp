/*
 * MtzManagerRefine.cpp
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>
#include "MtzManager.h"
#include "FileParser.h"



MtzManager *MtzManager::currentManager;

void MtzManager::applyUnrefinedPartiality()
{
    double wavelength = bestWavelength();
    double mosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
    double rlpSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
    double bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
    double exponent = FileParser::getKey("INITIAL_EXPONENT", INITIAL_EXPONENT);

    refreshPartialities(0, 0, 0, 0, 0, mosaicity, rlpSize, wavelength, bandwidth, exponent, cellDim[0], cellDim[1], cellDim[2]);
}

void MtzManager::setParams(double parameters[], int paramCount)
{
    wavelength = parameters[PARAM_WAVELENGTH];
	bandwidth = parameters[PARAM_BANDWIDTH];
	mosaicity = parameters[PARAM_MOS];
	spotSize = parameters[PARAM_SPOT_SIZE];
	hRot = parameters[PARAM_HROT];
	kRot = parameters[PARAM_KROT];
    aRot = parameters[PARAM_AROT];
    bRot = parameters[PARAM_BROT];
    cRot = parameters[PARAM_CROT];
    exponent = parameters[PARAM_EXPONENT];
    bFactor = parameters[PARAM_B_FACTOR];
    
    applyBFactor(bFactor);
    applyScaleFactor(parameters[PARAM_SCALE_FACTOR], 0, 0, true);
    
    cellDim[0] = parameters[PARAM_UNIT_CELL_A];
    cellDim[1] = parameters[PARAM_UNIT_CELL_B];
    cellDim[2] = parameters[PARAM_UNIT_CELL_C];
}

void MtzManager::getParams(double *parameters[], int paramCount)
{
	(*parameters)[PARAM_WAVELENGTH] = wavelength;
	(*parameters)[PARAM_BANDWIDTH] = bandwidth;
	(*parameters)[PARAM_MOS] = mosaicity;
	(*parameters)[PARAM_SPOT_SIZE] = spotSize;
	(*parameters)[PARAM_HROT] = hRot;
	(*parameters)[PARAM_KROT] = kRot;
    (*parameters)[PARAM_AROT] = aRot;
    (*parameters)[PARAM_BROT] = bRot;
    (*parameters)[PARAM_CROT] = cRot;
    (*parameters)[PARAM_EXPONENT] = exponent;
    (*parameters)[PARAM_B_FACTOR] = bFactor;
    (*parameters)[PARAM_SCALE_FACTOR] = scale;
    (*parameters)[PARAM_UNIT_CELL_A] = cellDim[0];
    (*parameters)[PARAM_UNIT_CELL_B] = cellDim[1];
    (*parameters)[PARAM_UNIT_CELL_C] = cellDim[2];
}

void MtzManager::getParamPointers(double ***parameters, int paramCount)
{
    (*parameters)[PARAM_WAVELENGTH] = &wavelength;
    (*parameters)[PARAM_BANDWIDTH] = &bandwidth;
    (*parameters)[PARAM_MOS] = &mosaicity;
    (*parameters)[PARAM_SPOT_SIZE] = &spotSize;
    (*parameters)[PARAM_HROT] = &hRot;
    (*parameters)[PARAM_KROT] = &kRot;
    (*parameters)[PARAM_AROT] = &aRot;
    (*parameters)[PARAM_BROT] = &bRot;
    (*parameters)[PARAM_CROT] = &cRot;
    (*parameters)[PARAM_EXPONENT] = &exponent;
    (*parameters)[PARAM_B_FACTOR] = &bFactor;
    (*parameters)[PARAM_SCALE_FACTOR] = &externalScale;
    (*parameters)[PARAM_UNIT_CELL_A] = &cellDim[0];
    (*parameters)[PARAM_UNIT_CELL_B] = &cellDim[1];
    (*parameters)[PARAM_UNIT_CELL_C] = &cellDim[2];
}

void MtzManager::getSteps(double *ranges[], int paramCount)
{
    (*ranges)[PARAM_HROT] = stepSizeOrientation;
    (*ranges)[PARAM_KROT] = stepSizeOrientation;
    (*ranges)[PARAM_AROT] = stepSizeOrientABC;
    (*ranges)[PARAM_BROT] = stepSizeOrientABC;
    (*ranges)[PARAM_CROT] = stepSizeOrientABC;
    (*ranges)[PARAM_MOS] = stepSizeMosaicity;
    (*ranges)[PARAM_SPOT_SIZE] = stepSizeRlpSize;
    (*ranges)[PARAM_WAVELENGTH] = stepSizeWavelength / 2;
    (*ranges)[PARAM_BANDWIDTH] = stepSizeBandwidth;
    (*ranges)[PARAM_B_FACTOR] = 0;
    (*ranges)[PARAM_SCALE_FACTOR] = 0;
    (*ranges)[PARAM_EXPONENT] = stepSizeExponent;
    (*ranges)[PARAM_UNIT_CELL_A] = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.5);
    (*ranges)[PARAM_UNIT_CELL_B] = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.5);
    (*ranges)[PARAM_UNIT_CELL_C] = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.5);
    
}

void MtzManager::refreshPartialities(double parameters[])
{
    refreshPartialities(parameters[PARAM_HROT],
                        parameters[PARAM_KROT],
                        parameters[PARAM_AROT],
                        parameters[PARAM_BROT],
                        parameters[PARAM_CROT],
                        parameters[PARAM_MOS],
                        parameters[PARAM_SPOT_SIZE],
                        parameters[PARAM_WAVELENGTH],
                        parameters[PARAM_BANDWIDTH],
                        parameters[PARAM_EXPONENT],
                        parameters[PARAM_UNIT_CELL_A],
                        parameters[PARAM_UNIT_CELL_B],
                        parameters[PARAM_UNIT_CELL_C]);
}

void MtzManager::refreshCurrentPartialities()
{
    if (externalScale != -1)
        this->applyScaleFactor(externalScale, 0, 0, true);
    
    applyBFactor(bFactor);
    refreshPartialities(this->hRot,
                        this->kRot,
                        this->aRot,
                        this->bRot,
                        this->cRot,
                        this->mosaicity,
                        this->spotSize,
                        this->wavelength,
                        this->bandwidth,
                        this->exponent,
                        this->cellDim[0],
                        this->cellDim[1],
                        this->cellDim[2]);
}

void MtzManager::refreshPartialities(double hRot, double kRot, double aRot, double bRot, double cRot, double mosaicity,
		double spotSize, double wavelength, double bandwidth, double exponent,
                                     double a, double b, double c)
{

    this->makeSuperGaussianLookupTable(exponent);

    if (!matrix)
        return;
    
    if (matrix->isComplex())
        this->matrix->changeOrientationMatrixDimensions(a, b, c, cellAngles[0], cellAngles[1], cellAngles[2]);
    
    MatrixPtr firstMatrix = MatrixPtr();
    MatrixPtr newMatrix = MatrixPtr();
    Miller::rotateMatrixABC(aRot, bRot, cRot, matrix, &firstMatrix);
    Miller::rotateMatrixHKL(hRot, kRot, 0, firstMatrix, &newMatrix);
    
	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			MillerPtr miller = reflections[i]->miller(j);
			miller->recalculatePartiality(newMatrix, mosaicity, spotSize,
					wavelength, bandwidth, exponent);
		}
	}
    
    logged << "Refreshed partialities with wavelength " << wavelength << ", spot size " << spotSize << " to generate " << accepted() << " accepted reflections." << std::endl;

    sendLog(LogLevelDebug);
}

static bool greaterThan(double num1, double num2)
{
    return (num1 > num2);
}

double MtzManager::medianWavelength(double lowRes, double highRes)
{
    vector<double> wavelengths;
    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>());
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            double wavelength = reflection(i)->miller(0)->getWavelength();
            
            if (wavelengthRange.size())
            {
                if (wavelength < wavelengthRange[0] || wavelength > wavelengthRange[1])
                    continue;
            }
            
            double isigi = reflection(i)->miller(j)->getRawestIntensity();
            
            if (isigi > refinementIntensityThreshold && isigi == isigi)
            {
                wavelengths.push_back(wavelength);
            }
        }
    }
    
    std::sort(wavelengths.begin(), wavelengths.end(), greaterThan);
    
    for (int i = 0; i < wavelengths.size(); i++)
    {
 //       std::cout << wavelengths[i] << std::endl;
    }
    
    if (wavelengths.size() < 2)
        return 1.29;
    
    return wavelengths[int(wavelengths.size() / 2)];
}

double MtzManager::bestWavelength(double lowRes, double highRes, bool usingReference)
{
    bool median = FileParser::getKey("MEDIAN_WAVELENGTH", false);
    
    if (median)
        return medianWavelength(lowRes, highRes);
    
    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>());
    
	double totalWavelength = 0;
	double count = 0;

    double low, high;
	StatisticsManager::convertResolutions(lowRes, highRes, &low, &high);

    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
    if (usingReference)
        refinementIntensityThreshold = 0.3;
    
	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			if (reflection(i)->getResolution() < low)
				continue;

			if (reflection(i)->getResolution() > high)
				continue;
            
            double wavelength = reflection(i)->miller(0)->getWavelength();
            
            if (wavelengthRange.size())
            {
                if (wavelength < wavelengthRange[0] || wavelength > wavelengthRange[1])
                {
                    continue;
                }
            }
            
			double weight = 1;
			double isigi = reflection(i)->miller(j)->getRawestIntensity();

            MtzManager *reference = MtzManager::getReferenceManager();
            Reflection *refReflection = NULL;
            
            if (reference != NULL)
            {
                reference->findReflectionWithId(reflection(i)->getReflId(), &refReflection);
                
                if (refReflection != NULL && refReflection->meanIntensity() < REFERENCE_WEAK_REFLECTION)
                    continue;
            }
            
            if (usingReference)
            {
                if (refReflection == NULL)
                    continue;
                
                double imgIntensity = reflection(i)->miller(j)->getRawestIntensity();
                double refIntensity = refReflection->meanIntensity();
                
                double proportion = imgIntensity / refIntensity;
                isigi = proportion;
            }
            
			if (isigi > refinementIntensityThreshold && isigi == isigi)
			{
				totalWavelength += wavelength
						* weight;
				count += weight;
			}
		}
	}

    if (count == 0)
    {
        this->setRejected(true);
    }

    double newWavelength = totalWavelength / count;
    
    logged << "Best wavelength " << newWavelength << " Å from " << count << " reflections." << std::endl;
    sendLog(LogLevelDetailed);
    return newWavelength;

}

double MtzManager::correlation(bool silent, double lowResolution,
		double highResolution)
{
	if (highResolution == -1)
		highResolution = maxResolutionAll;

	double correlation = StatisticsManager::cc_pearson(this,
			MtzManager::referenceManager, true, NULL,
			NULL, lowResolution, highResolution, false);

	return correlation;
}

double MtzManager::rSplitWithManager(MtzManager *otherManager, bool printHits,
		bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	StatisticsFunction *function = StatisticsManager::r_split;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::correlationWithManager(MtzManager *otherManager,
		bool printHits, bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	StatisticsFunction *function = StatisticsManager::cc_pearson;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
		StatisticsFunction *function, bool printHits, bool silent, double lowRes,
		double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	return statisticsWithManager(otherManager, function, NULL, RFactorNone,
			printHits, silent, lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::rFactorWithManager(RFactorType rFactor, bool printHits,
		bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	RFactorFunction *function = StatisticsManager::r_factor;
	return statisticsWithManager(NULL, NULL, function, rFactor, printHits,
			silent, lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
		StatisticsFunction *function, RFactorFunction *rFactorFunction,
		RFactorType rFactor, bool printHits, bool silent, double lowRes,
		double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	vector<double> shells;

	if (bins == 0)
		bins = 10;

	double maxResolution = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		if (reflection(i)->getResolution() > maxResolution
				&& reflection(i)->getResolution() < 1 / 1.2)
		{
            Reflection *bestReflection = reflection(i);
            
			maxResolution = bestReflection->getResolution();
		}
	}

	if (highRes == 0 || highRes < 1 / maxResolution)
		highRes = 1 / maxResolution;

	int hits = 0;
	double multiplicity = 0;

	double statistic = 0;

	if (rFactor == RFactorNone)
		statistic = function(this, otherManager, !printHits, &hits,
				&multiplicity, lowRes, highRes, shouldLog);
	else
		statistic = rFactorFunction(rFactor, this, &hits, &multiplicity, lowRes,
				highRes);

    std::cout << "N: " << "lowRes\thighRes\tValue\tHits\tMultiplicity" << std::endl;
    
	if (bins > 1 || !silent)
	{
		StatisticsManager::generateResolutionBins(lowRes, highRes, bins,
				&shells);

		for (int i = 0; i < shells.size() - 1; i++)
		{
			double statistic = 0;

			if (rFactor == RFactorNone)
				statistic = function(this, otherManager, 1, &hits,
						&multiplicity, shells[i], shells[i + 1], shouldLog);
			else
				statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
						shells[i], shells[i + 1]);

			if (!silent)
			{
				std::cout << "N: " << shells[i] << "\t" << shells[i + 1] << "\t"
						<< statistic << "\t" << hits << "\t" << multiplicity
						<< std::endl;
			}

			if (correlations != NULL)
			{
				boost::tuple<double, double, double, int> result = boost::make_tuple(
						shells[i], shells[i + 1], statistic, hits);
				correlations->push_back(result);
			}
		}
	}

	if (!silent)
	{
		double statistic = 0;

		if (rFactor == RFactorNone)
			statistic = function(this, otherManager, 1,  &hits,
					&multiplicity, lowRes, highRes, shouldLog);
		else
			statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
					lowRes, highRes);

		std::cout << "N: *** Overall ***" << std::endl;
		std::cout << "N: " << lowRes << "\t" << highRes << "\t" << statistic
				<< "\t" << hits << "\t" << multiplicity << std::endl;
	}

	return statistic;
}

double MtzManager::wavelengthStandardDeviation()
{
    double minResolution = maxResolutionRlpSize;
    double maxResolution = maxResolutionAll;
    int reflectionCount = 100;
    
    vector<Reflection *>refReflections, imageReflections;
    std::map<double, double> wavelengths;
    
    this->findCommonReflections(getReferenceManager(), imageReflections, refReflections);
    
    for (int i = 0; i < imageReflections.size(); i++)
    {
        double refIntensity = refReflections[i]->meanIntensity();
        if (refIntensity < 1000)
            continue;
        
        if (!imageReflections[i]->betweenResolutions(minResolution, maxResolution))
            continue;
        
        for (int j = 0; j < imageReflections[i]->millerCount(); j++)
        {
            double intensity = imageReflections[i]->miller(j)->getRawestIntensity();
            
            double wavelength = imageReflections[i]->miller(0)->getWavelength();
            
            if (wavelength < 1.25 || wavelength > 1.33)
                continue;
            
            double proportion = intensity / refIntensity;
            
            wavelengths[proportion] = wavelength;
        }
    }
    
    if (wavelengths.size() == 0)
        return FLT_MAX;
    
    int count = 0;
    vector<double> bigWavelengths;
    
    std::map<double, double>::iterator it = wavelengths.end();
    it--;
    
    for (; it != wavelengths.begin() && count < reflectionCount; --it)
    {
        bigWavelengths.push_back(it->second);
        count++;
    }
    
    
    double stDev = standard_deviation(&bigWavelengths);
    this->wavelength = median(&bigWavelengths);
    
//    std::cout << "St dev: " << stDev << ", median: " << wavelength << std::endl;
    
    return stDev;
}