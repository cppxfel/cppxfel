/*
 * IOMRefiner.cpp
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#include "IOMRefiner.h"
#include "Vector.h"
#include <cmath>
#include "gaussianfit.h"
#include "Holder.h"
#include "parameters.h"
#include <map>
#include <iostream>
#include "misc.h"
#include <fstream>
#include "Panel.h"
#include "FileParser.h"
#include "Miller.h"


#define BIG_BANDWIDTH 0.015
#define DISTANCE_TOLERANCE 0.01
#define WAVELENGTH_TOLERANCE 0.0001

#define ANGLE_TOLERANCE 0.0001
#define SPOT_DISTANCE_TOLERANCE 5

double IOMRefiner::intensityThreshold;
bool IOMRefiner::absoluteIntensity = false;
bool IOMRefiner::lowIntensityPenalty = false;

IOMRefiner::IOMRefiner(ImagePtr newImage, MatrixPtr matrix)
{
    int spgNum = FileParser::getKey("SPACE_GROUP", -1);

    Reflection::setSpaceGroup(spgNum);

    spaceGroup = ccp4spg_load_by_standard_num(spgNum);
    initialStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);

    testWavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    testDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    testBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
                                       OVER_PRED_BANDWIDTH) / 2;
    testSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
                                      OVER_PRED_SPOT_SIZE);
    maxResolution = FileParser::getKey(
                                       "MAX_INTEGRATED_RESOLUTION", MAX_INTEGRATED_RESOLUTION);;
    minResolution = FileParser::getKey(
                                       "MIN_INTEGRATED_RESOLUTION", 0.0);;
    searchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                    METROLOGY_SEARCH_SIZE);

    reference = NULL;
    lastMtz = MtzPtr();
    roughCalculation = FileParser::getKey("ROUGH_CALCULATION", false);
    hRot = 0;
    kRot = 0;
    lRot = 0;
    aRot = 0;
    bRot = 0;
    cRot = 0;
    bestHRot = 0;
    bestKRot = 0;
    bestLRot = 0;
    lastTotal = 0;
    lastStdev = 0;
    expectedSpots = FileParser::getKey("EXPECTED_SPOTS", 30);
    intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD", INTENSITY_THRESHOLD);
    refinement = RefinementTypeOrientationMatrixEarly;
    image = newImage;
    absoluteIntensity = FileParser::getKey("ABSOLUTE_INTENSITY", false);
    lowIntensityPenalty = FileParser::getKey("LOW_INTENSITY_PENALTY", false);
    this->matrix = matrix;
    refineA = FileParser::getKey("REFINE_UNIT_CELL_A", false);
    refineB = FileParser::getKey("REFINE_UNIT_CELL_B", false);;
    refineC = FileParser::getKey("REFINE_UNIT_CELL_C", false);;
    unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    orientationTolerance = FileParser::getKey("INDEXING_ORIENTATION_TOLERANCE", INDEXING_ORIENTATION_TOLERANCE);
    needsReintegrating = true;
    recalculateMillerPositions = false;

    int rotationModeInt = FileParser::getKey("ROTATION_MODE", 0);
    rotationMode = (RotationMode)rotationModeInt;

    complexUnitCell = false;
}

double IOMRefiner::getDetectorDistance()
{
    return getImage()->getDetectorDistance();
}

double IOMRefiner::getWavelength()
{
    return getImage()->getWavelength();
}

void IOMRefiner::setComplexMatrix()
{
    complexUnitCell = true;

    double *unitCellDouble = new double[3];
    matrix->unitCellLengths(&unitCellDouble);
    this->unitCell[0] = unitCellDouble[0];
    this->unitCell[1] = unitCellDouble[1];
    this->unitCell[2] = unitCellDouble[2];

    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());

    if (unitCell.size() == 0)
    {
        std::cout
        << "Please provide unit cell dimensions under keyword UNIT_CELL"
        << std::endl;
        exit(1);
    }

    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);

    if (fixUnitCell)
    {
        getImage()->setUnitCell(unitCell);
        setUnitCell(unitCell);

        matrix->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
    }

    delete [] unitCellDouble;
}

void IOMRefiner::dropMillers()
{
    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);
}

bool IOMRefiner::millerReachesThreshold(MillerPtr miller)
{
    double iSigI = miller->getRawIntensity() / miller->getCountingSigma();

    std::ostringstream logged;

    logged << "Absolute intensity is " << absoluteIntensity << ", iSigI is " << iSigI << ", raw intensity is " << miller->getRawIntensity() << std::endl;

    Logger::mainLogger->addStream(&logged, LogLevelDebug);

    if (absoluteIntensity)
    {
        return (miller->getRawIntensity() > intensityThreshold);
    }

    return (iSigI > intensityThreshold);
}

void IOMRefiner::getWavelengthHistogram(vector<double> &wavelengths,
                                     vector<int> &frequencies, LogLevel level, int whichAxis)
{
    wavelengths.clear();
    frequencies.clear();


    double wavelength = getImage()->getWavelength();
    vector<double> totals;

    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>(2, 0));

    double spread = testBandwidth * 2;
    double interval = (wavelength * spread * 2) / 20;

    double minLength = wavelength * (1 - spread);
    double maxLength = wavelength * (1 + spread);

    if (wavelengthRange[0] != 0 || wavelengthRange[1] != 0)
    {
        minLength = wavelengthRange[0];
        maxLength = wavelengthRange[1];
        interval = (maxLength - minLength) / 20;
    }

    std::ostringstream logged;
    logged << "Wavelength histogram for " << this->getImage()->getFilename() << std::endl;

    for (double i = minLength; i < maxLength; i += interval)
    {
        wavelengths.push_back(i);
        double frequency = 0;
        double total = 0;

        for (int j = 0; j < millers.size(); j++)
        {
            double ewald = millers[j]->getWavelength();

            if (ewald < i || ewald > i + interval)
                continue;

            bool strong = millerReachesThreshold(millers[j]);

            double weight = whichAxis == 0 ? 1 : millers[j]->getEwaldWeight(hRot, kRot, whichAxis == 1);

         //   if (resolutionWeights)
         //       weight *= millers[j]->getResolution();

            total += weight;
            if (strong)
            {
                frequency += weight;
            }
            else if (lowIntensityPenalty)
            {
                frequency -= weight / 4;
            }

            if (frequency < 0)
            {
                frequency = 0;
            }
        }

        logged << std::setprecision(4) << i << "\t";

        for (int i=0; i < frequency; i++)
        {
            logged << ".";
        }

        logged << std::endl;

        frequencies.push_back(frequency);
        totals.push_back(total);
    }

    int strong = getTotalReflections();

    logged << std::endl;
    logged << "Total strong reflections: " << strong << std::endl;


    Logger::mainLogger->addStream(&logged, level);
}

void IOMRefiner::calculateNearbyMillers(bool rough)
{
    MatrixPtr matrix = getMatrix();
    double wavelength = getImage()->getWavelength();

    double minBandwidth = wavelength * (1 - testBandwidth * 2);
    double maxBandwidth = wavelength * (1 + testBandwidth * 2);

    double sphereThickness = FileParser::getKey("SPHERE_THICKNESS", 0.01);

    double minSphere = 1 / wavelength * (1 - sphereThickness);
    double maxSphere = 1 / wavelength * (1 + sphereThickness);

    double minSphereSquared = pow(minSphere, 2);
    double maxSphereSquared = pow(maxSphere, 2);

    nearbyMillers.clear();

    int maxMillers[3];

    MatrixPtr newMatrix = matrix->copy();

    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;

    if (!(hRad == 0 && kRad == 0))
        newMatrix->rotate(hRad, kRad, 0);

    newMatrix->maxMillers(maxMillers, maxResolution);

    logged << "Integrating to maximum Miller indices: (" << maxMillers[0] << ", " << maxMillers[1] << ", " << maxMillers[2] << ")" << std::endl;

    int overRes = 0;
    int underRes = 0;
    double maxDistSquared = pow(1 / maxResolution, 2);
    double minResolutionSquared = pow(1 / minResolution, 2);

    for (int h = -maxMillers[0]; h < maxMillers[0]; h++)
    {
        for (int k = -maxMillers[1]; k < maxMillers[1]; k++)
        {
            for (int l = -maxMillers[2]; l < maxMillers[2]; l++)
            {
                if (ccp4spg_is_sysabs(spaceGroup, h, k, l))
                    continue;

                MillerPtr newMiller = MillerPtr(new Miller(NULL, h, k, l));
                newMiller->setImageAndIOMRefiner(getImage(), this);

                if (rough == false)
                {
                    nearbyMillers.push_back(newMiller);
                    continue;
                }

                vec hkl = new_vector(h, k, l);
                newMatrix->multiplyVector(&hkl);

                if (hkl.l > 0)
                    continue;

                vec beam = new_vector(0, 0, -1 / wavelength);
                vec beamToMiller = vector_between_vectors(beam, hkl);

                double sphereRadiusSquared = length_of_vector_squared(beamToMiller);

                if (sphereRadiusSquared < minSphereSquared || sphereRadiusSquared > maxSphereSquared)
                {
                    double ewaldSphere = getEwaldSphereNoMatrix(hkl);

                    if (ewaldSphere < minBandwidth || ewaldSphere > maxBandwidth)
                    {
                        continue;
                    }
                }

                double res = length_of_vector_squared(hkl);

                if (res > maxDistSquared)
                {
                    overRes++;
                    continue;
                }

                if (minResolution > 0 && res < minResolutionSquared)
                {
                    underRes++;
                    continue;
                }

                if (h == 0 && k == 0 && l == 0)
                    continue;

               nearbyMillers.push_back(newMiller);
            }
        }
    }

    logged << "Rejected " << overRes << " due to being over resolution edge of " << maxResolution << " Å." << std::endl;

    if (minResolution > 0)
        logged << "Rejected " << underRes << " due to being under resolution edge of " << minResolution << " Å." << std::endl;

    sendLog(LogLevelDetailed);
    needsReintegrating = true;
}

void IOMRefiner::lockUnitCellDimensions()
{
    double spgNum = spaceGroup->spg_num;
    if (spgNum >= 75 && spgNum <= 194)
    {
        unitCell[1] = unitCell[0];
    }
    if (spgNum >= 195)
    {
        unitCell[1] = unitCell[0];
        unitCell[2] = unitCell[0];
    }
}

void IOMRefiner::checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox, bool perfectCalculation)
{
    MatrixPtr matrix = getMatrix();

    if (complexUnitCell)
    {
   //     unitCell[1] = unitCell[0];

        lockUnitCellDimensions();
        matrix->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
    }
    double wavelength = getImage()->getWavelength();
    double maxD = 1 / maxResolution;
    if (maxResolution == 0)
        maxD = FLT_MAX;

    if (testDistance != 0)
        getImage()->setDetectorDistance(testDistance);

    if (testWavelength != 0)
    {
        getImage()->setWavelength(testWavelength);
        wavelength = testWavelength;
    }

    millers.clear();

    double averageEwald = 0;

    logged << "Testing " << nearbyMillers.size() << " reflections close to the Ewald sphere with wavelength " << wavelength << std::endl;

    int cutResolution = 0;
    int partialityTooLow = 0;
    int unacceptableIntensity = 0;

    MatrixPtr newMatrix = MatrixPtr();

    if (rotationMode == RotationModeHorizontalVertical)
    {
        logged << "Rotating by " << hRot << ", " << kRot << ", " << lRot << std::endl;
        sendLog(LogLevelDetailed);

        Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &newMatrix);
    }
    else if (rotationMode == RotationModeUnitCellABC)
    {
        Miller::rotateMatrixABC(aRot, bRot, cRot, matrix, &newMatrix);
    }

    lastRotatedMatrix = newMatrix;

    std::vector<MillerPtr> *chosenMillerArray = &nearbyMillers;

    if (!perfectCalculation && roughCalculation && !needsReintegrating && roughMillers.size())
    {
        chosenMillerArray = &roughMillers;
    }
    else if (needsReintegrating)
    {
        bandwidth *= 2;
        roughMillers.clear();
        std::vector<MillerPtr>().swap(roughMillers);
    }

    logged << "Checking " << chosenMillerArray->size() << " reflections." << std::endl;

    for (int i = 0; i < chosenMillerArray->size(); i++)
    {
        MillerPtr miller = (*chosenMillerArray)[i];

        vec hkl = new_vector(miller->getH(), miller->getK(), miller->getL());
        matrix->multiplyVector(&hkl);

    //    int roughX = 0;
    //    int roughY = 0;

        miller->setMatrix(matrix);

        double d = length_of_vector(hkl);
        if (d > maxD)
        {
            cutResolution++;
            continue;
        }

        miller->crossesBeamRoughly(newMatrix, 0.0, testSpotSize, wavelength, bandwidth);

      //  miller->recalculatePartiality(newMatrix, 0.0, testSpotSize,
       //                               wavelength, bandwidth, 1.5, true);

        if (i == 0)
        {
            logged << "Calculated partiality with parameters: hRot " << hRot << ", kRot " << kRot << ", spot size " << testSpotSize << ", wavelength " << wavelength << ", bandwidth " << bandwidth << std::endl;
            sendLog(LogLevelDebug);
        }

        if (miller->getPartiality() <= 0.05)
        {
        //    logged << "Rejected Miller partiality too low at\t" << roughX << "\t" << roughY << std::endl;
            //sendLog(LogLevelDebug);
            partialityTooLow++;
            continue;
        }

  //      miller->setPartialityModel(PartialityModelScaled);
        miller->getWavelength(newMatrix);

        if (complexShoebox)
        {
            double initialBandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
            double initialRlpSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
            double initialMosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);

            miller->makeComplexShoebox(wavelength, initialBandwidth, initialMosaicity, initialRlpSize);
        }

        if (needsReintegrating || !roughCalculation || recalculateMillerPositions)
        {
            miller->integrateIntensity(newMatrix);
        }

        double rawIntensity = miller->getRawIntensity();

        if (rawIntensity != rawIntensity)
        {
            unacceptableIntensity++;
            chosenMillerArray->erase(chosenMillerArray->begin() + i);
            i--;
            continue;
        }

        averageEwald += miller->getWavelength();

        this->millers.push_back(miller);

        if (needsReintegrating && millerReachesThreshold(miller))
        {
            roughMillers.push_back(miller);
        }
    }

    logged << "Using wavelength " << wavelength << " Å and distance " << getImage()->getDetectorDistance() << " mm" << std::endl;
    logged << "Beyond resolution cutoff: " << cutResolution << std::endl;
    logged << "Partiality equal to 0: " << partialityTooLow << std::endl;
    logged << "Image pixels were masked/flagged: " << unacceptableIntensity << std::endl;
    logged << "Reflections accepted: " << millers.size() << std::endl;

    sendLog(LogLevelDetailed);

    needsReintegrating = false;
}

double IOMRefiner::minimizeParameter(double *meanStep, double *param, int whichAxis)
{
    double param_trials[3];
    double param_scores[3];

    int j = 0;
    double param_min_score = FLT_MAX;
    int param_min_num = 1;

    double bestParam = *param;

    std::ostringstream logged;
    logged << "Scores for " << getImage()->getFilename() << ": ";

    for (double i = bestParam - *meanStep; j < 3; i += *meanStep)
    {
        *param = i;
        this->checkAllMillers(maxResolution, testBandwidth, false, false);
        param_scores[j] = score(whichAxis);
        logged << param_scores[j] << ", ";
        param_trials[j] = i;
        j++;
    }

    param_min_score = param_scores[1];

    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }

    *param = param_trials[param_min_num];
    this->checkAllMillers(maxResolution, testBandwidth, false, false);

    logged << "chosen no. " << param_min_num << std::endl;

    Logger::mainLogger->addStream(&logged, LogLevelDetailed);

    if (param_min_num == 1)
        *meanStep /= 2;

    return param_min_score;
}

void IOMRefiner::minimizeTwoParameters(double *meanStep1, double *meanStep2,
                                    double *param1, double *param2)
{
    double param_trials1[9];
    double param_trials2[9];
    double param_scores[9];

    int j = 0;
    double param_min_score = FLT_MAX;
    int param_min_num = 4;

    double bestParam1 = *param1;
    double bestParam2 = *param2;

    bool perfect = false;

    for (double i = bestParam1 - *meanStep1; j < 3; i += *meanStep1)
    {
        int l = 0;

        for (double k = bestParam2 - *meanStep2; l < 3; k += *meanStep2)
        {
            *param1 = i;
            *param2 = k;
            this->checkAllMillers(maxResolution, testBandwidth, false, perfect);
            param_scores[j * 3 + l] = score(0, false);
            param_trials1[j * 3 + l] = i;
            param_trials2[j * 3 + l] = k;
            l++;
        }
        j++;
    }

    param_min_score = param_scores[4];

    for (int i = 0; i < 9; i++)
    {
        if ((param_scores[i] < param_min_score))
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }

    }

    *param1 = param_trials1[param_min_num];
    *param2 = param_trials2[param_min_num];

    if (param_min_num == 4)
    {
        *meanStep1 /= 2;
        *meanStep2 /= 2;
    }
}

int IOMRefiner::getTotalReflections(double threshold)
{
    int count = 0;

    for (int i = 0; i < millers.size(); i++)
    {
        if (millers[i]->getRawIntensity() > threshold)
            count++;
    }

    return count;
}

int IOMRefiner::getTotalReflections()
{
    int count = 0;

    for (int i = 0; i < millers.size(); i++)
    {
        if (millerReachesThreshold(millers[i]))
            count++;
    }

    return count;
}

bool IOMRefiner::millerWithinBandwidth(MillerPtr miller)
{
    double minBandwidth = getImage()->getWavelength() * (1 - testBandwidth * 2);
    double maxBandwidth = getImage()->getWavelength() * (1 + testBandwidth * 2);

    double wavelength = miller->getWavelength();

    return (wavelength > minBandwidth && wavelength < maxBandwidth);
}

int IOMRefiner::getTotalReflectionsWithinBandwidth()
{
    int count = 0;

    for (int i = 0; i < millers.size(); i++)
    {
        if (!millerWithinBandwidth(millers[i]))
            continue;

        if (millerReachesThreshold(millers[i]))
            count++;
    }

    return count;
}

void IOMRefiner::findSpots()
{
    int tolerance = 60;

    //  getImage()->printBox(1180, 660, 20);

    for (int i = 0; i < getImage()->getXDim(); i += tolerance)
    {
        for (int j = 0; j < getImage()->getYDim(); j += tolerance)
        {
            Spot *spot = new Spot(getImage());

            double maxLift = 0;
            double maxX = 0;
            double maxY = 0;

            for (int tolX = 0; tolX < tolerance; tolX++)
            {
                for (int tolY = 0; tolY < tolerance; tolY++)
                {
                    double x = i + tolX;
                    double y = j + tolY;

                    double lift = spot->maximumLift(getImage(), x, y);

                    if (lift > maxLift)
                    {
                        maxLift = lift;
                        maxX = x;
                        maxY = y;
                    }
                }
            }

            if (maxLift > 0)
            {
                spot->setXY(maxX, maxY);
                spots.push_back(spot);

                getImage()->addSpotCover(maxX - 30, maxY - 30, maxX + 30, maxY + 30);
                //      getImage()->printBox(maxX, maxY, 8);
            }
        }
    }

    Spot::sortSpots(&spots);

    std::string name = "spots-" + getImage()->getFilename();
    int lastindex = (int)name.find_last_of(".");
    std::string rootName = name.substr(0, lastindex);
    std::string datName = rootName + ".dat";
    writeDatFromSpots(datName);

    std::ostringstream logged;
    logged << "Found " << spots.size() << " spots" << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelNormal);
}

void IOMRefiner::duplicateSpots(vector<ImagePtr> images)
{
    std::map<vector<int>, int> frequencies =
    std::map<vector<int>, int>();

    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->IOMRefinerCount(); j++)
        {
            IOMRefinerPtr IOMRefiner = images[i]->getIOMRefiner(j);

            vector<Spot *> spots = IOMRefiner->getSpots();

            std::cout << "Image: " << i << " Spot count: " << spots.size()
            << std::endl;

            for (int j = 0; j < spots.size(); j++)
            {
                vector<int> coord = vector<int>();
                coord.push_back(spots[j]->getX());
                coord.push_back(spots[j]->getY());

                if (frequencies.count(coord) == 0)
                {
                    frequencies[coord] = 1;
                }
                else
                {
                    frequencies[coord]++;
                }
            }
        }
    }

    int threshold = 0.2 * images.size();
    int spotsRemoved = 0;

    if (images.size() < 4)
        return;

    for (std::map<vector<int>, int>::iterator it = frequencies.begin();
         it != frequencies.end(); ++it)
    {
        std::cout << it->first[0] << "\t" << it->first[1] << "\t"
        << frequencies[it->first] << std::endl;

        if (frequencies[it->first] >= threshold)
        {
            int startX = it->first[0] - 2;
            int endX = it->first[0] + 2;
            int startY = it->first[1] - 2;
            int endY = it->first[1] + 2;

            Image::applyMaskToImages(images, startX, startY, endX, endY);
            spotsRemoved++;
        }
    }

    std::cout << "Spots removed: " << spotsRemoved << std::endl;
}

void IOMRefiner::writeDatFromSpots(std::string filename)
{
    std::ofstream dat;
    dat.open(filename);
    int count = 0;

    for (int i = 0; i < expectedSpots && i < spots.size(); i++)
    {
        //      if (spots[i]->isAcceptable(image))
        //      {
        dat << "0\t0\t0\t1\t1\t1\t" << spots[i]->getX() << "\t" << spots[i]->getY()
        << "\t1" << std::endl;

        std::cout << "Spot lift " << spots[i]->weight() << std::endl;

        count++;
        //      }
    }

    std::cout << count << " spots remain on image." << std::endl;

    dat.close();
}


double IOMRefiner::score(int whichAxis, bool silent)
{
    if (refinement == RefinementTypeDetectorWavelength)
    {
        double value = getTotalReflections();
        logged << "Total reflections within bandwidth: " << value << std::endl;
        sendLog();
        return 0 - value;
    }

    if (refinement == RefinementTypeOrientationMatrixReverse)
    {
        MatrixPtr invRotation = lastRotatedMatrix->getRotation()->inverse3DMatrix();
        MatrixPtr invTransform = lastRotatedMatrix->getUnitCell()->inverse3DMatrix();
        MatrixPtr newMatrix = MatrixPtr();
        Miller::rotateMatrixHKL(hRot, kRot, lRot, MatrixPtr(new Matrix()), &newMatrix);

        std::vector<double> ewaldWavelengths;

        for (int i = 0; i < getImage()->spotCount(); i++)
        {
            SpotPtr spot = getImage()->spot(i);

            vec estimatedVec = spot->estimatedVector();

            newMatrix->multiplyVector(&estimatedVec);

            double ewald = getEwaldSphereNoMatrix(estimatedVec);

            ewaldWavelengths.push_back(ewald);

          //  invRotation->printDescription();

            invTransform->multiplyVector(&estimatedVec);
            invRotation->multiplyVector(&estimatedVec);

        //    vec remainder = new_vector(fmod(estimatedVec.h, 1), fmod(estimatedVec.k, 1), fmod(estimatedVec.l, 1));
            logged << estimatedVec.h << ", " << estimatedVec.k << ", " << estimatedVec.l << std::endl;

         //   double addition = length_of_vector_squared(remainder);

         //   sumSqrDiff += addition;

         //   sendLog();
        }

        double stdev = standard_deviation(&ewaldWavelengths);

        logged << "Standard deviation of Ewald sphere wavelengths: " << stdev << std::endl;
        sendLog();

        return stdev;
    }

    if (refinement == RefinementTypeRefineLAxis)
    {
        double averageShift = 0;
        int count = 0;

        for (int i = 0; i < millers.size(); i++)
        {
            if (millerReachesThreshold(millers[i]))
            {
                std::pair<double, double> shift = millers[i]->getShift();
                double shiftDistance = sqrt(pow(shift.first, 2) + pow(shift.second, 2));

                averageShift += shiftDistance;
                count++;
            }
        }

        averageShift /= count;

        logged << "Average shift: " << averageShift << std::endl;
        sendLog(LogLevelDetailed);
        return averageShift;
    }

    if (refinement == RefinementTypeOrientationMatrixEarly || refinement == RefinementTypeOrientationMatrixEarlySeparated)
    {
        vector<double> wavelengths;
        vector<int> frequencies;

        switch (whichAxis)
        {
            case 0:
                Logger::mainLogger->addString("Optimising both axes", LogLevelDetailed);
                break;
            case 1:
                Logger::mainLogger->addString("Optimising H axis", LogLevelDetailed);
                break;
            case 2:
                Logger::mainLogger->addString("Optimising K axis", LogLevelDetailed);
                break;
        }

        LogLevel level = silent ? LogLevelDebug : LogLevelDetailed;
        if (whichAxis != 0 && silent == 0) level = LogLevelDebug;

        getWavelengthHistogram(wavelengths, frequencies, level, whichAxis);

        double mean = 0;
        double stdev = 0;
        histogram_gaussian(&wavelengths, &frequencies, mean, stdev);

        double total = getTotalReflectionsWithinBandwidth();
        //      double totalIntensity = getTotalIntegratedSignal();

        double totalWeight = 20;
        double stdevWeight = 15;

        double totalChange = total / lastTotal - 1;
        double totalStdev = stdev / lastStdev - 1;

        double score = (totalStdev * stdevWeight - totalChange * totalWeight);

        return score;
    }

    if (refinement == RefinementTypeOrientationMatrixEarlyWeighted)
    {
        std::vector<double> wavelengths, weights;

        for (int i = 0; i < millers.size(); i++)
        {
            MillerPtr miller = millers[i];

            if (!millerReachesThreshold(miller))
                continue;

            double weight = miller->getRawestIntensity();

            if (weight < 0) weight = 0;

            double wavelength = miller->getWavelength();

            if (!(wavelength == wavelength && weight == weight))
                continue;

            wavelengths.push_back(wavelength);
            weights.push_back(weight);
        }

        double stdev = standard_deviation(&wavelengths, &weights);

        return stdev;
    }

    if (refinement == RefinementTypeOrientationMatrixStdevOnly)
    {
        vector<double> wavelengths;
        vector<int> frequencies;
        getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed, 0);

        wavelengths.clear();
        vector<double>().swap(wavelengths);

        for (int i = 0; i < millers.size(); i++)
        {
            if (millerReachesThreshold(millers[i]) && millerWithinBandwidth(millers[i]))
            {
                wavelengths.push_back(millers[i]->getWavelength());
            }
        }

        double stdev = standard_deviation(&wavelengths);
        return stdev;
    }

    if (refinement == RefinementTypeOrientationMatrixHighestPeak)
    {
        vector<double> wavelengths;
        vector<int> frequencies;

        getWavelengthHistogram(wavelengths, frequencies);

        double mean = 0;
        double stdev = 0;
        histogram_gaussian(&wavelengths, &frequencies, mean, stdev);

        int *highest = new int[3];

        for (int i = 0; i < frequencies.size(); i++)
        {
            if (frequencies[i] > highest[0])
            {
                for (int j = 0; j < 2; j++)
                    highest[j + 1] = highest[j];
                highest[0] = frequencies[i];
            }
        }

        double sum = 0;
        for (int i = 0; i < 3; i++)
            sum += highest[i];

        return 1 / sum;
    }

    if (refinement == RefinementTypeOrientationMatrixRough)
    {
        vector<double> wavelengths;

        for (int i = 0; i < millers.size(); i++)
        {
            if (millerReachesThreshold(millers[i]))
            {
                wavelengths.push_back(millers[i]->getWavelength());
            }
        }

        double stdev = standard_deviation(&wavelengths);
        double num = wavelengths.size();

        double newScore = stdev / num;

        return newScore;
    }

    return 0;
}

void IOMRefiner::refineDetectorAndWavelength(MtzManager *reference)
{
    int oldSearch = getSearchSize();
    setSearchSize(1);
    this->reference = reference;
    this->calculateNearbyMillers(true);
    checkAllMillers(maxResolution, testBandwidth);
    refinement = RefinementTypeOrientationMatrixEarly;

    testDistance = getImage()->getDetectorDistance();
    testWavelength = getImage()->getWavelength();
    double oldDistance = testDistance;
    double oldWavelength = testWavelength;

    int count = 0;
    double distStep = 0.13 * 5;
    double waveStep = 0.01 * 5;

    double cStep = FileParser::getKey("STEP_UNIT_CELL_C", 0.2);;

    bool refinedC = false;


    bool refinedDist = false;
    bool refinedWave = false;

    vector<double> wavelengths;
    vector<int> frequencies;

    double mean = 0;
    double stdev = 0;
    getWavelengthHistogram(wavelengths, frequencies);
    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);

    double oldScore = score();
    double newScore = 0;

    while (!(refinedDist && refinedWave) && count < 50)
    {
     /*   if (!refinedDist)
            this->minimizeParameter(&distStep, &testDistance);

        if (!refinedWave)
            this->minimizeParameter(&waveStep, &testWavelength);
        */
        this->minimizeTwoParameters(&distStep, &waveStep, &testDistance, &testWavelength);

        if (!refinedC)
            minimizeParameter(&cStep, &unitCell[2]);

        newScore = score();

        sendLog(LogLevelNormal);

        if (cStep < 0.01)
            refinedC = true;

        if (distStep < DISTANCE_TOLERANCE)
            refinedDist = true;

        if (waveStep < WAVELENGTH_TOLERANCE)
            refinedWave = true;
    }

    getWavelengthHistogram(wavelengths, frequencies);
    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);

 //   testWavelength = mean;
 //   this->getImage()->setWavelength(mean);

    logged << "Distance refined for " << getImage()->getFilename() << ":\t" << oldDistance << "\t" << oldWavelength << "\t" << testDistance << "\t" << testWavelength << "\t" << oldScore << "\t" << newScore << std::endl;

    sendLog(LogLevelNormal);

    getImage()->setDetectorDistance(testDistance);

    setSearchSize(oldSearch);
}

bool compareScore(std::pair<vector<double>, double> a, std::pair<vector<double>, double> b)
{
    return (a.second < b.second);
}

void IOMRefiner::matchMatrixToSpots(RefinementType refinement)
{
    double wedge = FileParser::getKey("INDEXING_SLICE_ANGLE", 30.0); // degrees
    //  map<vector<double>, double> scores = map<vector<double>, double>();
    vector<std::pair<vector<double>, double> > scores;

    MatrixPtr copyMatrix = this->getMatrix()->copy();

    this->refinement = refinement;

    calculateNearbyMillers(false);

    for (double hr = 0; hr < 360; hr += wedge)
    {
        double hRad = hr * M_PI / 180;

        std::cout << hr / 360 * 100 << "% ..." << std::endl;

        for (double kr = 0; kr < 360; kr += wedge)
        {
            double kRad = kr * M_PI / 180;

            this->getMatrix()->rotate(hRad, kRad, 0);
            double theScore = score();

            vector<double> rotation = vector<double>();
            rotation.push_back(hRad);
            rotation.push_back(kRad);

            std::pair<vector<double>, double> pair = std::make_pair(rotation,
                                                             theScore);
            scores.push_back(pair);

            std::cout << hRad << "\t" << kRad << "\t" << theScore << std::endl;

            this->setMatrixCopy(copyMatrix);
        }
    }

    std::cout << "100%" << std::endl;

    vector<double> bestNum = vector<double>();

    std::sort(scores.begin(), scores.end(), compareScore);

    std::cout << std::endl << "Solutions to try: " << std::endl;

    for (int i = 0; i < 6; i++)
    {
        std::cout << scores[i].first[0] << "\t" << scores[i].first[1] << "\t"
        << scores[i].second << std::endl;

        solutions.push_back(scores[i].first);
    }

    this->getMatrix()->rotate(scores[0].first[0], scores[0].first[1], 0);
    std::cout << "Max score: " << score() << std::endl;

    this->setMatrix(copyMatrix);

    std::cout << std::endl;
}

double IOMRefiner::getRot(int rotNum)
{
    if (rotationMode == RotationModeHorizontalVertical)
    {
        switch (rotNum) {
            case 0:
                return hRot;
            case 1:
                return kRot;
            case 2:
                return lRot;
            default:
                break;
        }
    }
    else if (rotationMode == RotationModeUnitCellABC)
    {
        switch (rotNum) {
            case 0:
                return aRot;
            case 1:
                return bRot;
            case 2:
                return cRot;
            default:
                break;
        }
    }

    return 0;
}

void IOMRefiner::refineOrientationMatrix()
{
    int orientationScore = FileParser::getKey("ORIENTATION_SCORE", 0);
    RefinementType refinementType = (RefinementType)orientationScore;

    this->refineOrientationMatrix(refinementType);
}

void IOMRefiner::refineOrientationMatrix(RefinementType refinementType)
{
    refinement = refinementType;
    this->calculateNearbyMillers(true);

    testDistance = getImage()->getDetectorDistance();
    testWavelength = getImage()->getWavelength();

    vector<double> wavelengths;
    vector<int> frequencies;

    checkAllMillers(maxResolution, testBandwidth);
    Logger::mainLogger->addString("Wavelength histogram before refinement", LogLevelDetailed);
    sendLog(LogLevelDetailed);
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);

    double mean = 0;
    double stdev = 0;
    double theScore = 0;

    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);

    lastStdev = stdev;
    lastTotal = getTotalReflections();

    getImage()->setWavelength(mean);

    bool recalculated = false;

    for (int i = 0; i < 1; i++)
    {
        double hRotStep = initialStep;
        double kRotStep = initialStep;
        double lRotStep = initialStep / 4;
        double aRotStep = initialStep;
        double bRotStep = initialStep;
        double cRotStep = initialStep;

        double aStep = FileParser::getKey("STEP_UNIT_CELL_A", 0.2);
        double bStep = FileParser::getKey("STEP_UNIT_CELL_B", 0.2);
        double cStep = FileParser::getKey("STEP_UNIT_CELL_C", 0.2);

        double alphaStep = FileParser::getKey("STEP_UNIT_CELL_ALPHA", 0.5);
        double betaStep = FileParser::getKey("STEP_UNIT_CELL_BETA", 0.5);
        double gammaStep = FileParser::getKey("STEP_UNIT_CELL_GAMMA", 0.5);

        bool refinedH = (rotationMode == RotationModeHorizontalVertical) ? false : true;
        bool refinedK = (rotationMode == RotationModeHorizontalVertical) ? false : true;
        bool refinedL = (rotationMode == RotationModeHorizontalVertical) ? !FileParser::getKey("REFINE_IN_PLANE_OF_DETECTOR", false) : true;
        bool refinedA = !refineA;
        bool refinedB = !refineB;
        bool refinedC = !refineC;
        bool refinedAlpha = !FileParser::getKey("REFINE_UNIT_CELL_ALPHA", false);
        bool refinedBeta = !FileParser::getKey("REFINE_UNIT_CELL_BETA", false);
        bool refinedGamma = !FileParser::getKey("REFINE_UNIT_CELL_GAMMA", false);

        bool refinedARot = (rotationMode == RotationModeUnitCellABC) ? false : true;
        bool refinedBRot = (rotationMode == RotationModeUnitCellABC) ? false : true;
        bool refinedCRot = (rotationMode == RotationModeUnitCellABC) ? false : true;

        int count = 0;

        while (!(refinedH && refinedK && refinedL && (refinedA && refinedB && refinedC)) && count < 20)
        {
            if (!refinedL)
            {
                int oldSearchSize = searchSize;
                const int bigSize = 6;

                if (searchSize < bigSize)
                {
                    searchSize = bigSize;
                }

                recalculateMillerPositions = true;
                refinement = RefinementTypeRefineLAxis;

                this->minimizeParameter(&lRotStep, &lRot);
                refinement = refinementType;
                recalculateMillerPositions = false;

                searchSize = oldSearchSize;
            }
            if (!refinedH && !refinedK)
            {
                this->minimizeTwoParameters(&hRotStep, &kRotStep, &hRot, &kRot);
            }

            if (!refinedARot && !refinedBRot && !refinedCRot)
            {
                this->minimizeTwoParameters(&aRotStep, &bRotStep, &aRot, &bRot);
                this->minimizeTwoParameters(&bRotStep, &cRotStep, &bRot, &cRot);
                this->minimizeTwoParameters(&cRotStep, &aRotStep, &cRot, &aRot);
                this->minimizeTwoParameters(&aRotStep, &aStep, &aRot, &unitCell[0]);
                this->minimizeTwoParameters(&bRotStep, &bStep, &bRot, &unitCell[1]);
                this->minimizeTwoParameters(&cRotStep, &cStep, &cRot, &unitCell[2]);
            }

       //     refinementType = RefinementTypeOrientationMatrixPanelStdev;

            if (!refinedA)
                minimizeParameter(&aStep, &unitCell[0]);
            if (!refinedB)
                minimizeParameter(&bStep, &unitCell[1]);
            if (!refinedC)
                minimizeParameter(&cStep, &unitCell[2]);

            if (!refinedAlpha)
            {
                minimizeParameter(&alphaStep, &unitCell[3]);
            }

            if (!refinedBeta)
            {
                minimizeParameter(&betaStep, &unitCell[4]);
            }

            if (!refinedGamma)
            {
                minimizeParameter(&gammaStep, &unitCell[5]);
            }

       //     refinementType = RefinementTypeOrientationMatrixEarly;

            checkAllMillers(maxResolution, testBandwidth, false, false);

            getWavelengthHistogram(wavelengths, frequencies);
            histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
            lastStdev = stdev;
            lastTotal = getTotalReflections();

            double newScore = score();
            lastScore = newScore;

            logged << getRot(0) << "\t" << getRot(1) << "\t" << getRot(2) << "\t" << newScore << "\t" << aRotStep << std::endl;
            sendLog(LogLevelDetailed);

            if ((rotationMode == RotationModeHorizontalVertical && hRotStep < 0.25 && kRotStep < 0.25) ||
                (rotationMode == RotationModeUnitCellABC && aRotStep < 0.25 && bRotStep < 0.25 && cRotStep < 0.25))
            {
                if (!recalculated)
                {
                    recalculated = true;
                    this->calculateNearbyMillers(true);
                }

            //    refinement = RefinementTypeOrientationMatrixStdevOnly;
            }

            if (hRotStep < orientationTolerance)
                refinedH = true;
            if (kRotStep < orientationTolerance)
                refinedK = true;

            if (lRotStep < orientationTolerance)
                refinedL = true;

            if (aRotStep < orientationTolerance)
                refinedARot = true;
            if (bRotStep < orientationTolerance)
                refinedBRot = true;
            if (cRotStep < orientationTolerance)
                refinedCRot = true;

            if (aStep < 0.01)
                refinedA = true;
            if (bStep < 0.01)
                refinedB = true;
            if (cStep < 0.01)
                refinedC = true;

            if (alphaStep < 0.01)
                refinedAlpha = true;
            if (betaStep < 0.01)
                refinedBeta = true;
            if (gammaStep < 0.01)
                refinedGamma = true;

            count++;
        }

        count = 0;

        if (aStep < 0.01)
            refinedA = true;
        if (bStep < 0.01)
            refinedB = true;
        if (cStep < 0.01)
            refinedC = true;
    }

    refinement = RefinementTypeOrientationMatrixEarly;
    lastScore = score();

    logged << "Current wavelength: " << testWavelength << " Å." << std::endl;
    logged << "Rotation result:\t" << getImage()->getFilename() << "\t" << hRot
    << "\t" << kRot << "\t" << getTotalReflections() << "\t" << getLastScore() << std::endl;

    double hRad = getRot(0) * M_PI / 180;
    double kRad = getRot(1) * M_PI / 180;
    double lRad = getRot(2) * M_PI / 180;

    vector<double> originalUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    double *lengths = new double[3];
    getMatrix()->unitCellLengths(&lengths);

    double aRatio = originalUnitCell[0] / lengths[0];
    double bRatio = originalUnitCell[1] / lengths[1];
    double cRatio = originalUnitCell[2] / lengths[2];

    double aveRatio = (aRatio + bRatio + cRatio) / 3;

    lengths[0] *= aveRatio;
    lengths[1] *= aveRatio;
    lengths[2] *= aveRatio;

    delete [] lengths;

    if (rotationMode == RotationModeHorizontalVertical)
        getMatrix()->rotate(hRad, kRad, lRad);
    else
    {
        MatrixPtr oldMatrix = getMatrix()->copy();
        getMatrix()->rotateABC(oldMatrix, getRot(0), getRot(1), getRot(2));
    }

    bestHRot = getRot(0);
    bestKRot = getRot(1);
    bestLRot = getRot(2);

    hRot = 0;
    kRot = 0;
    lRot = 0;
    aRot = 0;
    bRot = 0;
    cRot = 0;

    needsReintegrating = true;
    checkAllMillers(maxResolution, testBandwidth);
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);

    sendLog(LogLevelNormal);
}

struct greater { template<class T> bool operator()(T const &a, T const &b) const { return a > b; } };

bool IOMRefiner::isBasicGoodSolution()
{
    vector<double> wavelengths;
    vector<int> frequencies;

    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);

    int maxFrequency = 0;
    int freqMaxNum = 0;
    int all = 0;

    for (int i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > maxFrequency)
        {
            maxFrequency = frequencies[i];
            freqMaxNum = i;
        }

        all += frequencies[i];
    }

    int totalSquashed = 0.33 * frequencies.size();
    int start = freqMaxNum - totalSquashed / 2;
    int end = freqMaxNum + totalSquashed / 2;

    if (start < 0) start = 0;
    if (end > frequencies.size()) end = frequencies.size();

    int highSum = 0;

    for (int i = start; i < end; i++)
    {
        highSum += frequencies[i];
    }

    double squashedProportion = (double)highSum / (double)all;

    logged << "(" << getImage()->getFilename() << ") squashed proportion: " << squashedProportion << std::endl;
    sendLog();

    return (squashedProportion > 0.6);
}

bool IOMRefiner::isGoodSolution()
{
    if (getImage()->getClass() == ImageClassDIALS)
    {
        logged << "Solution refinement has not been enabled for DIALS images yet - assuming good solution." << std::endl;
        sendLog();

        return true;
    }

    bool good = false;
    double goodSolutionStdev = FileParser::getKey("GOOD_SOLUTION_ST_DEV", 0.04);
    double badSolutionStdev = FileParser::getKey("BAD_SOLUTION_ST_DEV", 0.066);
    double goodSolutionSumRatio = FileParser::getKey("GOOD_SOLUTION_SUM_RATIO", 6.5);
    int goodSolutionHighestPeak = FileParser::getKey("GOOD_SOLUTION_HIGHEST_PEAK", 17);
    int badSolutionHighestPeak = FileParser::getKey("BAD_SOLUTION_HIGHEST_PEAK", 10);
    int minimumReflections = FileParser::getKey("MINIMUM_REFLECTION_CUTOFF", 30);
    std::ostringstream details;

    logged << "(" << getImage()->getFilename() << ") Standard deviation: " << lastStdev << std::endl;
    sendLog(LogLevelNormal);

    vector<double> wavelengths;
    vector<int> frequencies;

    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);


    double totalMean = 0;
    double totalStdev = 0;
    int maxFrequency = 0;

    histogram_gaussian(&wavelengths, &frequencies, totalMean, totalStdev);

    for (int i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > maxFrequency)
        {
            maxFrequency = frequencies[i];
        }
    }

    double diffSquared = 0;
    double worstSquared = 0;
    double maxHeight = super_gaussian(totalMean, totalMean, goodSolutionStdev * 0.8, 2) ;

    for (int i = 0; i < frequencies.size(); i++)
    {
        double length = wavelengths[i];
        double expected = super_gaussian(length, totalMean, goodSolutionStdev * 0.8, 2) * (double)maxFrequency / maxHeight;

    //    logged << "expected: " << expected << ", frequency: " << frequencies[i] << std::endl;

        diffSquared += pow(expected - frequencies[i], 2);
        worstSquared += pow(expected, 2);
    }

    diffSquared = sqrt(diffSquared);
    worstSquared = sqrt(worstSquared);

    std::sort(frequencies.begin(), frequencies.end(), greater());

    double highSum = 0;
    std::vector<double> lowOnly;
    int cutoff = 3;

    for (int i = 0; i < frequencies.size(); i++)
    {
        if (i < cutoff)
        {
            highSum += frequencies[i];
            continue;
        }

        lowOnly.push_back(frequencies[i]);
    }

    highSum /= cutoff;

    double stdevLow = standard_deviation(&lowOnly, NULL, 0);

    logged << "Stdev low: " << stdevLow << ", highAverage " << highSum << std::endl;
    sendLog();

    if (lastStdev < goodSolutionStdev)
    {
        good = true;
        details << "(" << getImage()->getFilename() << ") Standard deviation is sufficiently low (" << lastStdev << " vs " << goodSolutionStdev << ")" << std::endl;
    }

    if (highSum > stdevLow * goodSolutionSumRatio)
    {
        good = true;
        details << "(" << getImage()->getFilename() << ") Sum ratio is sufficiently high (" << highSum << " vs " << stdevLow << ")" << std::endl;
    }

  /*  if (highSum <= 5)
    {
        details << "(" << getImage()->getFilename() << ") However, high sum not high enough (" << highSum << ")" << std::endl;
        good = false;
    }
    */
    if (frequencies[0] > goodSolutionHighestPeak)
    {
        details << "(" << getImage()->getFilename() << ") Highest peak is high enough (" << frequencies[0] << " vs " << goodSolutionHighestPeak << ")" << std::endl;
        good = true;
    }

    if (frequencies[0] < badSolutionHighestPeak)
    {
        details << "(" << getImage()->getFilename() << ") Highest peak is too low (" << frequencies[0] << " vs " << badSolutionHighestPeak << ")" << std::endl;
        good = false;
    }

//    if (lastScore < 3)
 //       good = true;

    if (getTotalReflections() < minimumReflections)
    {
        details << "(" << getImage()->getFilename() << ") However, not enough reflections (" << getTotalReflections() << " vs " << minimumReflections << ")" << std::endl;
        good = false;
    }

    if (lastStdev > badSolutionStdev)
    {
        good = false;
        details << "(" << getImage()->getFilename() << ") However, standard deviation too high (" << lastStdev << " vs " << badSolutionStdev << ")" << std::endl;
    }

    details << "Decision: " << (good ? "keep." : "throw away.") << std::endl;

    Logger::mainLogger->addStream(&details, LogLevelNormal);

    return good;
}

void IOMRefiner::calculateOnce()
{
    calculateNearbyMillers(true);
    checkAllMillers(maxResolution, testBandwidth);

    vector<double> wavelengths;
    vector<int> frequencies;

    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed, 0);
}

void IOMRefiner::showHistogram(bool silent)
{
    vector<double> wavelengths;
    vector<int> frequencies;

    double mean = 0;
    double stdev = 0;
    double theScore = 0;

    getWavelengthHistogram(wavelengths, frequencies, silent ? LogLevelDebug : LogLevelNormal, 0);
    gaussian_fit(wavelengths, frequencies, (int)wavelengths.size(), &mean, &stdev,
                 &theScore, true);
}

MtzPtr IOMRefiner::newMtz(int index, bool silent)
{
    calculateNearbyMillers(true);
    checkAllMillers(maxResolution, testBandwidth);
    vector<double> wavelengths;
    vector<int> frequencies;


    double mean = 0;
    double stdev = 0;
    double theScore = 0;

    getWavelengthHistogram(wavelengths, frequencies, silent ? LogLevelDebug : LogLevelNormal, 0);
    gaussian_fit(wavelengths, frequencies, (int)wavelengths.size(), &mean, &stdev,
                 &theScore, true);

    bool complexShoebox = FileParser::getKey("COMPLEX_SHOEBOX", false);

    checkAllMillers(maxResolution, testBandwidth, complexShoebox);

    MatrixPtr newMat = getMatrix()->copy();
    /*
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    newMat->rotate(hRad, kRad, 0);*/

    MtzPtr mtz = MtzPtr(new MtzManager());
    mtz->setWavelength(mean);
    mtz->setFilename(getImage()->filenameRoot() + "_" + i_to_str(index) + ".mtz");
    mtz->setSpaceGroup(spaceGroup->spg_num);
    mtz->setUnitCell(unitCell);

    mtz->setMatrix(newMat);
    double distance = getImage()->getDetectorDistance();
    mtz->setDetectorDistance(distance);

    char *hallSymbol = ccp4spg_symbol_Hall(spaceGroup);

    space_group _spaceGroup = space_group(hallSymbol);
    space_group_type spgType = space_group_type(_spaceGroup);
    asu asymmetricUnit = asu(spgType);

    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        if (miller->getRawestIntensity() > 60000)
            continue;

        miller->incrementOverlapMask();
        miller->setMtzParent(&*mtz);

        int index = Reflection::indexForReflection(miller->getH(), miller->getK(), miller->getL(),
                                               mtz->getLowGroup(), false);

        Reflection *found = NULL;
        mtz->findReflectionWithId(index, &found);

        Panel::addMillerToPanelArray(miller);

        if (found != NULL)
        {
            found->addMiller(miller);
            miller->setParent(found);
        }
        else
        {
            Reflection *reflection = new Reflection();
            reflection->setSpaceGroup(spaceGroup->spg_num);
            reflection->addMiller(miller);
            reflection->setUnitCellDouble(&unitCell[0]);
            reflection->calculateResolution(&*mtz);
            miller->setParent(reflection);
            mtz->addReflection(reflection);
            mtz->sortLastReflection();
        }
    }

    this->sendLog(LogLevelDetailed);

    double cutoff = FileParser::getKey("SIGMA_RESOLUTION_CUTOFF", SIGMA_RESOLUTION_CUTOFF);

    if (cutoff != 0)
        mtz->cutToResolutionWithSigma(cutoff);

    std::string imgFilename = "img-" + getImage()->filenameRoot() + "_" + i_to_str(index) + ".mtz";
    mtz->writeToFile(imgFilename, true, true, false, true);
    mtz->writeToDat();

    nearbyMillers.clear();
    lastMtz = mtz;

    return mtz;
}

IOMRefiner::~IOMRefiner()
{
//    std::cout << "Deallocating IOMRefiner." << std::endl;

    lastMtz = MtzPtr();

    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);

    millers.clear();
    vector<MillerPtr>().swap(millers);

    if (spaceGroup != NULL)
        ccp4spg_free(&spaceGroup);
}

void IOMRefiner::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

std::string IOMRefiner::refinementSummaryHeader()
{
    return "Filename\tRefl num\tScore\tHoriz rot\tVert rot\tBeam rot\tStdev\tWavelength\tDistance\tUnitCellA\tUnitCellB\tUnitCellC\tAlpha\tBeta\tGamma";
}

std::string IOMRefiner::refinementSummary()
{
    std::string filename = getImage()->getFilename();
    int totalReflections = getTotalReflections();
    double lastScore = getLastScore();
    double wavelength = getWavelength();
    double distance = getDetectorDistance();
    MatrixPtr matrix = getMatrix();
    double *lengths = new double[3];

    for (int i = 0; i < 3; i++)
    {
        lengths[i] = 0;
    }

    matrix->unitCellLengths(&lengths);

    std::ostringstream summary;

    summary << filename << "\t" << totalReflections << "\t" << lastScore << "\t"
    << bestHRot << "\t" << bestKRot << "\t" << bestLRot << "\t" << lastStdev << "\t" << wavelength << "\t" << distance << "\t" << lengths[0] << "\t" << lengths[1] << "\t" << lengths[2] << "\t" << unitCell[3] << "\t" << unitCell[4] << "\t" << unitCell[5];

    delete [] lengths;

    return summary.str();
}
