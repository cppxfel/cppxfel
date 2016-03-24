//
//  SpotVector.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotVector.h"
#include "Matrix.h"
#include "misc.h"
#include "FileParser.h"

double SpotVector::distanceDifference(SpotVectorPtr standardVector)
{
    double standardDistance = standardVector->distance();
    double diff = fabs(standardDistance - distance());
    
    return diff;
}

double SpotVector::trustComparedToStandardVector(SpotVectorPtr standardVector)
{
    double standardDistance = standardVector->distance();
    double diff = fabs(standardDistance - distance());
    
    return 1 / diff;
}

SpotVector::SpotVector(vec transformedHKL, vec normalHKL)
{
    firstSpot = SpotPtr();
    secondSpot = SpotPtr();
    approxResolution = 0;
    minDistanceTolerance = 0;
  
    update = false;
    hkl = normalHKL;
    spotDiff = copy_vector(transformedHKL);
    cachedDistance = length_of_vector(spotDiff);
}

SpotVector::SpotVector(SpotPtr first, SpotPtr second)
{
    firstSpot = first;
    secondSpot = second;
    update = false;
    approxResolution = 0;
    minDistanceTolerance = 0;
    
    if (!first || !second)
        return;
    
    calculateDistance();
}

void SpotVector::calculateDistance()
{
    vec firstVector = firstSpot->estimatedVector();
    vec secondVector = secondSpot->estimatedVector();
    
    spotDiff = copy_vector(secondVector);
    take_vector_away_from_vector(firstVector, &spotDiff);
    cachedDistance = length_of_vector(spotDiff);
}

double SpotVector::distance()
{
    if (update)
    {
        calculateDistance();
        update = false;
    }
    
    return cachedDistance;
}

// in radians
double SpotVector::angleWithVector(SpotVectorPtr spotVector2)
{
    return angleBetweenVectors(spotVector2->spotDiff, spotDiff);
}

double SpotVector::cosineWithVector(SpotVectorPtr spotVector2)
{
    return cosineBetweenVectors(spotVector2->spotDiff, spotDiff);
}

double SpotVector::angleWithVector(SpotVectorPtr spotVector2, MatrixPtr mat)
{
    vec otherDiff = copy_vector(spotVector2->spotDiff);
    
    if (mat)
    {
        mat->multiplyVector(&otherDiff);
    }
    
    return angleBetweenVectors(spotDiff, otherDiff);
}

double SpotVector::angleWithVertical()
{
    vec vertical = new_vector(0, 1, 0);
    
    double angle = angleBetweenVectors(spotDiff, vertical);
    
    if (spotDiff.h < 0)
        angle = -angle;
    
    return angle;
}

void SpotVector::projectedXYDisplacement(double *x, double *y)
{
    double length = distance();
    
    double angle = angleWithVertical();
    
    *x = length * sin(angle);
    *y = length * cos(angle);
}

bool SpotVector::isCloseToSpotVector(SpotVectorPtr spotVector2, double maxDistance)
{
    vec spotDiff2 = spotVector2->getSpotDiff();
    
    if (fabs(spotDiff2.h - spotDiff.h) > maxDistance || fabs(spotDiff2.k - spotDiff.k) > maxDistance
        || fabs(spotDiff2.l - spotDiff.l) > maxDistance)
        return false;
    
    return true;
}

double SpotVector::similarityToSpotVector(SpotVectorPtr spotVector2)
{
    vec displace = copy_vector(spotDiff);
    take_vector_away_from_vector(spotVector2->getSpotDiff(), &displace);
    
    return length_of_vector(displace);
}

SpotVectorPtr SpotVector::copy()
{
    SpotVectorPtr newPtr = SpotVectorPtr(new SpotVector(firstSpot, secondSpot));
    newPtr->hkl = copy_vector(hkl);
    newPtr->spotDiff = copy_vector(spotDiff);
    newPtr->update = update;
    newPtr->sameLengthStandardVectors = sameLengthStandardVectors;
    
    return newPtr;
}

SpotVectorPtr SpotVector::vectorRotatedByMatrix(MatrixPtr mat)
{
    SpotVectorPtr newVec = copy();
    
    mat->multiplyVector(&newVec->hkl);
    mat->multiplyVector(&newVec->spotDiff);
    
    return newVec;
}

bool SpotVector::hasCommonSpotWithVector(SpotVectorPtr spotVector2)
{
    if (firstSpot == spotVector2->firstSpot || firstSpot == spotVector2->secondSpot)
    {
        return true;
    }
    
    if (secondSpot == spotVector2->firstSpot || secondSpot == spotVector2->secondSpot)
    {
        return true;
    }
    
    return false;
}

std::string SpotVector::description()
{
    return "(" + f_to_str(spotDiff.h) + ", " + f_to_str(spotDiff.k) + ", " + f_to_str(spotDiff.l) + ")";
}

void SpotVector::addSimilarLengthStandardVectors(std::vector<SpotVectorPtr> standardVectors, double tolerance)
{
    sameLengthStandardVectors.clear();
    tolerance = this->getMinDistanceTolerance();
    
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double trust = trustComparedToStandardVector(standardVectors[i]);
        
        if (trust > tolerance)
        {
            sameLengthStandardVectors.push_back(standardVectors[i]);
        }
    }
    
    std::ostringstream logged;
    logged << "Added " << sameLengthStandardVectors.size() << " similar standard lengths to vector (min tolerance " << getMinDistanceTolerance() << ")." << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

SpotVectorPtr SpotVector::differenceFromVector(SpotVectorPtr spotVec)
{
    vec spotDiffDiff = copy_vector(spotDiff);
    
    take_vector_away_from_vector(spotVec->spotDiff, &spotDiffDiff);
    
    SpotVectorPtr newVec = SpotVectorPtr(new SpotVector(spotDiffDiff, new_vector(0, 0, 0)));
    
    return newVec;
}

SpotVectorPtr SpotVector::vectorBetweenSpotsFromArray(std::vector<SpotVectorPtr> vectors, SpotPtr spot1, SpotPtr spot2)
{
    for (int i = 0; i < vectors.size(); i++)
    {
        SpotPtr firstSpot = vectors[i]->getFirstSpot();
        
        if (firstSpot == spot1 || firstSpot == spot2)
        {
            SpotPtr secondSpot = vectors[i]->getSecondSpot();
            
            if ((secondSpot == spot1 || secondSpot == spot2) && secondSpot != firstSpot)
            {
                return vectors[i];
            }
        }
    }
    
    return SpotVectorPtr();
}

double SpotVector::getResolution()
{
    if (approxResolution != 0)
        return approxResolution;
    
    double resol1 = firstSpot->resolution();
    double resol2 = secondSpot->resolution();
    
    double minResolution = std::max(resol1, resol2);
    
    approxResolution = minResolution;
    
    return approxResolution;
}

double SpotVector::getMinDistanceTolerance()
{
    if (minDistanceTolerance != 0)
        return minDistanceTolerance;
    
    double resolution = getResolution();
    
    double rlpSize = FileParser::getKey("INITIAL_RLP_SIZE", 0.0001);
    //double mosaicity = FileParser::getKey("INITIAL_MOSAICITY", 0.0);
    double bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", 0.0013);
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.);
    
    double minWavelength = wavelength * (1 - bandwidth * 2);
    double maxWavelength = wavelength * (1 + bandwidth * 2);
    
    // we set k = 0
    
    double minRadius = 1 / minWavelength;
    double minL = - pow(resolution, 2) / (2 * minRadius);
    double minH = sqrt(pow(resolution, 2) - pow(minL, 2));
    
    double maxRadius = 1 / maxWavelength;
    double maxL = - pow(resolution, 2) / (2 * maxRadius);
    double maxH = sqrt(pow(resolution, 2) - pow(maxL, 2));
    
    vec minVec = new_vector(minH, 0, minL);
    vec maxVec = new_vector(maxH, 0, maxL);
    
    take_vector_away_from_vector(minVec, &maxVec);
    
    minDistanceTolerance = length_of_vector(maxVec);
    minDistanceTolerance += rlpSize;
    
    minDistanceTolerance = 1 / minDistanceTolerance;
    
    return minDistanceTolerance;
}