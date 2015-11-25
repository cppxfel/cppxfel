//
//  SpotVector.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotVector.h"
#include "Matrix.h"

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
    
    hkl = normalHKL;
    spotDiff = copy_vector(transformedHKL);
}

SpotVector::SpotVector(SpotPtr first, SpotPtr second)
{
    firstSpot = first;
    secondSpot = second;
    
    if (!first || !second)
        return;
    
    vec firstVector = first->estimatedVector();
    vec secondVector = second->estimatedVector();
    
    spotDiff = copy_vector(secondVector);
    take_vector_away_from_vector(firstVector, &spotDiff);
}

double SpotVector::distance()
{
    return length_of_vector(spotDiff);
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
