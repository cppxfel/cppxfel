//
//  SpotVector.h
//  cppxfel
//
//  Created by Helen Ginn on 15/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SpotVector__
#define __cppxfel__SpotVector__

#include <stdio.h>
#include "Spot.h"
#include "parameters.h"

class SpotVector
{
private:
    SpotPtr firstSpot;
    SpotPtr secondSpot;
    bool update;
    double cachedDistance;
    double approxResolution;
    double minDistanceTolerance;

    std::vector<SpotVectorPtr> sameLengthStandardVectors;
    vec hkl;
    vec spotDiff;
public:

    static bool isGreaterThan(SpotVectorPtr spotVec1, SpotVectorPtr spotVec2)
    {
        return (spotVec2->distance() < spotVec1->distance());
    }

    SpotVector(SpotPtr first, SpotPtr second);
    SpotVector(vec transformedHKL, vec normalHKL);

    bool hasCommonSpotWithVector(SpotVectorPtr spotVector2);
    double distance();
    void calculateDistance();
    double getResolution();
    double getMinDistanceTolerance();
    double angleWithVertical();
    double angleWithVector(SpotVectorPtr spotVector2);
    double angleWithVector(SpotVectorPtr spotVector2, MatrixPtr mat);
    double similarityToSpotVector(SpotVectorPtr spotVector2);
    void projectedXYDisplacement(double *x, double *y);
    bool isCloseToSpotVector(SpotVectorPtr spotVector2, double maxDistance);
    double trustComparedToStandardVector(SpotVectorPtr standardVector);
    double distanceDifference(SpotVectorPtr standardVector);
    SpotVectorPtr copy();
    SpotVectorPtr vectorRotatedByMatrix(MatrixPtr mat);
    std::string description();
    void addSimilarLengthStandardVectors(std::vector<SpotVectorPtr> standardVectors, double tolerance);
    double cosineWithVector(SpotVectorPtr spotVector2);
    SpotVectorPtr differenceFromVector(SpotVectorPtr spotVec);
    static SpotVectorPtr vectorBetweenSpotsFromArray(std::vector<SpotVectorPtr> vectors, SpotPtr spot1, SpotPtr spot2);

    std::vector<SpotVectorPtr> standardVectorsOfSameDistance()
    {
        return sameLengthStandardVectors;
    }

    void setUpdate()
    {
        update = true;
    }

    vec getVector()
    {
        return spotDiff;
    }

    SpotPtr getFirstSpot()
    {
        return firstSpot;
    }

    SpotPtr getSecondSpot()
    {
        return secondSpot;
    }

    vec getSpotDiff()
    {
        return spotDiff;
    }

    double getH()
    {
        return hkl.h;
    }

    double getK()
    {
        return hkl.k;
    }

    double getL()
    {
        return hkl.l;
    }
};

#endif /* defined(__cppxfel__SpotVector__) */
