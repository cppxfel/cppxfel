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
    
    vec hkl;
    vec spotDiff;
public:
    
    SpotVector(SpotPtr first, SpotPtr second);
    SpotVector(vec transformedHKL, vec normalHKL);
    
    bool hasCommonSpotWithVector(SpotVectorPtr spotVector2);
    double distance();
    double angleWithVertical();
    double angleWithVector(SpotVectorPtr spotVector2, MatrixPtr mat = MatrixPtr());
    double similarityToSpotVector(SpotVectorPtr spotVector2);
    void projectedXYDisplacement(double *x, double *y);
    bool isCloseToSpotVector(SpotVectorPtr spotVector2, double maxDistance);
    double trustComparedToStandardVector(SpotVectorPtr standardVector);
    SpotVectorPtr copy();
    SpotVectorPtr vectorRotatedByMatrix(MatrixPtr mat);
    
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
