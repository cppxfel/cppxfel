//
//  FreeLattice.h
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__FreeLattice__
#define __cppxfel__FreeLattice__

#include <stdio.h>
#include <vector>
#include "parameters.h"

class FreeLattice
{
protected:
    std::vector<SpotVectorPtr> spotVectors;
    std::vector<SpotVectorPtr> expandedSpotVectors;
    
    void calculateExpandedVectors(bool originOnly);
    void startingAngles(double a, double b, double c, double alpha, double beta, double gamma);
public:
    FreeLattice();
    FreeLattice(double a, double b, double c, double alpha, double beta, double gamma);

    void addExpanded();
    void powderPattern(bool originOnly = false, std::string filename = "latticePowder.csv");
    void anglePattern(bool originOnly = false);
};

#endif /* defined(__cppxfel__FreeLattice__) */
