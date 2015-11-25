//
//  AmbiguityBreaker.h
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__AmbiguityBreaker__
#define __cppxfel__AmbiguityBreaker__

#include "MtzManager.h"
#include "parameters.h"
#include <vector>

class StatisticsManager;

class AmbiguityBreaker
{
private:
    vector<MtzPtr> mtzs;
    int ambiguityCount;
    StatisticsManager *statsManager;
    double gridCorrelation(int imageNumI, int imageNumJ);
    double evaluation();
    double gradientForImage(int imageNum, int axis);
    
    double distance(int vectorNum, int centreNum);
    double toggleValue(int slowCloud, int fastCloud);
    double evaluationCloudCluster();
    double gradientCloudCluster(int centre, int axis);
    
    double dotProduct(int imageNumI, int imageNumJ);
    
    scitbx::af::shared<double> x;
    scitbx::af::shared<double> clouds;
    
    void assignPartialities();
    void breakAmbiguity();
    void makeCorrelationGrid();
    void printResults();
    void split();
    void merge();
    MtzManager *merged;
    
public:
    void setMtzs(vector<MtzPtr> newMtzs);
    
    AmbiguityBreaker(vector<MtzPtr> newMtzs);
    void run();
    void overrideAmbiguity(int newAmbiguity);
    
    MtzManager *getMergedMtz()
    {
        return merged;
    }
};

#endif /* defined(__cppxfel__AmbiguityBreaker__) */
