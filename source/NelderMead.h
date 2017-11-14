//
//  NelderMead.h
//  cppxfel
//
//  Created by Helen Ginn on 07/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__NelderMead__
#define __cppxfel__NelderMead__

#include <stdio.h>
#include "parameters.h"
#include "Logger.h"

typedef std::pair<std::vector<double>, double> TestPoint;

class NelderMead
{
private:
    double alpha;
    double gamma;
    double rho;
    double sigma;
    bool unlimited;

    std::vector<double *> paramPtrs;
    std::vector<TestPoint> testPoints;
    std::ostringstream logged;
    double (*evaluationFunction)(void *object);
    void *evaluatingObject;

    void setWorstTestPoint(TestPoint &newPoint);
    TestPoint *worstTestPoint();
    void orderTestPoints();
    void evaluateTestPoint(int num);
    void evaluateTestPoint(TestPoint *testPoint);
    void setTestPointParameters(TestPoint *testPoint);
    std::vector<double> calculateCentroid();

    TestPoint reflectOrExpand(std::vector<double> centroid, double scale);
    TestPoint reflectedPoint(std::vector<double> centroid);
    TestPoint expandedPoint(std::vector<double> centroid);
    TestPoint contractedPoint(std::vector<double> centroid);
    void reduction();

    void addPoints(std::vector<double> *point, std::vector<double> pointToAdd);
    void scalePoint(std::vector<double> *point, double scale);
    void subtractPoints(std::vector<double> *point, std::vector<double> pointToSubtract);
public:
    NelderMead(std::vector<double *> newParamPtrs, std::vector<double> expectedRanges, void *object, double (*score)(void *object));
    void process();

    void setUnlimited(bool newLimited)
    {
        unlimited = newLimited;
    }

    int paramCount()
    {
        return (int)paramPtrs.size();
    }

    void sendLog(LogLevel logLevel = LogLevelDetailed);
};

#endif /* defined(__cppxfel__NelderMead__) */
