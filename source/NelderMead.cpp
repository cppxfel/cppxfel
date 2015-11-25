//
//  NelderMead.cpp
//  cppxfel
//
//  Created by Helen Ginn on 07/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "NelderMead.h"
#include "FileParser.h"

bool converged()
{
    return false;
}

void NelderMead::addPoints(std::vector<double> *point, std::vector<double> pointToAdd)
{
    for (int i = 0; i < point->size(); i++)
    {
        (*point)[i] += pointToAdd[i];
    }
}

void NelderMead::subtractPoints(std::vector<double> *point, std::vector<double> pointToSubtract)
{
    for (int i = 0; i < point->size(); i++)
    {
        (*point)[i] -= pointToSubtract[i];
    }
}

void NelderMead::scalePoint(std::vector<double> *point, double scale)
{
    for (int i = 0; i < point->size(); i++)
    {
        (*point)[i] *= scale;
    }
}

void NelderMead::setWorstTestPoint(TestPoint &newPoint)
{
    testPoints[testPoints.size() - 1] = newPoint;
}

TestPoint *NelderMead::worstTestPoint()
{
    orderTestPoints();
    return &testPoints[testPoints.size() - 1];
}

std::vector<double> NelderMead::calculateCentroid()
{
    std::vector<double> centroid;
    centroid.resize(paramCount());
    orderTestPoints();
    
    for (int i = 0; i < paramCount(); i++)
    {
        double total = 0;
        
        for (int j = 0; j < testPoints.size() - 1; j++)
        {
            total += testPoints[j].first[i];
        }
        
        total /= testPoints.size() - 1;
        
        centroid[i] = total;
    }
    
    logged << "Lowest test point: " << testPoints[0].second << std::endl;
    sendLog(LogLevelDebug);
    
    return centroid;
}

TestPoint NelderMead::reflectOrExpand(std::vector<double> centroid, double scale)
{
    TestPoint *maxPoint = worstTestPoint();
    
    std::vector<double> diffVec = centroid;
    subtractPoints(&diffVec, maxPoint->first);
    scalePoint(&diffVec, scale);
    std::vector<double> reflectedVec = centroid;
    addPoints(&reflectedVec, diffVec);
    
    TestPoint reflection = std::make_pair(reflectedVec, 0);
    evaluateTestPoint(&reflection);
    
    return reflection;
}

TestPoint NelderMead::reflectedPoint(std::vector<double> centroid)
{
    return reflectOrExpand(centroid, alpha);
}

TestPoint NelderMead::expandedPoint(std::vector<double> centroid)
{
    return reflectOrExpand(centroid, gamma);
}

TestPoint NelderMead::contractedPoint(std::vector<double> centroid)
{
    return reflectOrExpand(centroid, rho);
}

void NelderMead::reduction()
{
    TestPoint bestPoint = testPoints[0];
    
    for (int i = 1; i < testPoints.size(); i++)
    {
        TestPoint point = testPoints[i];
        
        std::vector<double> diffVec = point.first;
        subtractPoints(&diffVec, bestPoint.first);
        scalePoint(&diffVec, sigma);
        std::vector<double> finalVec = bestPoint.first;
        addPoints(&finalVec, diffVec);
        
        TestPoint contractPoint = std::make_pair(finalVec, 0);
        evaluateTestPoint(&contractPoint);
        
        testPoints[i] = contractPoint;
    }
}

static bool testPointWorseThanTestPoint(TestPoint one, TestPoint two)
{
    return one.second < two.second;
}

void NelderMead::orderTestPoints()
{
    std::sort(testPoints.begin(), testPoints.end(), testPointWorseThanTestPoint);
}

void NelderMead::evaluateTestPoint(int num)
{
    evaluateTestPoint(&testPoints[num]);
}

void NelderMead::evaluateTestPoint(TestPoint *testPoint)
{
    setTestPointParameters(testPoint);
    double eval = evaluationFunction(evaluatingObject);
    testPoint->second = eval;
}

void NelderMead::setTestPointParameters(TestPoint *testPoint)
{
    for (int i = 0; i < paramCount(); i++)
    {
        double *ptr = paramPtrs[i];
        
        if (ptr)
            *ptr = testPoint->first[i];
    }
}

void NelderMead::process()
{
    int count = 0;
    
    int maxCount = FileParser::getKey("NELDER_MEAD_CYCLES", 100);
    
    while ((!converged() && count < maxCount) || unlimited)
    {
        sendLog(LogLevelDebug);
        std::vector<double> centroid = calculateCentroid();
        count++;

        logged << "Evaluation of best point: " << testPoints[0].second << std::endl;
        sendLog(LogLevelDebug);

        TestPoint reflected = reflectedPoint(centroid);
        
        if (reflected.second < testPoints[1].second)
        {
            logged << "Reflecting" << std::endl;
            setWorstTestPoint(reflected);
            continue;
        }
        
        if (reflected.second < testPoints[0].second)
        {
            TestPoint expanded = expandedPoint(centroid);
            bool expandedBetter = (expanded.second < reflected.second);
            setWorstTestPoint(expandedBetter ? expanded : reflected);
            
            logged << (expandedBetter ? "Expanding" : "Reflecting") << std::endl;
            
            continue;
        }
        
        TestPoint contracted = contractedPoint(centroid);
        TestPoint *worstPoint = worstTestPoint();
        
        if (contracted.second < worstPoint->second)
        {
            logged << "Contracting" << std::endl;
            setWorstTestPoint(contracted);
            continue;
        }
        else
        {
            logged << "Reducing" << std::endl;
            reduction();
        }
    }
    
    orderTestPoints();
    setTestPointParameters(&testPoints[0]);
    
    logged << "Evaluation of best point: " << testPoints[0].second << std::endl;
    sendLog(LogLevelDetailed);

}

NelderMead::NelderMead(std::vector<double *> newParamPtrs, std::vector<double> expectedRanges, void *object, double (*score)(void *object))
{
    alpha = 1;
    gamma = 2;
    rho = -0.5;
    sigma = 0.5;
    unlimited = false;
    
    std::vector<double *>streamlinedPtrs;
    
    for (int i = 0; i < newParamPtrs.size(); i++)
    {
        if (newParamPtrs[i] != NULL)
            streamlinedPtrs.push_back(newParamPtrs[i]);
        else
        {
            expectedRanges.erase(expectedRanges.begin() + i);
            newParamPtrs.erase(newParamPtrs.begin() + i);
            i--;
        }
    }
    
    paramPtrs = streamlinedPtrs;
    evaluationFunction = score;
    evaluatingObject = object;
    
    int testPointCount = (int)paramCount() + 1;
    testPoints.resize(testPointCount);
    
    assert(paramCount() > 1);
    
    for (int i = 0; i < testPoints.size(); i++)
    {
        testPoints[i].second = 0;
        testPoints[i].first.resize(paramCount());
        
        logged << "Test point parameter " << i << ": ";
        
        for (int j = 0; j < paramCount(); j++)
        {
            if (i == 0)
            {
                testPoints[i].first[j] = *paramPtrs[j];
            }
            
            if (i > 0)
            {
                int minJ = i - 1;
                double scale = 2;
                
                testPoints[i].first[j] = testPoints[0].first[j] + (j == minJ) * scale * expectedRanges[j];
            }
            
            logged << testPoints[i].first[j] << ", " << std::endl;
        }
    }
   
    logged << "Test point evaluations: ";
    
    for (int i = 0; i < testPoints.size(); i++)
    {
        evaluateTestPoint(i);
        logged << testPoints[i].second << ", ";
    }
    
    logged << std::endl;
    
    std::vector<double> centroid = calculateCentroid();
    logged << "Starting centroid: " << std::endl;
    
    for (int i = 0; i < centroid.size(); i++)
    {
        logged << centroid[i] << ", ";
    }
    
    logged << std::endl;
    
    sendLog(LogLevelDebug);
    
    
}

void NelderMead::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

