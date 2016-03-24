//
//  FreeLattice.cpp
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "FreeLattice.h"
#include "parameters.h"
#include <math.h>
#include "SpotVector.h"
#include "Vector.h"
#include "FileParser.h"
#include "CSV.h"

FreeLattice::FreeLattice()
{
    
}

void FreeLattice::addExpanded()
{
    double reciprocalTolerance = FileParser::getKey("RECIPROCAL_TOLERANCE", 0.0015);
    calculateExpandedVectors(false);
    
    for (int i = 0; i < expandedSpotVectors.size() - 1; i++)
    {
        for (int j = i + 1; j < expandedSpotVectors.size(); j++)
        {
            SpotVectorPtr iSpot = expandedSpotVectors[i];
            SpotVectorPtr jSpot = expandedSpotVectors[j];
            
            if (iSpot->isCloseToSpotVector(jSpot, reciprocalTolerance))
            {
                expandedSpotVectors.erase(expandedSpotVectors.begin() + j);
                j--;
            }
        }
    }
    
    spotVectors = expandedSpotVectors;
    expandedSpotVectors.clear();
}

void FreeLattice::calculateExpandedVectors(bool originOnly)
{
    double maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.075) / 2;
    expandedSpotVectors.clear();
    
    std::vector<SpotVectorPtr> origin;
    origin.push_back(SpotVectorPtr(new SpotVector(new_vector(0, 0, 0), new_vector(0, 0, 0))));
    
    std::vector<SpotVectorPtr> *firstVecs = &origin;
    int arrayCount = (originOnly ? 1 : 2);
    
    for (int array = 0; array < arrayCount; array++)
    {
        if (array == 1)
        {
            firstVecs = &spotVectors;
        }
        
        for (int i = 0; i < firstVecs->size(); i++)
        {
            if ((*firstVecs)[i]->distance() > maxDistance)
                continue;
            
            for (int j = 0; j < spotVectors.size(); j++)
            {
                if (spotVectors[j]->distance() > maxDistance)
                    continue;

                SpotVectorPtr firstVec = (*firstVecs)[i];
                SpotVectorPtr secondVec = spotVectors[j];
                
                SpotVectorPtr difference = firstVec->differenceFromVector(secondVec);
                
                expandedSpotVectors.push_back(difference);
            }
        }
    }
}

void FreeLattice::powderPattern(bool originOnly, std::string filename)
{
    double powderPatternStep = FileParser::getKey("POWDER_PATTERN_STEP", 0.001);
    std::vector<double> allDistances;
    
    calculateExpandedVectors(originOnly);
    
    for (int i = 0; i < expandedSpotVectors.size(); i++)
    {
        allDistances.push_back(expandedSpotVectors[i]->distance());
    }
    
    std::map<double, int> map = histogram(allDistances, powderPatternStep);
    CSV csv = CSV(2, "distance", "frequency");
    csv.histogram(map);
    csv.plotColumns(0, 1);
    csv.writeToFile(filename);
}

void FreeLattice::anglePattern(bool originOnly)
{
    std::vector<double> probeDistances = FileParser::getKey("PROBE_DISTANCES", std::vector<double>());
    double distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 1000.);
    double angleStep = FileParser::getKey("POWDER_PATTERN_STEP_ANGLE", 2.0);
    
    if (probeDistances.size() < 2)
        return;
    
    double distance1 = probeDistances[0];
    double distance2 = probeDistances[1];
    
    std::vector<double> allAngles;
    
    calculateExpandedVectors(originOnly);
    
    for (int i = 0; i < expandedSpotVectors.size(); i++)
    {
        double iDiff = fabs(distance1 - expandedSpotVectors[i]->distance());
        if (1 / iDiff < distanceTolerance)
            continue;

        for (int j = 0; j < expandedSpotVectors.size(); j++)
        {
            double jDiff = fabs(distance2 - expandedSpotVectors[j]->distance());
            if (1 / jDiff < distanceTolerance)
                continue;
            
            if (i == j) continue;
            
            double angle = expandedSpotVectors[i]->angleWithVector(expandedSpotVectors[j]) * 180 / M_PI;
            
            if (angle != angle)
                continue;
            
            double otherAngle = 180 - angle;
            

            std::cout << angle << ", " << otherAngle << std::endl;
            
            allAngles.push_back(angle);
            allAngles.push_back(otherAngle);
        }
    }
    
    std::map<double, int> map = histogram(allAngles, angleStep);
    CSV csv = CSV(2, "angle", "frequency");
    csv.histogram(map);
    csv.writeToFile("freeLatticeAngle.csv");

}

FreeLattice::FreeLattice(double a, double b, double c, double alpha, double beta, double gamma)
{
    startingAngles(a, b, c, alpha, beta, gamma);
}

// in degrees
void FreeLattice::startingAngles(double a, double b, double c, double alpha, double beta, double gamma)
{
    double radAlpha = alpha * M_PI / 180;
    double radBeta = beta * M_PI / 180;
    double radGamma = gamma * M_PI / 180;
    
    double a0, a1, a2;
    double b0, b1, b2;
    double c0, c1, c2;
    
    a0 = a;
    a1 = 0;
    a2 = 0;
    b2 = 0;
    
    b0 = b * cos(radGamma);
    c0 = c * cos(radBeta);
    b1 = b * sin(radGamma);
    
    c1 = b * c * pow(cos(radAlpha), 2);
    c1 -= b * cos(radGamma) * c * cos(radBeta);
    c1 /= b * sin(radGamma);
    
    c2 = pow(c, 2);
    c2 -= pow(c1, 2);
    c2 -= pow(c0, 2);
    c2 = sqrt(c2);
    
    vec origin = new_vector(0, 0, 0);
    vec spotA1 = new_vector(a0, a1, a2);
    vec spotA2 = reverseVector(spotA1);

    vec spotB1 = new_vector(b0, b1, b2);
    vec spotB2 = reverseVector(spotB1);

    vec spotC1 = new_vector(c0, c1, c2);
    vec spotC2 = reverseVector(spotC1);

    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotA1, origin)));
    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotA2, origin)));
    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotB1, origin)));
    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotB2, origin)));
    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotC1, origin)));
    spotVectors.push_back(SpotVectorPtr(new SpotVector(spotC2, origin)));
}
