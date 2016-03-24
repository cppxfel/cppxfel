//
//  GetterSetterMap.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__GetterSetterMap__
#define __cppxfel__GetterSetterMap__

#include <stdio.h>
#include "parameters.h"

typedef enum
{
    GetterSetterStepSearch,
    GetterSetterNelderMead
} GetterSetterRefinementType;

class GetterSetterMap
{
private:
    Getter evaluationFunction;
    int maxCycles;
    void *evaluateObject;
    
    std::vector<void *> objects;
    std::vector<Getter> getters;
    std::vector<Setter> setters;
    std::vector<double> stepSizes;
    std::vector<double> stepConvergences;
    double minimizeParameter(int i);
    
public:
    GetterSetterMap()
    {
        evaluationFunction = NULL;
        maxCycles = 30;
    };
    
    void refine(GetterSetterRefinementType type);
    
    void addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence);
    void setEvaluationFunction(Getter function, void *evaluatedObject)
    {
        evaluationFunction = function;
        evaluateObject = evaluatedObject;
    }
    
    void setCycles(int num)
    {
        maxCycles = num;
    }
};

#endif /* defined(__cppxfel__GetterSetterMap__) */
