//
//  Scaler.h
//  cppxfel
//
//  Created by Helen Ginn on 01/04/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Scaler__
#define __cppxfel__Scaler__

#include <stdio.h>
#include <vector>
#include "MtzManager.h"


typedef std::map<MtzPtr, vector<Reflection *> > MtzDataMap;
typedef std::map<MtzPtr, vector<double *> > ParameterMap;
typedef std::map<MtzPtr, vector<double> > StepMap;

class Scaler
{
private:
    MtzDataMap mtzData;
    ParameterMap parameterMap;
    StepMap stepMap;
    MtzPtr groupedMtz;
    double evaluate();
    double gradientForImageParameter(MtzPtr mtz, int paramNum);
    double hessianGradientForParamType(int paramNum);
    bool isRefiningParameter(int paramNum);
    double stepForParam(int paramNum);
    MtzPtr pointerForMtz(MtzManager *pointer);
    
    double xStepNorm(double step, int paramNum);
    double gNorm(int paramNum);
    double hessianGradientForParam(MtzPtr mtz, int paramNum);
    void calculateDiagonals();
    void calculateGradients();
    void loadParametersFromMtzs();
    void loadParametersIntoMtzs();
    
    int paramsPerImage();
    int parameterCount();
    scitbx::af::shared<double> g;
    scitbx::af::shared<double> diag;
    scitbx::af::shared<double> x;
    static Scaler *publicScaler;
public:
    
    static Scaler *getScaler(vector<MtzPtr> mtzs = vector<MtzPtr>(), MtzManager **grouped = NULL);
    
    Scaler(vector<MtzPtr> mtzs, MtzManager **grouped);
    void minimizeRMerge();
    void minimizeRMergeLBFGS();
    void minimizeRMergeNelderMead();
    
    bool mtzIsBeneficial(MtzPtr mtz);
    double evaluateForImage(MtzManager *mtz);
    double evaluateForImage(MtzPtr mtz);
    static double evaluateStatic(void *object);
};

#endif /* defined(__cppxfel__Scaler__) */
