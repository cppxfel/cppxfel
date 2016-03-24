//
//  UnitCellLattice.h
//  cppxfel
//
//  Created by Helen Ginn on 10/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__UnitCellLattice__
#define __cppxfel__UnitCellLattice__

#include <stdio.h>
#include "FreeLattice.h"
#include "csymlib.h"
#include "Matrix.h"
#include "parameters.h"
#include <scitbx/vec3.h>

using scitbx::vec3;

class UnitCellLattice : public FreeLattice
{
private:
    CSym::CCP4SPG *spaceGroup;
    MatrixPtr unitCellMatrix;
    MatrixPtr unitCellOnly;
    MatrixPtr unitCellMatrixInverse;
    static std::vector<MatrixPtr> symOperators;
    int maxMillerIndexTrial;
    double maxDistance;
    double minDistance;
    std::vector<vec3<int> > integerVectors;
    void getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax);

public:
    void setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution = 0);
    
    int standardVectorCount()
    {
        return (int)spotVectors.size();
    }
    
    vec3<int> intVector(int i)
    {
        return (vec3<int>)integerVectors[i];
    }
    
    SpotVectorPtr standardVector(int i)
    {
        return spotVectors[i];
    }
    
    int symOperatorCount()
    {
        return (int)symOperators.size();
    }
    
    MatrixPtr symOperator(int i)
    {
        return symOperators[i];
    }
    
    double getMaxDistance()
    {
        return maxDistance;
    }
    
    MatrixPtr getUnitCellOnly()
    {
        return unitCellOnly;
    }
    
    double getMinDistance()
    {
        return minDistance;
    }
    
    std::vector<SpotVectorPtr> getStandardVectors()
    {
        return spotVectors;
    }
    
    UnitCellLattice(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution = 0) : FreeLattice()
    {
        setup(a, b, c, alpha, beta, gamma, spaceGroupNum, resolution);
    }
};

#endif /* defined(__cppxfel__UnitCellLattice__) */
