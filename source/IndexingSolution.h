//
//  IndexingSolution.h
//  cppxfel
//
//  Created by Helen Ginn on 18/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexingSolution__
#define __cppxfel__IndexingSolution__

#include <stdio.h>
#include "parameters.h"
#include "SpotVector.h"
#include "Matrix.h"
#include "csymlib.h"
#include "LoggableObject.h"
#include "Holder.h"
#include <mutex>
#include "UnitCellLattice.h"

typedef std::map<SpotVectorPtr, SpotVectorPtr> SpotVectorMap;
typedef std::map<SpotVectorPtr, MatrixPtr> SpotVectorMatrixMap;
typedef std::map<SpotVectorPtr, SpotVectorMatrixMap > SpotVectorMatrixMap2D;

class IndexingSolution : LoggableObject
{
private:
    static Reflection *newReflection;
    SpotVectorMap spotVectors;
    SpotVectorMatrixMap2D matrices;
    static UnitCellLatticePtr lattice;
    static std::mutex setupMutex;

    bool vectorAgreesWithExistingVectors(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    static bool vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr *firstMatch, SpotVectorPtr *secondMatch);
    MatrixPtr createSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstStandard = SpotVectorPtr());
    bool vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2);
    void addVectorToList(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    void addMatrix(SpotVectorPtr observedVector1, SpotVectorPtr observedVector2, MatrixPtr solution);
    bool vectorSolutionsAreCompatible(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    static bool spotVectorHasAnAppropriateDistance(SpotVectorPtr observedVector);

    static double distanceTolerance;
    static double distanceToleranceReciprocal;
    static double angleTolerance;
    static double solutionAngleSpread;
    static double approximateCosineDelta;
    static bool checkingCommonSpots;
    static int spaceGroupNum;
    static CSym::CCP4SPG *spaceGroup;
    static std::vector<double> unitCell;
    static int maxMillerIndexTrial;
    static MatrixPtr unitCellOnly;
    static MatrixPtr unitCellMatrix;
    static MatrixPtr unitCellMatrixInverse;
    static bool notSetup;
    static bool finishedSetup;

    double averageTheta;
    double averagePhi;
    double averagePsi;
    int matrixCount;

    IndexingSolution();
    IndexingSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstMatch, SpotVectorPtr secondMatch);
    IndexingSolution(SpotVectorMap firstMap, SpotVectorMap secondMap, SpotVectorMatrixMap2D matrixMap1, SpotVectorMatrixMap2D matrixMap2, MatrixPtr symOperator);

public:
    static std::vector<IndexingSolutionPtr> startingSolutionsForVectors(SpotVectorPtr firstVector, SpotVectorPtr secondVector);
    int extendFromSpotVectors(std::vector<SpotVectorPtr> *possibleVectors, int limit = 0);
    MatrixPtr createSolution();
    static bool matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2, bool force = false);
    static void setupStandardVectors();
    std::string getNetworkPDB();
    std::string printNetwork();
    static void pruneSpotVectors(std::vector<SpotVectorPtr> *spotVectors);
    void removeSpotVectors(std::vector<SpotVectorPtr> *spotVectors);
    static void calculateSimilarStandardVectorsForImageVectors(std::vector<SpotVectorPtr> vectors);
    std::vector<double> totalDistances();
    std::vector<double> totalAngles();
    std::vector<double> totalDistanceTrusts();
    bool spotsAreNotTooClose(SpotVectorPtr observedVector);
    static IndexingSolutionPtr startingSolutionsForThreeSpots(std::vector<SpotPtr> *spots, std::vector<SpotVectorPtr> *spotVectors);
    static void reset();

    IndexingSolutionPtr copy();
    ~IndexingSolution();

    static int standardVectorCount()
    {
        return lattice->standardVectorCount();
    }

    static SpotVectorPtr standardVector(int i)
    {
        return lattice->standardVector(i);
    }

    static int symOperatorCount()
    {
        return lattice->symOperatorCount();
    }

    static MatrixPtr symOperator(int i)
    {
        return lattice->symOperator(i);
    }

    int spotVectorCount()
    {
        return (int)spotVectors.size();
    }

    bool wasSuccessful()
    {
        return spotVectors.size();
    }

    static double getMinDistance();
};

#endif /* defined(__cppxfel__IndexingSolution__) */
