//
//  UnitCellLattice.cpp
//  cppxfel
//
//  Created by Helen Ginn on 10/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "UnitCellLattice.h"
#include "csymlib.h"
#include "FileParser.h"
#include "SpotVector.h"
#include "Logger.h"

std::vector<MatrixPtr> UnitCellLattice::symOperators;

void UnitCellLattice::getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax)
{
    double *lengths = new double[3];

    unitCellMatrix->unitCellLengths(&lengths);

    *hMax = lengths[0] / resolution;
    *kMax = lengths[1] / resolution;
    *lMax = lengths[2] / resolution;

    delete [] lengths;
}

void UnitCellLattice::setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution)
{
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);

    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup, a, b, c, alpha, beta, gamma);
    /*
    Logger::mainLogger->addString("\nSymmetry operators:");

    for (int i = 0; i < symOperatorCount(); i++)
    {
        symOperator(i)->printDescription();
    }

    Logger::mainLogger->addString("\n");
    */
    unitCellOnly = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);

    MatrixPtr rotationMat = MatrixPtr(new Matrix());

    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);

    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();

    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 0);
    maxDistance = 0;

    int maxMillerIndexTrialH, maxMillerIndexTrialK, maxMillerIndexTrialL;

    if (resolution == 0 && maxMillerIndexTrial == 0)
    {
        resolution = 1 / FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    }

    if (resolution == 0)
    {
        maxMillerIndexTrialH = maxMillerIndexTrial;
        maxMillerIndexTrialK = maxMillerIndexTrial;
        maxMillerIndexTrialL = maxMillerIndexTrial;
    }
    else
    {
        getMaxMillerIndicesForResolution(resolution, &maxMillerIndexTrialH, &maxMillerIndexTrialK, &maxMillerIndexTrialL);
    }

    int count = 0;

    for (int i = -maxMillerIndexTrialH; i <= maxMillerIndexTrialH; i++)
    {
        for (int j = -maxMillerIndexTrialK; j <= maxMillerIndexTrialK; j++)
        {
            for (int k = -maxMillerIndexTrialL; k <= maxMillerIndexTrialL; k++)
            {
                if (spaceGroupNum != 19 && spaceGroupNum != 178)
                    if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                        continue;

                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);

                if (resolution == 0 && length_of_vector(hkl) > maxMillerIndexTrial)
                    continue;

                unitCellMatrix->multiplyVector(&hkl_transformed);

                double distance = length_of_vector(hkl_transformed);

                if (distance > maxDistance)
                    maxDistance = distance;

                SpotVectorPtr newStandardVector = SpotVectorPtr(new SpotVector(hkl_transformed, hkl));

                vec3<int> integer = vec3<int>(i, j, k);
                integerVectors.push_back(integer);

                spotVectors.push_back(newStandardVector);
                count++;
            }
        }
    }

    std::ostringstream logged;
    logged << "Added " << count << " test vectors from space group / unit cell." << std::endl;
    Logger::mainLogger->addStream(&logged);

    minDistance = FLT_MAX;

    for (int i = 0; i < 3; i++)
    {
        vec hkl = new_vector((i == 0), (i == 1), (i == 2));
        unitCellMatrix->multiplyVector(&hkl);

        if (length_of_vector(hkl) < minDistance)
            minDistance = length_of_vector(hkl);
    }
}
