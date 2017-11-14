//
//  IndexingSolution.cpp
//  cppxfel
//
//  Created by Helen Ginn on 18/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//


#include "IndexingSolution.h"
#include "FileParser.h"
#include "Logger.h"
#include <string>
#include "misc.h"
#include <iomanip>
#include "UnitCellLattice.h"
#include "Holder.h"

using cctbx::sgtbx::space_group;
using cctbx::sgtbx::rt_mx;
using cctbx::sgtbx::rot_mx;
using cctbx::uctbx::unit_cell;

int IndexingSolution::spaceGroupNum = 0;
CSym::CCP4SPG *IndexingSolution::spaceGroup = NULL;
std::vector<double> IndexingSolution::unitCell;
int IndexingSolution::maxMillerIndexTrial = 4;
MatrixPtr IndexingSolution::unitCellOnly;
MatrixPtr IndexingSolution::unitCellMatrix;
MatrixPtr IndexingSolution::unitCellMatrixInverse;
double IndexingSolution::distanceTolerance = 4000;
double IndexingSolution::distanceToleranceReciprocal = 0.0025;
double IndexingSolution::angleTolerance = 1.0 * M_PI / 180;
double IndexingSolution::solutionAngleSpread = 8.0 * M_PI / 180;
double IndexingSolution::approximateCosineDelta = 0.012;
UnitCellLatticePtr IndexingSolution::lattice;
Reflection *IndexingSolution::newReflection;
bool IndexingSolution::notSetup = true;
bool IndexingSolution::finishedSetup = false;
bool IndexingSolution::checkingCommonSpots = true;
std::mutex IndexingSolution::setupMutex;

void IndexingSolution::reset()
{
    notSetup = true;
    setupStandardVectors();
}

void IndexingSolution::setupStandardVectors()
{
    if (finishedSetup) return;

    setupMutex.lock();

    if (finishedSetup)
    {
        setupMutex.unlock();
        return;
    }

    notSetup = false;

    distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
    distanceToleranceReciprocal = 1 / distanceTolerance;
    angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
    solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 8.0) * M_PI / 180;
    checkingCommonSpots = FileParser::getKey("CHECKING_COMMON_SPOTS", true);

    double fortyFiveAngle = M_PI / 4;
    double slightlySmallerAngle = fortyFiveAngle - angleTolerance;

    approximateCosineDelta = fabs(cos(fortyFiveAngle) - cos(slightlySmallerAngle));

    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);

    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);

    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());

    std::vector<double> testCell = FileParser::getKey("RECIPROCAL_UNIT_CELL", std::vector<double>());

    if (testCell.size() == 6)
    {
        unitCell = Matrix::unitCellFromReciprocalUnitCell(testCell[0], testCell[1], testCell[2],
                                                          testCell[3], testCell[4], testCell[5]);
        std::ostringstream stream;

        stream << "Unit cell from reciprocal: ";

        for (int i = 0; i < 6; i++)
        {
            stream << unitCell[i] << " ";
        }

        stream << std::endl;
        Logger::mainLogger->addStream(&stream);
    }

    if (unitCell.size() < 6 && testCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL or RECIPROCAL_UNIT_CELL." << std::endl;
        exit(1);
    }

    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 4);

    lattice = UnitCellLatticePtr(new UnitCellLattice(unitCell[0], unitCell[1], unitCell[2],
                                                     unitCell[3], unitCell[4], unitCell[5], spaceGroupNum));

    newReflection = new Reflection();
    newReflection->setUnitCellDouble(&unitCell[0]);
    newReflection->setSpaceGroup(spaceGroupNum);

    finishedSetup = true;

    setupMutex.unlock();
}

void IndexingSolution::calculateSimilarStandardVectorsForImageVectors(std::vector<SpotVectorPtr> vectors)
{
    for (int i = 0; i < vectors.size(); i++)
    {
        vectors[i]->addSimilarLengthStandardVectors(lattice->getStandardVectors(), distanceTolerance);
    }
}

bool IndexingSolution::matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2, bool force)
{
    double minTrace = FLT_MAX;

    std::ostringstream logged;


    for (int k = 0; k < symOperatorCount(); k++)
    {
        MatrixPtr symOp = symOperator(k);

        for (int l = 0; l < newReflection->ambiguityCount(); l++)
        {
            MatrixPtr mat3 = mat2->copy();
            MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
            mat3->getRotation()->preMultiply(*symOp);
            mat3->getRotation()->preMultiply(*ambiguity);

            MatrixPtr subtractedMat = mat1->getRotation()->copy();
            subtractedMat->subtract(mat3->getRotation());
            MatrixPtr transposedMat = subtractedMat->copy()->transpose();
            transposedMat->multiply(*subtractedMat);

            double trace = transposedMat->trace();

            if (trace < minTrace)
                minTrace = trace;
        }
    }

    double badMinTrace = sqrt(4 * (1 - cos(solutionAngleSpread)));

    return (minTrace < badMinTrace);
}

bool IndexingSolution::vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2)
{
    double standardCos = standard1->cosineWithVector(standard2);
    double realCos = firstObserved->cosineWithVector(secondObserved);

    double difference = fabs(realCos - standardCos);

    return (difference < approximateCosineDelta);
    /*
    double standardAngle = standard1->angleWithVector(standard2);
    double realAngle = firstObserved->angleWithVector(secondObserved);

    double difference = fabs(realAngle - standardAngle);

    return (difference < angleTolerance);*/
}

bool IndexingSolution::vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr *firstMatch, SpotVectorPtr *secondMatch)
{
    for (int i = 0; i < standardVectorCount(); i++)
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(standardVector(i));
        double firstTolerance = firstVector->getMinDistanceTolerance();

        if (firstVectorTrust > firstTolerance)
        {
            for (int j = 0; j < standardVectorCount(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVector(j));
                double secondTolerance = secondVector->getMinDistanceTolerance();

                if (secondVectorTrust > secondTolerance)
                {
                    double realAngle = firstVector->angleWithVector(secondVector);
                    double expectedAngle = standardVector(i)->angleWithVector(standardVector(j));

                    double difference = fabs(realAngle - expectedAngle);

                    if (difference < angleTolerance)
                    {
                        *firstMatch = standardVector(i);
                        *secondMatch = standardVector(j);

                        return true;
                    }
                }
            }
        }
    }

    return false;
}

bool match_greater_than_match(std::pair<MatrixPtr, double> a, std::pair<MatrixPtr, double> b)
{
    return a.second > b.second;
}

MatrixPtr IndexingSolution::createSolution()
{
    std::vector<SpotVectorPtr> allSpotVectors;
    std::vector<MatrixPtr> matrices;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        allSpotVectors.push_back(it->first);
    }

    for (int i = 0; i < allSpotVectors.size() && i < 1; i++)
    {
        for (int j = i + 1; j < allSpotVectors.size(); j++)
        {
            MatrixPtr mat = createSolution(allSpotVectors[i], allSpotVectors[j]);
            matrices.push_back(mat);

            double theta, phi, psi;
            mat->eulerAngles(&theta, &phi, &psi);
        }
    }

    std::vector<std::pair<MatrixPtr, double> > scoredSolutions;
    for (int j = 0; j < matrices.size(); j++)
    {
        MatrixPtr aMat = matrices[j];
        double score = 0;
        int count = 0;

        for (int k = 0; k < matrices.size(); k++)
        {
            MatrixPtr bMat = matrices[k]->copy();

            double angle = bMat->similarityToRotationMatrix(aMat, solutionAngleSpread);

            if (angle < solutionAngleSpread && angle != -1)
            {
                double addition = pow(solutionAngleSpread - angle, 2);
                score += addition;
                count++;
            }
        }

#ifndef __OPTIMIZE__
        logged << matrices[j]->summary() << "\t" << score << "\t" << count << std::endl;
#endif

        scoredSolutions.push_back(std::make_pair(aMat, score));
    }

#ifndef __OPTIMIZE__
    sendLog(LogLevelDebug);
#endif

    std::sort(scoredSolutions.begin(), scoredSolutions.end(), match_greater_than_match);


    MatrixPtr chosenMat = MatrixPtr(new Matrix());

    if (scoredSolutions.size())
    {
        /*     MatrixPtr rotationMat = Matrix::matrixFromEulerAngles(averageTheta, averagePhi, averagePsi);

         double theta, phi, psi;
         rotationMat->eulerAngles(&theta, &phi, &psi);

         MatrixPtr unitCellMat = unitCellOnly->copy();

         chosenMat->setComplexMatrix(unitCellMat, rotationMat);*/

        chosenMat = scoredSolutions[0].first;
#ifndef __OPTIMIZE__
        logged << "Chosen solution:\t" << chosenMat->summary() << "\tfrom " << spotVectors.size() << " vectors"<< std::endl;
        sendLog(LogLevelDebug);

        logged << chosenMat->description() << std::endl;
        sendLog(LogLevelDebug);
#endif
    }
    return chosenMat;
}

MatrixPtr IndexingSolution::createSolution(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr firstStandard)
{
    if (!firstStandard)
    {
        firstStandard = spotVectors[firstObserved];
    }

    SpotVectorPtr secondStandard = spotVectors[secondObserved];

    // Get our important variables.
    vec observedVec1 = firstObserved->getVector();
    vec observedVec2 = secondObserved->getVector();

    vec simulatedVec1 = firstStandard->getVector();
    vec simulatedVec2 = secondStandard->getVector();

    // Rotate reciprocal space so that the first observed vector lines up with the simulated vector.
    MatrixPtr rotateSpotDiffMatrix = rotation_between_vectors(observedVec1, simulatedVec1);

    // The first simulated vector shall become the new axis to twirl around
    vec firstAxisUnit = copy_vector(simulatedVec1); // checked

    // Make this a unit vector
    scale_vector_to_distance(&firstAxisUnit, 1);

    // Make a copy of the second observed vector
    vec rotatedObservedVec2 = copy_vector(observedVec2); // checked

    // Rotate this vector by first matrix so that we have treated the two observed vectors equally
    rotateSpotDiffMatrix->multiplyVector(&rotatedObservedVec2); // checked

    double resultantAngle = 0;

    // Now we twirl around the firstAxisUnit until the rotated observed vector matches the second simulated vector
    // as closely as possible.
    MatrixPtr secondTwizzleMatrix = closest_rotation_matrix(rotatedObservedVec2, simulatedVec2, firstAxisUnit, &resultantAngle);

    //   logged << "Resultant angle: " << resultantAngle << std::endl;

    // We want to apply the first matrix and then the second matrix, so we multiply these.
    MatrixPtr combinedMatrix = rotateSpotDiffMatrix->copy();
    combinedMatrix->multiply(*secondTwizzleMatrix);

    // But we actually need the inverse rotation to go from "true" coordinates to "crystal" coordinates.
    MatrixPtr rotateFinalMatrix = combinedMatrix->inverse3DMatrix();

    // Rotate because of funny business from cctbx.xfel (rotated matrix is 90º from DIALS).
    rotateFinalMatrix->rotate(0, 0, -M_PI/2);

    // Create the goods.
    MatrixPtr fullMat = MatrixPtr(new Matrix());
    fullMat->setComplexMatrix(lattice->getUnitCellOnly()->copy(), rotateFinalMatrix);

    // Send back the goods.
    return fullMat;
}

bool IndexingSolution::spotVectorHasAnAppropriateDistance(SpotVectorPtr observedVector)
{
    double myDistance = observedVector->distance();
    double myTolerance = observedVector->getMinDistanceTolerance();

    for (int i = 0; i < standardVectorCount(); i++)
    {
        double distance = standardVector(i)->distance();

        if (fabs(distance - myDistance) < 1 / myTolerance)
        {
            return true;
        }
    }

    return false;
}

void IndexingSolution::pruneSpotVectors(std::vector<SpotVectorPtr> *spotVectors)
{
    int count = 0;

    for (int i = 0; i < spotVectors->size(); i++)
    {
        bool appropriate = spotVectorHasAnAppropriateDistance((*spotVectors)[i]);

        if (!appropriate)
        {
            count++;
            spotVectors->erase(spotVectors->begin() + i);
            i--;
        }
    }
}

void IndexingSolution::removeSpotVectors(std::vector<SpotVectorPtr> *spotVectors)
{
    int count = 0;

    for (int i = 0; i < spotVectors->size(); i++)
    {
        bool exists = this->spotVectors.count((*spotVectors)[i]);

        if (exists)
        {
            count++;
            spotVectors->erase(spotVectors->begin() + i);
            i--;
        }
    }

    logged << "Removed " << count << " spot vectors which led to a bad solution." << std::endl;
    sendLog();
}

bool IndexingSolution::spotsAreNotTooClose(SpotVectorPtr observedVector)
{
    SpotPtr spot1 = observedVector->getFirstSpot();
    SpotPtr spot2 = observedVector->getSecondSpot();

    double squareMinDistance = pow(lattice->getMinDistance(), 2);

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotPtr spot3 = vector->getFirstSpot();
        SpotPtr spot4 = vector->getSecondSpot();

        if (spot3->closeToSecondSpot(spot1, squareMinDistance) || spot3->closeToSecondSpot(spot2, squareMinDistance)
            || spot4->closeToSecondSpot(spot1, squareMinDistance) || spot4->closeToSecondSpot(spot2, squareMinDistance))
            return false;
    }

    return true;
}

bool IndexingSolution::vectorAgreesWithExistingVectors(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr myVector = it->first;

        bool similar = vectorPairLooksLikePair(it->first, observedVector, it->second, standardVector);

        if (!similar)
            return false;
    }

    return true;
}

bool IndexingSolution::vectorSolutionsAreCompatible(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    int count = 0;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr myVector = it->first;

        MatrixPtr newSolution = createSolution(observedVector, myVector, standardVector);

        double theta, phi, psi;

        newSolution->eulerAngles(&theta, &phi, &psi);

        double thetaDiff = fabs(averageTheta - theta);
        double phiDiff = fabs(averagePhi - phi);
        double psiDiff = fabs(averagePsi - psi);

        bool similar = ((thetaDiff < solutionAngleSpread) && (phiDiff < solutionAngleSpread) && (psiDiff < solutionAngleSpread));

        //    logged << solutionAngleSpread << ": " << thetaDiff << ", " << phiDiff << ", " << psiDiff << " = " << similar << std::endl;
        //    sendLog();

        if (!similar)
            return false;

        count++;
    }

    return true;
}

void IndexingSolution::addMatrix(SpotVectorPtr observedVector1, SpotVectorPtr observedVector2, MatrixPtr solution)
{
    matrices[observedVector1][observedVector2] = solution;
    matrices[observedVector2][observedVector1] = solution;

    double theta, phi, psi;
    solution->eulerAngles(&theta, &phi, &psi);

    double allThetas = averageTheta * matrixCount;
    double allPhis = averagePhi * matrixCount;
    double allPsis = averagePsi * matrixCount;

    matrixCount++;
    allThetas += theta;
    allPhis += phi;
    allPsis += psi;

    allThetas /= matrixCount;
    allPhis /= matrixCount;
    allPsis /= matrixCount;

    averageTheta = allThetas;
    averagePhi = allPhis;
    averagePsi = allPsis;
}

void IndexingSolution::addVectorToList(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    spotVectors[observedVector] = standardVector;

    for (SpotVectorMap::iterator i = spotVectors.begin(); i != spotVectors.end(); i++)
    {
        if (i->first == observedVector)
            continue;

        MatrixPtr newSolution = createSolution(observedVector, i->first);
        addMatrix(i->first, observedVector, newSolution);
    }
}


int IndexingSolution::extendFromSpotVectors(std::vector<SpotVectorPtr> *possibleVectors, int limit)
{
    int added = 0;

    for (int i = 0; i < possibleVectors->size() && i < 3000; i++)
    {
        SpotVectorPtr possibleVector = (*possibleVectors)[i];
        bool commonSpots = false;
        bool duplicates = false;

        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr myVector = it->first;

            if (checkingCommonSpots && myVector->hasCommonSpotWithVector(possibleVector))
            {
                commonSpots = true;
            }

            if (myVector == possibleVector)
            {
                duplicates = true;
            }
        }

        if (checkingCommonSpots && !commonSpots)
            continue;

        if (duplicates)
            continue;

        // Find a standard vector which works with all the other vectors
        std::vector<SpotVectorPtr> standardSelection = possibleVector->standardVectorsOfSameDistance();

        for (int j = 0; j < standardSelection.size(); j++)
        {
            bool agrees = vectorAgreesWithExistingVectors(possibleVector, standardSelection[j]);

            if (!agrees)
                continue;

            agrees = vectorSolutionsAreCompatible(possibleVector, standardSelection[j]);

            if (!agrees)
            {
                continue;
            }

            agrees = spotsAreNotTooClose(possibleVector);

            if (!agrees)
            {
                logged << "Too close!" << std::endl;
                sendLog(LogLevelDetailed);

                continue;
            }

            if (agrees)
            {
                added++;

                addVectorToList(possibleVector, standardSelection[j]);

                possibleVectors->erase(possibleVectors->begin() + i);
                i--;

                if (added == limit && limit > 0)
                    return added;
            }
        }
    }

    return added;
}

IndexingSolution::IndexingSolution(SpotVectorMap firstMap, SpotVectorMap secondMap, SpotVectorMatrixMap2D matrixMap1, SpotVectorMatrixMap2D matrixMap2, MatrixPtr symOperator)
{
    spotVectors = firstMap;

    for (SpotVectorMap::iterator it = secondMap.begin(); it != secondMap.end(); it++)
    {
        SpotVectorPtr myVector = it->first;
        spotVectors[myVector] = it->second->vectorRotatedByMatrix(symOperator);
    }

    matrices = matrixMap1;
}

std::string IndexingSolution::getNetworkPDB()
{
    std::ostringstream pdbLog;
    int count = 0;

    pdbLog << "HETATM";
    pdbLog << std::fixed;
    pdbLog << std::setw(5) << count << "                   ";
    pdbLog << std::setprecision(2) << std::setw(8)  << 0;
    pdbLog << std::setprecision(2) << std::setw(8) << 0;
    pdbLog << std::setprecision(2) << std::setw(8) << 0;
    pdbLog << "                       N" << std::endl;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        count++;
        SpotVectorPtr standard = it->second;

        pdbLog << "HETATM";
        pdbLog << std::fixed;
        pdbLog << std::setw(5) << count << "                   ";
        pdbLog << std::setprecision(2) << std::setw(8)  << standard->getH();
        pdbLog << std::setprecision(2) << std::setw(8) << standard->getK();
        pdbLog << std::setprecision(2) << std::setw(8) << standard->getL();
        pdbLog << "                       O" << std::endl;
    }

    return pdbLog.str();
}

std::string IndexingSolution::printNetwork()
{
    std::ostringstream printed;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotVectorPtr standard = it->second;

        double x1 = vector->getFirstSpot()->getX();
        double y1 = vector->getFirstSpot()->getY();

        double x2 = vector->getSecondSpot()->getX();
        double y2 = vector->getSecondSpot()->getY();

        double h = standard->getH();
        double k = standard->getK();
        double l = standard->getL();

        printed << "vec\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << h << "\t" << k << "\t" << l << std::endl;
    }

    return printed.str();
}

IndexingSolutionPtr IndexingSolution::copy()
{
    IndexingSolutionPtr newPtr = IndexingSolutionPtr(new IndexingSolution());

    newPtr->distanceTolerance = distanceTolerance;
    newPtr->angleTolerance = angleTolerance;
    newPtr->spotVectors = spotVectors;
    newPtr->matrices = matrices;
    newPtr->averageTheta = averageTheta;
    newPtr->averagePhi = averagePhi;
    newPtr->averagePsi = averagePsi;

    return newPtr;
}

IndexingSolution::IndexingSolution()
{

}

IndexingSolution::IndexingSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstMatch, SpotVectorPtr secondMatch)
{
    averageTheta = 0;
    averagePhi = 0;
    averagePsi = 0;
    matrixCount = 0;

    addVectorToList(firstVector, firstMatch);
    addVectorToList(secondVector, secondMatch);
}

std::vector<IndexingSolutionPtr> IndexingSolution::startingSolutionsForVectors(SpotVectorPtr firstVector, SpotVectorPtr secondVector)
{
    if (notSetup)
    {
        setupStandardVectors();
    }

    if (!firstVector->hasCommonSpotWithVector(secondVector))
        return std::vector<IndexingSolutionPtr>();

    std::vector<IndexingSolutionPtr> solutions;

    SpotVectorPtr firstMatch = SpotVectorPtr();
    SpotVectorPtr secondMatch = SpotVectorPtr();

    vectorMatchesVector(firstVector, secondVector, &firstMatch, &secondMatch);

    if (!(firstMatch && secondMatch))
        return solutions;

    IndexingSolutionPtr newSolution = IndexingSolutionPtr(new IndexingSolution(firstVector, secondVector, firstMatch, secondMatch));
    solutions.push_back(newSolution);

    return solutions;
}

IndexingSolutionPtr IndexingSolution::startingSolutionsForThreeSpots(std::vector<SpotPtr> *spots, std::vector<SpotVectorPtr> *spotVectors)
{
    if (notSetup)
    {
        setupStandardVectors();
    }

    for (int i = 0; i < spots->size() - 2; i++)
    {
        for (int j = i + 1; j < spots->size() - 1; j++)
        {
            for (int k = j + 1; k < spots->size(); k++)
            {
                SpotVectorPtr vector1 = SpotVector::vectorBetweenSpotsFromArray(*spotVectors, spots->at(i), spots->at(j));

                if (!vector1)
                    continue;

                SpotVectorPtr vector2 = SpotVector::vectorBetweenSpotsFromArray(*spotVectors, spots->at(j), spots->at(k));

                if (!vector2)
                    continue;

                SpotVectorPtr vector3 = SpotVector::vectorBetweenSpotsFromArray(*spotVectors, spots->at(k), spots->at(i));

            //    if (!vector3)
            //        continue;

                SpotVectorPtr match12_1 = SpotVectorPtr();
                SpotVectorPtr match12_2 = SpotVectorPtr();
                SpotVectorPtr match23_1 = SpotVectorPtr();
                SpotVectorPtr match23_2 = SpotVectorPtr();
                SpotVectorPtr match31_1 = SpotVectorPtr();
                SpotVectorPtr match31_2 = SpotVectorPtr();

                if ((vectorMatchesVector(vector1, vector2, &match12_1, &match12_2))/* &&
                    (vectorMatchesVector(vector2, vector3, &match23_1, &match23_2)) &&
                    (vectorMatchesVector(vector1, vector2, &match31_1, &match31_2))*/)
                {
                    // there is still some way in which this will still not match

                    IndexingSolutionPtr newSolution = IndexingSolutionPtr(new IndexingSolution(vector1, vector2, match12_1, match12_2));

                   // if (newSolution->vectorSolutionsAreCompatible(vector3, match23_2))
                   // {
                        spots->erase(spots->begin() + i);
                        i--;
                        spots->erase(spots->begin() + j);
                        j--;
                  //      spots->erase(spots->begin() + k);
                  //      k--;

                        return newSolution;
                   // }
                }
            }
        }
    }

    return IndexingSolutionPtr();
}

IndexingSolution::~IndexingSolution()
{
    spotVectors.clear();
    SpotVectorMap().swap(spotVectors);
}

std::vector<double> IndexingSolution::totalDistances()
{
    std::vector<double> distances;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        //   SpotVectorPtr standard = it->second;

        distances.push_back(vector->distance());
    }

    return distances;
}

std::vector<double> IndexingSolution::totalAngles()
{
    std::vector<double> angles;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector1 = it->first;
        SpotVectorPtr standard1 = it->second;

        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr vector2 = it->first;
            SpotVectorPtr standard2 = it->second;

            double angle1 = standard2->angleWithVector(standard1);
            double angle2 = vector2->angleWithVector(vector1);

            double angleDiff = fabs(angle2 - angle1);

            angles.push_back(angleDiff * 180 / M_PI);
        }
    }

    return angles;
}

std::vector<double> IndexingSolution::totalDistanceTrusts()
{
    std::vector<double> trusts;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotVectorPtr standard = it->second;

        double distanceDiff = fabs(vector->distance() - standard->distance());
        double trust = 1 / distanceDiff;

        if (trust > distanceTolerance * 10)
            trust = distanceTolerance * 10;

        trusts.push_back(trust);
    }

    return trusts;
}

double IndexingSolution::getMinDistance()
{
    setupStandardVectors();
    return lattice->getMinDistance();
}
