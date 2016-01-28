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

std::vector<SpotVectorPtr> IndexingSolution::standardVectors;
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
std::vector<MatrixPtr> IndexingSolution::symOperators;
Reflection *IndexingSolution::newReflection;
bool IndexingSolution::notSetup = true;
bool IndexingSolution::finishedSetup = false;


void IndexingSolution::setupStandardVectors()
{
    notSetup = false;
    
    distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
    distanceToleranceReciprocal = 1 / distanceTolerance;
    angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
    solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 8.0) * M_PI / 180;
    
    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);
    
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup);
    
    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2],
                                              unitCell[3], unitCell[4], unitCell[5]);
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }
    
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 4);
    double maxDistance = 0;
    
    
    for (int i = -maxMillerIndexTrial; i <= maxMillerIndexTrial; i++)
    {
        for (int j = -maxMillerIndexTrial; j <= maxMillerIndexTrial; j++)
        {
            for (int k = -maxMillerIndexTrial; k <= maxMillerIndexTrial; k++)
            {
                if (spaceGroupNum != 19 && spaceGroupNum != 178)
                    if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                        continue;
                
                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);
                
                if (length_of_vector(hkl) > maxMillerIndexTrial)
                    continue;
                
                unitCellMatrix->multiplyVector(&hkl_transformed);
                
                double distance = length_of_vector(hkl_transformed);
                
                if (distance > maxDistance)
                    maxDistance = distance;
                
                SpotVectorPtr newStandardVector = SpotVectorPtr(new SpotVector(hkl_transformed, hkl));
                
                standardVectors.push_back(newStandardVector);
            }
        }
    }
    
    newReflection = new Reflection();
    newReflection->setUnitCellDouble(&unitCell[0]);
    newReflection->setSpaceGroup(spaceGroupNum);
    
    finishedSetup = true;
}

bool IndexingSolution::matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2, bool force)
{
    std::ostringstream logged;
    
    for (int k = 0; k < symOperators.size(); k++)
    {
        MatrixPtr symOperator = symOperators[k];
        
        for (int l = 0; l < newReflection->ambiguityCount(); l++)
        {
            
            MatrixPtr mat3 = mat2->copy();
            MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
            mat3->getRotation()->preMultiply(*symOperator);
            mat3->getRotation()->preMultiply(*ambiguity);
            
            double angle = mat1->similarityToRotationMatrix(mat3, solutionAngleSpread * 2.5, force);
            
            double theta, phi, psi;
            mat3->eulerAngles(&theta, &phi, &psi);
            
            double theta2, phi2, psi2;
            mat1->eulerAngles(&theta2, &phi2, &psi2);
            
            logged << "Similarity angle: " << angle << " for angles (" << theta << ", " << phi << ", " << psi << ") and (" << theta2 << ", " << phi2 << ", " << psi2 << ")" << std::endl;
            Logger::mainLogger->addStream(&logged, LogLevelDebug);
            
            if (angle < solutionAngleSpread * 2.5 && angle != -1)
            {
                return true;
            }
        }
    }
    
    return false;
}

bool IndexingSolution::vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2)
{
    double firstVectorTrust = firstObserved->trustComparedToStandardVector(standard1);
    double secondVectorTrust = secondObserved->trustComparedToStandardVector(standard2);
    
    if (firstVectorTrust < distanceTolerance)
        return false;
    
    if (secondVectorTrust < distanceTolerance)
        return false;
    
    double standardAngle = standard1->angleWithVector(standard2);
    double realAngle = firstObserved->angleWithVector(secondObserved);
    
    double difference = fabs(realAngle - standardAngle);
    
    if (difference > angleTolerance)
        return false;
    
    return true;
}

bool IndexingSolution::allVectorMatches(SpotVectorPtr firstVector, SpotVectorPtr secondVector, std::vector<SpotVectorPtr> *firstMatches, std::vector<SpotVectorPtr> *secondMatches)
{
    bool good = false;
    
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(standardVectors[i]);
        
        if (firstVectorTrust > distanceTolerance)
        {
            for (int j = 0; j < standardVectors.size(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVectors[j]);
                
                if (secondVectorTrust > distanceTolerance)
                {
                    double realAngle = firstVector->angleWithVector(secondVector);
                    double expectedAngle = standardVectors[i]->angleWithVector(standardVectors[j]);
                    
                    double difference = fabs(realAngle - expectedAngle);
                    
                    if (difference < angleTolerance)
                    {
                        firstMatches->push_back(standardVectors[i]);
                        secondMatches->push_back(standardVectors[j]);
                        
                        good = true;
                    }
                }
            }
        }
    }
    
    return false;
}

bool IndexingSolution::vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr *firstMatch, SpotVectorPtr *secondMatch)
{
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(standardVectors[i]);
        
        if (firstVectorTrust > distanceTolerance)
        {
            for (int j = 0; j < standardVectors.size(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVectors[j]);
                
                if (secondVectorTrust > distanceTolerance)
                {
                    double realAngle = firstVector->angleWithVector(secondVector);
                    double expectedAngle = standardVectors[i]->angleWithVector(standardVectors[j]);
                    
                    double difference = fabs(realAngle - expectedAngle);
                    
                    if (difference < angleTolerance)
                    {
                        *firstMatch = standardVectors[i];
                        *secondMatch = standardVectors[j];
                        
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
    
    for (int i = 0; i < allSpotVectors.size() - 1; i++)
    {
        for (int j = i + 1; j < allSpotVectors.size(); j++)
        {
            MatrixPtr mat = createSolution(allSpotVectors[i], allSpotVectors[j]);
            matrices.push_back(mat);
            
            double theta, phi, psi;
            mat->eulerAngles(&theta, &phi, &psi);
        }
    }
    
    sendLog(LogLevelDetailed);

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
        
        logged << matrices[j]->summary() << "\t" << score << "\t" << count << std::endl;
        
        scoredSolutions.push_back(std::make_pair(aMat, score));
    }

    sendLog(LogLevelDetailed);

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
        logged << "Chosen solution:\t" << chosenMat->summary() << "\tfrom " << spotVectors.size() << " vectors"<< std::endl;
        sendLog(LogLevelDebug);
        
        logged << chosenMat->description() << std::endl;
        sendLog(LogLevelNormal);
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
    fullMat->setComplexMatrix(unitCellOnly->copy(), rotateFinalMatrix);
    
    // Send back the goods.
    return fullMat;
}

bool IndexingSolution::solutionCompatibleForMerge(IndexingSolutionPtr otherSolution)
{
    if (spotVectors.size() >= 3 && otherSolution->spotVectors.size() >= 3)
        return true;
    
    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        for (SpotVectorMap::iterator it2 = otherSolution->spotVectors.begin(); it2 != otherSolution->spotVectors.end(); it2++)
        {
            SpotVectorPtr myVector = it->first;
            SpotVectorPtr herVector = it2->first;
            
            if (myVector->hasCommonSpotWithVector(herVector))
            {
                return true;
            }
        }
    }
    
    return false;
}

bool IndexingSolution::spotVectorHasAnAppropriateDistance(SpotVectorPtr observedVector)
{
    double myDistance = observedVector->distance();
    
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double distance = standardVectors[i]->distance();
        
        if (fabs(distance - myDistance) < 1 / distanceTolerance)
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
        SpotVectorPtr herVector = (*possibleVectors)[i];
        bool commonSpots = false;
        bool checkingCommonSpots = FileParser::getKey("CHECKING_COMMON_SPOTS", true);
        bool duplicates = false;
        
        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr myVector = it->first;
            
            if (checkingCommonSpots && myVector->hasCommonSpotWithVector(herVector))
            {
                commonSpots = true;
            }
            
            if (myVector == herVector)
            {
                duplicates = true;
            }
        }
        
        if (checkingCommonSpots && !commonSpots)
            continue;
        
        if (duplicates)
            continue;
        
        sendLog(LogLevelDebug);
        
        // Find a standard vector which works with all the other vectors
        for (int j = 0; j < standardVectors.size(); j++)
        {
            bool agrees = vectorAgreesWithExistingVectors(herVector, standardVectors[j]);
            
            if (!agrees)
                continue;
            
            agrees = vectorSolutionsAreCompatible(herVector, standardVectors[j]);
            
            if (agrees)
            {
                added++;
               
                addVectorToList(herVector, standardVectors[j]);
                
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
  //  distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
  //  angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
  //  solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 8.0) * M_PI / 180;
    
    spotVectors = firstMap;
    
    for (SpotVectorMap::iterator it = secondMap.begin(); it != secondMap.end(); it++)
    {
        SpotVectorPtr myVector = it->first;
        spotVectors[myVector] = it->second->vectorRotatedByMatrix(symOperator);
    }
    
    matrices = matrixMap1;
    /*
    for (SpotVectorMatrixMap2D::iterator it = matrixMap2.begin(); it != matrixMap2.end(); it++)
    {
        SpotVectorPtr firstVec = it->first;
        SpotVectorMatrixMap internalMap = it->second;
        
        for (SpotVectorMatrixMap::iterator it2 = internalMap.begin(); it2 != internalMap.end(); it2++)
        {
            SpotVectorPtr secondVec = it2->first;
            MatrixPtr matrix = it2->second;
            
            matrices[firstVec][secondVec] = matrix;
        }
    }*/
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

IndexingSolution::~IndexingSolution()
{
    spotVectors.clear();
    SpotVectorMap().swap(spotVectors);
}
