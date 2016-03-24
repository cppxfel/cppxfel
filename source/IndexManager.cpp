//
//  IndexManager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "IndexManager.h"
#include "FileParser.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "Logger.h"
#include <algorithm>
#include "parameters.h"
#include <fstream>
#include "SpotVector.h"
#include "IndexingSolution.h"
#include "FreeLattice.h"
#include "CSV.h"
#include "Holder.h"
#include "Miller.h"

IndexManager::IndexManager(std::vector<ImagePtr> newImages)
{
    images = newImages;
    
    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);
    
    if (spaceGroupNum == 0)
    {
        std::cout << "Please provide space group number in SPACE_GROUP" << std::endl;
        exit(1);
    }
    
    spaceGroup = ccp4spg_load_by_ccp4_num(spaceGroupNum);

    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }
    
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
    
    newReflection = new Reflection();
    newReflection->setUnitCellDouble(&unitCell[0]);
    newReflection->setSpaceGroup(spaceGroupNum);
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2],
                                                    unitCell[3], unitCell[4], unitCell[5]);
    
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    minimumTrustAngle = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0);
    minimumTrustDistance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);

    solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 10.0);
    
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 4);
    maxDistance = 0;
    
    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup, unitCell[0], unitCell[1], unitCell[2],
                                           unitCell[3], unitCell[4], unitCell[5]);
    
    logged << "Calculating distances of unit cell " << spaceGroupNum << std::endl;
    sendLog();
    
    lattice = UnitCellLatticePtr(new UnitCellLattice(unitCell[0], unitCell[1], unitCell[2],
                                                     unitCell[3], unitCell[4], unitCell[5], spaceGroupNum));
    
    for (int i = -maxMillerIndexTrial; i <= maxMillerIndexTrial; i++)
    {
        for (int j = -maxMillerIndexTrial; j <= maxMillerIndexTrial; j++)
        {
            for (int k = -maxMillerIndexTrial; k <= maxMillerIndexTrial; k++)
            {
                if (spaceGroupNum != 19 && spaceGroupNum != 178 && spaceGroupNum != 5)
                {
                    if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                    {
                        continue;
                    }
                    
                    if (spaceGroupNum == 5 && (i + j % 2 != 0))
                        continue;
                    
              //      if (spaceGroupNum == 146 && ((-i + j + k) % 3 != 0))
              //          continue;

                }
                
             //   if (k != 0) continue;
                
                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);
                
                if (length_of_vector(hkl) > maxMillerIndexTrial)
                    continue;
                
                unitCellMatrix->multiplyVector(&hkl_transformed);
                
                logged << hkl.h << "\t" << hkl.k << "\t" << hkl.l << "\t" << hkl_transformed.h << "\t" << hkl_transformed.k << "\t" << hkl_transformed.l << std::endl;
                sendLog(LogLevelDebug);
                
                
                double distance = length_of_vector(hkl_transformed);
                
                if (distance > maxDistance)
                    maxDistance = distance;
                
                vectorDistances.push_back(std::make_pair(hkl, distance));
            }
        }
    }

    maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    
    smallestDistance = FLT_MAX;
    
    if (1 / unitCell[0] < smallestDistance)
        smallestDistance = 1 / unitCell[0];
    if (1 / unitCell[1] < smallestDistance)
        smallestDistance = 1 / unitCell[1];
    if (1 / unitCell[2] < smallestDistance)
        smallestDistance = 1 / unitCell[2];
    
    smallestDistance *= 0.85;
    
    for (int i = 0; i < vectorDistances.size(); i++)
    {
        vec transformed = copy_vector(vectorDistances[i].first);
        unitCellMatrix->multiplyVector(&transformed);
        
        logged << vectorDistances[i].first.h << "\t"
         << vectorDistances[i].first.k << "\t"
         << vectorDistances[i].first.l << "\t"
         << transformed.h << "\t"
         << transformed.k << "\t"
         << transformed.l << "\t"
        << vectorDistances[i].second << std::endl;
    }
    
    sendLog(LogLevelDebug);
}

bool match_greater_than_match(Match a, Match b)
{
    return a.second > b.second;
}

bool greater_than_scored_matrix(std::pair<MatrixPtr, std::pair<double, double> > a, std::pair<MatrixPtr, std::pair<double, double> > b)
{
    return a.second.first > b.second.first;
}

bool IndexManager::matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2)
{
    std::ostringstream logged;
    MatrixPtr mat3 = mat2->copy();
    
    for (int k = 0; k < symOperators.size(); k++)
    {
        MatrixPtr symOperator = symOperators[k];
        
        for (int l = 0; l < newReflection->ambiguityCount(); l++)
        {
            MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
            mat3->getRotation()->preMultiply(*symOperator);
            mat3->getRotation()->preMultiply(*ambiguity);
            
            double radianSpread = solutionAngleSpread * M_PI / 180;
            double angle = mat1->similarityToRotationMatrix(mat3, radianSpread, true);
            
            double theta, phi, psi;
            mat3->eulerAngles(&theta, &phi, &psi);
            
            double theta2, phi2, psi2;
            mat1->eulerAngles(&theta2, &phi2, &psi2);
            
         //   logged << "Similarity angle: " << angle << " for angles (" << theta << ", " << phi << ", " << psi << ") and (" << theta2 << ", " << phi2 << ", " << psi2 << ")" << std::endl;
         //   Logger::mainLogger->addStream(&logged);

            
            if (angle < radianSpread && angle != -1)
            {
                return true;
            }
        }
    }
    
    return false;
}

int IndexManager::indexOneImage(ImagePtr image, std::vector<MtzPtr> *mtzSubset)
{
    int successes = 0;
    int spotNum = image->spotCount();
    double angleTolerance = minimumTrustAngle * M_PI / 180;
    double trustTolerance = minimumTrustDistance;
    double finalTolerance = solutionAngleSpread * M_PI / 180;
    bool alwaysAccept = FileParser::getKey("ACCEPT_ALL_SOLUTIONS", false);
    int spotsPerLattice = FileParser::getKey("SPOTS_PER_LATTICE", 100);
    std::ostringstream logged;
    bool refineOrientations = FileParser::getKey("REFINE_ORIENTATIONS", true);
    
    std::vector<Match> possibleMatches;
    std::vector<MatrixPtr> possibleSolutions;
    
    for (int j = 0; j < image->spotVectorCount(); j++)
    {
        SpotVectorPtr dataPoint = image->spotVector(j);
        
        for (int k = 0; k < vectorDistances.size(); k++)
        {
            VectorDistance trial;
            
            trial.first = copy_vector(vectorDistances[k].first);
            trial.second = vectorDistances[k].second;
            
            double trust = 1 / (fabs(dataPoint->distance() - trial.second));
            
            if (trust > trustTolerance)
            {
                std::pair<SpotVectorPtr, VectorDistance> vectorPair = std::make_pair(dataPoint, trial);
                
                possibleMatches.push_back(std::make_pair(vectorPair, trust));
            }
        }
    }
    
    std::sort(possibleMatches.begin(), possibleMatches.end(), match_greater_than_match);
    
    for (int j = 0; j < possibleMatches.size(); j++)
    {
        double trust = possibleMatches[j].second;
        SpotVectorPtr dataPoint = possibleMatches[j].first.first;
        VectorDistance trial = possibleMatches[j].first.second;
        
        logged << trust << "\t" << trial.second << "\t" << dataPoint->distance() << "\t" << trial.first.h << "\t" <<
        trial.first.k << "\t" <<
        trial.first.l << "\t" << std::endl;
    }
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    
    logged << "Total spots for " << image->getFilename() << ": " << image->spotCount() << std::endl;
    logged << "Total vectors: " << image->spotVectorCount() << std::endl;
    Logger::mainLogger->addStream(&logged); logged.str("");
    logged << "Total possible vector matches: " << possibleMatches.size() << std::endl;
    Logger::mainLogger->addStream(&logged); logged.str("");
    
    bool thorough_searching = FileParser::getKey("THOROUGH_SOLUTION_SEARCHING", false);
    
    int maxSearchNumberMatches = FileParser::getKey("MAX_SEARCH_NUMBER_MATCHES", 5000);
    int maxSearchNumberSolutions = FileParser::getKey("MAX_SEARCH_NUMBER_SOLUTIONS", 8000);
    
    
    for (int j = 0; j < possibleMatches.size() && j < maxSearchNumberMatches; j++)
    {
        for (int k = j - 1; k >= 0 && j < possibleMatches.size(); k--)
        {
            std::pair<SpotVectorPtr, VectorDistance> vectorPair1 = possibleMatches[j].first;
            std::pair<SpotVectorPtr, VectorDistance> vectorPair2 = possibleMatches[k].first;
            
            SpotVectorPtr observed1 = vectorPair1.first;
            SpotVectorPtr observed2 = vectorPair2.first;
            
            vec observedVec1 = observed1->getVector();
            vec observedVec2 = observed2->getVector();
            
            if (!vectors_are_equal(observedVec1, observedVec2))
            {
                vec simulatedVec1 = possibleMatches[j].first.second.first;
                vec simulatedVec2 = possibleMatches[k].first.second.first;
                
                double angle1 = fabs(angleBetweenVectors(observedVec1, observedVec2));
                double angle2 = fabs(angleBetweenVectors(simulatedVec1, simulatedVec2));
                
                double angleDiff = fabs(angle1 - angle2);
                
                if (angleDiff < angleTolerance)
                {
                    // rotation to get from observedvec1 to observedvec2.
                    MatrixPtr rotateSpotDiffMatrix = rotation_between_vectors(observedVec1, simulatedVec1);
                    
                    vec firstSpotVec = observed1->getFirstSpot()->estimatedVector();
                    vec reverseSpotVec = reverseVector(firstSpotVec);
                    rotateSpotDiffMatrix->multiplyVector(&reverseSpotVec);
                    //        vec rereversed = reverseVector(reverseSpotVec);
                    
                    // we have now found one axis responsible for the first observation
                    // now we must rotate around that axis (rereversed) until
                    // observedVec2 matches simulatedVec2
                    // first we make a unit vector to rotate around
                    
                    vec firstAxisUnit = copy_vector(simulatedVec1);
                    scale_vector_to_distance(&firstAxisUnit, 1);
                    
                    // we need to rotate observedVec2 by the first rotation matrix
                    vec rotatedObservedVec2 = copy_vector(observedVec2);
                    rotateSpotDiffMatrix->multiplyVector(&rotatedObservedVec2);
                    
                    // and now we can find the appropriate rotation matrix given our rotation axis firstAxisUnit
                    
                    double resultantAngle = 0;
                    
                    MatrixPtr secondTwizzleMatrix = closest_rotation_matrix(rotatedObservedVec2, simulatedVec2, firstAxisUnit, &resultantAngle);
                    
                    if (resultantAngle > angleTolerance)
                        continue;
                    
                    MatrixPtr combinedMatrix = rotateSpotDiffMatrix->copy();
                    combinedMatrix->multiply(*secondTwizzleMatrix);
                    
                    vec possibleMiller = copy_vector(firstSpotVec);
                    combinedMatrix->multiplyVector(&possibleMiller);
                    unitCellMatrixInverse->multiplyVector(&possibleMiller);
                    
                    if (Logger::getPriorityLevel() == LogLevelDetailed || Logger::getPriorityLevel() == LogLevelDebug)
                    {
                        logged << "Relative Millers (" << simulatedVec1.h << ", " << simulatedVec1.k << ", " << simulatedVec1.l << ") and (" << simulatedVec2.h << ", " << simulatedVec2.k << ", " << simulatedVec2.l << ") from trusts " << possibleMatches[j].second << " and " << possibleMatches[k].second << " with angles " << angleDiff * 180 / M_PI << ", " << resultantAngle * 180 / M_PI << std::endl;
                        
                        Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
                    }
                    
                    SpotPtr pair1Spot1 = observed1->getFirstSpot();
                    SpotPtr pair1Spot2 = observed1->getSecondSpot();
                    
                    SpotPtr pair2Spot1 = observed2->getFirstSpot();
                    SpotPtr pair2Spot2 = observed2->getSecondSpot();
                    
                    MatrixPtr rotateFinalMatrix = combinedMatrix->inverse3DMatrix();
                    rotateFinalMatrix->rotate(0, 0, -M_PI/2);
                    
                    MatrixPtr fullMat = MatrixPtr(new Matrix());
                    fullMat->setComplexMatrix(unitCellOnly, rotateFinalMatrix);
                    
                    possibleSolutions.push_back(fullMat);
                    
                    if (!thorough_searching)
                        j++;
                }
            }
        }
    }
    
    // dummy Reflection to give us indexing ambiguity info
    
    MillerPtr miller = MillerPtr(new Miller(NULL, 0, 0, 0));
    
    logged << possibleSolutions.size() << " total solutions for image " << image->getFilename() << std::endl;
    
    if (possibleSolutions.size() == 0)
        return 0;
    
    std::ostringstream dummyVecStr;
    
    std::vector<std::pair<MatrixPtr, std::pair<double, double> > > scoredSolutions;
    
    Logger::mainLogger->addStream(&logged); logged.str("");
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    double countSum = 0;
    int solutionsSearched = 0;
    
    for (int j = 0; j < possibleSolutions.size() && j < maxSearchNumberSolutions - 1; j++)
    {
        MatrixPtr aMat = possibleSolutions[j];
      //  vec aDummyVec = new_vector(1, 0, 0);
      //  aMat->multiplyVector(&aDummyVec);
        double score = 0;
        int count = 0;
        solutionsSearched ++;
        
        for (int k = 0; k < possibleSolutions.size() && k < maxSearchNumberSolutions; k++)
        {
            MatrixPtr bMat = possibleSolutions[k]->copy();
            
            double angle = bMat->similarityToRotationMatrix(aMat, finalTolerance);
            //   double angle = angleBetweenVectors(aMat, bMat);
            
            if (angle < finalTolerance && angle != -1)
            {
                double addition = pow(finalTolerance - angle, 2);
                score += addition;
                count++;
                countSum++;
            }
        }
        
        logged << j << "\t" << aMat->summary() << "\t" << score << "\t" << count << std::endl;
        
        scoredSolutions.push_back(std::make_pair(aMat, std::make_pair(score, count)));
    }
    
    countSum /= solutionsSearched;
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    Logger::mainLogger->addStream(&dummyVecStr, LogLevelDetailed); dummyVecStr.str("");
    
    std::vector<MatrixPtr> chosenSolutions;
    double expectedLattices = (double)spotNum / (double)spotsPerLattice;
    int trialCount = FileParser::getKey("SOLUTION_ATTEMPTS", 12);
    
    std::sort(scoredSolutions.begin(), scoredSolutions.end(), greater_than_scored_matrix);
    
    for (int j = 0; j < trialCount && j < scoredSolutions.size(); j++)
    {
        MatrixPtr solution = scoredSolutions[j].first;
    //    vec aDummyVec = new_vector(1, 0, 0);
     //   solution->multiplyVector(&aDummyVec);
        
        for (int k = j + 1; k < scoredSolutions.size() && k < j + trialCount; k++)
        {
            MatrixPtr secondSol = scoredSolutions[k].first;
            
    //        logged << solution->description() << std::endl;
    //        logged << secondSol->description() << std::endl;
            
            bool erase = matrixSimilarToMatrix(solution, secondSol);
            
       //     logged << erase << std::endl;
            
            if (erase)
            {
                scoredSolutions.erase(scoredSolutions.begin() + k);
                j--;
                break;
            }
        }
    }
    
    logged << "Scored solution size: " << scoredSolutions.size() << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    
    std::sort(scoredSolutions.begin(), scoredSolutions.end(), greater_than_scored_matrix);
    
    for (int j = 0; j < trialCount && j < scoredSolutions.size(); j++)
    {
        vec aDummyVec = new_vector(1, 0, 0);
        scoredSolutions[j].first->multiplyVector(&aDummyVec);
        double score = scoredSolutions[j].second.first;
        logged << "High score solution:\t" << score << "\t" << aDummyVec.h << "\t" << aDummyVec.k << "\t" << aDummyVec.l << std::endl;
        chosenSolutions.push_back(scoredSolutions[j].first);
    }
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    int minNeighbours = FileParser::getKey("MINIMUM_NEIGHBOURS", 0);
    
    if (scoredSolutions.size() > 0)
    {
        double score = scoredSolutions[0].second.first;
        double count = scoredSolutions[0].second.second;
        logged << "Top solution score for " << image->getFilename() << ": " << score << " from " << count << " solutions. Average solution count: " << countSum << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelNormal); logged.str("");
        
        if (count < minNeighbours)
            return 0;
    }
    
    // weed out duplicates from previous solutions
    
    for (int j = 0; j < image->IOMRefinerCount(); j++)
    {
        MatrixPtr previousMat = image->getIOMRefiner(j)->getMatrix();
        
        for (int m = 0; m < chosenSolutions.size(); m++)
        {
            bool erase = matrixSimilarToMatrix(chosenSolutions[m], previousMat);
            
            if (erase == true)
            {
                chosenSolutions.erase(chosenSolutions.begin() + m);
                j--;
                break;
            }
        }
    }
    
    // weed out duplicates from current possible solutions
    
    for (int j = 0; j < chosenSolutions.size(); j++)
    {
        MatrixPtr chosenSolution = chosenSolutions[j];
        vec aDummyVec = new_vector(1, 0, 0);
        chosenSolution->multiplyVector(&aDummyVec);
        
        for (int m = j + 1; m < chosenSolutions.size(); m++)
        {
            MatrixPtr otherSolution = chosenSolutions[m];
            
            bool erase = matrixSimilarToMatrix(chosenSolution, otherSolution);
            
            if (erase == true)
            {
                chosenSolutions.erase(chosenSolutions.begin() + m);
                j--;
                break;
            }
        }
    }
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
    
    for (int j = 0; j < chosenSolutions.size(); j++)
    {
        MatrixPtr solution = chosenSolutions[j];
        logged << solution->summary() << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
        
    }
    
    logged << "Chosen " << chosenSolutions.size() << " matrices to integrate." << std::endl;
    Logger::mainLogger->addStream(&logged); logged.str("");
    
    Logger::mainLogger->addStream(&logged); logged.str("");
    
    for (int j = 0; j < chosenSolutions.size(); j++)
    {
        image->setUpIOMRefiner(chosenSolutions[j]);
        int lastRefiner = image->IOMRefinerCount() - 1;
        IOMRefinerPtr refiner = image->getIOMRefiner(lastRefiner);
        if (refineOrientations)
        {
            refiner->refineOrientationMatrix();
        }
        else
        {
            refiner->calculateOnce();
        }
        
        bool successfulImage = refiner->isGoodSolution();
        
        if (successfulImage || alwaysAccept)
        {
            logged << "Successful crystal " << j + 1 << "/" << chosenSolutions.size() << " for " << image->getFilename() << std::endl;
            mtzSubset->push_back(refiner->newMtz(lastRefiner));
            Logger::mainLogger->addStream(&logged); logged.str("");
            successes++;
        }
        else
        {
            logged << "Unsuccessful crystal " << j + 1 << "/" << chosenSolutions.size() << " for " << image->getFilename() << std::endl;
            image->removeRefiner(lastRefiner);
            Logger::mainLogger->addStream(&logged); logged.str("");
        }
    }
    
    int removed = image->throwAwayIntegratedSpots(*mtzSubset);
    int spotCount = image->spotCount();
    
    logged << "Removed " << removed << " spots after finding " << successes << " crystals, leaving " << spotCount << " spots." << std::endl;
    
    Logger::mainLogger->addStream(&logged); logged.str("");
    
    return successes;
}

void IndexManager::indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset)
{
    bool newMethod = FileParser::getKey("NEW_INDEXING_METHOD", true);
    int maxThreads = FileParser::getMaxThreads();
    std::ostringstream logged;
    
    if (newMethod)
    {
        while (true)
        {
            ImagePtr image = indexer->getNextImage();
            
            if (!image)
            {
                logged << "Finishing thread " << offset << std::endl;
                Logger::mainLogger->addStream(&logged); logged.str("");
                return;
            }
            
            logged << "Starting image " << image->getFilename() << " on thread " << offset << std::endl;
            Logger::mainLogger->addStream(&logged); logged.str("");
            
            image->findIndexingSolutions();
            
            std::vector<MtzPtr> mtzs = image->getLastMtzs();
            
            mtzSubset->reserve(mtzSubset->size() + mtzs.size());
            mtzSubset->insert(mtzSubset->begin(), mtzs.begin(), mtzs.end());
        }
        
        return;
    }
    
    bool alwaysAccept = FileParser::getKey("ACCEPT_ALL_SOLUTIONS", false);
    bool oneCycleOnly = FileParser::getKey("ONE_INDEXING_CYCLE_ONLY", false);
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);
    
    for (int i = offset; i < indexer->images.size(); i += maxThreads)
    {
        ImagePtr image = indexer->images[i];
        logged << "Starting image " << i << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");

        bool finished = false;
        
        while (!finished)
        {
            int extraSolutions = 1;
            
            while (extraSolutions > 0)
            {
                image->compileDistancesFromSpots(indexer->maxDistance, indexer->smallestDistance, alwaysFilterSpots);
                extraSolutions = indexer->indexOneImage(image, mtzSubset);
                
                if (alwaysAccept)
                    break;
                
                if (oneCycleOnly)
                    break;
            }
            
            if (!alwaysAccept && !oneCycleOnly && !alwaysFilterSpots)
            {
                image->compileDistancesFromSpots(indexer->maxDistance, indexer->smallestDistance, true);
                extraSolutions += indexer->indexOneImage(image, mtzSubset);
            }
            
            if (extraSolutions == 0 || alwaysAccept || oneCycleOnly)
                finished = true;
        }
        
        logged << "N: Finished image " << image->getFilename() << " on " << image->IOMRefinerCount() << " crystals and " << image->spotCount() << " spots." << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");

        image->writeSpotsList();
        image->dropImage();
    }
}

void IndexManager::powderPattern()
{
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);
    
    std::ostringstream pdbLog;
    
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->compileDistancesFromSpots(maxDistance, smallestDistance, alwaysFilterSpots);
    }
    
    refineMetrology();
    
    PowderHistogram frequencies = generatePowderHistogram();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);
            
            vec spotDiff = copy_vector(spotVec->getSpotDiff());
            spotDiff.h *= 106 * 20;
            spotDiff.k *= 106 * 20;
            spotDiff.l *= 106 * 20;

            if (i == 0)
            {
                pdbLog << "HETATM";
                pdbLog << std::fixed;
                pdbLog << std::setw(5) << j << "                   ";
                pdbLog << std::setprecision(2) << std::setw(8)  << spotDiff.h;
                pdbLog << std::setprecision(2) << std::setw(8) << spotDiff.k;
                pdbLog << std::setprecision(2) << std::setw(8) << spotDiff.l;
                pdbLog << "                       O" << std::endl;
            }
        }
    }
    
    if (images.size() == 0)
    {
        logged << "No images specified." << std::endl;
        sendLog();
        return;
    }
    
    pdbLog << "HETATM";
    pdbLog << std::fixed;
    pdbLog << std::setw(5) << images[0]->spotVectorCount() + 1 << "                   ";
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << "                       N" << std::endl;
    
    std::ofstream pdbLog2;
    pdbLog2.open("projection.pdb");
    
    for (int j = 0; j < images[0]->spotCount(); j++)
    {
        vec estimatedVector = images[0]->spot(j)->estimatedVector();
        
        pdbLog2 << "HETATM";
        pdbLog2 << std::fixed;
        pdbLog2 << std::setw(5) << j << "                   ";
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.h * 106;
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.k * 106;
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.l * 106;
        pdbLog2 << "                       O" << std::endl;
    }
    
    pdbLog2 << "HETATM";
    pdbLog2 << std::fixed;
    pdbLog2 << std::setw(5) << images[0]->spotCount() + 1 << "                   ";
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << "                       N" << std::endl;
    
    pdbLog2.close();
    
    logged << "******* DISTANCE FREQUENCY *******" << std::endl;
    
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);

    CSV powder(3, "Distance", "Frequency", "Perfect frequency");
    
    for (PowderHistogram::iterator it = frequencies.begin(); it != frequencies.end(); it++)
    {
        double distance = it->first * step;
        double freq = it->second.first;
        double perfect = it->second.second;
        
        powder.addEntry(0, distance, freq, perfect);
        
//        powderLog << distance << "," << freq << "," << perfect << std::endl;
    }
    
    powder.writeToFile("powder.csv");
    logged << powder.plotColumns(0, 1) << std::endl;
    sendLog();
    
    lattice->powderPattern(false, "unitCellLatticePowder.csv");
    
  //  powderLog.close();
    
    /// angles
    
    std::vector<double> probeDistances = FileParser::getKey("PROBE_DISTANCES", std::vector<double>());
    
    if (probeDistances.size() < 2)
        return;
    
    double distance1 = probeDistances[0]; double distance2 = probeDistances[1];
    
    double angleStep = FileParser::getKey("POWDER_PATTERN_STEP_ANGLE", 2.0) * M_PI / 180;
    double distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 500.);
    std::map<int, int> angleHistogram;
    std::map<int, int> perfectAngleHistogram;
    
    for (double i = 0; i < M_PI; i += angleStep)
    {
        int angleCategory = i / angleStep;
        angleHistogram[angleCategory] = 0;
        perfectAngleHistogram[angleCategory] = 0;
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        std::vector<double> newAngles = images[i]->anglesBetweenVectorDistances(distance1, distance2, distanceTolerance);
        
        for (int j = 0; j < newAngles.size(); j++)
        {
            int angleCategory = newAngles[j] / angleStep;
            angleHistogram[angleCategory]++;
        }
    }
    
    std::vector<VectorDistance> perfectProbe1, perfectProbe2;
    
    for (int i = 0; i < vectorDistances.size(); i++)
    {
        double distance = vectorDistances[i].second;
        bool chosen = false;
        
        if (1 / fabs(probeDistances[0] - distance) > distanceTolerance)
        {
            perfectProbe1.push_back(vectorDistances[i]);
            chosen = true;
        }

        if (1 / fabs(probeDistances[1] - distance) > distanceTolerance)
        {
            perfectProbe2.push_back(vectorDistances[i]);
            chosen = true;
        }
        
        if (chosen)
        {
            logged << "Matching vector: " << vectorDistances[i].first.h << "\t" << vectorDistances[i].first.k << "\t" << vectorDistances[i].first.l << std::endl;
            sendLog();
        }

    }
    
    logged << "Comparing " << perfectProbe1.size() << " against " << perfectProbe2.size() << std::endl;
    sendLog();
    
    for (int i = 0; i < perfectProbe1.size(); i++)
    {
        for (int j = 0; j < perfectProbe2.size(); j++)
        {
            vec probeTransformed1 = copy_vector(perfectProbe1[i].first);
            vec probeTransformed2 = copy_vector(perfectProbe2[j].first);
            
            if (perfectProbe1[i].first.h == perfectProbe2[i].first.h &&
                perfectProbe1[i].first.k == perfectProbe2[i].first.k &&
                perfectProbe1[i].first.l == perfectProbe2[i].first.l)
            {
         //       continue;
            }
            
            unitCellMatrix->multiplyVector(&probeTransformed1);
            unitCellMatrix->multiplyVector(&probeTransformed2);
            
            double angle = angleBetweenVectors(probeTransformed1, probeTransformed2);
            double mirror = M_PI - angle;
            
            int angleCategory = angle / angleStep;
            perfectAngleHistogram[angleCategory]++;
            
            angleCategory = mirror / angleStep;
            perfectAngleHistogram[angleCategory]++;
        }
    }
    
    
    CSV csv(3, "Angle", "Frequency", "Perfect frequency");
    
    for (std::map<int, int>::iterator it = angleHistogram.begin(); it != angleHistogram.end(); it++)
    {
        double angle = it->first * angleStep * 180 / M_PI;
        double freq = it->second;
        double perfectFreq = 0;
        
        if (perfectAngleHistogram.count(it->first))
            perfectFreq = perfectAngleHistogram[it->first];
        
        csv.addEntry(0, angle, freq, perfectFreq);
        
 //       angleLog << angle << "," << freq << "," << perfectFreq << "," << std::endl;
    }
    
    csv.writeToFile("angle.csv");
    csv.plotColumns(0, 1);
    
   // angleLog.close();
}

ImagePtr IndexManager::getNextImage()
{
    std::lock_guard<std::mutex> lock(indexMutex);
    
    nextImage++;
    
    if (nextImage >= images.size())
    {
        return ImagePtr();
    }
    else
    {
        return images[nextImage];
    }
}

bool IndexManager::modifyParameters()
{
    int minimumSolutionNetworkCount = FileParser::getKey("MINIMUM_SOLUTION_NETWORK_COUNT", 20);
    double distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.);
    double maxReciprocalDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    bool finish = false;
    
    int goodSolutions = 0;
    int failedSolutions = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        goodSolutions += images[i]->IOMRefinerCount();
        failedSolutions += images[i]->failedRefinerCount();
    }
    
    double percentIndexed = (double)goodSolutions / (double)images.size();
    double percentRubbish = (double)failedSolutions / (double)goodSolutions;
    
    logged << "Percentage indexed correctly: " << percentIndexed * 100 << "%" << std::endl;
    logged << "Incorrect results: " << percentRubbish * 100 << "%" << std::endl;
    
    sendLog();
    
    if (goodSolutions == 0)
        percentRubbish = 0;
    
    int maxThreads = FileParser::getMaxThreads();
    double timePerImage = lastTime * (double)maxThreads / (double)images.size();
    
    // now we need to decide if we want to lower one of these...
    
    if (percentIndexed < 0.8)
    {
        logged << "Indexing percentage is not yet above 80%. Deciding what to do about it..." << std::endl;
        
        if (percentRubbish > 0.3)
        {
            logged << "The percentage of incorrect results is high. Multiplying the MINIMUM_SOLUTION_VECTOR_NETWORK by 1.2." << std::endl;
            minimumSolutionNetworkCount *= 0.8;
        }
        else if (percentRubbish < 0.1)
        {
            logged << "The percentage of incorrect results is low and so might be missing good solutions as well. Multiplying the MINIMUM_SOLUTION_VECTOR_NETWORK by 0.8." << std::endl;
            minimumSolutionNetworkCount *= 0.8;
        }
        
        if (percentRubbish < 0.3)
        {
            logged << "The percentage of incorrect results is low and perhaps the distance tolerance is too stringent. Multiplying the distance tolerance by 0.9." << std::endl;
            distanceTolerance *= 0.9;
        }
        else if (percentRubbish > 0.6)
        {
            logged << "The percentage of incorrect results is low and perhaps the distance tolerance is not stringent enough. Multiplying the distance tolerance by 1.1." << std::endl;
            distanceTolerance *= 1.1;
        }
        
        if (timePerImage < 45)
        {
            logged << "There are not many solutions but the time taken per image is very low. Perhaps we increase the distance of considered vectors in the lattice. Multiplying MAX_RECIPROCAL_DISTANCE by 1.2 and reducing MINIMUM_SOLUTION_VECTOR_NETWORK by 0.9." << std::endl;
            minimumSolutionNetworkCount *= 0.9;
            maxReciprocalDistance *= 1.2;
        }
    }
    else if (percentIndexed >= 0.8)
    {
        logged << "Indexing percentage is good." << std::endl;
        
        if (timePerImage > 180)
        {
            logged << "Indexing takes longer than 3 mins per image." << std::endl;
            
            if (percentRubbish > 0.3)
            {
                logged << "A lot of time is spent wading through incorrect results. Multiplying the MINIMUM_SOLUTION_VECTOR_NETWORK value by " << 1 + percentRubbish << std::endl;
                minimumSolutionNetworkCount *= 1 + percentRubbish / 2;
            }
            else
            {
                logged << "Not much time is wasted on incorrect results." << std::endl;
                finish = true;
            }
        }
        else
        {
            finish = true;

        }
    }
    
    FileParser::setKey("MINIMUM_SOLUTION_NETWORK_COUNT", minimumSolutionNetworkCount);
    FileParser::setKey("MINIMUM_TRUST_DISTANCE", distanceTolerance);
    FileParser::setKey("MAX_RECIPROCAL_DISTANCE", maxReciprocalDistance);
    
    logged << std::endl << " ************************ " << std::endl;
    logged << " ****** New values ****** " << std::endl;
    logged << " ************************ " << std::endl;
    logged << std::endl << "MINIMUM_SOLUTION_VECTOR_NETWORK " << minimumSolutionNetworkCount << std::endl;
    logged << "MINIMUM_TRUST_DISTANCE " << distanceTolerance << std::endl;
    logged << "MAX_RECIPROCAL_DISTANCE " << maxReciprocalDistance << std::endl << std::endl;
 
    if (finish)
    {
        logged << "Please copy and paste these parameters into your index.txt file." << std::endl;
        logged << "When you wish to stop teaching, switch LEARNING_TO_INDEX to OFF." << std::endl;
    }
    
    sendLog();
    
    return finish;
}

void IndexManager::index()
{
    bool learningToIndex = FileParser::getKey("LEARNING_TO_INDEX", false);
    
    if (learningToIndex)
    {
        for (int i = 0; i < images.size(); i++)
        {
            images[i]->setLearningToIndex(true);
        }
    }
    
    int maxThreads = FileParser::getMaxThreads();
    IndexingSolution::setupStandardVectors();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    nextImage = -1;
    int maxLearnCycles = 50;
    
    for (int num = 0; num < (learningToIndex ? maxLearnCycles : 1); num++)
    {
        time_t startcputime;
        time(&startcputime);

        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(indexThread, this, &managerSubsets[i], i);
            threads.add_thread(thr);
        }

        threads.join_all();
    
        time_t endcputime;
        time(&endcputime);
        
        clock_t difference = endcputime - startcputime;
        double seconds = difference;
        lastTime = seconds;
        
        if (learningToIndex && num != maxLearnCycles - 1)
        {
            IndexingSolution::reset();
            bool finish = modifyParameters();

            if (finish)
                break;
            
            for (int i = 0; i < images.size(); i++)
            {
                images[i]->reset();
            }
        }

        nextImage = -1;
    }
    
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    mtzs.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzs.insert(mtzs.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
}

void IndexManager::indexFromScratch()
{
    std::vector<double> startingData = FileParser::getKey("INDEX_STARTING_DATA", std::vector<double>());
    
    if (startingData.size() < 6)
    {
        std::cout << "INDEX_STARTING_DATA should be set to dist1, dist2, dist3, angle1, angle2, angle3" << std::endl;
    }
    
    FreeLattice lattice = FreeLattice(startingData[0], startingData[1], startingData[2], startingData[3], startingData[4], startingData[5]);
    lattice.addExpanded();
//    lattice.addExpanded();
    lattice.powderPattern();
    lattice.anglePattern(false);
}

PowderHistogram IndexManager::generatePowderHistogram()
{
    PowderHistogram frequencies;
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    
    double maxCatDistance = 0.05;
    int maxCategory = maxCatDistance / step;
    
    for (int i = 0; i < maxCategory; i++)
    {
        frequencies[i] = std::make_pair(0, 0);
    }
    
    
    for (int i = 0; i < vectorDistances.size(); i++)
    {
        double distance = vectorDistances[i].second;
        int categoryNum = distance / step;
        
        frequencies[categoryNum].second++;
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);
            
            double distance = spotVec->distance();
            int categoryNum = distance / step;
            
            frequencies[categoryNum].first++;
        }
    }
    
    return frequencies;
}


double IndexManager::metrologyTarget(void *object)
{
    static_cast<IndexManager *>(object)->updateAllSpots();
    PowderHistogram histogram = static_cast<IndexManager *>(object)->generatePowderHistogram();
    double highPercentage = 0.0524;
    std::vector<double> frequencies;
    
    for (PowderHistogram::iterator it = histogram.begin(); it != histogram.end(); it++)
    {
        frequencies.push_back(it->second.first);
    }
    
    std::sort(frequencies.begin(), frequencies.end(), std::greater<double>());
    
    int maxCap = highPercentage * frequencies.size();
    int maxSum = 0;
    int minSum = 0;
    int maxCount = 0;
    int minCount = 0;
    
    for (int i = 0; i < frequencies.size(); i++)
    {
        if (i <= maxCap)
        {
            maxSum += frequencies[i];
            maxCount++;
        }
        else
        {
            minSum += frequencies[i];
            minCount++;
        }
    }
    
    double result = -maxSum;
    
    static_cast<IndexManager *>(object)->logged << "Metrology target result: " << result << std::endl;
    static_cast<IndexManager *>(object)->sendLog();
    
    return result;
}

void IndexManager::updateAllSpots()
{
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->updateAllSpots();
    }
}

void IndexManager::refineMetrology()
{
    return;
    
    Panel::setUsePanelInfo(true);
    
    std::map<int, std::pair<double, double> > steps;
    
    for (int i = 0; i < Panel::panelCount(); i++)
    {
        steps[i] = std::make_pair(2, 2);
    }
    
    logged << "Panel count: " << Panel::panelCount() << std::endl;
    sendLog();
    
    bool converged = false;
    double convergeValue = 0.5;
    
    while (!converged)
    {
        converged = true;
        
        for (int i = 0; i < Panel::panelCount(); i++)
        {
            PanelPtr panel = Panel::getPanel(i);
            
            double *xShiftPtr = panel->pointerToBestShiftX();
            double *yShiftPtr = panel->pointerToBestShiftY();
            
            minimizeParameter(steps[i].first, xShiftPtr, metrologyTarget, this);
            minimizeParameter(steps[i].second, yShiftPtr, metrologyTarget, this);
            
            if (steps[i].first > convergeValue || steps[i].second > convergeValue)
            {
                converged = false;
            }
        }
        
    }
    
    logged << Panel::printAll() << std::endl;
    sendLog();
}

void IndexManager::indexingParameterAnalysis()
{
    std::vector<double> goodNetworkCount, badNetworkCount;
    std::vector<double> goodDistanceTrusts, badDistanceTrusts;
    std::vector<double> goodAngleTrusts, badAngleTrusts;
    std::vector<double> goodDistances, badDistances;
    
    logged << "Indexing parameter analysis." << std::endl;
    sendLog();
    
    for (int i = 0; i < images.size(); i++)
    {
        logged << "Processing image " << images[i]->getFilename() << std::endl;
        sendLog();
        
        for (int good = 0; good < 2; good++)
        {
            std::vector<double> *networkCount = good ? &goodNetworkCount : &badNetworkCount;
            std::vector<double> *distanceTrusts = good ? &goodDistanceTrusts : &badDistanceTrusts;
            std::vector<double> *angles = good ? &goodAngleTrusts : &badAngleTrusts;
            std::vector<double> *distances = good ? &goodDistances : &badDistances;
            
            for (int j = 0; j < images[i]->goodOrBadSolutionCount(good); j++)
            {
                IndexingSolutionPtr solution = images[i]->getGoodOrBadSolution(j, good);
                
                networkCount->push_back(solution->spotVectorCount());
                
                std::vector<double> additionalTrusts = solution->totalDistanceTrusts();
                distanceTrusts->reserve(distanceTrusts->size() + additionalTrusts.size());
                distanceTrusts->insert(distanceTrusts->begin(), additionalTrusts.begin(), additionalTrusts.end());
                
                std::vector<double> additionalDistances = solution->totalDistances();
                distances->reserve(distances->size() + additionalDistances.size());
                distances->insert(distances->begin(), additionalDistances.begin(), additionalDistances.end());
                
                std::vector<double> additionalAngles = solution->totalAngles();
                angles->reserve(angles->size() + additionalAngles.size());
                angles->insert(angles->begin(), additionalAngles.begin(), additionalAngles.end());
            }
        }
    }
    /*
    std::vector<double> goodNetworkCount, badNetworkCount;
    std::vector<double> goodDistanceTrusts, badDistanceTrusts;
    std::vector<double> goodAngleTrusts, badAngleTrusts;
    std::vector<double> goodDistances, badDistances;
    */
    
    std::map<double, int> goodNetworkCountHistogram = histogram(goodNetworkCount, 2);
    std::map<double, int> badNetworkCountHistogram = histogram(badNetworkCount, 2);
    histogramCSV("networkCountAnalysis.csv", goodNetworkCountHistogram, badNetworkCountHistogram);

    std::map<double, int> goodDistanceTrustsHistogram = histogram(goodDistanceTrusts, 100);
    std::map<double, int> badDistanceTrustsHistogram = histogram(badDistanceTrusts, 100);
    histogramCSV("distanceTrustAnalysis.csv", goodDistanceTrustsHistogram, badDistanceTrustsHistogram);

    double powderPatternStep = FileParser::getKey("POWDER_PATTERN_STEP", 0.001);
    std::map<double, int> goodDistancesHistogram = histogram(goodDistances, powderPatternStep);
    std::map<double, int> badDistancesHistogram = histogram(badDistances, powderPatternStep);
    histogramCSV("distanceAnalysis.csv", goodDistancesHistogram, badDistancesHistogram);

    std::map<double, int> goodAnglesHistogram = histogram(goodAngleTrusts, 0.05);
    std::map<double, int> badAnglesHistogram = histogram(badAngleTrusts, 0.05);
    histogramCSV("angleTrustAnalysis.csv", goodAnglesHistogram, badAnglesHistogram);
}

std::vector<IOMRefinerPtr> IndexManager::consolidateOrientations(ImagePtr image1, ImagePtr image2, int *oneHand, int *otherHand, int *both)
{
    std::vector<IOMRefinerPtr> refiners;
    
    for (int i = 0; i < image1->IOMRefinerCount(); i++)
    {
        bool found = false;
        IOMRefinerPtr refiner1 = image1->getIOMRefiner(i);
        
        for (int j = 0; j < image2->IOMRefinerCount(); j++)
        {
            IOMRefinerPtr refiner2 = image2->getIOMRefiner(j);
            
            MatrixPtr mat1 = refiner1->getMatrix();
            MatrixPtr mat2 = refiner2->getMatrix();
            
            if (IndexingSolution::matrixSimilarToMatrix(mat1, mat2))
            {
                found = true;
                refiners.push_back(refiner1);
                (*both)++;

                // ignore refiner2 as it is duplicate
            }
        }
        
        if (found == false)
        {
            refiners.push_back(refiner1);
            (*oneHand)++;
        }
    }
    
    for (int i = 0; i < image2->IOMRefinerCount(); i++)
    {
        IOMRefinerPtr refiner1 = image2->getIOMRefiner(i);
        bool found = false;
        
        for (int j = 0; j < image1->IOMRefinerCount(); j++)
        {
            IOMRefinerPtr refiner2 = image1->getIOMRefiner(j);
            
            MatrixPtr mat1 = refiner1->getMatrix();
            MatrixPtr mat2 = refiner2->getMatrix();
            
            if (IndexingSolution::matrixSimilarToMatrix(mat1, mat2))
            {
                found = true;
                // ignore both as it would have been done in the previous loop
                break;
            }
        }
        
        if (!found)
        {
            refiners.push_back(refiner1);
            (*otherHand)++;
        }
    }
    
    return refiners;
}

void IndexManager::combineLists()
{
    IndexingSolution::setupStandardVectors();
    
    std::vector<ImagePtr> mergedImages;
    logged << "Combining lists from " << images.size() << " images and " << mergeImages.size() << " images." << std::endl;
    
    int oneHand = 0;
    int otherHand = 0;
    int both = 0;
    
    // first we iterate through array 1
    for (int i = 0; i < images.size(); i++)
    {
        ImagePtr image1 = images[i];
        bool found = false;
        
        for (int j = 0; j < mergeImages.size(); j++)
        {
            ImagePtr image2 = mergeImages[j];
            
            // Are image1 and image2 the same
            
            if (image1->getFilename() == image2->getFilename())
            {
                found = true;
                // we have a duplicate! check individual orientations
                
                std::vector<IOMRefinerPtr> newRefiners = consolidateOrientations(image1, image2, &oneHand, &otherHand, &both);
                image1->clearIOMRefiners();
                image1->setIOMRefiners(newRefiners);
                mergedImages.push_back(image1);
                
                break;
            }
        }
        
        if (!found)
        {
            mergedImages.push_back(image1);
            oneHand += image1->IOMRefinerCount();
        }
    }
    
    for (int i = 0; i < mergeImages.size(); i++)
    {
        ImagePtr image1 = mergeImages[i];
        bool found = false;
        
        for (int j = 0; j < images.size(); j++)
        {
            ImagePtr image2 = images[j];
            // Are image1 and image2 the same
            
            if (image1->getFilename() == image2->getFilename())
            {
                // we do not consider this image because it was covered in the last loop
                found = true;
                break;
            }
        }
        
        if (!found)
        {
            mergedImages.push_back(image1);
            otherHand += image1->IOMRefinerCount();
        }
    }
    
    logged << "Unique to list 1: " << oneHand << std::endl;
    logged << "Unique to list 2: " << otherHand << std::endl;
    logged << "Present in both lists: " << both << std::endl;

    sendLog();
    
    images = mergedImages;
    mergeImages.clear();
    std::vector<ImagePtr>().swap(mergeImages);
}