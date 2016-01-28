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
    
    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup);
    
    logged << "Calculating distances of unit cell " << spaceGroupNum << std::endl;
    sendLog();
    
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
                sendLog();
                
                
                double distance = length_of_vector(hkl_transformed);
                
                if (distance > maxDistance)
                    maxDistance = distance;
                
                vectorDistances.push_back(std::make_pair(hkl, distance));
            }
        }
    }

    maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.02);
    
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
        logged << vectorDistances[i].first.h << "\t"
         << vectorDistances[i].first.k << "\t"
         << vectorDistances[i].first.l << "\t"
        << vectorDistances[i].second << std::endl;
    }
    
    sendLog();
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
    bool newMethod = FileParser::getKey("NEW_INDEXING_METHOD", false);
    int maxThreads = FileParser::getMaxThreads();
    std::ostringstream logged;
    
    if (newMethod)
    {
        for (int i = offset; i < indexer->images.size(); i += maxThreads)
        {
            ImagePtr image = indexer->images[i];
            logged << "Starting image " << i << std::endl;
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
            
            double angle = spotVec->angleWithVertical();
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
    
    std::ofstream powderLog;
    powderLog.open("powder.csv");
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);

    for (PowderHistogram::iterator it = frequencies.begin(); it != frequencies.end(); it++)
    {
        double distance = it->first * step;
        double freq = it->second.first;
        double perfect = it->second.second;
        
        powderLog << distance << "," << freq << "," << perfect << std::endl;
    }
    
    powderLog.close();
    
    /// angles
    
    std::vector<double> probeDistances = FileParser::getKey("PROBE_DISTANCES", std::vector<double>());
    
    if (probeDistances.size() < 2)
        return;
    
    double distance1 = probeDistances[0]; double distance2 = probeDistances[1];
    
    double angleStep = FileParser::getKey("POWDER_PATTERN_STEP_ANGLE", 2.0) * M_PI / 180;
    double distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 500.);
    std::map<int, int> angleHistogram;
    
    for (double i = 0; i < M_PI; i += angleStep)
    {
        int angleCategory = i / angleStep;
        angleHistogram[angleCategory] = 0;
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
    
    std::ofstream angleLog;
    angleLog.open("angle.csv");
    
    for (std::map<int, int>::iterator it = angleHistogram.begin(); it != angleHistogram.end(); it++)
    {
        double angle = it->first * angleStep * 180 / M_PI;
        double freq = it->second;
        
        angleLog << angle << "," << freq << "," << std::endl;
    }
    
    angleLog.close();
}

void IndexManager::index()
{
    int maxThreads = FileParser::getMaxThreads();
    IndexingSolution::setupStandardVectors();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(indexThread, this, &managerSubsets[i], i);
        threads.add_thread(thr);
    }
    
    
    threads.join_all();
    
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
    
    double minAve = minSum / minCount;
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
