//
//  DIALSImage.cpp
//  cppxfel
//
//  Created by Helen Ginn on 11/03/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "DIALSImage.h"
#include <scitbx/vec2.h>
#include "Matrix.h"
#include "IndexingSolution.h"
#include <dxtbx/model/panel.h>
#include "FileParser.h"
#include "InputFileParser.h"

using coord_type = dxtbx::model::Detector::coord_type;
using scitbx::vec2;
using iterator =  dxtbx::model::Detector::iterator;

using namespace cppxfel;

InputFileParser *DIALSImage::personalParser = new InputFileParser("");

void DIALSImage::init()
{
    Logger::mainLogger = LoggerPtr(new Logger());
    boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);
}

void DIALSImage::setMaxReciprocalDistance(double newMaxReciprocalDistance)
{
    FileParser::setKey("MAX_RECIPROCAL_DISTANCE", newMaxReciprocalDistance);
}

void DIALSImage::setInverseDomainSize(double inverseDomainSize)
{
    FileParser::setKey("INITIAL_RLP_SIZE", inverseDomainSize);
}

void DIALSImage::setCrystalParameters(int spaceGroupNum, double params[6])
{
    FileParser::setKey("SPACE_GROUP", spaceGroupNum);

    std::vector<double> unitCell;

    for (int i = 0; i < 6; i++)
    {
        unitCell.push_back(params[i]);
    }

    FileParser::setKey("UNIT_CELL", unitCell);
}

// convert from one format into another and also change origin
vec3<double> DIALSImage::cppxfelVecToScitbxVec3(vec hkl)
{
    double radius = 1 / getWavelength();

    hkl.l += radius;

    return vec3<double>(hkl.h, hkl.k, hkl.l);
}

vec DIALSImage::scitbxVec3ToCppxfelVec(vec3<double> hkl3)
{
    double radius = 1 / getWavelength();

    hkl3[2] -= radius;

    return new_vector(hkl3[0], hkl3[1], hkl3[2]);
}

std::pair<double, double> DIALSImage::reciprocalCoordinatesToPixels(vec hkl)
{
    vec3<double> scitbxVec3 = cppxfelVecToScitbxVec3(hkl);

    coord_type ray = dxtbxDetector->get_ray_intersection(scitbxVec3);

    // first of ray pair is panel number
    // second of ray pair is currently mm

    int panel_num = ray.first;

    // find true pixel intersection off panel num

    vec2<double> pix = (*dxtbxDetector)[panel_num].get_ray_intersection_px(scitbxVec3);

    return std::make_pair(pix[0], pix[1]);
}

void DIALSImage::addSpots(std::vector<vec3<double> > rays)
{
    double radius = 1 / getWavelength();
    rays.push_back(vec3<double>(0, 0, radius));

    for (int i = 0; i < rays.size(); i++)
    {
        vec3<double> ray = rays[i];
        vec hkl = scitbxVec3ToCppxfelVec(ray);

        SpotPtr spot = SpotPtr(new Spot(shared_from_this()));
        spot->setXYFromEstimatedVector(hkl);

        spots.push_back(spot);
    }
}

IndexingSolutionStatus DIALSImage::tryIndexingSolution(IndexingSolutionPtr solutionPtr)
{
    if (solutionPtr->spotVectorCount() < minimumSolutionNetworkCount)
        return IndexingSolutionTrialFailure;

    MatrixPtr solutionMatrix = solutionPtr->createSolution();
    bool similar = checkIndexingSolutionDuplicates(solutionMatrix);

    if (!similar)
    {
        solutions.push_back(solutionMatrix);
        return IndexingSolutionTrialSuccess;
    }
    else
    {
        return IndexingSolutionTrialDuplicate;
    }
}

bool DIALSImage::checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast)
{
    for (int i = 0; i < solutions.size() - excludeLast; i++)
    {
        bool similar = IndexingSolution::matrixSimilarToMatrix(newSolution, solutions[i], true);

        if (similar)
            return true;
    }

    return false;
}

void DIALSImage::getSolution(int i, double *matrixParams[9])
{
    MatrixPtr solution = solutions[i];

    solution->sensibleComponents(matrixParams);
}
