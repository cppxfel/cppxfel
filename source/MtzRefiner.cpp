/*
 * MtzRefiner.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "MtzRefiner.h"
#include <vector>
#include "misc.h"
#include "Vector.h"
#include <boost/thread/thread.hpp>
#include "GraphDrawer.h"
#include "Image.h"
#include "Miller.h"
#include <fstream>
#include "AmbiguityBreaker.h"
#include "IndexManager.h"
#include "FreeMillerLibrary.h"

#include "FileParser.h"
#include "parameters.h"
#include "FileReader.h"
#include "PanelParser.h"
#include "Panel.h"
#include "Logger.h"
#include "XManager.h"
#include "MtzGrouper.h"

bool MtzRefiner::hasPanelParser;
int MtzRefiner::imageLimit;
int MtzRefiner::cycleNum;



MtzRefiner::MtzRefiner()
{
    // TODO Auto-generated constructor stub
    reference = NULL;
    panelParser = NULL;

    hasPanelParser = false;
    int logInt = FileParser::getKey("VERBOSITY_LEVEL", 0);
    Logger::mainLogger->changePriorityLevel((LogLevel)logInt);

    imageLimit = FileParser::getKey("IMAGE_LIMIT", 0);

    hasRefined = false;
    isPython = false;

    indexManager = NULL;
}

void print_parameters()
{
    std::cout << "N: ===== Default parameters =====" << std::endl;

    std::cout << "N: Initial parameters for MTZ grid search:" << std::endl;
    std::cout << "N: Bandwidth: " << INITIAL_BANDWIDTH << std::endl;
    std::cout << "N: Mosaicity: " << INITIAL_MOSAICITY << std::endl;
    std::cout << "N: Spot size: " << INITIAL_SPOT_SIZE << std::endl;
    std::cout << "N: Exponent: " << INITIAL_EXPONENT << std::endl;
    std::cout << std::endl;

    std::cout << "N: Parameters refined during grid search: ";
    std::cout << (OPTIMISED_ROT ? "rotations; " : "");
    std::cout << (OPTIMISED_BANDWIDTH ? "bandwidth; " : "");
    std::cout << (OPTIMISED_WAVELENGTH ? "wavelength; " : "");
    std::cout << (OPTIMISED_EXPONENT ? "exponent; " : "");
    std::cout << (OPTIMISED_MOSAICITY ? "mosaicity; " : "");
    std::cout << (OPTIMISED_SPOT_SIZE ? "spot size; " : "");
    std::cout << std::endl << std::endl;

    std::cout << "N: Grid search steps: ";
    std::cout << "N: Bandwidth: " << BANDWIDTH_STEP << std::endl;
    std::cout << "N: Rotation: " << ROT_STEP << std::endl;
    std::cout << "N: Spot size: " << SPOT_STEP << std::endl;
    std::cout << "N: Exponent: " << EXPONENT_STEP << std::endl;
    std::cout << "N: Wavelength: " << MEAN_STEP << std::endl;
    std::cout << std::endl;

    std::cout << "N: Grid search tolerances: ";
    std::cout << "N: Bandwidth: " << BANDWIDTH_TOLERANCE << std::endl;
    std::cout << "N: Rotation: " << ROT_TOLERANCE << std::endl;
    std::cout << "N: Spot size: " << SPOT_STEP << std::endl;
    std::cout << "N: Exponent: " << SPOT_SIZE_TOLERANCE << std::endl;
    std::cout << "N: Wavelength: " << MEAN_TOLERANCE << std::endl;
    std::cout << std::endl;

    std::cout << "N: Maximum resolution used for most grid searches: "
    << MAX_OPTIMISATION_RESOLUTION << std::endl;
    std::cout << "N: Maximum resolution used for spot size: "
    << MAX_SPOT_SIZE_OPT_RESOLUTION << std::endl;
    std::cout << "N: Further optimisation is "
    << (FURTHER_OPTIMISATION ? "" : "not ") << "used." << std::endl;

    std::cout << "N: Default detector distance: " << DEFAULT_DETECTOR_DISTANCE
    << std::endl;
    std::cout << "N: Default wavelength: " << DEFAULT_WAVELENGTH << std::endl;
    std::cout << "N: Intensity threshold: " << INTENSITY_THRESHOLD << std::endl;

    std::cout << std::endl << "N: Polarisation correction: "
    << (POLARISATION_CORRECTION ? "on" : "off") << std::endl;
    std::cout << "N: Polarisation factor (horizontal): "
    << HORIZONTAL_POLARISATION_FACTOR << std::endl;
    std::cout << std::endl << "N: Minimum miller count for rejection: "
    << MIN_MILLER_COUNT << std::endl;
    std::cout << "N: Rejecting miller indices: "
    << (REJECTING_MILLERS ? "on" : "off") << std::endl;
}

// MARK: Refinement

void MtzRefiner::cycleThreadWrapper(MtzRefiner *object, int offset)
{
    object->cycleThread(offset);
}

void MtzRefiner::cycleThread(int offset)
{
    int img_num = (int)mtzManagers.size();
    int j = 0;

    bool partialitySpectrumRefinement = FileParser::getKey("REFINE_ENERGY_SPECTRUM", false);

    std::vector<int> targets = FileParser::getKey("TARGET_FUNCTIONS", std::vector<int>());

    int maxThreads = FileParser::getMaxThreads();

    for (int i = offset; i < img_num; i += maxThreads)
    {
        j++;

        std::ostringstream logged;

        MtzPtr image = mtzManagers[i];

        if (!image->isRejected())
        {
            logged << "Refining image " << i << " " << image->getFilename() << std::endl;
            Logger::mainLogger->addStream(&logged);

            bool silent = (targets.size() > 0);

            if (partialitySpectrumRefinement)
            {
                image->refinePartialities();
                image->replaceBeamWithSpectrum();
            }
            else
            {
                image->gridSearch(silent);
            }

            if (targets.size() > 0)
            {
                ScoreType firstScore = image->getScoreType();
                for (int i = 0; i < targets.size(); i++)
                {
                    silent = (i < targets.size() - 1);
                    image->setDefaultScoreType((ScoreType)targets[i]);

                    if (partialitySpectrumRefinement)
                    {
                        image->refinePartialities();
                    }
                    else
                    {
                        image->gridSearch(silent);
                    }
                }
                image->setDefaultScoreType(firstScore);
            }

        //    image->writeToDat();
        }
    }
}

void MtzRefiner::cycle()
{
    MtzManager::setReference(reference);

    time_t startcputime;
    time(&startcputime);

    boost::thread_group threads;

    int maxThreads = FileParser::getMaxThreads();

    std::ostringstream logged;
    logged << "Filename\tScore type\t\tCorrel\tRfactor\tPart correl\tRmerge\tB factor\tHits" << std::endl;
    Logger::mainLogger->addStream(&logged);

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(cycleThreadWrapper, this, i);
        threads.add_thread(thr);
    }

    threads.join_all();

    time_t endcputime;
    time(&endcputime);

    clock_t difference = endcputime - startcputime;
    double seconds = difference;

    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;

    std::cout << "N: Refinement cycle: " << minutes << " minutes, "
    << finalSeconds << " seconds." << std::endl;
}

void MtzRefiner::initialMerge()
{
    MtzManager *originalMerge = NULL;
    /*
     Lbfgs_Cluster *lbfgs = new Lbfgs_Cluster();
     lbfgs->initialise_cluster_lbfgs(mtzManagers, &originalMerge);
     reference = originalMerge;
     delete lbfgs;

     reference->writeToFile("initialMerge.mtz");
     */

    AmbiguityBreaker breaker = AmbiguityBreaker(mtzManagers);
    breaker.run();
    originalMerge = breaker.getMergedMtz();
    reference = originalMerge;

    reference->writeToFile("initialMerge.mtz");
}

void MtzRefiner::setupFreeMillers()
{
    std::string freeMiller = FileParser::getKey("FREE_MILLER_LIST", std::string(""));
    double freeMillerProportion = FileParser::getKey("FREE_MILLER_PROPORTION", 0.0);

    if (freeMiller.length() && freeMillerProportion > 0)
    {
        FreeMillerLibrary::setup();
    }
}

void MtzRefiner::refine()
{
    setupFreeMillers();

    MtzManager *originalMerge = NULL;

    bool initialExists = loadInitialMtz();
    readMatricesAndMtzs();

    if (!initialExists)
    {
        initialMerge();
    }

    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", false);

    if (fixUnitCell)
    {
        for (int i = 0; i < mtzManagers.size(); i++)
        {
            vector<double> givenUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());

            mtzManagers[i]->getMatrix()->changeOrientationMatrixDimensions(givenUnitCell[0], givenUnitCell[1], givenUnitCell[2], givenUnitCell[3], givenUnitCell[4], givenUnitCell[5]);
            mtzManagers[i]->setUnitCell(givenUnitCell);
        }
    }

    originalMerge = reference;

    std::cout << "N: Total images loaded: " << mtzManagers.size() << std::endl;

    originalMerge->writeToFile("originalMerge.mtz");

    MtzManager::currentManager = originalMerge;
    MtzManager::setReference(reference);
    double correl = originalMerge->correlation(true);
    std::cout << "Merged correlation = " << correl << std::endl;

    MtzManager::setLowRes(0);
    MtzManager::setHighRes(0);

    refineCycle();

    hasRefined = true;
}

void MtzRefiner::refineCycle(bool once)
{
    int i = 0;
    bool finished = false;

    int minimumCycles = FileParser::getKey("MINIMUM_CYCLES", 6);
    int maximumCycles = FileParser::getKey("MAXIMUM_CYCLES", 50);

    bool stop = FileParser::getKey("STOP_REFINEMENT", true);
    double correlationThreshold = FileParser::getKey("CORRELATION_THRESHOLD",
                                                     CORRELATION_THRESHOLD);
    double initialCorrelationThreshold = FileParser::getKey(
                                                            "INITIAL_CORRELATION_THRESHOLD", INITIAL_CORRELATION_THRESHOLD);
    int thresholdSwap = FileParser::getKey("THRESHOLD_SWAP",
                                           THRESHOLD_SWAP);
    bool exclusion = FileParser::getKey("EXCLUDE_OWN_REFLECTIONS", false);
    bool replaceReference = FileParser::getKey("REPLACE_REFERENCE", true);

    double resolution = FileParser::getKey("MAX_RESOLUTION_ALL",
                                           MAX_OPTIMISATION_RESOLUTION);
    int scalingInt = FileParser::getKey("SCALING_STRATEGY",
                                        (int) SCALING_STRATEGY);
    ScalingType scaling = (ScalingType) scalingInt;

    while (!finished)
    {
        cycleNum = i;
        cycle();

        std::cout << "Grouping final MTZs" << std::endl;
        MtzGrouper *grouper = new MtzGrouper();
        MtzManager::setReference(reference);
        if (i >= 0)
            grouper->setScalingType(scaling);
        if (once)
            grouper->setScalingType(ScalingTypeAverage);
        grouper->setWeighting(WeightTypePartialitySigma);
        grouper->setExpectedResolution(resolution);
        grouper->setMtzManagers(mtzManagers);
        grouper->setExcludeWorst(true);
        if (i < thresholdSwap)
            grouper->setCorrelationThreshold(initialCorrelationThreshold);
        else
            grouper->setCorrelationThreshold(correlationThreshold);

        MtzManager *mergedMtz = NULL;
        MtzManager *unmergedMtz = NULL;
        grouper->merge(&mergedMtz, &unmergedMtz, i);

        MtzManager::currentManager = mergedMtz;
        MtzManager::setReference(reference);

        std::cout << "Reflections: " << mergedMtz->reflectionCount() << std::endl;
        if (!once)
        {
            double correl = mergedMtz->correlation(true);
            std::cout << "N: Merged correlation = " << correl << std::endl;
        }
        mergedMtz->description();

        double scale = 1000 / mergedMtz->averageIntensity();
        mergedMtz->applyScaleFactor(scale);

        std::cout << "Here" << std::endl;

        std::string filename = "allMerge" + i_to_str(i) + ".mtz";
        mergedMtz->writeToFile(filename.c_str(), true);

        MtzManager *reloadMerged = new MtzManager();

        reloadMerged->setFilename(filename.c_str());
        reloadMerged->loadReflections(1);
        reloadMerged->description();

        MtzManager::currentManager = reloadMerged;
        double reloaded = reloadMerged->correlation(true);
        std::cout << "Reloaded correlation = " << reloaded << std::endl;

        if (reloaded > 0.999 && i >= minimumCycles && stop)
            finished = true;

        if (once)
            finished = true;

        if (i == maximumCycles - 1 && stop)
            finished = true;

        delete grouper;
        if (replaceReference)
        {
            if (!once)
                delete reference;
            reference = exclusion ? unmergedMtz : mergedMtz;
            MtzManager::setReference(reference);
        }
        i++;
    }
}


// MARK: Symmetry-related reflection refinement

void MtzRefiner::refineSymmetry()
{
    readMatricesAndMtzs();

    std::cout << "N: Total images loaded: " << mtzManagers.size() << std::endl;

    std::cout << "Refining images using symmetry: " << std::endl;

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        int symmCount = mtzManagers[i]->symmetryRelatedReflectionCount();
        int paramCount = mtzManagers[i]->refinedParameterCount();

        if (symmCount < paramCount * 8 && paramCount > 30)
        {
            mtzManagers[i]->setRejected(true);
            continue;
        }

        std::cout << mtzManagers[i]->getFilename() << "\t" << symmCount << "\t" << paramCount << std::endl;

        mtzManagers[i]->setDefaultScoreType(ScoreTypeSymmetry);
    }

    refineCycle(true);

    MtzManager::setReference(reference);

    int defaultScoreInt = FileParser::getKey("DEFAULT_TARGET_FUNCTION",
                                             (int) DEFAULT_SCORE_TYPE);
    ScoreType desiredScoreType = (ScoreType) defaultScoreInt;

    double spotSum = 0;
    double mosaicitySum = 0;
    int count = 0;

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        double spotSize = fabs(mtzManagers[i]->getSpotSize());
        double mosaicity = fabs(mtzManagers[i]->getMosaicity());

        if (mtzManagers[i]->isRejected() == false)
        {
            spotSum += spotSize;
            mosaicitySum += mosaicity;
            count++;
        }
    }

    spotSum /= count;
    mosaicitySum /= count;

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if (mtzManagers[i]->isRejected())
        {
            mtzManagers[i]->setSpotSize(spotSum);
            mtzManagers[i]->setMosaicity(mosaicitySum);
            mtzManagers[i]->setRejected(false);
        }

        mtzManagers[i]->setOptimisingOrientation(true);
        mtzManagers[i]->setOptimisingWavelength(true);

        mtzManagers[i]->setDefaultScoreType(desiredScoreType);
    }
}

// MARK: Loading data

bool MtzRefiner::loadInitialMtz(bool force)
{
    bool hasInitialMtz = FileParser::hasKey("INITIAL_MTZ");

    if (reference != NULL && !force)
        return true;

    std::ostringstream logged;

    logged << "Initial MTZ has "
    << (hasInitialMtz ? "" : "not ") << "been provided." << std::endl;

    Logger::mainLogger->addStream(&logged);

    if (hasInitialMtz)
    {

        std::string referenceFile = FileParser::getKey("INITIAL_MTZ",
                                                       std::string(""));

        reference = new MtzManager();

        reference->setFilename(referenceFile.c_str());
        reference->loadReflections(1);
        reference->setSigmaToUnity();
    }

    MtzManager::setReference(reference);

    return hasInitialMtz;
}


int MtzRefiner::imageMax(size_t lineCount)
{
    int skip = imageSkip(lineCount);

    int end = (int)lineCount;

    if (imageLimit != 0)
        end = imageLimit < lineCount ? imageLimit : (int)lineCount;

    end += skip;

    if (end > lineCount)
        end = (int)lineCount;

    return end;
}

void MtzRefiner::applyParametersToImages()
{
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());

    if (unitCell.size() == 0)
    {
        std::cout
        << "Please provide unit cell dimensions under keyword UNIT_CELL"
        << std::endl;
        exit(1);
    }

    int spg_num = FileParser::getKey("SPACE_GROUP", -1);

    double distance = 0;
    double wavelength = 0;

    if (spg_num == -1)
    {
        std::cout << "Please set space group number under keyword SPACE_GROUP"
        << std::endl;
    }

    if (isFromPython())
    {
        distance = FileParser::getKey("DETECTOR_DISTANCE", 0.);
        wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.);

        if (distance != 0 || wavelength != 0)
        {
            std::cout << "Overriding image parameters with distance and/or wavelength specified in input file." << std::endl;
        }
    }

    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);

    for (int i = 0; i < images.size(); i++)
    {
        ImagePtr newImage = images[i];

        newImage->setPinPoint(true);

        CCP4SPG *spg = ccp4spg_load_by_standard_num(spg_num);

        newImage->setSpaceGroup(spg);

        if (distance != 0)
        {
            newImage->setDetectorDistance(distance);
        }

        if (wavelength != 0)
        {
            newImage->setWavelength(wavelength);
        }

        if (fixUnitCell)
            newImage->setUnitCell(unitCell);

        if (fixUnitCell)
            for (int j = 0; j < newImage->IOMRefinerCount(); j++)
            {
                newImage->getIOMRefiner(j)->getMatrix()->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
            }
    }
}

void MtzRefiner::readSingleImageV2(std::string *filename, vector<ImagePtr> *newImages, vector<MtzPtr> *newMtzs, int offset)
{
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    double tolerance = FileParser::getKey("ACCEPTABLE_UNIT_CELL_TOLERANCE", 0.0);
    vector<double> givenUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());

    bool checkingUnitCell = false;

    if (givenUnitCell.size() > 0 && tolerance > 0.0)
        checkingUnitCell = true;

    if (wavelength == 0 && newImages != NULL)
    {
        std::cout << "No integration wavelength provided. If this is not deliberate, please provide initial wavelength for integration under keyword INTEGRATION_WAVELENGTH" << std::endl;
    }

    if (detectorDistance == 0 && newImages != NULL)
    {
        std::cout << "No detector distance provided. If this is not deliberate, please provide detector distance for integration under keyword DETECTOR_DISTANCE" << std::endl;

    }

    bool hasBeamCentre = FileParser::hasKey("BEAM_CENTRE");

    bool setSigmaToUnity = FileParser::getKey("SET_SIGMA_TO_UNITY", true);

    bool ignoreMissing = FileParser::getKey("IGNORE_MISSING_IMAGES", false);
    const std::string contents = FileReader::get_file_contents(
                                                               filename->c_str());

    vector<std::string> imageList = FileReader::split(contents, "\nimage ");

    int maxThreads = FileParser::getMaxThreads();

    int skip = imageSkip(imageList.size());
    int end = imageMax(imageList.size());

    if (skip > 0)
    {
        std::ostringstream logged;
        logged << "Skipping " << skip << " lines" << std::endl;
        Logger::mainLogger->addStream(&logged);
    }

    for (int i = offset + skip; i < end; i += maxThreads)
    {
        if (imageList[i].length() == 0)
            continue;

        vector<std::string> lines = FileReader::split(imageList[i], '\n');

        vector<std::string> components = FileReader::split(lines[0], ' ');

        if (components.size() <= 1)
            continue;

        std::string imgName = components[1];

        imgName.erase(std::remove(imgName.begin(), imgName.end(), '\r'), imgName.end());
        imgName.erase(std::remove(imgName.begin(), imgName.end(), '\n'), imgName.end());

        if (newImages)
            imgName += ".img";
        else if (newMtzs)
            imgName += ".mtz";

        std::ostringstream logged;

        if (!FileReader::exists(imgName) && !ignoreMissing)
        {
            logged << "Skipping image " << imgName << std::endl;
            Logger::mainLogger->addStream(&logged);
            continue;
        }

        logged << "Loading image " << i << " (" << imgName << ")"
        << std::endl;
        Logger::mainLogger->addStream(&logged);

        double usedWavelength = wavelength;
        double usedDistance = detectorDistance;

        bool fromDials = FileParser::getKey("FROM_DIALS", true);

        if (components.size() >= 3)
        {
            usedWavelength = atof(components[1].c_str());
            usedDistance = atof(components[2].c_str());
        }

        MatrixPtr unitCell;
        MatrixPtr newMatrix;
        std::string paramsLine = "";

        ImagePtr newImage = ImagePtr(new Image(imgName, wavelength,
                                    detectorDistance));
        bool hasSpots = false;

        for (int i = 1; i < lines.size(); i++)
        {
            vector<std::string> components = FileReader::split(lines[i], ' ');

            if (components.size() == 0)
                continue;

            if (components[0] == "spots")
            {
                std::string spotsFile = components[1];
                newImage->setSpotsFile(spotsFile);
                hasSpots = true;
            }

            if (components[0] == "matrix")
            {
                // individual matrices
                double matrix[9];
                readMatrix(matrix, lines[i]);

                newMatrix = MatrixPtr(new Matrix(matrix));

                if (fromDials)
                {
                    newMatrix->rotate(0, 0, M_PI / 2);
                }

                if (newImages)
                {
                    newImage->setUpIOMRefiner(newMatrix);
                }
            }

            if (components[0] == "wavelength" && wavelength == 0)
            {
                double newWavelength = wavelength;

                if (components.size() >= 2)
                    newWavelength = atof(components[1].c_str());

                newImage->setWavelength(newWavelength);
            }

            if (components[0] == "distance" && detectorDistance == 0)
            {
                double newDistance = detectorDistance;

                if (components.size() >= 2)
                    newDistance = atof(components[1].c_str());

                newImage->setDetectorDistance(newDistance);
            }


            if (components[0] == "centre" && !hasBeamCentre)
            {
                if (components.size() < 3)
                    continue;

                double centreX = atof(components[1].c_str());
                double centreY = atof(components[2].c_str());

                newImage->setBeamX(centreX);
                newImage->setBeamY(centreY);
            }

            if (components[0] == "unitcell")
            {
                // individual matrices
                double matrix[9];
                readMatrix(matrix, lines[i]);

                unitCell = MatrixPtr(new Matrix(matrix));
            }

            if (components[0] == "rotation")
            {
                // individual matrices
                double matrix[9];
                readMatrix(matrix, lines[i]);

                MatrixPtr rotation = MatrixPtr(new Matrix(matrix));

                if (unitCell)
                {
                    if (newImages)
                    {
                        vector<double> correction = FileParser::getKey("ORIENTATION_CORRECTION", vector<double>());

                        if (fromDials)
                        {
                            rotation->rotate(0, 0, M_PI / 2);
                        }

                        if (correction.size() >= 2)
                        {
                            double rightRot = correction[0] * M_PI / 180;
                            double upRot = correction[1] * M_PI / 180;
                            double swivelRot = 0;

                            if (correction.size() > 2)
                                swivelRot = correction[2] * M_PI / 180;

                            rotation->rotate(rightRot, upRot, swivelRot);
                        }
                    }

                    newMatrix = MatrixPtr(new Matrix);
                    newMatrix->setComplexMatrix(unitCell, rotation);

                    if (newImages)
                    {
                        newImage->setUpIOMRefiner(newMatrix);
                    }

                    unitCell = MatrixPtr();
                }
                else
                {
                    std::cout << "Warning, unmatched unitcell / rotation pair?" << std::endl;
                }
            }

            if (components[0] == "params" && newMtzs)
            {
                paramsLine = lines[i];
            }
        }
        if (newImages)
        {
            if (hasSpots)
            {
                newImage->processSpotList();
            }

            if (checkingUnitCell && newImage->checkUnitCell(givenUnitCell[0], givenUnitCell[1], givenUnitCell[2], tolerance))
            {
                newImages->push_back(newImage);

            }
            else if (!checkingUnitCell)
                newImages->push_back(newImage);
        }

        if (newMtzs)
        {
            MtzPtr newManager = MtzPtr(new MtzManager());

            newManager->setFilename(imgName.c_str());
            newManager->setMatrix(newMatrix);
            newManager->loadReflections(PartialityModelScaled, true);

            if (setSigmaToUnity)
                newManager->setSigmaToUnity();
            newManager->loadParametersMap();
            newManager->setParamLine(paramsLine);

            if (newManager->reflectionCount() > 0)
            {
                if (checkingUnitCell && newManager->checkUnitCell(givenUnitCell[0], givenUnitCell[1], givenUnitCell[2], tolerance))
                {
                    newMtzs->push_back(newManager);
                }
                else if (!checkingUnitCell)
                {
                    newMtzs->push_back(newManager);
                }
                else
                {
                    Logger::mainLogger->addString("Skipping file " + *filename + " due to poor unit cell");
                }
            }
        }
    }

}

void MtzRefiner::readMatricesAndImages(std::string *filename, bool areImages, std::vector<ImagePtr> *targetImages)
{
    if (targetImages == NULL && images.size() > 0)
        return;

    std::string aFilename = "";

    if (filename == NULL)
    {
        aFilename = FileParser::getKey("ORIENTATION_MATRIX_LIST", std::string(""));
        filename = &aFilename;

        if (filename->length() == 0)
        {
            std::cout << "No orientation matrix list provided. Exiting now." << std::endl;
            exit(1);
        }
    }

    double version = FileParser::getKey("MATRIX_LIST_VERSION", 2.0);

    // thought: turn the vector concatenation into a templated function

    boost::thread_group threads;

    int maxThreads = FileParser::getMaxThreads();

    vector<vector<ImagePtr> > imageSubsets;
    vector<vector<MtzPtr> > mtzSubsets;
    imageSubsets.resize(maxThreads);
    mtzSubsets.resize(maxThreads);

    for (int i = 0; i < maxThreads; i++)
    {
        if (version == 1.0)
        {
            boost::thread *thr = new boost::thread(singleLoadImages, filename,
                                                   &imageSubsets[i], i);
            threads.add_thread(thr);
        }
        else if (version == 2.0)
        {
            vector<MtzPtr> *chosenMtzs = areImages ? NULL : &mtzSubsets[i];
            vector<ImagePtr> *chosenImages = areImages ? &imageSubsets[i] : NULL;
            boost::thread *thr = new boost::thread(readSingleImageV2, filename,
                                                   chosenImages, chosenMtzs, i);
            threads.add_thread(thr);
        }
    }

    threads.join_all();


    int total = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        if (areImages)
            total += imageSubsets[i].size();

        if (!areImages)
            total += mtzSubsets[i].size();
    }

    if (targetImages == NULL)
    {
        targetImages = &images;
    }

    mtzManagers.reserve(total);
    targetImages->reserve(total);
    int lastPos = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        if (areImages)
        {
            targetImages->insert(targetImages->begin() + lastPos,
                          imageSubsets[i].begin(), imageSubsets[i].end());
            lastPos += imageSubsets[i].size();
        }
        else
        {
            mtzManagers.insert(mtzManagers.begin() + lastPos,
                               mtzSubsets[i].begin(), mtzSubsets[i].end());
            lastPos += mtzSubsets[i].size();
        }
    }

 //   if (version == 2.0 && areImages)
 //       applyParametersToImages();
}

void MtzRefiner::singleLoadImages(std::string *filename, vector<ImagePtr> *newImages, int offset)
{
    const std::string contents = FileReader::get_file_contents(
                                                               filename->c_str());

    vector<std::string> lines = FileReader::split(contents, '\n');

    int spg_num = FileParser::getKey("SPACE_GROUP", -1);

    if (spg_num == -1)
    {
        std::cout << "Please set space group number under keyword SPACE_GROUP"
        << std::endl;
    }

    double overPredSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
                                                 OVER_PRED_SPOT_SIZE);
    double overPredBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
                                                  OVER_PRED_BANDWIDTH);
    overPredBandwidth /= 2;

    double orientationStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);

    // @TODO intensityThreshold should be flexible
    double intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD",
                                                   INTENSITY_THRESHOLD);

    double metrologySearchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                                    METROLOGY_SEARCH_SIZE);

    double maxIntegratedResolution = FileParser::getKey(
                                                        "MAX_INTEGRATED_RESOLUTION", MAX_INTEGRATED_RESOLUTION);

    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);

    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);

    if (wavelength == 0)
    {
        std::cout
        << "Please provide initial wavelength for integration under keyword INTEGRATION_WAVELENGTH"
        << std::endl;
        exit(1);
    }

    if (detectorDistance == 0)
    {
        std::cout
        << "Please provide detector distance for integration under keyword DETECTOR_DISTANCE"
        << std::endl;
        exit(1);
    }

    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());

    if (unitCell.size() == 0)
    {
        std::cout
        << "Please provide unit cell dimensions under keyword UNIT_CELL"
        << std::endl;
        exit(1);
    }

    int maxThreads = FileParser::getMaxThreads();

    int end = imageMax(lines.size());

    int skip = imageSkip(lines.size());

    if (skip > 0)
    {
        std::cout << "Skipping " << skip << " lines" << std::endl;
    }

    for (int i = offset + skip; i < end + skip; i += maxThreads)
    {
        vector<std::string> components = FileReader::split(lines[i], ' ');

        if (components.size() == 0)
            continue;

        std::string imgName = components[0] + ".img";
        if (!FileReader::exists(imgName))
        {
            continue;
        }

        std::cout << "Loading image " << i << " (" << imgName << ")"
        << std::endl;

        if (components.size() > 11)
        {
            wavelength = atof(components[10].c_str());
            detectorDistance = atof(components[11].c_str());
        }


        ImagePtr newImage = ImagePtr(new Image(imgName, wavelength,
                                    detectorDistance));

        MatrixPtr newMat = MatrixPtr();

        if (components.size() >= 10)
        {
            double matrix[9];
            readMatrix(matrix, lines[i]);
            newMat = MatrixPtr(new Matrix(matrix));
            if (fixUnitCell)
                newMat->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
        }
        else
        {
            newMat = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
        }

        newImage->setUpIOMRefiner(newMat);
        newImage->setPinPoint(true);

        CCP4SPG *spg = ccp4spg_load_by_standard_num(spg_num);

        newImage->setSpaceGroup(spg);
        newImage->setMaxResolution(maxIntegratedResolution);
        newImage->setSearchSize(metrologySearchSize);
        newImage->setIntensityThreshold(intensityThreshold);
        newImage->setUnitCell(unitCell);
        newImage->setInitialStep(orientationStep);

        // @TODO: masks should be addable through input file
        /*
         newImage->addMask(0, 840, 1765, 920);
         newImage->addMask(446, 1046, 1324, 1765);
         newImage->addMask(820, 450, 839, 473);
         */

        newImage->setTestSpotSize(overPredSpotSize);
        newImage->setTestBandwidth(overPredBandwidth);

        newImages->push_back(newImage);
    }
}

void MtzRefiner::readMatricesAndMtzs()
{
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              std::string(""));

    double version = FileParser::getKey("MATRIX_LIST_VERSION", 2.0);

    if (version == 2.0)
    {
        readMatricesAndImages(NULL, false);

        return;
    }


    if (filename == "")
    {
        std::cout
        << "Orientation matrix list has not been provided for refinement. Please provide under keyword ORIENTATION_MATRIX_LIST."
        << std::endl;
        exit(1);
    }

    std::ostringstream logged;

    if (mtzManagers.size() > 0)
    {
        logged << "Mtzs already present; not reloading from list" << std::endl;
        Logger::mainLogger->addStream(&logged);
        return;
    }

    logged << "Loading MTZs from list" << std::endl;

    const std::string contents = FileReader::get_file_contents(
                                                               filename.c_str());

    vector<std::string> lines = FileReader::split(contents, '\n');

    int maxThreads = FileParser::getMaxThreads();

    boost::thread_group threads;

    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(singleThreadRead, lines,
                                               &managerSubsets[i], i);
        threads.add_thread(thr);
    }

    threads.join_all();

    int total = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }

    mtzManagers.reserve(total);
    int lastPos = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        mtzManagers.insert(mtzManagers.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }

    logged << "Mtz count: " << mtzManagers.size() << std::endl;
    Logger::mainLogger->addStream(&logged);
}

void MtzRefiner::singleThreadRead(vector<std::string> lines,
                                  vector<MtzPtr> *mtzManagers, int offset)
{
    int maxThreads = FileParser::getMaxThreads();

    int end = imageMax(lines.size());

    int skip = FileParser::getKey("IMAGE_SKIP", 0);
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    double tolerance = FileParser::getKey("ACCEPTABLE_UNIT_CELL_TOLERANCE", 0.0);

    bool checkingUnitCell = false;

    if (unitCell.size() > 0 && tolerance > 0.0)
        checkingUnitCell = true;

    if (skip > lines.size())
    {
        std::cout << "Mtz skip beyond Mtz count" << std::endl;
        exit(1);
    }

    if (skip > 0)
    {
        std::cout << "Skipping " << skip << " lines" << std::endl;
    }

    for (int i = skip + offset; i < end; i += maxThreads)
    {
        std::ostringstream log;

        vector<std::string> components = FileReader::split(lines[i], ' ');

        if (components.size() < 10)
            continue;

        std::string mtzName = components[0] + ".mtz";
        if (!FileReader::exists(mtzName))
        {
            log << "Skipping file " << mtzName << std::endl;
            continue;
        }

        log << "Loading file " << mtzName << std::endl;

        double matrix[9];

        readMatrix(matrix, lines[i]);

        MtzPtr newManager = MtzPtr(new MtzManager());

        newManager->setFilename(mtzName.c_str());
        newManager->setMatrix(matrix);
        newManager->loadReflections(1);
        newManager->setSigmaToUnity();
        newManager->loadParametersMap();

        if (newManager->reflectionCount() > 0)
        {
            if (checkingUnitCell && newManager->checkUnitCell(unitCell[0], unitCell[1], unitCell[2], tolerance))
            {
                mtzManagers->push_back(newManager);
            }
            else if (!checkingUnitCell)
            {
                mtzManagers->push_back(newManager);
            }
            else
            {
                log << "Skipping file " << mtzName << " due to poor unit cell" << std::endl;
            }
        }

        Logger::mainLogger->addStream(&log);
    }
}

void MtzRefiner::readMatrix(double (&matrix)[9], std::string line)
{
    vector<std::string> components = FileReader::split(line, ' ');

    for (int j = 1; j <= 9; j++)
    {
        std::string component = components[j];

        /* Locate the substring to replace. */
        int index = (int)component.find(std::string("*^"), 0);
        if (index != std::string::npos)
        {
            component.replace(index, 2, std::string("e"));
        }

        double matVar = atof(component.c_str());
        matrix[j - 1] = matVar;
    }
}

// MARK: Merging

void MtzRefiner::correlationAndInverse(bool shouldFlip)
{
    if (MtzManager::getReferenceManager() == NULL)
    {
        MtzManager::setReference(reference);
    }


    for (int i = 0; i < mtzManagers.size(); i++)
    {
        double correl = mtzManagers[i]->correlation(true);
        double invCorrel = correl;

        if (MtzManager::getReferenceManager()->ambiguityCount() > 1)
        {
            mtzManagers[i]->setActiveAmbiguity(1);
            invCorrel = mtzManagers[i]->correlation(true);
            mtzManagers[i]->setActiveAmbiguity(0);

            if (invCorrel > correl && shouldFlip)
                mtzManagers[i]->setActiveAmbiguity(1);
        }
        double newCorrel = mtzManagers[i]->correlation(true);
        mtzManagers[i]->setRefCorrelation(newCorrel);

        std::cout << mtzManagers[i]->getFilename() << "\t" << correl << "\t"
        << invCorrel << std::endl;
    }
}

void MtzRefiner::merge()
{
    setupFreeMillers();

    loadInitialMtz();
    MtzManager::setReference(this->reference);
    readMatricesAndMtzs();

    bool partialityRejection = FileParser::getKey("PARTIALITY_REJECTION", false);

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->excludeFromLogCorrelation();
        if (partialityRejection)
            mtzManagers[i]->excludePartialityOutliers();
    }

    correlationAndInverse(true);

    double correlationThreshold = FileParser::getKey("CORRELATION_THRESHOLD",
                                                     CORRELATION_THRESHOLD);

    vector<MtzPtr> idxOnly;

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if ((mtzManagers[i]->getActiveAmbiguity() == 0))
            idxOnly.push_back(mtzManagers[i]);
    }


    bool mergeAnomalous = FileParser::getKey("MERGE_ANOMALOUS", false);
    int scalingInt = FileParser::getKey("SCALING_STRATEGY",
                                        (int) SCALING_STRATEGY);
    ScalingType scaling = (ScalingType) scalingInt;


    MtzGrouper *grouper = new MtzGrouper();
    grouper->setScalingType(scaling);
    grouper->setWeighting(WeightTypePartialitySigma);
    grouper->setMtzManagers(mtzManagers);

    grouper->setCutResolution(false);

    grouper->setCorrelationThreshold(correlationThreshold);

    MtzManager *mergedMtz = NULL;
    MtzManager *unmergedMtz = NULL;
    grouper->merge(&mergedMtz, NULL, -1, mergeAnomalous);
    mergedMtz->writeToFile("remerged.mtz");

    delete unmergedMtz;
    delete mergedMtz;
    delete grouper;
}

// MARK: Integrating images

void MtzRefiner::integrateImagesWrapper(MtzRefiner *object,
                                        vector<MtzPtr> *&mtzSubset, int offset, bool orientation)
{
    object->integrateImages(mtzSubset, offset, orientation);
}

void MtzRefiner::integrateImages(vector<MtzPtr> *&mtzSubset,
                                 int offset, bool orientation)
{
    int maxThreads = FileParser::getMaxThreads();

    bool refineDistances = FileParser::getKey("REFINE_DISTANCES", false);

    for (int i = offset; i < images.size(); i += maxThreads)
    {
        std::ostringstream logged;
        logged << "Integrating image " << i << std::endl;
        Logger::mainLogger->addStream(&logged);

        if (refineDistances)
            images[i]->refineDistances();

        if (orientation)
            images[i]->refineOrientations();

        vector<MtzPtr> mtzs = images[i]->currentMtzs();
        mtzSubset->insert(mtzSubset->end(), mtzs.begin(), mtzs.end());
        images[i]->dropImage();
    }
}

void MtzRefiner::loadImageFiles()
{
    loadPanels();

    if ((!isFromPython()) || (isFromPython() && images.size() == 0))
    {
        std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                                  std::string(""));
        if (filename == "")
        {
            std::cout
            << "Filename list has not been provided for integration. Please provide under keyword ORIENTATION_MATRIX_LIST."
            << std::endl;
            exit(1);
        }

        readMatricesAndImages(&filename);
    }
    else
    {
        int num = (int)images.size();
        std::cout << "There have been " << num << " images loaded before starting integration command." << std::endl;

        if (num == 0)
        {
            std::cout << "WARNING! You have not loaded any images from Python! Please load at least one image from the python command line. Aborting integration." << std::endl;
            return;
        }

        std::cout << "Applying parameters to images." << std::endl;

   //     applyParametersToImages();
    }
}

void MtzRefiner::integrate()
{
    bool orientation = FileParser::getKey("REFINE_ORIENTATIONS", false);

    loadImageFiles();

    int crystals = 0;

    for (int i = 0; i < images.size(); i++)
    {
        crystals += images[i]->IOMRefinerCount();
    }

    std::cout << images.size() << " images with " << crystals << " crystal orientations." << std::endl;

    mtzManagers.clear();

    int maxThreads = FileParser::getMaxThreads();

    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(integrateImagesWrapper, this,
                                               &managerSubsets[i], i, orientation);
        threads.add_thread(thr);
    }

    threads.join_all();

    Panel::finaliseMillerArrays();

    int total = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }

    std::cout << "N: Total images loaded: " << total << std::endl;

    mtzManagers.reserve(total);
    int lastPos = 0;

    for (int i = 0; i < maxThreads; i++)
    {
        mtzManagers.insert(mtzManagers.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }

    writeNewOrientations(false, true);

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->writeToDat();
    }

    integrationSummary();
}


void MtzRefiner::integrationSummary()
{
    std::ostringstream refineSummary;

    refineSummary << IOMRefiner::refinementSummaryHeader() << std::endl;

    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->IOMRefinerCount(); j++)
        {
            refineSummary << images[i]->getIOMRefiner(j)->refinementSummary() << std::endl;
        }
    }

    std::string summaryString = refineSummary.str();
    Logger::mainLogger->addStream(&refineSummary);

    std::replace(summaryString.begin(), summaryString.end(), '\t', ',');

    std::ofstream summaryCSV;
    summaryCSV.open("integration.csv");
    summaryCSV << summaryString;
    summaryCSV.close();

    Logger::mainLogger->addString("Written integration summary to integration.csv");

}

void MtzRefiner::loadPanels()
{
    std::string panelList = FileParser::getKey("PANEL_LIST", std::string(""));

    panelParser = new PanelParser(panelList);

    if (panelList != "" && Panel::panelCount() == 0)
    {
        std::cout << "Loading panels" << std::endl;
        panelParser->parse();
        hasPanelParser = true;
    }
    else if (panelList == "")
    {
        std::cout << "No panel list provided. Continuing regardless..."
        << std::endl;
    }
    else if (Panel::panelCount() > 0)
    {
        std::cout << "Panels already present" << std::endl;
    }
}

void MtzRefiner::findSteps()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        std::string mtzFile = mtzManagers[i]->getFilename();
        std::string baseName = getBaseFilename(mtzFile);
        std::string csvHK = "csv-hk-" + baseName + ".csv";
        std::string csvHW = "csv-hw-" + baseName + ".csv";
        std::string csvKW = "csv-kw-" + baseName + ".csv";

        mtzManagers[i]->gridSearch();

        mtzManagers[i]->findSteps(PARAM_HROT, PARAM_KROT, csvHK);
        mtzManagers[i]->findSteps(PARAM_HROT, PARAM_WAVELENGTH, csvHW);
        mtzManagers[i]->findSteps(PARAM_KROT, PARAM_WAVELENGTH, csvKW);
    }
}

void MtzRefiner::refineDetectorGeometry()
{
    if (hasPanelParser)
    {
        for (int i = 0; i < 5; i++)
        {
            Logger::mainLogger->addString("Refining beam center and detector distance.");
            Panel::expectedBeamCentre();
            Panel::refineDetectorDistance();

            Panel::clearAllMillers();
            integrate();
        }
    }
    else
    {
        Logger::mainLogger->addString("Individual panel information has not been set up.");
    }

    Panel::plotAll(PlotTypeAbsolute);
}

// MARK: Dry integrating

/*
// ARK: Find spots

void MtzRefiner::findSpotsWrapper(MtzRefiner *object, int offset, bool orientation)
{
    object->findSpotsThread(offset);
}

void MtzRefiner::findSpotsThread(offset)
{

}

void MtzRefiner::findSpots()
{
    this->readMatricesAndImages();
    loadPanels();
    std::cout << "N: Total images loaded: " << images.size() << std::endl;


}
*/

// MARK: indexing

void MtzRefiner::writeNewOrientations(bool includeRots, bool detailed)
{
    std::ofstream mergeMats;
    std::ofstream refineMats;
    std::ofstream integrateMats;

    std::string filename = FileParser::getKey("NEW_MATRIX_LIST", std::string("new_orientations.dat"));

    mergeMats.open("merge-" + filename);
    refineMats.open("refine-" + filename);
    integrateMats.open("integrate-" + filename);

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzPtr manager = mtzManagers[i];

        // write out matrices etc.
        std::string imgFilename = manager->filenameRoot();

        MatrixPtr matrix = manager->getMatrix()->copy();

        if (includeRots)
        {
            double hRad = manager->getHRot() * M_PI / 180;
            double kRad = manager->getKRot() * M_PI / 180;

            matrix->rotate(hRad, kRad, 0);
        }

        std::string prefix = "img-";

        if (imgFilename.substr(0, 3) == "img")
            prefix = "";


        if (!detailed)
        {
            refineMats << prefix << imgFilename << " ";
            std::string description = matrix->description();
            refineMats << description << std::endl;
        }
        else
        {
            refineMats << "image " << prefix << imgFilename << std::endl;
            mergeMats << "image ref-" << prefix << imgFilename << std::endl;
            std::string description = matrix->description(true);
            refineMats << description << std::endl;
            mergeMats << description << std::endl;

            if (hasRefined)
            {
                refineMats << manager->getParamLine() << std::endl;
            }
        }
    }

    for (int i = 0; i < images.size(); i++)
    {
        ImagePtr image = images[i];

        // write out matrices etc.
        std::string imgFilename = image->filenameRoot();

        integrateMats << "image " << imgFilename << std::endl;

        for (int j = 0; j < image->IOMRefinerCount(); j++)
        {
            MatrixPtr matrix = image->getIOMRefiner(j)->getMatrix()->copy();

            matrix->rotate(0, 0, -M_PI / 2);
            std::string desc90 = matrix->description(true);
            integrateMats << desc90 << std::endl;
        }

        if (image->getSpotsFile().length() > 0)
        {
            integrateMats << "spots " << image->getSpotsFile() << std::endl;
        }
    }

    Logger::mainLogger->addString("Written to filename set: " + filename);

    refineMats.close();
    mergeMats.close();
    integrateMats.close();
}

void MtzRefiner::index()
{
    loadPanels();
    this->readMatricesAndImages();
    std::cout << "N: Total images loaded: " << images.size() << std::endl;

    if (!indexManager)
        indexManager = new IndexManager(images);

    indexManager->index();

    mtzManagers = indexManager->getMtzs();

    writeNewOrientations();
    integrationSummary();
}

void MtzRefiner::combineLists()
{
    loadPanels();
    this->readMatricesAndImages();
    std::cout << "N: Total images loaded in file 1: " << images.size() << std::endl;

    std::string secondList = FileParser::getKey("SECOND_MATRIX_LIST", std::string(""));
    std::vector<ImagePtr> secondImages;

    this->readMatricesAndImages(&secondList, true, &secondImages);
    std::cout << "N: Total images loaded in file 2 (" << secondList << "): " << secondImages.size() << std::endl;

    if (!indexManager)
        indexManager = new IndexManager(images);

    indexManager->setMergeImages(secondImages);
    indexManager->combineLists();

    writeNewOrientations();
    integrationSummary();
}

void MtzRefiner::indexFromScratch()
{
    loadPanels();
    this->readMatricesAndImages();
    std::cout << "N: Total images loaded: " << images.size() << std::endl;

    if (!indexManager)
        indexManager = new IndexManager(images);

    indexManager->indexFromScratch();

}

void MtzRefiner::powderPattern()
{
    loadPanels();
    this->readMatricesAndImages();

    if (!indexManager)
        indexManager = new IndexManager(images);

    indexManager->powderPattern();
}

void MtzRefiner::indexingParameterAnalysis()
{
    if (!indexManager)
        indexManager = new IndexManager(images);

    indexManager->indexingParameterAnalysis();
}

// MARK: Miscellaneous

void MtzRefiner::refineMetrology()
{
    Panel::printToFile("new_panels.txt");
    Panel::plotAll(PlotTypeAbsolute);
}

void MtzRefiner::xFiles()
{
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              std::string(""));

    if (filename == "")
    {
        std::cout << "Orientation matrix list has not been provided.";
        exit(1);
    }

    readXFiles(filename);

    loadInitialMtz();
}

void MtzRefiner::readXFiles(std::string filename)
{
    if (mtzManagers.size() > 0)
    {
        return;
    }

    const std::string contents = FileReader::get_file_contents(filename.c_str());

    vector<std::string> lines = FileReader::split(contents, '\n');

    for (int i = 0; i < lines.size(); i++)
    {
        vector<std::string> filenames = FileReader::split(lines[i], ' ');

        XManager *xManager = new XManager();
        xManager->setFilenames(filenames);
        xManager->loadReflections(PartialityModelFixed);
        xManager->loadParametersMap();

        MtzPtr mtz = MtzPtr();
        mtz.reset(xManager);

        mtzManagers.push_back(mtz);
    }

    std::ostringstream logged;
    logged << "Mtz count: " << mtzManagers.size() << std::endl;
    Logger::mainLogger->addStream(&logged);
}

void MtzRefiner::polarisationGraph()
{
    GraphDrawer drawer = GraphDrawer(reference);

#ifdef MAC
    drawer.plotPolarisation(mtzManagers);
#endif
}

MtzRefiner::~MtzRefiner()
{
//    std::cout << "Deallocating MtzRefiner." << std::endl;

    delete reference;
    delete panelParser;
    mtzManagers.clear();

    images.clear();

    // TODO Auto-generated destructor stub
}

void MtzRefiner::removeSigmaValues()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->setSigmaToUnity();
    }

    if (MtzManager::getReferenceManager() != NULL)
        MtzManager::getReferenceManager()->setSigmaToUnity();
}


void MtzRefiner::displayIndexingHands()
{
    std::ostringstream idxLogged, invLogged;

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        std::ostringstream &which = mtzManagers[i]->getActiveAmbiguity() == 0 ? idxLogged : invLogged;

        which << mtzManagers[i]->getFilename() << std::endl;
    }

    Logger::mainLogger->addString("**** Indexing hands ****\n - First indexing hand");
    Logger::mainLogger->addStream(&idxLogged);
    Logger::mainLogger->addString("- Second indexing hand");
    Logger::mainLogger->addStream(&invLogged);

}

void MtzRefiner::applyUnrefinedPartiality()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->applyUnrefinedPartiality();
    }
}

void MtzRefiner::refineDistances()
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->IOMRefinerCount(); j++)
        {
            images[i]->getIOMRefiner(j)->refineDetectorAndWavelength(reference);
        }
    }
}

void MtzRefiner::orientationPlot()
{
    GraphDrawer drawer = GraphDrawer(reference);
    drawer.plotOrientationStats(mtzManagers);
}

void MtzRefiner::addMatrixToLastImage(scitbx::mat3<double> unit_cell, scitbx::mat3<double> rotation)
{
    ImagePtr lastImage = images.back();
    MatrixPtr newMat = MatrixPtr(new Matrix(unit_cell, rotation));
    lastImage->setUpIOMRefiner(newMat);

    double *lengths = new double[3];
    newMat->unitCellLengths(&lengths);

    std::cout << "Added crystal of unit cell dimensions " << lengths[0] << ", " << lengths[1] << ", " << lengths[2] << " Å." << std::endl;

    delete [] lengths;
}

void MtzRefiner::loadDxtbxImage(std::string imageName, vector<int> imageData, double distance, double wavelength)
{
    ImagePtr newImage = ImagePtr(new Image(imageName, wavelength, distance));
    newImage->setImageData(imageData);


    images.push_back(newImage);
    std::cout << "Loaded image " << imageName << std::endl;
}

int MtzRefiner::imageSkip(size_t totalCount)
{
    int skip = FileParser::getKey("IMAGE_SKIP", 0);

    if (skip > totalCount)
    {
        std::cout << "Image skip beyond image count" << std::endl;
        exit(1);
    }

    return skip;
}
