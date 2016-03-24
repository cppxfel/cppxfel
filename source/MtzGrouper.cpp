/*
 * MtzGrouper.cpp
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#include "MtzGrouper.h"

#include "Scaler.h"
#include <cmath>
#include "StatisticsManager.h"
#include "misc.h"
#include "csymlib.h"
#include <boost/thread/thread.hpp>
#include "Miller.h"
#include "Holder.h"
#include "lbfgs_scaling.h"
#include "ccp4_general.h"
#include "ccp4_parser.h"
#include <fstream>

#include "FileParser.h"
#include "GraphDrawer.h"
#include "FreeMillerLibrary.h"


MtzGrouper::MtzGrouper()
{
	correlationThreshold = 0;
	excludeWorst = true;
	weighting = WeightTypePartialitySigma;
	scalingType = ScalingTypeAverage;
	acceptableResolution = 1;
	cutResolution = false;
    expectedResolution = FileParser::getKey("MAX_RESOLUTION_ALL", 1.6);
    usingNewRefinement = FileParser::getKey("MASS_SCALING", false);
    
    exclusionByCCHalf = FileParser::getKey("EXCLUSION_BY_CC_HALF", false);
        
}

MtzGrouper::~MtzGrouper()
{

}

bool MtzGrouper::isMtzAccepted(MtzPtr mtz)
{
    if (!excludeWorst)
        return true;
    
    int minimumReflectionCutoff = FileParser::getKey(
                                                     "MINIMUM_REFLECTION_CUTOFF",
                                                     MINIMUM_REFLECTION_CUTOFF);

    double refPartCorrelThreshold = FileParser::getKey(
                                                     "PARTIALITY_CORRELATION_THRESHOLD",
                                                     0.0);

    if (refPartCorrelThreshold > 0)
    {
        if (mtz->getRefPartCorrel() < refPartCorrelThreshold)
        {
            Logger::mainLogger->addString("Rejecting due to low partiality correlation", LogLevelDetailed);
            return false;
        }
    }
    
    if (mtz->accepted() < minimumReflectionCutoff)
    {
        Logger::mainLogger->addString("Rejecting due to not reaching minimum number of reflections", LogLevelDetailed);
        return false;
    }
    
    double refCorrelation = mtz->getRefCorrelation();
    
    if (refCorrelation < 0 || refCorrelation == 1)
    {
        Logger::mainLogger->addString("Rejecting due to suspicious correlation with reference", LogLevelDetailed);
        
        return false;
    }
    
    double minimumRSplit = FileParser::getKey("R_FACTOR_THRESHOLD", 0.0);
    
    bool needsRSplit = mtz->getReferenceManager() != NULL;
    double rSplit = 0;
    
    if (needsRSplit)
        rSplit = mtz->rSplit(0, 0);
    
    if ((refCorrelation < correlationThreshold && excludeWorst
         && !mtz->isFreePass()))
    {
        Logger::mainLogger->addString("Rejecting due to poor correlation with reference", LogLevelDetailed);
        
        return false;
    }
    
    if (mtz->isRejected())
    {
        Logger::mainLogger->addString("Rejecting due to true rejection flag", LogLevelDetailed);
        
        return false;
    }
    
    if (needsRSplit && minimumRSplit > 0 && rSplit > minimumRSplit)
    {
        Logger::mainLogger->addString("Rejecting due to R split being too high", LogLevelDetailed);
        
        return false;
    }
    
    return true;
}

void MtzGrouper::setMtzManagers(const vector<MtzPtr>& mtzManagers)
{
	this->mtzManagers = mtzManagers;

	double averageCorrelation = 0;

	for (int i = 0; i < mtzManagers.size(); i++)
	{
		averageCorrelation += mtzManagers[i]->getRefCorrelation();
	}

	averageCorrelation /= mtzManagers.size();

	this->setCorrelationThreshold(averageCorrelation - 0.05);
}

void MtzGrouper::merge(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		int cycle, bool anom)
{
	logged << "N: ==== Merge cycle " << cycle << " ====" << std::endl;
	logged << "N: Scaling type: ";

	switch (scalingType)
	{
	case ScalingTypeAverage:
		logged << "Average merge" << std::endl;
		break;
	case ScalingTypeReference:
		logged << "Against reference (gradient)" << std::endl;
		break;
	case ScalingTypeReferenceLeastSquares:
		logged << "Against reference (least squares)" << std::endl;
		break;
	case ScalingTypeMinimizeRMerge:
		logged << "Minimize R merge" << std::endl;
		break;
	case ScalingTypeBFactor:
		logged << "B factor and scale" << std::endl;
		break;
	case ScalingTypeResolutionShells:
		logged << "Gradients per resolution shell" << std::endl;
		break;
	default:
		logged << "Unsupported option" << std::endl;
		break;
	}

	double averageCorrelation = 0;
	double averageAboveCutoff = 0;
	int aboveCutoffNum = 0;
	double rotationCorrection = 0;
    
    int rotMode = FileParser::getKey("ROTATION_MODE", 0);
    RotationMode mode = (RotationMode)rotMode;

    
    logged << "Filename\tCorrel\tRsplit\tPartcorrel\tRefcount\tMosaicity\tWavelength\tBandwidth\t";
    
    logged << (mode == RotationModeHorizontalVertical ? "hRot\tkRot\t" : "aRot\tbRot\tcRot\t");
    
    logged << "rlpSize\texp\tcellA\tcellB\tcellC" << std::endl;

	for (int i = 0; i < mtzManagers.size(); i++)
	{
		double correl = mtzManagers[i]->getRefCorrelation();
        if (correl != -1)
            averageCorrelation += correl;

		double hRot = mtzManagers[i]->getHRot();
		double kRot = mtzManagers[i]->getKRot();
        
        double aRot = mtzManagers[i]->getARot();
        double bRot = mtzManagers[i]->getBRot();
        double cRot = mtzManagers[i]->getCRot();

		double correction = sqrt(hRot * hRot + kRot * kRot);
        
		rotationCorrection += correction;

		if (correl > correlationThreshold)
		{
			averageAboveCutoff += correl;
			aboveCutoffNum++;
		}

        double a, b, c;
        mtzManagers[i]->getMatrix()->orientationMatrixUnitCell(&a, &b, &c);
        
        double rSplit = 0;
        if (MtzManager::getReferenceManager() != NULL)
            rSplit = mtzManagers[i]->rSplit(0, expectedResolution, true, true);
        
        double partCorrel = mtzManagers[i]->getRefPartCorrel();
        
        double *cellDims = new double[3];
        mtzManagers[i]->getMatrix()->unitCellLengths(&cellDims);
        
        
        
		logged << mtzManagers[i]->getFilename() << "\t" << correl << "\t" << rSplit << "\t" << partCorrel << "\t"
				<< mtzManagers[i]->accepted() << "\t"
				<< mtzManagers[i]->getMosaicity() << "\t"
				<< mtzManagers[i]->getWavelength() << "\t"
				<< mtzManagers[i]->getBandwidth() << "\t";
        if (mode == RotationModeHorizontalVertical)
            logged << hRot << "\t" << kRot << "\t";
        else if (mode == RotationModeUnitCellABC)
            logged << aRot << "\t" << bRot << "\t" << cRot << "\t";
        
        logged << mtzManagers[i]->getSpotSize() << "\t"
            << mtzManagers[i]->getExponent() << "\t" << cellDims[0] << "\t" << cellDims[1] << "\t" << cellDims[2] << std::endl;
        
        delete [] cellDims;
	}
    
    std::string tabbedParams = logged.str();
    std::replace(tabbedParams.begin(), tabbedParams.end(), '\t', ',');
    
    std::ofstream paramLog;
    std::string paramLogName = "params_cycle_" + i_to_str(cycle) + ".csv";
    paramLog.open(paramLogName);
    paramLog << tabbedParams << std::endl;
    paramLog.close();
    
    logged << "Written parameter table to " << paramLogName << "." << std::endl;

	averageCorrelation /= mtzManagers.size();
	rotationCorrection /= mtzManagers.size();
	averageAboveCutoff /= aboveCutoffNum;

	logged << "N: Average correlation per image: " << averageCorrelation
			<< std::endl;
	logged << "N: Average correlation for those above threshold: "
			<< averageAboveCutoff << std::endl;
	logged << "N: Average rotation correction in degrees: "
    << rotationCorrection << std::endl;
    
    sendLog();
    
    if (MtzManager::getReferenceManager() != NULL)
    {
        double refScale = 1000 / MtzManager::getReferenceManager()->averageIntensity();
        MtzManager::getReferenceManager()->applyScaleFactor(refScale);
    }
    
	for (int i = 0; i < mtzManagers.size(); i++)
	{
        double scale = 1;

        if (scalingType == ScalingTypeAverage)
		{
			scale = 1000 / mtzManagers[i]->averageIntensity();
		}
		else if (scalingType == ScalingTypeReference
                 || scalingType == ScalingTypeMinimizeRMerge)
		{
			scale = mtzManagers[i]->gradientAgainstManager(
					MtzManager::getReferenceManager(), false);
       //     std::cout << mtzManagers[i]->getFilename() << " " << scale << std::endl;
		}
		else if (scalingType == ScalingTypeBFactor)
		{
			double newScale = 1;
			double bFactor = 0;

			mtzManagers[i]->bFactorAndScale(&newScale, &bFactor);
			mtzManagers[i]->applyScaleFactor(newScale, bFactor);
		}
		else if (scalingType == ScalingTypeResolutionShells)
		{
			mtzManagers[i]->applyScaleFactorsForBins();
		}

		mtzManagers[i]->applyScaleFactor(scale);
	}
/*
#ifdef MAC
	GraphDrawer drawer = GraphDrawer(MtzManager::getReferenceManager());

	drawer.resolutionStatsPlot(mtzManagers, "resolution");
	drawer.resolutionStatsPlot(mtzManagers, "intensity_ref", GraphMap(), true, false);
	drawer.resolutionStatsPlot(mtzManagers, "intensity_img", GraphMap(), true, true);
#endif*/

	if (scalingType == ScalingTypeMinimizeRMerge)
	{
        vector<MtzPtr> acceptedMtzs;
        
        for (int i = 0; i < mtzManagers.size(); i++)
        {
            if (isMtzAccepted(mtzManagers[i]))
                acceptedMtzs.push_back(mtzManagers[i]);
        }
        
		Lbfgs_Scaling scaling = Lbfgs_Scaling(acceptedMtzs);
		scaling.run();
	}

	logged << "Altered scales." << std::endl;

	MtzManager *idxMerge = NULL;
	MtzManager **unmerged = NULL;
	MtzManager *invMerge = NULL;

    std::string idxName = std::string("half1Merge.mtz");
    std::string invName = std::string("half2Merge.mtz");
    std::string unmergedName = std::string("unmerged.mtz");
    
    if (cycle >= 0)
    {
        unmergedName = std::string("unmerged") + i_to_str(cycle) + std::string(".mtz");
        idxName = std::string("half1Merge") + i_to_str(cycle) + std::string(".mtz");
        invName = std::string("half2Merge") + i_to_str(cycle) + std::string(".mtz");
    }
    
	if (anom == false)
	{
        time_t startcputime;
        time(&startcputime);
        
        logged << "**** MERGING ALL DATA ****" << std::endl;
        sendLog();
        merge(mergeMtz, unmergedMtz, false, true, &unmergedName);
		
        time_t endcputime;
        time(&endcputime);
        
        time_t difference = endcputime - startcputime;
        double seconds = difference;
        
        int finalSeconds = (int) seconds % 60;
        int minutes = seconds / 60;
        
        logged << "N: Clock time " << minutes << " minutes, " << finalSeconds << " seconds to merge full data set" << std::endl;
        
        logged << "**** MERGING HALF DATA (1) ****" << std::endl;
        sendLog();
        merge(&idxMerge, unmerged, true, false);

        logged << "**** MERGING HALF DATA (2) ****" << std::endl;
        sendLog();
        merge(&invMerge, unmerged, false, false);
	}
	else
	{
		mergeAnomalous(mergeMtz, unmergedMtz, false, true);

		mergeAnomalous(&idxMerge, unmerged, true, false);
		mergeAnomalous(&invMerge, unmerged, false, false);
	}

	idxMerge->writeToFile(idxName, true);
	invMerge->writeToFile(invName, true);
    
	logged << "N: === R split ===" << std::endl;
    sendLog();
	idxMerge->rSplitWithManager(invMerge, false, false, 0, expectedResolution, 20, NULL, true);
	logged << "N: === CC half ===" << std::endl;
    sendLog();
	idxMerge->correlationWithManager(invMerge, false, false, 0,
			expectedResolution, 20, NULL, true);
    
    if (FreeMillerLibrary::active())
    {
        logged << "N: === R split (free only) ===" << std::endl;
        sendLog();
        idxMerge->rSplitWithManager(invMerge, false, false, 0, expectedResolution, 20, NULL, true, true);
        logged << "N: === CC half (free only) ===" << std::endl;
        sendLog();
        idxMerge->correlationWithManager(invMerge, false, false, 0,
                                         expectedResolution, 20, NULL, false, true);
    }
        
    sendLog();

	delete idxMerge;
	delete invMerge;
}

void MtzGrouper::mergeWrapper(void *object, MtzManager **mergeMtz,
		MtzManager **unmergedMtz, bool firstHalf, bool all, bool anom)
{
	if (!anom)
	{
		static_cast<MtzGrouper *>(object)->merge(mergeMtz, unmergedMtz,
				firstHalf, all);
	}
	else
	{
		static_cast<MtzGrouper *>(object)->mergeAnomalous(mergeMtz, unmergedMtz,
				firstHalf, all);
	}
}

void MtzGrouper::checkCCHalf(vector<MtzPtr> *managers, int offset, int *total)
{
    Scaler *scaler = Scaler::getScaler();
    int accepted = 0;
    int maxThreads = FileParser::getMaxThreads();
    
    if (scaler != NULL)
    {
        for (int i = offset; i < managers->size(); i += maxThreads)
        {
            accepted += scaler->mtzIsBeneficial((*managers)[i]);
        }
    }
    
    *total += accepted;
}

void MtzGrouper::merge(MtzManager **mergeMtz, MtzManager **unmergedMtz,
                       bool firstHalf, bool all, std::string *unmergedName)
{
	*mergeMtz = new MtzManager();
    
    std::string filename = std::string("merged_") + (all ? std::string("all_") : std::string("half_")) + (firstHalf ? std::string("0") : std::string("1"));
    
	(*mergeMtz)->setFilename(filename);
	(*mergeMtz)->copySymmetryInformationFromManager(mtzManagers[0]);

	int start = firstHalf ? 0 : (int)mtzManagers.size() / 2;
	int end = firstHalf ? (int)mtzManagers.size() / 2 : (int)mtzManagers.size();

	if (all)
	{
		start = 0;
		end = (int)mtzManagers.size();
	}

	int mtzCount = groupMillers(mergeMtz, unmergedMtz, start, end);

    Scaler *scaler = Scaler::getScaler(mtzManagers, mergeMtz);
    
    if (usingNewRefinement)
    {
        scaler->minimizeRMerge();
    }
    
    int total = 0;
    
    if (exclusionByCCHalf && all)
    {
        boost::thread_group threads;
        
        int maxThreads = FileParser::getMaxThreads();
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(checkCCHalf, &mtzManagers, i, &total);
            threads.add_thread(thr);
                    }
        
        threads.join_all();
    }

 //   std::cout << "N: Accepted " << total << " due to increase in CC half" << std::endl;

    
	if (unmergedMtz != NULL)
	{
		(*unmergedMtz) = &*(*mergeMtz)->copy();
	}
    
    if (unmergedName != NULL)
    {
        (*mergeMtz)->writeToFile(*unmergedName);
    }

	logged << "Merging miller indices" << std::endl;

	bool outlier_rejection = FileParser::getKey("OUTLIER_REJECTION",
	REJECTING_MILLERS);

	mergeMillers(mergeMtz, all && outlier_rejection, mtzCount);

	unflipMtzs();
}

void MtzGrouper::mergeAnomalous(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		bool firstHalf, bool all, std::string filename)
{
	MtzManager *positive = new MtzManager();
	MtzManager *negative = new MtzManager();
	positive->setFilename("positive_pair");
	negative->setFilename("negative_pair");
	positive->copySymmetryInformationFromManager(mtzManagers[0]);
	negative->copySymmetryInformationFromManager(mtzManagers[0]);

	*mergeMtz = new MtzManager();
	(*mergeMtz)->setFilename(filename);
	(*mergeMtz)->copySymmetryInformationFromManager(mtzManagers[0]);

	int start = firstHalf ? 0 : (int)mtzManagers.size() / 2;
	int end = firstHalf ? (int)mtzManagers.size() / 2 : (int)mtzManagers.size();

	if (all)
	{
		start = 0;
		end = (int)mtzManagers.size();
	}

	int mtzCount = groupMillersWithAnomalous(&positive, &negative, start, end);

	bool outlier_rejection = FileParser::getKey("OUTLIER_REJECTION", REJECTING_MILLERS);

	mergeMillers(&positive, outlier_rejection && all, mtzCount);
	mergeMillers(&negative, outlier_rejection && all, mtzCount);

	std::cout << "Merging miller indices" << std::endl;

	differenceBetweenMtzs(mergeMtz, &positive, &negative);

	if (all)
	{
		writeAnomalousMtz(&positive, &negative, "anomalous_diff.mtz");
	}

	delete negative;
	delete positive;

	(*mergeMtz)->description();

	unflipMtzs();
}

int MtzGrouper::groupMillers(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		int start, int end)
{
	int mtzCount = 0;
    bool fastMerge = FileParser::getKey("FAST_MERGE", false);
    
    std::map<int, int> flipCounts;
    
    for (int i = 0; i < mtzManagers[0]->ambiguityCount(); i++)
    {
        flipCounts[i] = 0;
    }
    
	for (int i = start; i < end; i++)
	{
        if (!exclusionByCCHalf && !isMtzAccepted(mtzManagers[i]))
        {
            mtzManagers[i]->incrementFailedCount();
            continue;
        }
        
        if (!exclusionByCCHalf)
            mtzManagers[i]->resetFailedCount();
        
        mtzManagers[i]->flipToActiveAmbiguity();
        int ambiguity = mtzManagers[i]->getActiveAmbiguity();
        flipCounts[ambiguity]++;

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->reflectionCount(); j++)
		{
            if (fastMerge && mtzManagers[i]->reflection(j)->acceptedCount() == 0)
                continue;
            
            if (mtzManagers[i]->reflection(j)->getResolution() > cutoffRes)
				continue;
            
            
			long unsigned int refl_id = mtzManagers[i]->reflection(j)->getReflId();

			Reflection *reflection = NULL;
			(*mergeMtz)->findReflectionWithId(refl_id, &reflection);

			if (reflection == NULL)
			{
				Reflection *newReflection = mtzManagers[i]->reflection(j)->copy(false);
				(*mergeMtz)->addReflection(newReflection);
			}
			else
			{
				for (int k = 0; k < mtzManagers[i]->reflection(j)->millerCount();
						k++)
				{
                    MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(k);
                    
                    if (fastMerge && !newMiller->accepted())
                    {
                        continue;
                    }
					reflection->addMiller(newMiller);
                }
			}
		}
	}

	long unsigned int last_refl_id = 0;

	for (int i = 0; i < (*mergeMtz)->reflectionCount(); i++)
	{
		Reflection *reflection = (*mergeMtz)->reflection(i);

		if (reflection->getReflId() == last_refl_id)
		{
			std::cout << "Same" << std::endl;
		}
		if (reflection->getReflId() < last_refl_id)
		{
			std::cout << "Less" << std::endl;
		}

		last_refl_id = reflection->getReflId();
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Flip ratios: ";
    
    for (std::map<int, int>::iterator it = flipCounts.begin(); it != flipCounts.end(); it++)
    {
        std::cout << flipCounts[it->first] << " ";
    }
    
    std::cout << std::endl;
    
	std::cout << "N: Reflections used: " << (*mergeMtz)->reflectionCount() << std::endl;

	return mtzCount;
}

int MtzGrouper::groupMillersWithAnomalous(MtzManager **positive,
		MtzManager **negative, int start, int end)
{
	int flipCount = 0;
	int mtzCount = 0;

    logged << "Grouping Millers for anomalous merge." << std::endl;
    sendLog();
    
	for (int i = start; i < end; i++)
	{
        if (!isMtzAccepted(mtzManagers[i]))
			continue;

        mtzManagers[i]->flipToActiveAmbiguity();

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->reflectionCount(); j++)
		{
			for (int k = 0; k < mtzManagers[i]->reflection(j)->millerCount(); k++)
			{
                bool friedel = false;
				mtzManagers[i]->reflection(j)->miller(k)->positiveFriedel(&friedel);
                
                logged << "Acquired Friedel " << friedel << std::endl;
                sendLog(LogLevelDebug);
                
				MtzManager *friedelMtz = (friedel ? *positive : *negative);

				if (mtzManagers[i]->reflection(j)->getResolution() > cutoffRes)
					continue;

				Reflection *reflection = NULL;
				friedelMtz->findReflectionWithId(
						mtzManagers[i]->reflection(j)->getReflId(), &reflection);

				// there is a bug here which would mis-sort friedel pairs if
				// there were opposing ones of the same symmetry in the same image

				if (reflection == NULL)
				{
					Reflection *newReflection = mtzManagers[i]->reflection(j)->copy(true);
					newReflection->clearMillers();
					MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(k);
					newReflection->addMiller(newMiller);
					friedelMtz->addReflection(newReflection);
					friedelMtz->sortLastReflection();
				}
				else
				{
					{
						MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(
								k);
						reflection->addMiller(newMiller);
					}
				}
			}
		}
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Total flipped: " << flipCount << std::endl;

	return mtzCount;
}

void MtzGrouper::mergeMillers(MtzManager **mergeMtz, bool reject, int mtzCount)
{
	int reflectionCount = 0;
	int millerCount = 0;
	int rejectCount = 0;
	double aveStdErr = 0;
	int aveStdErrCount = 0;
    
    bool recalculateSigma = FileParser::getKey("RECALCULATE_SIGMA", false);
    bool minimumMultiplicity = FileParser::getKey("MINIMUM_MULTIPLICITY", 0);
    
    for (int i = 0; i < (*mergeMtz)->reflectionCount(); i++)
	{
		Reflection *reflection = (*mergeMtz)->reflection(i);

        int accepted = reflection->acceptedCount();
        
        if (accepted <= minimumMultiplicity || accepted == 0)
        {
            (*mergeMtz)->removeReflection(i);
            i--;
            continue;
        }
        sendLog();

		double totalIntensity = 0;
		double totalStdev = 0;

		reflection->merge(weighting, &totalIntensity, &totalStdev, reject);

		double error = totalStdev / totalIntensity;

		if (error == error && totalStdev != 100)
		{
			aveStdErr += fabs(error);
			aveStdErrCount++;

		}

        double totalSigma = reflection->meanSigma() / reflection->meanPartiality();
        
        if (recalculateSigma)
        {
            totalSigma = reflection->mergeSigma();
            
            if (totalSigma == 0)
            {
                totalIntensity = std::nan(" ");
                totalSigma = std::nan(" ");
            }
        }
        
		millerCount += reflection->acceptedCount();
		rejectCount += reflection->rejectCount();
		reflectionCount++;

		int h, k, l = 0;
		MillerPtr firstMiller = (*mergeMtz)->reflection(i)->miller(0);
		ccp4spg_put_in_asu((*mergeMtz)->getLowGroup(), firstMiller->getH(),
				firstMiller->getK(), firstMiller->getL(), &h, &k, &l);

		MillerPtr newMiller = MillerPtr(new Miller((*mergeMtz), h, k, l));

        if (totalSigma == 0)
        {
            std::cout << "Error!" << std::endl;
        }
        
		newMiller->setData(totalIntensity, totalSigma, 1, 0);
		newMiller->setParent((*mergeMtz)->reflection(i));
		(*mergeMtz)->reflection(i)->calculateResolution(*mergeMtz);

		(*mergeMtz)->reflection(i)->clearMillers();
		(*mergeMtz)->reflection(i)->addMiller(newMiller);
	}

	aveStdErr /= aveStdErrCount;

	double multiplicity = (double) millerCount / (double) reflectionCount;
	double aveRejection = (double) rejectCount / (double) mtzCount;

	std::cout << "N: Total MTZs: " << mtzManagers.size() << std::endl;
	std::cout << "N: Multiplicity before merge: " << multiplicity << std::endl;
	std::cout << "N: Rejects per image: " << aveRejection << std::endl;
	std::cout << "N: Average error per reflection: " << aveStdErr << std::endl;

	(*mergeMtz)->insertionSortReflections();
}

void MtzGrouper::writeAnomalousMtz(MtzManager **positive, MtzManager **negative,
		std::string filename)
{
    double doubleCell[6];
	float cell[6], wavelength, fdata[9];

	/* variables for symmetry */
	CCP4SPG *mtzspg = (*positive)->getLowGroup();
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	MTZ *mtzout;
	MTZXTAL *xtal;
	MTZSET *set;
	MTZCOL *colout[9];

	/*  Removed: General CCP4 initializations e.g. HKLOUT on command line */

	(*positive)->getUnitCell(&doubleCell[0], &doubleCell[1], &doubleCell[2], &doubleCell[3], &doubleCell[4],
			&doubleCell[5]);
    
    for (int i = 0; i < 6; i++)
    {
        cell[i] = (float)doubleCell[i];
    }

	wavelength = (*positive)->getWavelength();

	mtzout = MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Anomalous dataset ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = MtzOpenForWrite(filename.c_str());

// then add symm headers...
	for (int i = 0; i < mtzspg->nsymop; ++i)
		CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
	strncpy(ltypex, mtzspg->symbol_old, 1);
	ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
			mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);

// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "XFEL crystal", "XFEL project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "I(+)", "K");
	colout[4] = MtzAddColumn(mtzout, set, "SIGI(+)", "M");
	colout[5] = MtzAddColumn(mtzout, set, "I(-)", "K");
	colout[6] = MtzAddColumn(mtzout, set, "SIGI(-)", "M");
	colout[7] = MtzAddColumn(mtzout, set, "IMEAN", "J");
	colout[8] = MtzAddColumn(mtzout, set, "SIGIMEAN", "Q");

	int num = 0;

	int hits = 0;
	vector<Reflection *> posReflections;
	vector<Reflection *> negReflections;

	(*positive)->findCommonReflections(*negative, posReflections, negReflections, &hits);

	for (int i = 0; i < posReflections.size(); i++)
	{
		double intensityPlus = posReflections[i]->meanIntensity();
		double sigmaPlus = posReflections[i]->meanSigma();

		double intensityMinus = negReflections[i]->meanIntensity();
		double sigmaMinus = negReflections[i]->meanSigma();

        double intensity = (intensityPlus + intensityMinus) / 2;
		double sigma = (sigmaPlus + sigmaMinus) / 2;
        
        if (intensityPlus != intensityPlus && intensityMinus != intensityMinus)
            continue;
        
        if (sigmaPlus < 0.0001)
            sigmaPlus = nan(" ");

        if (sigmaMinus < 0.0001)
            sigmaMinus = nan(" ");

        if (intensityPlus != intensityPlus)
        {
            intensity = intensityMinus;
            sigma = sigmaMinus;
        }
        
        if (intensityMinus != intensityMinus)
        {
            intensity = intensityPlus;
            sigma = sigmaPlus;
        }

		num++;

		int h = negReflections[i]->miller(0)->getH();
		int k = negReflections[i]->miller(0)->getK();
		int l = negReflections[i]->miller(0)->getL();
		int _h, _k, _l;
		ccp4spg_put_in_asu(mtzspg, h, k, l, &_h, &_k, &_l);

		fdata[0] = _h;
		fdata[1] = _k;
		fdata[2] = _l;
		fdata[3] = intensityMinus;
		fdata[4] = sigmaMinus;
		fdata[5] = intensityPlus;
		fdata[6] = sigmaPlus;
		fdata[7] = intensity;
		fdata[8] = sigma;
		ccp4_lwrefl(mtzout, fdata, colout, 9, num);
	}

// print header information, just for info
//	ccp4_lhprt(mtzout, 1);
	MtzPut(mtzout, " ");
	MtzFree(mtzout);
}

void MtzGrouper::differenceBetweenMtzs(MtzManager **mergeMtz,
		MtzManager **positive, MtzManager **negative)
{
	int hits = 0;
	vector<Reflection *> posReflections;
	vector<Reflection *> negReflections;

	(*positive)->findCommonReflections(*negative, posReflections, negReflections, &hits);

	for (int i = 0; i < posReflections.size(); i++)
	{
		Reflection *newReflection = posReflections[i]->copy(true);
		double posInt = posReflections[i]->meanIntensity();
		double negInt = negReflections[i]->meanIntensity();

		double posSigma = posReflections[i]->meanSigma();
		double negSigma = negReflections[i]->meanSigma();

		newReflection->miller(0)->setRawIntensity(posInt - negInt);
		newReflection->miller(0)->setPartiality(1);
		newReflection->miller(0)->setSigma(posSigma + negSigma);

		(*mergeMtz)->addReflection(newReflection);
	}

	(*mergeMtz)->insertionSortReflections();
}

void MtzGrouper::unflipMtzs()
{
	for (int i = 0; i < mtzManagers.size(); i++)
	{
        mtzManagers[i]->resetFlip();
	}
}

void MtzGrouper::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

