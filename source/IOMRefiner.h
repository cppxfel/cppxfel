/*
 * IOMRefiner.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IOMRefiner_H_
#define IOMRefiner_H_

#include "Image.h"
#include "Miller.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "csymlib.h"
#include "Spot.h"
#include "parameters.h"
#include "Logger.h"

class Image;
class Miller;

using namespace CSym;

typedef enum
{
    RefinementTypeOrientationMatrixEarly = 0,
	RefinementTypeDetectorWavelength = 1,
    RefinementTypeOrientationMatrixVeryEarly = 2,
    RefinementTypeOrientationMatrixLate = 3,
	RefinementTypeOrientationMatrixSpots = 4,
    RefinementTypeOrientationMatrixExactSpots = 5,
	RefinementTypeOrientationMatrixRough = 6,
	RefinementTypeOrientationMatrixMedian = 7,
    RefinementTypeOrientationMatrixTotalSignal = 8,
    RefinementTypeOrientationMatrixHighestPeak = 9,
    RefinementTypeOrientationMatrixEarlySeparated = 10,
    RefinementTypeOrientationMatrixPanelStdev = 11,
    RefinementTypeOrientationMatrixStdevOnly = 12,
} RefinementType;

class IOMRefiner
{
private:
	Image *image;
	vector<MillerPtr> millers;
	vector<MillerPtr> nearbyMillers;
    vector<MillerPtr> roughMillers;
	vector<Spot *>spots;
    MatrixPtr matrix;
    std::vector<Match> indexingMatches;

    double minResolution;
    bool roughCalculation;
    bool needsReintegrating;
    CCP4SPG *spaceGroup;
    bool complexUnitCell;
    RotationMode rotationMode;
	vector<double> unitCell;
    bool refineA;
    bool refineB;
    bool refineC;
	double hRot;
	double kRot;
    double lRot;
    double aRot;
    double bRot;
    double cRot;
    double bestHRot;
    double bestKRot;
    double bestLRot;
	int search;
	double testDistance;
	double testWavelength;
	double testSpotSize;
	double initialStep;
    double getTotalStandardDeviation();
	int expectedSpots;
	int searchSize;
	static double intensityThreshold;
    static bool absoluteIntensity;
	double maxResolution;
	MtzManager *reference;
	void calculateNearbyMillers(bool rough);
	RefinementType refinement;
	double lastTotal;
	double lastStdev;
	double testBandwidth;
    double lastScore;
	vector<vector<double> > solutions;
    
    double orientationTolerance;
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);

public:
	IOMRefiner(Image *newImage = NULL, MatrixPtr matrix = MatrixPtr());
    void setComplexMatrix();
    virtual ~IOMRefiner();

    void calculateOnce();
	void checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox = false, bool perfectCalculation = true);
	MtzPtr newMtz(int i);
	void getWavelengthHistogram(vector<double> &wavelengths,
			vector<int> &frequencies, LogLevel level = LogLevelDetailed, int whichAxis = 0);
	double score(int whichAxis = 0, bool silent = false);

	static bool millerReachesThreshold(MillerPtr miller);
	void findSpots();
	static void duplicateSpots(vector<Image *>images);
	void writeDatFromSpots(std::string filename);
	static void scatterSpots(vector<Image *> images);

	void matchMatrixToSpots();
	void matchMatrixToSpots(RefinementType refinement);
	double minimizeParameter(double *meanStep, double *param, int whichAxis = 0);
	void minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double *param1, double *param2);

	void refineDetectorAndWavelength(MtzManager *reference = NULL);
	void refineOrientationMatrix();
	void refineOrientationMatrix(RefinementType refinementType);
	double refineRoundBeamAxis(double start, double end, double wedge, bool allSolutions);
	void refineRoundBeamAxis();

    bool millerWithinBandwidth(MillerPtr miller);
	int getTotalReflections();
	int getTotalReflections(double threshold);
    int getTotalReflectionsWithinBandwidth();
    double medianIntensity();
	int identicalSpotsAndMillers();
	double getTotalIntegratedSignal();

    double getRot(int rotNum);
	void dropMillers();

    bool isGoodSolution();
    double getDetectorDistance();
    double getWavelength();
    std::string refinementSummary();
    static std::string refinementSummaryHeader();
    
    void getBestRots(double *rot1, double *rot2, double *rot3)
    {
        *rot1 = bestHRot;
        *rot2 = bestKRot;
        *rot3 = bestLRot;
    }
    
    void setMatch(Match newMatch1, Match newMatch2)
    {
        indexingMatches.push_back(newMatch1);
        indexingMatches.push_back(newMatch2);
    }
    
    std::vector<Match> getMatch()
    {
        return indexingMatches;
    }
    
    MatrixPtr getMatrix()
    {
        return matrix;
    }
    
    void setMatrix(MatrixPtr matrix)
    {
        this->matrix = matrix;
    }
    
    void setMatrixCopy(MatrixPtr matrix)
    {
        this->matrix = matrix->copy();
    }
    
    double getLastScore()
    {
        return lastScore;
    }
    
	Image*& getImage()
	{
		return image;
	}

	void setImage(Image*& image)
	{
		this->image = image;
	}

	int getSearch() const
	{
		return search;
	}

	void setSearch(int search)
	{
		this->search = search;
	}

	double getHRot() const
	{
		return hRot;
	}

	void setHRot(double rot)
	{
		hRot = rot;
	}

	double getKRot() const
	{
		return kRot;
	}

	void setKRot(double rot)
	{
		kRot = rot;
	}

	double getTestBandwidth() const
	{
		return testBandwidth;
	}

	void setTestBandwidth(double testBandwidth)
	{
		this->testBandwidth = testBandwidth;
	}

	double getTestSpotSize() const
	{
		return testSpotSize;
	}

	void setTestSpotSize(double testSpotSize)
	{
		this->testSpotSize = testSpotSize;
	}

	double getMaxResolution() const
	{
		return maxResolution;
	}

	void setMaxResolution(double maxResolution)
	{
		this->maxResolution = maxResolution;
	}

	int getSearchSize() const
	{
		return searchSize;
	}

	void setSearchSize(int searchSize)
	{
		this->searchSize = searchSize;
	}

	double getInitialStep() const
	{
		return initialStep;
	}

	void setInitialStep(double initialStep = 1)
	{
		this->initialStep = initialStep;
	}

	CCP4SPG*& getSpaceGroup()
	{
		return spaceGroup;
	}

	void setSpaceGroup(CCP4SPG*& spaceGroup)
	{
        if (this->spaceGroup != NULL)
            ccp4spg_free(&this->spaceGroup);
        
		this->spaceGroup = spaceGroup;
	}

	const vector<Spot *>& getSpots() const
	{
		return spots;
	}

	int getExpectedSpots() const
	{
		return expectedSpots;
	}

	void setExpectedSpots(int expectedSpots)
	{
		this->expectedSpots = expectedSpots;
	}

	const vector<vector<double> >& getSolutions() const
	{
		return solutions;
	}

	void setSolutions(const vector<vector<double> >& solutions)
	{
		this->solutions = solutions;
	}

	double getIntensityThreshold() const
	{
		return intensityThreshold;
	}

	void setIntensityThreshold(double intensityThreshold)
	{
		this->intensityThreshold = intensityThreshold;
	}

	vector<double>& getUnitCell()
	{
		return unitCell;
	}

	void setUnitCell(vector<double>& unitCell)
	{
		this->unitCell = unitCell;
	}
    
    void setOrientationTolerance(double newTolerance)
    {
        orientationTolerance = newTolerance;
    }
    
    double getBestLRot()
    {
        return bestLRot;
    }
    
    bool isCalculatingRough()
    {
        return roughCalculation;
    }
    
    void setCalculatingRough(bool rough)
    {
        roughCalculation = rough;
    }
};

#endif /* IOMRefiner_H_ */
