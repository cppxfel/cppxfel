/*
 * IOMRefiner.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IOMRefiner_H_
#define IOMRefiner_H_

#include "Image.h"
#include "parameters.h"
#include "IndexManager.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "csymlib.h"
#include "Spot.h"
#include "Logger.h"

class Miller;

using namespace CSym;

typedef enum
{
    RefinementTypeOrientationMatrixEarly = 0,
	RefinementTypeDetectorWavelength = 1,
    RefinementTypeOrientationMatrixVeryEarly = 2,
    RefinementTypeOrientationMatrixRough = 6,
	RefinementTypeOrientationMatrixHighestPeak = 9,
    RefinementTypeOrientationMatrixEarlySeparated = 10,
    RefinementTypeOrientationMatrixStdevOnly = 12,
    RefinementTypeRefineLAxis = 13,
    RefinementTypeOrientationMatrixEarlyWeighted = 14,
    RefinementTypeOrientationMatrixReverse = 15,
    
} RefinementType;

class IOMRefiner : public boost::enable_shared_from_this<IOMRefiner>
{
private:
	ImageWeakPtr image;
	vector<MillerPtr> millers;
	vector<MillerPtr> nearbyMillers;
    vector<MillerPtr> roughMillers;
	vector<Spot *>spots;
    MatrixPtr matrix;
    MatrixPtr lastRotatedMatrix;
    std::vector<Match> indexingMatches;
    MtzPtr lastMtz;

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
	double testDistance;
	double testWavelength;
	double testSpotSize;
	double initialStep;
    double getTotalStandardDeviation();
	int expectedSpots;
	int searchSize;
	static double intensityThreshold;
    static bool absoluteIntensity;
    static bool lowIntensityPenalty;
	double maxResolution;
	MtzManager *reference;
	void calculateNearbyMillers(bool rough);
	RefinementType refinement;
	double lastTotal;
	double lastStdev;
	double testBandwidth;
    double lastScore;
	vector<vector<double> > solutions;
    bool recalculateMillerPositions;
    
    double orientationTolerance;
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);

public:
	IOMRefiner(ImagePtr newImage = ImagePtr(), MatrixPtr matrix = MatrixPtr());
    void setComplexMatrix();
    virtual ~IOMRefiner();

    void lockUnitCellDimensions();
    void calculateOnce();
	void checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox = false, bool perfectCalculation = true);
	MtzPtr newMtz(int i, bool silent = false);
	void getWavelengthHistogram(vector<double> &wavelengths,
			vector<int> &frequencies, LogLevel level = LogLevelDetailed, int whichAxis = 0);
	double score(int whichAxis = 0, bool silent = false);

	static bool millerReachesThreshold(MillerPtr miller);
	void findSpots();
	static void duplicateSpots(vector<ImagePtr>images);
	void writeDatFromSpots(std::string filename);
	
	void matchMatrixToSpots();
	void matchMatrixToSpots(RefinementType refinement);
	double minimizeParameter(double *meanStep, double *param, int whichAxis = 0);
	void minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double *param1, double *param2);

	void refineDetectorAndWavelength(MtzManager *reference = NULL);
	void refineOrientationMatrix();
	void refineOrientationMatrix(RefinementType refinementType);
	
    void showHistogram(bool silent);
    bool millerWithinBandwidth(MillerPtr miller);
	int getTotalReflections();
	int getTotalReflections(double threshold);
    int getTotalReflectionsWithinBandwidth();
    
    double getRot(int rotNum);
	void dropMillers();

    bool isGoodSolution();
    bool isBasicGoodSolution();

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
    
	ImagePtr getImage()
	{
		return image.lock();
	}

	void setImage(ImagePtr image)
	{
		this->image = image;
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
    
    MtzPtr getLastMtz()
    {
        return lastMtz;
    }
};

#endif /* IOMRefiner_H_ */
