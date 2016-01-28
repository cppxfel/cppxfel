/*
 * Image.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IMAGE_H_
#define IMAGE_H_

#include "Matrix.h"
#include "IOMRefiner.h"
#include "parameters.h"
#include "Logger.h"
#include "csymlib.h"
#include "SpotVector.h"
#include "LoggableObject.h"

class IOMRefiner;
class ImageCluster;

typedef enum
{
    IndexingSolutionTrialSuccess,
    IndexingSolutionTrialFailure,
    IndexingSolutionTrialDuplicate,
    IndexingSolutionBranchFailure,
} IndexingSolutionStatus;

class Image : LoggableObject, public boost::enable_shared_from_this<Image>
{
private:
    int pixelCountCutoff;
	std::string filename;
	vector<int> data;
    vector<unsigned char> overlapMask;
	void loadImage();
    void findSpots();
    vector<IOMRefinerPtr> indexers;
    bool shouldMaskValue;
    bool maskedValue;
    bool fitBackgroundAsPlane;
    std::string spotsFile;
    IndexingSolutionStatus extendIndexingSolution(IndexingSolutionPtr solutionPtr, std::vector<SpotVectorPtr> existingVectors, int *failures = NULL, int added = 0);
    
	/* Shoebox must be n by n where n is an odd number */
	int shoebox[7][7];

	int xDim;
	int yDim;

	int beamX;
	int beamY;
	double mmPerPixel;
    bool noCircles;
    double detectorGain;

	double detectorDistance; // mm
	double wavelength;
	bool pinPoint;

    int indexingFailureCount;
    int minimumSolutionNetworkCount;
    
    std::vector<SpotPtr> spots;
    std::vector<SpotVectorPtr> spotVectors;
    double commonCircleThreshold;
    bool _hasSeeded;
    std::map<ImageCluster *, bool> unexpectedMatches;
    
	vector<vector<int> > masks;
	vector<vector<int> > spotCovers;

	int shoeboxLength();
	Mask flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    double integrateFitBackgroundPlane(int x, int y, ShoeboxPtr shoebox, double *error);
    double integrateSimpleSummation(int x, int y, ShoeboxPtr shoebox, double *error);
	double integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, double *error);
	bool checkShoebox(ShoeboxPtr shoebox, int x, int y);
    double weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    bool checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast = false);
public:
    void incrementOverlapMask(int x, int y, ShoeboxPtr shoebox);
    void incrementOverlapMask(int x, int y);
    void processSpotList();
    void writeSpotsList(std::string spotFile = "");
    
    unsigned char overlapAt(int x, int y);
    unsigned char maximumOverlapMask(int x, int y, ShoeboxPtr shoebox);
	Image(std::string filename = "", double wavelength = 0,
			double distance = 0);
	void focusOnSpot(int *x, int *y, int tolerance1, int tolerance2);
	void focusOnAverageMax(int *x, int *y, int tolerance1, int tolerance2 = 1, bool even = false);
    void focusOnMaximum(int *x, int *y, int tolerance = 0, double shiftX = 0, double shiftY = 0);
	void dropImage();
	virtual ~Image();
	void setUpIOMRefiner(MatrixPtr matrix);
    void setUpIOMRefiner(MatrixPtr unitcell, MatrixPtr rotation);
	std::string filenameRoot();
	void printBox(int x, int y, int tolerance);
	void addMask(int startX, int startY, int endX, int endY);
	void addSpotCover(int startX, int startY, int endX, int endY);
	bool coveredBySpot(int x, int y);
	static void applyMaskToImages(vector<ImagePtr> images, int startX,
			int startY, int endX, int endY);
    void refineDistances();
    IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr);
    std::vector<double> anglesBetweenVectorDistances(double distance1, double distance2, double tolerance);
    
    void rotatedSpotPositions(MatrixPtr rotationMatrix, std::vector<vec> *spotPositions, std::vector<std::string> *spotElements);

	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
	}
    
    std::string getBasename()
    {
        int fullStopIndex = (int)filename.rfind(".");
        if (fullStopIndex == std::string::npos)
            return filename;
        
        std::string basename = filename.substr(0, fullStopIndex);
        
        return basename;
    }
    
    std::string getSpotsFile()
    {
        return spotsFile;
    }
    
    void setSpotsFile(std::string newFile)
    {
        spotsFile = newFile;
    }

	int valueAt(int x, int y);
	bool accepted(int x, int y);
	double intensityAt(int x, int y, ShoeboxPtr shoebox, double *error, int tolerance = 0);

	void index();
	void refineIndexing(MtzManager *reference);
	void refineOrientations();
	vector<MtzPtr> currentMtzs();
    std::vector<MtzPtr> getLastMtzs();
	bool isLoaded();
    
    void setSpaceGroup(CSym::CCP4SPG *spg);
    void setMaxResolution(double res);
    void setSearchSize(int searchSize);
    void setIntensityThreshold(double threshold);
    void setUnitCell(vector<double> dims);
    void setInitialStep(double step);
    void setTestSpotSize(double spotSize);
    void setTestBandwidth(double bandwidth);
    void setOrientationTolerance(double newTolerance);
    
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
    void findIndexingSolutions();
    void compileDistancesFromSpots(double maxReciprocalDistance = 0, double tooCloseDistance = 0, bool filter = false);
    void filterSpotVectors();
    int throwAwayIntegratedSpots(std::vector<MtzPtr> mtzs);
    void updateAllSpots();
    
    void removeRefiner(int j)
    {
        indexers.erase(indexers.begin() + j);
    }
    
    int spotVectorCount()
    {
        return (int)spotVectors.size();
    }
    
    SpotPtr spot(int i)
    {
        return spots[i];
    }
    
    int spotCount()
    {
        return (int)spots.size();
    }
    SpotVectorPtr spotVector(int i)
    {
        return spotVectors[i];
    }
    
    bool hasSeeded()
    {
        return _hasSeeded;
    }
    
    void setSeeded(bool seed = true)
    {
        _hasSeeded = seed;
    }
    
    int IOMRefinerCount()
    {
        return (int)indexers.size();
    }
    
    IOMRefinerPtr getIOMRefiner(int i)
    {
        return indexers[i];
    }
    
    void addIOMRefiner(IOMRefinerPtr newRefiner)
    {
        indexers.push_back(newRefiner);
    }
    
    void clearIOMRefiners()
    {
        indexers.clear();
        std::vector<IOMRefinerPtr>().swap(indexers);
    }
    
	int getXDim() const
	{
		return xDim;
	}

	void setXDim(int dim)
	{
		xDim = dim;
	}

	int getYDim() const
	{
		return yDim;
	}

	void setYDim(int dim)
	{
		yDim = dim;
	}

	double getDetectorDistance() const
	{
		return detectorDistance;
	}

	void setDetectorDistance(double detectorDistance)
	{
		this->detectorDistance = detectorDistance;
        
        
	}

	double getWavelength() const
	{
		return wavelength;
	}

	void setWavelength(double wavelength)
	{
		this->wavelength = wavelength;
	}

	int getBeamX() const
	{
		return beamX;
	}

	void setBeamX(int beamX)
	{
		this->beamX = beamX;
	}

	int getBeamY() const
	{
		return beamY;
	}

	void setBeamY(int beamY)
	{
		this->beamY = beamY;
	}

	double getMmPerPixel() const
	{
		return mmPerPixel;
	}

	void setMmPerPixel(double mmPerPixel)
	{
		this->mmPerPixel = mmPerPixel;
	}

	bool isPinPoint() const
	{
		return pinPoint;
	}

	void setPinPoint(bool pinPoint)
	{
		this->pinPoint = pinPoint;
	}
   
    void setImageData(vector<int> newData);
    
    double getDetectorGain()
    {
        return detectorGain;
    }
    
    void setDetectorGain(double newGain)
    {
        detectorGain = newGain;
    }
};

#endif /* IMAGE_H_ */
