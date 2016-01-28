#ifndef mtz_manager
#define mtz_manager

#include <string>
#include <iostream>
#include "cmtzlib.h"
#include "csymlib.h"
#include <vector>
#include "Matrix.h"

#include "Holder.h"
#include "Miller.h"
//#include <tr1/memory>
#include "definitions.h"
#include "parameters.h"
#include <tuple>
#include "Logger.h"

using namespace CMtz;
using namespace CSym;

typedef enum
{
	ScoreTypeMinimizeRSplit = 0,
	ScoreTypeCorrelation = 1,
	ScoreTypePartialityCorrelation = 2,
	ScoreTypeMinimizeRSplitLog = 3,
	ScoreTypeCorrelationLog = 4,
	ScoreTypeAngle = 5,
	ScoreTypePartialityLeastSquares = 6,
	ScoreTypePartialityGradient = 7,
    ScoreTypeSymmetry = 8,
    ScoreTypeStandardDeviation = 9,
    ScoreTypeMinimizeRMeas = 10,
    ScoreTypeMaximiseArea = 11,
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;

class MtzManager
{

protected:
    std::string filename;
	CCP4SPG *low_group;
	static bool reflection_comparison(Reflection *i, Reflection *j);

	double extreme_index(MTZ *mtz, int max);
	void hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k, int *l,
			int *multiplier, int *offset);
	int index_for_reflection(int h, int k, int l, bool inverted);
	void findMultiplier(MTZ *mtz, int *multiplier, int *offset);


	// minimisation stuff

	vector<Reflection *> reflections;
	vector<Reflection *> refReflections;
	vector<Reflection *> matchReflections;
    MtzManager *previousReference;
    int previousAmbiguity;
    bool allowTrust;
    bool setInitialValues;
    double penaltyWeight;
    double penaltyResolution;
    
    RotationMode rotationMode;
	double hRot;
	double kRot;
    double aRot;
    double bRot;
    double cRot;
	double mosaicity;
	double spotSize;
	double wavelength;
	double bandwidth;
	double exponent;
	double refCorrelation;
    double refPartCorrel;
	double scale;
	double cellDim[3];
	double cellAngles[3];
    int activeAmbiguity;
    
	bool usingFixedWavelength;

	bool optimisingWavelength;
	bool optimisingBandwidth;
	bool optimisingMosaicity;
	bool optimisingRlpSize;
	bool optimisingOrientation;
	bool optimisingExponent;

	double stepSizeWavelength;
	double stepSizeBandwidth;
	double stepSizeMosaicity;
	double stepSizeRlpSize;
	double stepSizeOrientation;
    double stepSizeOrientABC;
	double stepSizeExponent;

	double toleranceWavelength;
	double toleranceBandwidth;
	double toleranceMosaicity;
	double toleranceRlpSize;
	double toleranceOrientation;
	double toleranceExponent;

	ScoreType defaultScoreType;
	bool usePartialityFunction;
	double maxResolutionAll;
	double maxResolutionRlpSize;

    vector<double> trialUnitCell;
    double superGaussianScale;
    double lastExponent;
    double *params;
    double detectorDistance;
    
	bool finalised;
	bool inverse;
	bool flipped;
    int failedCount;
	static double lowRes;
	static double highRes;
	bool freePass;
	bool rejected;

	TrustLevel trust;
	ScoreType scoreType;

    MatrixPtr baseMatrix;
	static MtzManager *referenceManager;
	MtzManager *lastReference;

    static double unitCellScore(void *object);
    double wavelengthStandardDeviation();
    std::ostringstream logged;
public:
    vector<double> superGaussianTable;
    double bFactor;
    double externalScale;
    int removeStrongSpots(std::vector<SpotPtr> *spots);
    int rejectOverlaps();

	MtzManager(void);
	virtual ~MtzManager(void);
	MtzPtr copy();
	void loadParametersMap();

    void addReflections(vector<Reflection *>reflections);
	void clearReflections();
	void addReflection(Reflection *reflection);
	void removeReflection(int i);
	void excludeFromLogCorrelation();
	void excludePartialityOutliers();

    MatrixPtr matrix;
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
	void setFilename(std::string name);
    std::string getFilename(void);
	void description(void);
    void incrementFailedCount();
    void resetFailedCount();
    void resetDefaultParameters();
    void chooseAppropriateTarget();
    
    void setDefaultMatrix();
	void setMatrix(double *components);
	void setMatrix(MatrixPtr newMat);
	void insertionSortReflections(void);
	void sortLastReflection(void);
    void applyUnrefinedPartiality();
    void incrementActiveAmbiguity();
    void cutToResolutionWithSigma(double acceptableSigma);
    double maxResolution();

    std::string filenameRoot();
	void setSpaceGroup(int spgnum);
	virtual void loadReflections(PartialityModel model, bool special = false);
	void loadReflections(int partiality);
	static void setReference(MtzManager *reference);
	int findReflectionWithId(long unsigned int refl_id, Reflection **reflection, bool insertionPoint = false);
	void findCommonReflections(MtzManager *other,
			vector<Reflection *> &reflectionVector1, vector<Reflection *> &reflectionVector2,
			int *num = NULL, bool force = false);
	double meanCorrectedIntensity(Reflection *reflection);
	double gradientAgainstManager(MtzManager *otherManager, bool withCutoff = true, double lowRes = 0, double highRes = 0);
	double minimizeGradient(MtzManager *otherManager, bool leastSquares);
    void bFactorAndScale(double *scale, double *bFactor, double exponent = 1, vector<std::pair<double, double> > *dataPoints = NULL);
	double minimizeRFactor(MtzManager *otherManager);
    void applyBFactor(double bFactor);
	void applyScaleFactor(double scaleFactor, double lowRes = 0, double highRes = 0, bool absolute = false);
	void applyScaleFactorsForBins(int binCount = 20);
	void clearScaleFactor();
	void makeScalesPermanent();
	double averageIntensity(void);
	void setSigmaToUnity();
	void setPartialityToUnity();
	double partialityRatio(Reflection *imgReflection, Reflection *refReflection);
	void getRefReflections(vector<Reflection *> *refPointer,
			vector<Reflection *> *matchPointer);
	void reallowPartialityOutliers();

	void setUnitCell(double a, double b, double c, double alpha, double beta,
			double gamma);
	void setUnitCell(vector<double> unitCell);
	void getUnitCell(double *a, double *b, double *c, double *alpha, double *beta,
			double *gamma);
	void copySymmetryInformationFromManager(MtzPtr toCopy);
	void applyPolarisation(void);

	void writeToFile(std::string newFilename, bool announce = false, bool shifts = false, bool includeAmbiguity = false);
	void writeToDat();
    void sendLog(LogLevel priority = LogLevelNormal);

	double correlationWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double rSplitWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double statisticsWithManager(MtzManager *otherManager, StatisticsFunction *function,
			bool printHits, bool silent, double lowRes, double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double statisticsWithManager(MtzManager *otherManager,
			StatisticsFunction *function, RFactorFunction *rFactorFunction,
			RFactorType rFactor, bool printHits, bool silent, double lowRes,
			double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double rFactorWithManager(RFactorType rFactor, bool printHits = false, bool silent = true, double lowRes = 0,
			double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL,
			bool shouldLog = false);

    int symmetryRelatedReflectionCount();
    
// minimisation stuff

	static MtzManager *currentManager;

    int refinedParameterCount();
    
	static int paramMult, spotMult;
	void getParams(double *parameters[], int paramCount = PARAM_NUM);
    void getParamPointers(double ***parameters, int paramCount = PARAM_NUM);
    void getSteps(double *parameters[], int paramCount = PARAM_NUM);
    void setParams(double parameters[], int paramCount = PARAM_NUM);

    double medianWavelength(double lowRes, double highRes);
	double bestWavelength(double lowRes = 0.0, double highRes = 0, bool usingReference = false);
	double weightedBestWavelength(double lowRes, double highRes);
	int accepted(void);
	static double exclusionScoreWrapper(void *object, double lowRes = 0,
			double highRes = 0);
    static double bFactorScoreWrapper(void *object);
    static double scoreNelderMead(void *object);
	double exclusionScore(double lowRes, double highRes, ScoreType scoreType);
	double leastSquaresPartiality(double low, double high, ScoreType typeOfScore = ScoreTypePartialityCorrelation);
	double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high, bool withCutoff = false, bool set = false);
    double belowPartialityPenalty(double low, double high);
	std::string describeScoreType();
    
    void refreshCurrentPartialities();
	void refreshPartialities(double hRot, double kRot, double aRot, double bRot, double cRot, double mosaicity,
                             double spotSize, double wavelength, double bandwidth, double exponent,
                             double a, double b, double c);
	void refreshPartialities(double parameters[]);
	double minimizeParameter(double *meanStep, double **params, int paramNum,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes);
	double partialHits(vector<double> *partials, vector<double> *percentages);
	double partialityGradient();
    double maximisePartialityArea(double low, double high);

// more grid search

    void setParamLine(std::string line);
    std::string getParamLine();
    void findSteps(int param1, int param2, std::string csvName);
	void gridSearch(bool silent = false);
	double leastSquaresPartiality(ScoreType typeOfScore);
    double minimize();
	double minimize(double (*score)(void *object, double lowRes, double highRes),
			void *object);
	double minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double **params, int paramNum1, int paramNum2,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes, double low);

    void denormaliseMillers();
    void makeSuperGaussianLookupTable(double exponent);
    
    int ambiguityCount();
    
    void flipToActiveAmbiguity();
    void resetFlip();
    void setAdditionalWeight(double weight);
    
    double getDetectorDistance()
    {
        return detectorDistance;
    }
    
    void setDetectorDistance(double distance)
    {
        detectorDistance = distance;
    }
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity);

    virtual int reflectionCount()
    {
        return (int)reflections.size();
    }
    
    virtual Reflection *reflection(int i)
    {
        return reflections[i];
    }
    
    void setDefaultScoreType(ScoreType scoreType)
    {
        defaultScoreType = scoreType;
    }
    
    vector<double> getTrialUnitCell()
    {
        return trialUnitCell;
    }
    
    double getSuperGaussianScale() const
    {
        return superGaussianScale;
    }
    
    double getLastExponent() const
    {
        return lastExponent;
    }
    
	double getBandwidth() const
	{
		return bandwidth;
	}

	void setBandwidth(double bandwidth)
	{
		this->bandwidth = bandwidth;
	}

	double getMosaicity() const
	{
		return mosaicity;
	}

	void setMosaicity(double mosaicity)
	{
		this->mosaicity = mosaicity;
	}

	double getSpotSize() const
	{
		return spotSize;
	}

	void setSpotSize(double spotSize)
	{
		this->spotSize = spotSize;
	}

	double getWavelength() const
	{
		return wavelength;
	}

	void setWavelength(double wavelength)
	{
		this->wavelength = wavelength;
	}

	double getHRot() const
	{
		return hRot;
	}

	void setHRot(double hRot)
	{
		this->hRot = hRot;
	}

	double getKRot() const
	{
		return kRot;
	}

	void setKRot(double kRot)
	{
		this->kRot = kRot;
	}
    
    double getARot() const
    {
        return aRot;
    }
    
    double getBRot() const
    {
        return bRot;
    }
    
    double getCRot() const
    {
        return cRot;
    }

	double getRefCorrelation() const
	{
		return refCorrelation;
	}

	void setRefCorrelation(double refCorrelation)
	{
		this->refCorrelation = refCorrelation;
	}
    
    double getRefPartCorrel() const
    {
        return refPartCorrel;
    }
    
    void setRefPartCorrel(double refPartCorrel)
    {
        this->refPartCorrel = refPartCorrel;
    }

	double getExponent() const
	{
		return exponent;
	}

	void setExponent(double exponent)
	{
		this->exponent = exponent;
	}

	ScoreType getScoreType()
	{
		return scoreType;
	}

	void setScoreType(ScoreType newScoreType)
	{
		scoreType = newScoreType;
	}

	static double getHighRes()
	{
		return highRes;
	}

	static void setHighRes(double _highRes)
	{
		highRes = _highRes;
	}

	static double getLowRes()
	{
		return lowRes;
	}

	static void setLowRes(double _lowRes)
	{
		lowRes = _lowRes;
	}

	CCP4SPG*& getLowGroup()
	{
		return low_group;
	}

	void setLowGroup(CCP4SPG*& lowGroup)
	{
		low_group = lowGroup;
	}

	TrustLevel getTrust() const
	{
		return trust;
	}

	void setTrust(TrustLevel trust)
	{
		this->trust = trust;
	}

	bool isFinalised() const
	{
		return finalised;
	}

	void setFinalised(bool finalised)
	{
		this->finalised = finalised;
	}

	bool isFreePass() const
	{
		return freePass;
	}

	void setFreePass(bool freePass)
	{
		this->freePass = freePass;
	}

	bool isFlipped() const
	{
		return flipped;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	static MtzManager*& getReferenceManager()
	{
        return referenceManager;
	}

	bool isUsingFixedWavelength() const
	{
		return usingFixedWavelength;
	}

	void setUsingFixedWavelength(bool usingFixedWavelength)
	{
		this->usingFixedWavelength = usingFixedWavelength;
	}

	bool isOptimisingBandwidth() const
	{
		return optimisingBandwidth;
	}

	void setOptimisingBandwidth(bool optimisingBandwidth)
	{
		this->optimisingBandwidth = optimisingBandwidth;
	}

	bool isOptimisingExponent() const
	{
		return optimisingExponent;
	}

	void setOptimisingExponent(bool optimisingExponent)
	{
		this->optimisingExponent = optimisingExponent;
	}

	bool isOptimisingMosaicity() const
	{
		return optimisingMosaicity;
	}

	void setOptimisingMosaicity(bool optimisingMosaicity)
	{
		this->optimisingMosaicity = optimisingMosaicity;
	}

	bool isOptimisingOrientation() const
	{
		return optimisingOrientation;
	}

	void setOptimisingOrientation(bool optimisingOrientation)
	{
		this->optimisingOrientation = optimisingOrientation;
	}

	bool isOptimisingRlpSize() const
	{
		return optimisingRlpSize;
	}

	void setOptimisingRlpSize(bool optimisingRlpSize)
	{
		this->optimisingRlpSize = optimisingRlpSize;
	}

	bool isOptimisingWavelength() const
	{
		return optimisingWavelength;
	}

	void setOptimisingWavelength(bool optimisingWavelength)
	{
		this->optimisingWavelength = optimisingWavelength;
	}

	double getStepSizeBandwidth() const
	{
		return stepSizeBandwidth;
	}

	void setStepSizeBandwidth(double stepSizeBandwidth)
	{
		this->stepSizeBandwidth = stepSizeBandwidth;
	}

	double getStepSizeExponent() const
	{
		return stepSizeExponent;
	}

	void setStepSizeExponent(double stepSizeExponent)
	{
		this->stepSizeExponent = stepSizeExponent;
	}

	double getStepSizeMosaicity() const
	{
		return stepSizeMosaicity;
	}

	void setStepSizeMosaicity(double stepSizeMosaicity)
	{
		this->stepSizeMosaicity = stepSizeMosaicity;
	}

	double getStepSizeOrientation() const
	{
		return stepSizeOrientation;
	}

	void setStepSizeOrientation(double stepSizeOrientation)
	{
		this->stepSizeOrientation = stepSizeOrientation;
	}

	double getStepSizeRlpSize() const
	{
		return stepSizeRlpSize;
	}

	void setStepSizeRlpSize(double stepSizeRlpSize)
	{
		this->stepSizeRlpSize = stepSizeRlpSize;
	}

	double getStepSizeWavelength() const
	{
		return stepSizeWavelength;
	}

	void setStepSizeWavelength(double stepSizeWavelength)
	{
		this->stepSizeWavelength = stepSizeWavelength;
	}

	bool isRejected() const
	{
		return rejected;
	}

	void setRejected(bool rejected)
	{
		this->rejected = rejected;
	}

	double getScale() const
	{
		return scale;
	}

	void setScale(double scale)
	{
		this->scale = scale;
	}
};

#endif

