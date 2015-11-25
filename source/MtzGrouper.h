/*
 * MtzGrouper.h
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#ifndef MTZGROUPER_H_
#define MTZGROUPER_H_

#include <vector>

#include "MtzManager.h"
#include "parameters.h"
#include "Logger.h"

typedef enum
{
	ScalingTypeAverage = 0,
	ScalingTypeReference = 1,
	ScalingTypeReferenceLeastSquares = 2,
	ScalingTypeMinimizeRMerge = 3,
	ScalingTypeBFactor = 4,
	ScalingTypeResolutionShells = 5
} ScalingType;

class MtzGrouper
{
private:
    bool usingNewRefinement;
    std::ostringstream logged;
	double correlationThreshold;
	ScalingType scalingType;
	bool excludeWorst;
	WeightType weighting;
	double acceptableResolution;
	bool cutResolution;
	double expectedResolution;
    bool isMtzAccepted(MtzPtr mtz);
    bool exclusionByCCHalf;

    static void checkCCHalf(vector<MtzPtr> *managers, int offset, int *total);
	void merge(MtzManager **mergeMtz, MtzManager **unmergedMtz, bool firstHalf,
			bool all, std::string *unmergedName = NULL);
public:
	MtzGrouper();
	virtual ~MtzGrouper();

	vector<MtzPtr> mtzManagers;

    void merge(MtzManager **mergeMtz, MtzManager **unmergedMtz = NULL, int cycle = -1, bool anom = false);

	void mergeAnomalous(MtzManager **mergeMtz, MtzManager **unmergedMtz,
			bool firstHalf, bool all, std::string filename = "anomalous");
	void differenceBetweenMtzs(MtzManager **mergeMtz, MtzManager **positive, MtzManager **negative);

	static void mergeWrapper(void *object, MtzManager **mergeMtz,
			MtzManager **unmergedMtz, bool firstHalf, bool all, bool anom);
	int groupMillers(MtzManager **mergeMtz, MtzManager **unmergedMtz, int start,
			int end);
	void mergeMillers(MtzManager **mergeMtz, bool reject, int mtzCount);
	void unflipMtzs();

	int groupMillersWithAnomalous(MtzManager **positive, MtzManager **negative, int start,
			int end);
	void writeAnomalousMtz(MtzManager **positive, MtzManager **negative,
			std::string filename);
    
    void sendLog(LogLevel priority = LogLevelNormal);

	double getCorrelationThreshold() const
	{
		return correlationThreshold;
	}

	void setCorrelationThreshold(double correlationThreshold)
	{
		this->correlationThreshold = correlationThreshold;
	}

	const vector<MtzPtr>& getMtzManagers() const
	{
		return mtzManagers;
	}

	void setMtzManagers(const vector<MtzPtr>& mtzManagers);

	bool isExcludeWorst() const
	{
		return excludeWorst;
	}

	void setExcludeWorst(bool excludeWorst)
	{
		this->excludeWorst = excludeWorst;
	}

	WeightType getWeighting() const
	{
		return weighting;
	}

	void setWeighting(WeightType weighting)
	{
		this->weighting = weighting;
	}

	ScalingType getScalingType() const
	{
		return scalingType;
	}

	void setScalingType(ScalingType scalingType)
	{
		this->scalingType = scalingType;
	}

	double getAcceptableResolution() const
	{
		return acceptableResolution;
	}

	bool isCutResolution() const
	{
		return cutResolution;
	}

	void setCutResolution(bool cutResolution)
	{
		this->cutResolution = cutResolution;
	}

	double getExpectedResolution() const
	{
		return expectedResolution;
	}

	void setExpectedResolution(double expectedResolution)
	{
		this->expectedResolution = expectedResolution;
	}
};

#endif /* MTZGROUPER_H_ */
