/*
 * GraphDrawer.h
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

#ifndef GRAPHDRAWER_H_
#define GRAPHDRAWER_H_

#include "parameters.h"
#include "StatisticsManager.h"
#include "MtzManager.h"
#include <boost/variant.hpp>
#include <map>

typedef std::map<std::string, boost::variant<double, std::string> > GraphMap;

class GraphDrawer
{
private:
	MtzManager *mtz;
//	vector<MtzManager *>mtzs;

public:
	GraphDrawer(MtzManager *mtz);
	virtual ~GraphDrawer();

	static std::string generateFilename(std::string stem);
	static std::string generateFilename(std::string stem, std::string ext);


	std::string plot(std::string filename, GraphMap properties,
			vector<double> x, vector<double> y,
			vector<double> x2, vector<double> y2);

	std::string plot(std::string filename, GraphMap properties,
			vector<vector<double> > x, vector<vector<double> > y,
			vector<double> x2, vector<double> y2);

	std::string plot(std::string filename,
			GraphMap properties,
			vector<double> x, vector<double> y);

	std::string plot(std::string filename, GraphMap properties,
			vector<vector<double> > xs, vector<vector<double> > ys);

	void resolutionStatsPlot(vector<MtzManager *>& managers, std::string filename = "resolution_stats",
			GraphMap properties = GraphMap(), bool intensityBins = false, bool image = false);

    void plotPolarisation(vector<MtzPtr> mtzs);
	void correlationPlot(std::string filename, double xMax = 0, double yMax = 0);
	void partialityPlot(std::string filename, GraphMap properties = GraphMap(), double maxRes = 1.6);
	void bFactorPlot(vector<MtzManager *>& managers,
                     std::string filename = "all_gradients", GraphMap properties = GraphMap());
    
    void plotReflectionFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l);
    void plotOrientationStats(vector<MtzPtr> mtzs);
    void plotPartialityStats(int h = 0, int k = 0, int l = 0);
    MtzManager*& getMtz()
	{
		return mtz;
	}

	void setMtz(MtzManager*& mtz)
	{
		this->mtz = mtz;
	}
};

#endif /* GRAPHDRAWER_H_ */
