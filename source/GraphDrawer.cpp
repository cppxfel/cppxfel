/*
 * GraphDrawer.cpp
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

#include "GraphDrawer.h"
#include "misc.h"
#include <vector>
#include<iostream>
#include <algorithm>
#include <cmath>
#include <boost/variant.hpp>
#include <map>
#include <sstream>
#include <tuple>
#include "Panel.h"
#include <fstream>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/space_group.h>



GraphDrawer::GraphDrawer(MtzManager *mtz)
{
	// TODO Auto-generated constructor stub
	this->mtz = mtz;
}

bool sortByWavelength(Partial x, Partial y)
{
	return (x.wavelength < y.wavelength);
}

bool sortByPercentage(Partial x, Partial y)
{
	return (x.percentage < y.percentage);
}



std::string GraphDrawer::generateFilename(std::string stem)
{
	return generateFilename(stem, "png");
}

std::string GraphDrawer::generateFilename(std::string stem, std::string ext)
{
	time_t cputime;
	time(&cputime);

	std::ostringstream timeString;
	timeString << stem;

	timeString << "." << ext;

	return timeString.str();
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<vector<double> > xs,
		vector<vector<double> > ys)
{
	vector<double> x2, y2;

	return plot(filename, properties, xs, ys, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<double> x, vector<double> y)
{
	vector<double> x2, y2;

	return plot(filename, properties, x, y, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<double> x, vector<double> y, vector<double> x2,
		vector<double> y2)
{
	vector<vector<double> > xs;
	vector<vector<double> > ys;

    xs.push_back(x); xs.push_back(x2);
    ys.push_back(y); ys.push_back(y2);

	return plot(filename, properties, xs, ys, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<vector<double> > xs,
		vector<vector<double> > ys, vector<double> x2,
		vector<double> y2)
{
#ifdef MAC
/*	std::string plotType = "line";
	std::string plotType2 = "line";
	std::string title = "Title";
	std::string xTitle = "x axis";
	std::string yTitle = "y axis";

	if (properties.count("title") > 0)
		title = boost::get<std::string>(properties["title"]);

	if (properties.count("xTitle") > 0)
		xTitle = boost::get<std::string>(properties["xTitle"]);

	if (properties.count("yTitle") > 0)
		yTitle = boost::get<std::string>(properties["yTitle"]);

	if (properties.count("plotType") > 0)
		plotType = boost::get<std::string>(properties["plotType"]);

	if (properties.count("plotType2") > 0)
		plotType2 = boost::get<std::string>(properties["plotType2"]);

	plsdev("png");

	char geometry[] = "-geometry";
	char dims[] = "1600x1250";

	plsetopt(geometry, dims);

	std::string fullName = generateFilename(filename);

	plsfnam(fullName.c_str());
	plscolbg(255, 255, 255);
	plscol0(1, 0, 0, 0);

	double minX = FLT_MAX;
	double maxX = -FLT_MAX;
	double minY = FLT_MAX;
	double maxY = -FLT_MAX;

	for (int j = 0; j < xs.size(); j++)
	{
		for (int i = 0; i < xs[j].size(); i++)
		{
			if (xs[j][i] < minX)
				minX = xs[j][i];
			if (xs[j][i] > maxX)
				maxX = xs[j][i];
			if (ys[j][i] < minY)
				minY = ys[j][i];
			if (ys[j][i] > maxY)
				maxY = ys[j][i];
		}
	}

	if (properties.count("xMin") > 0)
		minX = boost::get<double>(properties["xMin"]);

	if (properties.count("xMax") > 0)
		maxX = boost::get<double>(properties["xMax"]);

	if (properties.count("yMin") > 0)
		minY = boost::get<double>(properties["yMin"]);

	if (properties.count("yMax") > 0)
		maxY = boost::get<double>(properties["yMax"]);

	// Initialize plplot
	plinit();

	// Create a labelled box to hold the plot.
	plenv(minX, maxX, minY, maxY, 0, 0);
	pllab(xTitle.c_str(), yTitle.c_str(), title.c_str());

	for (int i = 0; i < xs.size(); i++)
	{
		vector<double> x = xs[i];
		vector<double> y = ys[i];

		std::ostringstream plotTypeStream;
		plotTypeStream << "plotType_" << i;
		std::string plotTypeName = plotTypeStream.str();

		plssym(0, 0.5);

		std::ostringstream colourStream;
		colourStream << "colour_" << i;
		std::string colourString = colourStream.str();

		if (properties.count(colourString) > 0)
		{
			int colour = boost::get<double>(properties[colourString]);
			plcol0(colour);
		}

		std::string usedPlotType = plotType;
		if (properties.count(plotTypeName) > 0)
			usedPlotType = boost::get<std::string>(properties[plotTypeName]);

		if (usedPlotType == "fill")
		{
			double firstX = xs[i][0];
			double lastX = xs[i][xs[i].size() - 1];

			x.insert(x.begin(), firstX);
			y.insert(y.begin(), 0);
			x.push_back(lastX);
			y.push_back(0);
		}

		int nsize = (int)x.size();

		double *allX = &x[0];
		double *allY = &y[0];

		// Plot the data that was prepared above.

		if (usedPlotType == "point")
			plpoin(nsize, allX, allY, 4);

		if (usedPlotType == "line")
			plline(nsize, allX, allY);

		if (usedPlotType == "fill")
		{
            plscol0(3, 230, 230, 240);
			plcol0(3);
			plfill(nsize, allX, allY);
			plcol0(1);
			plline(nsize, allX, allY);
		}
	}

	if (x2.size() > 0)
	{
		double yMin2 = 0;
		if (properties.count("yMin2") > 0)
			yMin2 = boost::get<double>(properties["yMin2"]);

		double yMax2 = 1;
		if (properties.count("yMax2") > 0)
			yMax2 = boost::get<double>(properties["yMax2"]);

		plscol0(2, 200, 0, 0);
		plcol0(2);
		plwind(minX, maxX, yMin2, yMax2);
		plbox("", 0, 0, "cmstv", 0, 0);
		plwidth(2);

		double *allX = &x2[0];
		double *allY = &y2[0];
		int nsize = (int)x2.size();

		if (plotType2 == "point")
			plpoin(nsize, allX, allY, 4);

		if (plotType2 == "line")
		{
			plline(nsize, allX, allY);
			plpoin(nsize, allX, allY, 4);
		}
	}

	// Close PLplot library
	plend();

    std::ostringstream logged;
	logged << "Plotted " << fullName << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelDetailed);
*/
#endif
    
    std::string csvName = generateFilename(filename, "csv");
    
    std::ofstream csv;
    csv.open(csvName);
    
    int biggest = 0;
    
    for (int i = 0; i < xs.size(); i++)
    {
        if (xs[i].size() > biggest)
            biggest = (int)xs[i].size();
    }
    
    for (int i = 0; i < biggest; i++)
    {
        for (int j = 0; j < xs.size(); j++)
        {
            if (i < xs[j].size())
                csv << xs[j][i] << "," << ys[j][i] << ",";
            else csv << ",,";
            
        }
        
        csv << std::endl;
    }
    
    csv.close();
    
    return csvName;
}

void GraphDrawer::correlationPlot(std::string filename, double xMax, double yMax)
{
	GraphMap properties;

	properties["title"] = "Correlation plot";
	properties["xTitle"] = "Intensity (image)";
	properties["yTitle"] = "Intensity (reference MTZ)";
	properties["plotType"] = "point";

	if (xMax > 0)
		properties["xMax"] = xMax;

	if (yMax > 0)
		properties["yMax"] = yMax;

	vector<double> x;
	vector<double> y;

	MtzManager *shot1 = this->mtz;
	MtzManager *shot2 = mtz->getReferenceManager();
//	double scale = mtz->gradientAgainstManager(*shot2);
//	mtz->applyScaleFactor(scale);

	vector<Reflection *> reflections1;
	vector<Reflection *> reflections2;
	int num = 0;

	shot1->findCommonReflections(shot2, reflections1, reflections2, &num);

	for (int i = 0; i < num; i++)
	{
		double int1 = reflections1[i]->meanIntensity();
		std::string filename = shot1->getFilename();
		double int2 = reflections2[i]->meanIntensity();

		if (int1 != int1 || int2 != int2)
			continue;

		x.push_back(int1);
		y.push_back(int2);
	}

	std::string fullName = this->plot(filename, properties, x, y);
	std::cout << "N: plot " << fullName << std::endl;
}

void GraphDrawer::resolutionStatsPlot(vector<MtzManager *>& managers,
		std::string filename, GraphMap properties, bool intensityBins, bool image)
{
	if (!intensityBins)
		properties["title"] = "Intensity agreement vs resolution";
	else
		properties["title"] = "Intensity agreement vs strength of signal";

	if (!intensityBins)
	{
		properties["xTitle"] = "Resolution (1 / d^2)";
	}
	else
	{
		if (image)
		{
			properties["xTitle"] = "Log(image intensity)";
		}
		else
		{
			properties["xTitle"] = "Reference intensity";
		}
	}

	properties["yTitle"] =
			"Percentage of merged data set (individual reflections)";
	properties["plotType"] = "point";

	properties["xMin"] = 0;
//	properties["xMax"] = 0.5;
	properties["yMin"] = 0;
	properties["yMax"] = 4;

	properties["colour_0"] = 10;
	properties["colour_1"] = 4;

	vector<double> bins;
	StatisticsManager::generateResolutionBins(50, 1.6, 18, &bins);

	MtzManager *referenceManager = MtzManager::getReferenceManager();
	double grad = 1000 / referenceManager->averageIntensity();
	referenceManager->applyScaleFactor(grad);

	vector<double> blueX, blueY, redX, redY;
    
    double percentSum = 0;
    double sum = 0;
    
    for (int i = 0; i < managers.size(); i++)
	{
		mtz = managers[i];

	//	mtz->excludeFromLogCorrelation();
        
        double correl = mtz->correlation(true);
        double invCorrel = correl;
        
        if (correl < 0.85)
            continue;
        
        if (MtzManager::getReferenceManager()->ambiguityCount() == 2)
        {
            mtz->setActiveAmbiguity(1);
            invCorrel = mtz->correlation(true);
            mtz->setActiveAmbiguity(0);
            
            if (invCorrel > correl)
                mtz->setActiveAmbiguity(1);
        }
        double newCorrel = mtz->correlation(true);
        mtz->setRefCorrelation(newCorrel);
        
        std::cout << mtz->getFilename() << "\t" << correl << "\t"
        << invCorrel << std::endl;

		double scale = mtz->gradientAgainstManager(referenceManager);
        std::cout << "Scale: " << scale << std::endl;
        mtz->applyScaleFactor(scale);
     
		vector<Reflection *> refReflections, imgReflections;

		mtz->findCommonReflections(referenceManager, imgReflections, refReflections,
		NULL);

		for (int j = 0; j < refReflections.size(); j++)
		{
			if (!imgReflections[j]->anyAccepted())
				continue;

			double imgIntensity = imgReflections[j]->meanIntensity();
			double refIntensity = refReflections[j]->meanIntensity();
      //      double meanPartiality = imgReflections[j]->meanPartiality();
            
			double percent = imgIntensity / refIntensity;

			vector<double> *chosenX = &redX;
			vector<double> *chosenY = &redY;

			if (log(imgIntensity) > 7)
			{
				chosenX = &blueX;
				chosenY = &blueY;
			}
            
			if (!intensityBins)
			{
				double resolution = refReflections[j]->getResolution();
				double res_squared = pow(resolution, 2);
				chosenX->push_back(res_squared);
			}
			else
			{
				double logIntensity = 0;
				if (image)
					logIntensity = log(imgIntensity);
				else
					logIntensity = log(refIntensity);

				chosenX->push_back(logIntensity);
			}
			chosenY->push_back(percent);

        }
    }
    
    double aveSum = percentSum / (double)sum;
    
    std::cout << "Average sum = " << aveSum << std::endl;
    
	if (redX.size() == 0)
		properties["colour_0"] = 1;

	vector<vector<double> > xs, ys;
	xs.push_back(blueX);
	xs.push_back(redX);
	ys.push_back(blueY);
	ys.push_back(redY);

	plot(filename, properties, xs, ys);
}

void GraphDrawer::bFactorPlot(vector<MtzManager *>& managers, std::string filename,
		GraphMap properties)
{
	vector<vector<std::string> > filenames;
	time_t cputime;
	time(&cputime);

	for (int i = 0; i < managers.size(); i++)
	{
		filenames.push_back(vector<std::string>());

		for (int j = 0; j < 3; j++)
		{
			std::ostringstream stream;
			stream << cputime << "_bfactor_" << i << "_" << j;
			filenames[i].push_back(stream.str());
		}

		cputime++;
	}

	std::ostringstream results;

	properties["title"] = "Scale and B factor plot";
	properties["xTitle"] = "Resolution (1 / d^2)";
	properties["yTitle"] = "Ratio of intensities (image over reference)";
	properties["plotType"] = "line";

	properties["xMin"] = 0;
//	properties["xMax"] = 0.3;
	properties["yMin"] = 0;
	properties["yMax"] = 4;

	vector<double> bins;
	StatisticsManager::generateResolutionBins(50, 2.1, 18, &bins);

	vector<vector<double> > xs;
	vector<vector<double> > ys;

	MtzManager *referenceManager = MtzManager::getReferenceManager();
	/*
	double grad = 1000 / referenceManager->averageIntensity();
	referenceManager->applyScaleFactor(grad);
*/
	for (int i = 0; i < managers.size(); i++)
	{
		mtz = managers[i];
		mtz->excludeFromLogCorrelation();

		double correl = managers[i]->correlation(true);
		if (correl < 0.85)
			continue;

		 double gradient = managers[i]->gradientAgainstManager(
		 MtzManager::getReferenceManager());

	//	 double bFactor = 0;

	//	 managers[i]->bFactorAndScale(&gradient, &bFactor);
		 //	gradient = 1000 / managers[i]->averageIntensity();

		 	managers[i]->applyScaleFactor(gradient);


		vector<double> bins;
		StatisticsManager::generateResolutionBins(50, 1.6, 10, &bins);

		vector<double> x, y;

		for (int shell = 0; shell < bins.size() - 3; shell++)
		{
			double low = bins[shell];
			double high = bins[shell + 3];

			vector<Reflection *> refReflections, imgReflections;

			mtz->findCommonReflections(referenceManager, imgReflections, refReflections,
			NULL);

			double weights = 0;
			double refMean = 0;
			double imgMean = 0;
			int count = 0;

			for (int i = 0; i < imgReflections.size(); i++)
			{
				if (!imgReflections[i]->anyAccepted())
					continue;

				if (imgReflections[i]->betweenResolutions(low, high))
				{
					weights += imgReflections[i]->meanPartiality();
					refMean += refReflections[i]->meanIntensity()
							* imgReflections[i]->meanPartiality();
					imgMean += imgReflections[i]->meanIntensity()
							* imgReflections[i]->meanPartiality();
					count++;
				}
			}

			refMean /= weights;
			imgMean /= weights;

			double ratio = refMean / imgMean;
			ratio = 1 / ratio;

			if (ratio != ratio)
				continue;

			x.push_back(1 / pow(low, 2));
			y.push_back(ratio);
		}

		xs.push_back(x);
		ys.push_back(y);
	}

	this->plot(filename, properties, xs, ys);

	xs.clear();
	ys.clear();
/*
	return;

	for (int i = 0; i < managers.size(); i++)
	{
		vector<vector<double> > xs;
		vector<vector<double> > ys;

		mtz = managers[i];
		mtz->excludeFromLogCorrelation();

		double preRSplit = 0;
		double preCorrel = 0;

		double gradient = managers[i]->gradientAgainstManager(
				*MtzManager::getReferenceManager(), false, bins[0], bins[5]);
		managers[i]->applyScaleFactor(gradient);

		this->mtz = managers[i];
		this->correlationPlot(filenames[i][0], 12000, 12000);

		preRSplit = mtz->rSplitWithManager(MtzManager::getReferenceManager());
		preCorrel = mtz->correlationWithManager(
				MtzManager::getReferenceManager());

		vector<double> x;
		vector<double> y;

		bool include = true;

		for (int shell = 0; shell < bins.size() - 1; shell++)
		{
			double minD = 0;
			double maxD = 0;
			StatisticsManager::convertResolutions(bins[shell], bins[shell + 1],
					&minD, &maxD);

			double d = StatisticsManager::midPointBetweenResolutions(bins[shell],
					bins[shell + 1]);

			double one_over_d_squared = pow(1 / d, 2);

			double gradient = managers[i]->gradientAgainstManager(
					*MtzManager::getReferenceManager(), false, bins[shell],
					bins[shell + 1]);

			if (gradient > 3)
				include = false;

			x.push_back(one_over_d_squared);
			y.push_back(1 / gradient);
		}

		if (include && false)
		{
			xs.push_back(x);
			ys.push_back(y);
		}

		properties["plotType_0"] = "point";
		properties["plotType_1"] = "line";

		double bFactor = 0;
		double k = 1;
		double exponent_exponent = 1.0;

		vector<pair<double, double> > dataPoints;

		mtz->bFactorAndScale(&k, &bFactor, exponent_exponent, &dataPoints);
		mtz->applyScaleFactor(k, bFactor);

		vector<double> smoothBins;
		StatisticsManager::generateResolutionBins(50, 1.5, 20, &smoothBins);

		vector<double> newx, x2;
		vector<double> newy, y2;

		for (int j = 0; j < dataPoints.size(); j++)
		{
			x2.push_back(dataPoints[j].first);
			y2.push_back(dataPoints[j].second);
		}

		xs.push_back(x2);
		ys.push_back(y2);

		newx.push_back(0);
		newy.push_back(1 / k);

		for (int shell = 0; shell < smoothBins.size() - 1; shell++)
		{
			double midResolution = StatisticsManager::midPointBetweenResolutions(
					smoothBins[shell], smoothBins[shell + 1]);

			double scale = 1
					/ Miller::scaleForScaleAndBFactor(k, bFactor, midResolution,
							exponent_exponent);

			double one_over_d_squared = pow(1 / midResolution, 2);

			newx.push_back(one_over_d_squared);
			newy.push_back(scale);
		}

		xs.push_back(newx);
		ys.push_back(newy);

		plot(filenames[i][2], properties, xs, ys);

		xs.clear();
		ys.clear();

		this->correlationPlot(filenames[i][1], 12000, 12000);

		double rsplit = mtz->rSplitWithManager(
				MtzManager::getReferenceManager());

		double correl = mtz->correlationWithManager(
				MtzManager::getReferenceManager());

		results << "Results\t" << k << "\t" << bFactor << "\t" << preRSplit
				<< "\t" << rsplit << "\t" << preCorrel << "\t" << correl
				<< std::endl;
	}

	std::cout << results.str();*/
}

void GraphDrawer::plotPolarisation(vector<MtzPtr> mtzs)
{
    int count = 36;
    int divide = 360 / count;
    
    std::map<int, double> histogram;
    std::map<int, int> counts;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        for (int j = 0; j < mtzs[i]->reflectionCount(); j++)
        {
            if (!mtzs[i]->reflection(j)->betweenResolutions(1.8, 1.4))
                continue;
            
            for (int k = 0; k < mtzs[i]->reflection(j)->millerCount(); k++)
            {
                MillerPtr miller = mtzs[i]->reflection(j)->miller(k);
                
                if (miller->getRawIntensity() < 0)
                    continue;
                
                vec hkl = new_vector(miller->getH(), miller->getK(), miller->getL());
                mtzs[i]->getMatrix()->multiplyVector(&hkl);
                
                double angle = cartesian_to_angle(hkl.h, hkl.k);
                double degrees = angle * 180 / M_PI;
                
                int category = (int)(degrees / divide) * divide;
                
                if (histogram.count(category) == 0)
                {
                    histogram[category] = 0;
                    counts[category] = 0;
                }
                
                histogram[category] += miller->getRawIntensity();
                counts[category]++;
            }
        }
    }
    
    vector<double> xs, ys;
    GraphMap graphMap;
    graphMap["yMin"] = 0;
    graphMap["yMax"] = 200;
    graphMap["title"] = "Average raw intensity vs angle on detector";
    graphMap["xTitle"] = "Angle on detector";
    graphMap["yTitle"] = "Average raw intensity";
    
    int num = 0;
    
    for (int i = 0; i <= 0; i++)
    {
        for (std::map<int, double>::iterator it = histogram.begin(); it != histogram.end(); it++)
        {
            if (i == 0)
                histogram[it->first] /= counts[it->first];
            
            xs.push_back(it->first + (i * 360));
            ys.push_back(histogram[it->first]);
            num++;
        }
    }
    
    std::cout << "Number of X values: " << num << std::endl;
    
    plot("polarisation", graphMap, xs, ys);
}

void GraphDrawer::partialityPlot(std::string filename, GraphMap properties, double maxRes)
{
    double correl = mtz->correlation();
    mtz->setActiveAmbiguity(1);
    double invCorrel = mtz->correlation();
    
    std::cout << "Ambiguity 0: " << correl << ", ambiguity 1: " << invCorrel << std::endl;
    
    if (correl > invCorrel)
        mtz->setActiveAmbiguity(0);
    
    double gradient = mtz->gradientAgainstManager(MtzManager::getReferenceManager());
    mtz->applyScaleFactor(gradient);
    
    properties["title"] = "Partiality plot " + mtz->getFilename();
	properties["xTitle"] = "Wavelength (Å)";
	properties["yTitle"] = "Fraction of merged data set (%)";
	properties["plotType"] = "fill";

	vector<double> resolutions;
	StatisticsManager::generateResolutionBins(0, maxRes, 4, &resolutions);

	vector<std::string> files;

	vector<Partial> partials;
	StatisticsManager::twoImagePartialityStatsWritten(&partials,
			&MtzManager::getReferenceManager(), &mtz);

	std::cout << "Total number of reflections in MTZ: " << partials.size() << std::endl;

	double maxPercentage = 1000;
	std::sort(partials.begin(), partials.end(), sortByPercentage);

	maxPercentage = partials[partials.size() - 10].percentage;

    if (maxPercentage < 500 || !std::isfinite(maxPercentage))
		maxPercentage = 500;

	maxPercentage = 500;

	std::sort(partials.begin(), partials.end(), sortByWavelength);

	double minX = partials[50].wavelength;
	double maxX = partials[partials.size() - 50].wavelength;
    double middle = partials[partials.size() / 2].wavelength;

    double wantedWidth = 0.08;
    
  //  std::cout << "minX: " << minX << ", maxX: " << maxX << ", middle: " << middle << std::endl;
    
	if (maxX - minX > wantedWidth)
	{
        maxX = middle + wantedWidth / 2;
		minX = middle - wantedWidth / 2;
	}

	int resCount = (int)resolutions.size();

	properties["xMin"] = minX;
	properties["xMax"] = maxX;
	properties["yMin2"] = 0;
	properties["yMax2"] = 5.0;
	properties["yMax"] = maxPercentage;

	vector<double> xs, xs2, xs3, xs4, ys, ys2, ys3, ys4, scatterX, scatterY;

    std::cout << std::setw(12) << "Low res" << std::setw(12) << "High res" << std::setw(12) << "Num refl." << std::endl;
    
    
	for (int shell = 0; shell < resCount - 1; shell++)
	{
		xs.clear();
		ys.clear();
        xs2.clear();
		ys2.clear();
        xs3.clear();
        ys3.clear();

		double lowRes = 1 / resolutions[shell];
		double highRes = 1 / resolutions[shell + 1];

        std::ofstream partLog;
        partLog.open("partiality_" + i_to_str(shell) + ".csv");
        partLog << "h,k,l,wavelength,partiality,percentage,intensity,resolution" << std::endl;
        
        int count = 0;
        
		for (int i = 0; i < partials.size(); i++)
		{
			if (partials[i].resolution > lowRes
					&& partials[i].resolution < highRes)
			{
                if (partials[i].wavelength != partials[i].wavelength
                    || partials[i].percentage != partials[i].percentage
                    || partials[i].partiality != partials[i].partiality)
                    continue;
                
                if (partials[i].percentage > maxPercentage)
                    continue;

                double h = partials[i].miller->getH();
                double k = partials[i].miller->getK();
                double l = partials[i].miller->getL();
                double wavelength = partials[i].wavelength;
                double partiality = partials[i].partiality;
                double percentage = partials[i].percentage;
                double intensity = partials[i].miller->getRawestIntensity();
                double resolution = partials[i].resolution;

                partLog << h << "," << k << "," << l << "," << wavelength << "," << partiality << "," <<
                percentage << "," << intensity << "," << resolution << std::endl;
                count++;
            }
		}
        
        partLog.close();

		std::cout << std::setw(12) << 1 / lowRes << std::setw(12) <<  1 / highRes << std::setw(12) << count
				<< std::endl;
	}
/*
	int milliseconds = 1200 / resCount;

	std::string animate = "convert -delay " + i_to_str(milliseconds) + " -loop 0 ";
	std::string rmv = "rm ";

	for (int i = 0; i < resCount - 1; i++)
	{
		animate += files[i] + " ";
		rmv += files[i] + " ";
	}

	rmv += "\n";

	std::string fullName = generateFilename(filename, "gif");

	animate += fullName;
	animate += "\n";
//	std::cout << animate << std::endl;

	std::string all = animate + rmv;

	system(all.c_str());

	std::cout << "N: plot " << fullName << std::endl;

	properties["title"] = "Partiality scatter";
	properties["xTitle"] = "Predicted partiality";
	properties["yTitle"] = "Fraction of merged data set (%)";
	properties["plotType"] = "point";

	properties["xMin"] = 0;
	properties["xMax"] = 1.2;
	properties["yMin"] = 0;
    properties["yMax"] = maxPercentage;
*/
//	plot(filename + "_scatter", properties, scatterX, scatterY);
}

void GraphDrawer::plotPartialityStats()
{
 //   double threshold = 100;
    GraphMap map;
    map["xMin"] = 0.0;
    map["xMax"] = 1;
    map["yMax"] = 3;
    map["xTitle"] = "Calculated partiality";
    map["yTitle"] = "Proportion of merged value";
    map["plotType"] = "point";
    
    vector<double> xs, ys;
    
    for (int i = 0; i < mtz->reflectionCount(); i++)
    {
        bool okay = false;
        
        for (int j = 0; j < mtz->reflection(i)->millerCount(); j++)
        {
            MillerPtr miller = mtz->reflection(i)->miller(j);
            
            if (miller->getH() == 1 && miller->getK() == 0 && miller->getL() == 4)
            {
                okay = true;
            }
        }
        
        if (okay == false)
            continue;
    /*
        if (mtz->reflection(i)->millerCount() < 2)
            continue;
        
        if (!mtz->reflection(i)->anyAccepted())
            continue;
        
        if (mtz->reflection(i)->meanIntensity() < threshold)
            continue;
        
        if (!mtz->reflection(i)->betweenResolutions(0, 2.0))
            continue;
        */
     //   double max_intensity = mtz->reflection(i)->meanIntensity();
        
        for (int j = 0; j < mtz->reflection(i)->millerCount(); j++)
        {
            double partiality = mtz->reflection(i)->miller(j)->getPartiality();
/*
            if (partiality < 0.2 || partiality > 1)
            	continue;*/
           // double wavelength = mtz->reflection(i)->miller(j)->getWavelength();
            
            double percentage = mtz->reflection(i)->miller(j)->getRawIntensity();
/*            if (percentage < 0)
                percentage = 0;*/
            
        //    if (rand() % 20 == 0)
            {
                xs.push_back(partiality);
                ys.push_back(percentage);
            }
        }
    }
    
    plot("partialityStats", map, xs, ys);
}

void GraphDrawer::plotOrientationStats(vector<MtzPtr> mtzs)
{
    cctbx::miller::index<> genericIndex = cctbx::miller::index<>(0, 0, 1);
    cctbx::sgtbx::space_group *spaceGroup = mtzs[0]->reflection(0)->getSpaceGroup();
    
    vector<double> xs, ys, zs;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MatrixPtr matrix = mtzs[i]->getMatrix();
        
        cctbx::miller::sym_equiv_indices indices = cctbx::miller::sym_equiv_indices(*spaceGroup, genericIndex);
        
        cctbx::miller::index<double> position = matrix->multiplyIndex(&genericIndex);
        
        vec pos = new_vector(position[0], position[1], position[2]);
        scale_vector_to_distance(&pos, 1);
        
        std::cout << mtzs[i]->getFilename() << "\t" << pos.h << "\t" << pos.k << "\t" << pos.l << std::endl;
    }

    for (int i = 0; i < mtzs.size(); i++)
    {
        MatrixPtr matrix = mtzs[i]->getMatrix();
        
        cctbx::miller::sym_equiv_indices indices = cctbx::miller::sym_equiv_indices(*spaceGroup, genericIndex);
        
   //     unsigned long size = indices.indices().size();
        
   //     cctbx::miller::index<double> position = matrix->multiplyIndex(&genericIndex);
        
     //   double h = position[0];
     //   double k = position[1];
     //   double l = position[2];
        
        /*
        for (int i = 0; i < size; i++)
        {
            cctbx::miller::index<> asymIndex = indices.indices()[i].h();
            cctbx::miller::index<double> position = matrix->multiplyIndex(&asymIndex);
            
            xs.push_back(position[0]);
            ys.push_back(position[1]);
            zs.push_back(position[2]);
            
            std::cout << position[0] << "\t" << position[1] << "\t" << position[2] << std::endl;
            std::cout << -position[0] << "\t" << -position[1] << "\t" << -position[2] << std::endl;
        }*/
    }
}

void GraphDrawer::plotReflectionFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l)
{
    int index = Reflection::reflectionIdForCoordinates(h, k, l);
    
    std::cout << "filename,intensity,+-,resolution" << std::endl;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MtzPtr mtz = mtzs[i];
        
        for (int j = 0; j < mtz->reflectionCount(); j++)
        {
            Reflection *refl = mtz->reflection(j);
            
            long unsigned int thisID = refl->getReflId();
            
            if (thisID != index)
            {
                continue;
            }
            
            std::cout << mtz->getFilename().substr(0, 4) << "\t" << refl->meanIntensity() << "\t" << refl->meanSigma() << "\t" << refl->miller(0)->getPhase() << "\t" << 1 / refl->getResolution() << std::endl;
        }
    }
}

GraphDrawer::~GraphDrawer()
{
// TODO Auto-generated destructor stub
}

