/*
 * Panel.h
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#ifndef PANEL_H_
#define PANEL_H_

#include "parameters.h"
#include <sstream>
#include <iostream>
#include <vector>
#include "LoggableObject.h"
#include <boost/thread/thread.hpp>

typedef enum
{
    PlotTypeAbsolute,
    PlotTypeRelative
} PlotType;

typedef std::pair<double, double> Coord;

class Panel : public LoggableObject
{
private:
	Coord topLeft;
	Coord bottomRight;
	Coord bestShift;
    Coord originalShift;
	Coord tilt;
    double tiltWindowSize;
    double averageIntensity;
    bool tiltHorizontalAxis;
	bool beingWrittenTo;
	int allowedSearchSpace;
    PanelTag tag;
    double swivel;
    double gainScale;
    double height();
    double width();
    static double distanceMultiplier;
    static Coord beamCentre;
    Coord midPoint();

    Coord getSwivelShift(Miller *miller);
    Coord getSwivelCoords(Miller *miller);
    Coord getTiltShift(Miller *miller);
	Coord getTotalShift(Miller *miller);
    
    void fractionalCoordinates(Coord coord, Coord *frac);
    void fractionalCoordinates(Miller *miller, Coord *frac);

	static bool usePanelInfo;
    static vector<PanelPtr> panels;
    static vector<PanelPtr> badPanels;
    bool addMiller(MillerPtr miller);

    void refineAllParameters(double windowSize);
    void centreWindowShift();
    void findAllParameters();
    void findShift(double windowSize, double step, double x = 0, double y = 0);
    void findAxisDependence(double windowSize);

    double tiltShiftScore(double stdev = true);
    static double tiltShiftScoreWrapper(void *object);
    static double swivelShiftScoreWrapper(void *object);

    static double scoreDetectorDistance(void *object);
    
    std::map<boost::thread::id, vector<MillerPtr> > tempMillers;
    
    static void calculateMetrology(PanelPtr thisPanel);
    static std::string printAllThreaded();
    double detectorGain(double *error);

    double stdevScore(double minRes, double maxRes);
    vector<MillerPtr> millers;
	int defaultShift;

public:
	Panel(vector<double> dimensions, PanelTag tag = PanelTagNormal);
    Panel(double x1, double y1, double x2, double y2, PanelTag newTag);
    void init(vector<double> dimensions, PanelTag newTag);
    virtual ~Panel();

    static double scoreBetweenResolutions(double minRes, double maxRes);
    bool isCoordInPanel(Coord coord, Coord *topLeft = NULL, Coord *bottomRight = NULL);
	static void addMillerToPanelArray(MillerPtr miller);
	static PanelPtr panelForMiller(Miller *miller);
    static PanelPtr panelForSpot(Spot *spot);
    static PanelPtr panelForCoord(Coord coord);
	static void setupPanel(PanelPtr panel);
    static void removePanel(PanelPtr panel);
	void plotVectors(int i, PlotType plotType);
	static void plotAll(PlotType plotType);
    static void expectedBeamCentre();
    static void refineDetectorDistance();
    
    void finaliseMillerArray();
    static void finaliseMillerArrays();
    static Coord shiftForMiller(Miller *miller);
    static Coord shiftForSpot(Spot *spot);
    static double scaleForMiller(Miller *miller);
	void print(std::ostringstream *stream);
	static void printToFile(std::string filename);
	static std::string printAll();
    
    static bool hasMillers();
    static int panelCount();
    static void clearAllMillers();
    void clearMillers();
    static void checkDetectorGains();

	static bool shouldUsePanelInfo()
	{
		return usePanelInfo;
	}

	static void setUsePanelInfo(bool useInfo)
	{
		usePanelInfo = useInfo;
	}

	const Coord& getBestShift() const
	{
		return bestShift;
	}
    
    double *pointerToBestShiftX()
    {
        return &(bestShift.first);
    }
    
    double *pointerToBestShiftY()
    {
        return &(bestShift.second);
    }

	void setBestShift(const Coord& bestShift)
	{
		this->bestShift = bestShift;
	}
    
    double getGainScale()
    {
        return gainScale;
    }
    
    static PanelPtr getPanel(int i)
    {
        return panels[i];
    }
    
    PanelTag getTag()
    {
        return tag;
    }
    
    void setTag(PanelTag newBad)
    {
        tag = newBad;
    }
    
};

#endif /* PANEL_H_ */
