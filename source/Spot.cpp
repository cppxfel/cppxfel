/*
 * Spot.cpp
 *
 *  Created on: 30 Dec 2014
 *      Author: helenginn
 */

#include "Spot.h"
#include <iostream>
#include <cmath>
#include <cfloat>
#include "Image.h"
#include <algorithm>
#include "Vector.h"
#include <fstream>
#include "FileParser.h"

double Spot::maxResolution = 0;
double Spot::minIntensity = 0;
double Spot::minCorrelation = 0;

Spot::Spot(ImagePtr image)
{
	// TODO Auto-generated constructor stub
	probe = vector<vector<double> >();
	parentImage = image;
    angleDetectorPlane = 0;
    setAngle = false;
    checked = false;
    successfulCommonLines = 0;
    correctedX = -1;
    correctedY = -1;
    rejected = false;
    x = 0;
    y = 0;
    height = FileParser::getKey("IMAGE_SPOT_PROBE_HEIGHT", 100);
    background = FileParser::getKey("IMAGE_SPOT_PROBE_BACKGROUND", 10);
    length = FileParser::getKey("IMAGE_SPOT_PROBE_PADDING", 1) * 2 + 1;
    
    if (minCorrelation == 0)
    {
        maxResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 3.0);
        minIntensity = FileParser::getKey("IMAGE_MIN_SPOT_INTENSITY", 600.);
        minCorrelation = FileParser::getKey("IMAGE_MIN_CORRELATION", 0.7);
    }
    
	makeProbe(height, background, length);
}

Spot::~Spot()
{

}

bool Spot::isAcceptable(ImagePtr image)
{
	int length = (int)probe.size();
	int tolerance = (length - 1) / 2;

	for (int i = x - tolerance; i < x + tolerance; i++)
	{
		for (int j = y - tolerance; j < y + tolerance; j++)
		{
			if (!image->accepted(i, j))
				return false;
		}
	}

	return true;
}

double Spot::weight()
{
	return maximumLift(getParentImage(), x, y, true);
}

double Spot::maximumLift(ImagePtr image, int x, int y)
{
	return maximumLift(getParentImage(), x, y, false);
}

double Spot::maximumLift(ImagePtr image, int x, int y, bool ignoreCovers)
{
	int length = (int)probe.size();
	int tolerance = (length - 1) / 2;

	double minDifference = FLT_MAX;
	double penultimate = FLT_MAX;

	for (int i = x - tolerance; i < x + tolerance; i++)
	{
		for (int j = y - tolerance; j < y + tolerance; j++)
		{
			double imageValue = image->valueAt(i, j);

			if (image->coveredBySpot(i, j) && !ignoreCovers)
			{
				imageValue = 0;
			}

			if (imageValue == 0)
				continue;

			int probeX = i - x + tolerance;
			int probeY = j - y + tolerance;
			double probeValue = probe[probeX][probeY];

			double difference = (imageValue - probeValue);

			if (difference < minDifference)
			{
				penultimate = minDifference;
				minDifference = difference;
			}
		}
	}

	return (penultimate > 0 && penultimate != FLT_MAX ? penultimate : 0);
}

bool Spot::focusOnNearbySpot(double maxShift, double trialX, double trialY, int round)
{
  //  logged << "Finding new spot for round " << round << std::endl;
  //  sendLog(round > 1 ? LogLevelDetailed : LogLevelDebug);
        int focusedX = trialX;
    int focusedY = trialY;
    this->getParentImage()->focusOnAverageMax(&focusedX, &focusedY, maxShift);
    setXY(focusedX, focusedY);
  //  logged << "Original position (" << trialX << ", " << trialY << ") focusing on " << focusedX << ", " << focusedY << std::endl;
  //  sendLog(LogLevelDetailed);
    
    if (this->getParentImage()->valueAt(focusedX, focusedY) < minIntensity)
        return false;
    
    double resol = this->resolution();
    
    if (resol > (1. / maxResolution)) return false;
    if (resol < 0) return false;
    
    
    std::vector<double> probeIntensities, realIntensities;
    
    int padding = (length - 1) / 2;
    
    for (int i = -padding; i < padding + 1; i++)
    {
        for (int j = -padding; j < padding + 1; j++)
        {
            if (!(this->getParentImage()->accepted(focusedX + i, focusedY + j)))
            {
                logged << "Unacceptable pixel for round " << round << std::endl;
                return false;
            }
            
            int probeX = i + padding;
            int probeY = j + padding;
            
            double probeIntensity = probe[probeX][probeY];
            double realIntensity = this->getParentImage()->valueAt(focusedX + i, focusedY + j);
            
            if (!(probeIntensity == probeIntensity && realIntensity == realIntensity))
                continue;
            
            probeIntensities.push_back(probeIntensity);
            realIntensities.push_back(realIntensity);
        }
    }
    
    double correlation = correlation_between_vectors(&probeIntensities, &realIntensities);
    
    logged << "Correlation: " << correlation << " for round " << round << std::endl;
    if (round > 1)
        sendLog(LogLevelDetailed);

    if (correlation < minCorrelation)
        return false;
    
    sendLog(LogLevelDetailed);

    return true;
}

void Spot::makeProbe(int height, int background, int length)
{
	for (int i = 0; i < probe.size(); i++)
		probe[i].clear();

	probe.clear();

	if (length % 2 == 0)
		length++;

	int size = (length - 1) / 2;

	int centre = size;

	for (int i = 0; i < length; i++)
	{
		probe.push_back(vector<double>());

		for (int j = 0; j < length; j++)
		{
			probe[i].push_back(0);
		}
	}

	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length; j++)
		{
			int from_centre_i = i - centre;
			int from_centre_j = j - centre;

			double distance_from_centre = sqrt(
					from_centre_i * from_centre_i
							+ from_centre_j * from_centre_j);

			double fraction = distance_from_centre / (double) (size + 1);

            // further out the fraction, the lower the value
			if (fraction > 1)
				probe[i][j] = background;
			else
				probe[i][j] = (1 - fraction) * height + background;

		}

	}
}

void Spot::setXY(int x, int y)
{
	this->x = x;
	this->y = y;
}

double Spot::scatteringAngle(ImagePtr image)
{
    if (image == NULL)
        image = this->getParentImage();
    
	double beamX = image->getBeamX();
	double beamY = image->getBeamY();

	double distance_from_centre = sqrt(pow(getX() - beamX, 2) + pow(getY() - beamY, 2));
	double distance_pixels = distance_from_centre * image->getMmPerPixel();
	double detector_distance = image->getDetectorDistance();

	double sinTwoTheta = sin(distance_pixels / detector_distance);

	return sinTwoTheta;
}

double Spot::resolution()
{
    vec spotVec = this->estimatedVector();
    
    return length_of_vector(spotVec);
}

double Spot::angleFromSpotToCentre(double centreX, double centreY)
{
    double beamX = getParentImage()->getBeamX();
    double beamY = getParentImage()->getBeamY();
    
    vec beamCentre = new_vector(beamX, beamY, 0);
    vec circleCentre = new_vector(centreX, centreY, 0);
    vec circleCentreToBeam = vector_between_vectors(circleCentre, beamCentre);
    
    return angleInPlaneOfDetector(centreX, centreY, circleCentreToBeam);
}

double Spot::angleInPlaneOfDetector(double centreX, double centreY, vec upBeam)
{
    if (centreX == 0 && centreY == 0)
    {
        centreX = getParentImage()->getBeamX();
        centreY = getParentImage()->getBeamY();
    }
    
    vec centre = new_vector(centreX, centreY, 0);
    vec spotVec = new_vector(getX(), getY(), 0);
    vec spotVecFromCentre = vector_between_vectors(centre, spotVec);
    /*
    if (centreX == 679.5 && centreY == 843.75)
    {
        std::ostringstream logged;
        logged << "spotXY:\t" << getX() << "\t" << getY() << std::endl;
        logged << "beam:\t" << getParentImage()->getBeamX() << "\t" << getParentImage()->getBeamY() << std::endl;
        logged << "upBeam:\t" << upBeam.h << "\t" << upBeam.k << std::endl;
        logged << "centreXY:\t" << centreX << "\t" << centreY << std::endl;
        logged << "spotFromCentre:\t" << spotVecFromCentre.h << "\t" << spotVecFromCentre.k << std::endl;
        Logger::mainLogger->addStream(&logged);
    }*/
    
    angleDetectorPlane = angleBetweenVectors(upBeam, spotVecFromCentre);
    
    return angleDetectorPlane;
}

void Spot::sortSpots(vector<Spot *> *spots)
{
	std::cout << "Sorting spots" << std::endl;
	std::sort(spots->begin(), spots->end(), spotComparison);
}

bool Spot::spotComparison(Spot *a, Spot *b)
{
	return (a->weight() > b->weight());
}

void Spot::setUpdate()
{
    getX(true);
    getY(true);
}

Coord Spot::getXY()
{
    Coord shift = Panel::shiftForSpot(this);
    
 //   logged << x << "," << shift.first << "," << y << "," << shift.second << std::endl;
    
    if (shift.first == FLT_MAX)
    {
        shift.first = 0;
        shift.second = 0;
    }

    return std::make_pair(x - shift.first, y - shift.second);
}

double Spot::getX(bool update)
{
    if (update || correctedX == -1)
    {
        Coord shift = getXY();
        correctedX = shift.first;
    }
    
    return correctedX;
}

double Spot::getY(bool update)
{
    if (update || correctedY == -1)
    {
        Coord shift = getXY();
        correctedY = shift.second;
    }
    
    return correctedY;
}

Coord Spot::getRawXY()
{
    return std::make_pair(x, y);
}

bool Spot::isOnSameLineAsSpot(SpotPtr otherSpot, double toleranceDegrees)
{
    double tolerance = toleranceDegrees * M_PI / 180;
    
    double thisAngle = angleInPlaneOfDetector();
    double otherAngle = otherSpot->angleInPlaneOfDetector();
    
    if (thisAngle < 0)
        thisAngle += M_PI;
    
    if (otherAngle < 0)
        otherAngle += M_PI;
    
    return (thisAngle > otherAngle - tolerance && thisAngle < otherAngle + tolerance);
}

void Spot::writeDatFromSpots(std::string filename, std::vector<SpotPtr> spots)
{
    std::ofstream dat;
    dat.open(filename);
    int count = 0;
    
    for (int i = 0; i < spots.size(); i++)
    {
        dat << "0\t0\t0\t" << spots[i]->angleInPlaneOfDetector() << "\t1\t1\t"
        << spots[i]->x << "\t" << spots[i]->y
        << "\t1" << std::endl;
        
        count++;
    }
    
    dat.close();
}

vec Spot::estimatedVector()
{
 /*   if (lastEstimatedVector.h != 0 && lastEstimatedVector.k != 0 && lastEstimatedVector.l != 0)
    {
        return lastEstimatedVector;
    }
    */
    double beamX = getParentImage()->getBeamX() * getParentImage()->getMmPerPixel();
    double beamY = getParentImage()->getBeamY() * getParentImage()->getMmPerPixel();
    
    double wavelength = getParentImage()->getWavelength();
    
    double height = getParentImage()->getYDim();
    
    double mmX = getX() * getParentImage()->getMmPerPixel();
 //   double mmY = (height - getY()) * getParentImage()->getMmPerPixel();
    double mmY = getY() * getParentImage()->getMmPerPixel();
    
    double detector_distance = getParentImage()->getDetectorDistance();
    
    vec crystalVec = new_vector(beamX, beamY, 0 - detector_distance);
    vec spotVec = new_vector(mmX, mmY, 0);
    vec reciprocalCrystalVec = new_vector(0, 0, 0 - 1 / wavelength);
    
    vec crystalToSpot = vector_between_vectors(crystalVec, spotVec);
    scale_vector_to_distance(&crystalToSpot, 1 / wavelength);
    add_vector_to_vector(&reciprocalCrystalVec, crystalToSpot);
    
    reciprocalCrystalVec.k *= -1;
    
    return reciprocalCrystalVec;
}

std::string Spot::spotLine()
{
    std::ostringstream line;
    
    line << x << "\t" << y << std::endl;
    
    return line.str();
}

bool Spot::isSameAs(SpotPtr spot2)
{
    return (spot2->getX() == getX() && spot2->getY() == getY());
}


