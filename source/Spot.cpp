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

Spot::Spot(Image *image)
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
    
	makeProbe(500, 1);
}

Spot::~Spot()
{

}

bool Spot::isAcceptable(Image *image)
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
	return maximumLift(parentImage, x, y, true);
}

double Spot::maximumLift(Image *image, int x, int y)
{
	return maximumLift(image, x, y, false);
}

double Spot::maximumLift(Image *image, int x, int y, bool ignoreCovers)
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

void Spot::makeProbe(int height, int length)
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

			if (fraction > 1)
				probe[i][j] = 0;
			else
				probe[i][j] = (1 - fraction) * height; //(double)height * y;

		}

	}
}

void Spot::setXY(int x, int y)
{
	this->x = x;
	this->y = y;
}

double Spot::scatteringAngle(Image *image)
{
    if (image == NULL)
        image = this->parentImage;
    
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
    double twoTheta = asin(scatteringAngle());
    double theta = twoTheta / 2;
    double wavelength = parentImage->getWavelength();
    
    double d = wavelength / (2 * sin(theta));
    
    return 1 / d;
}

double Spot::angleFromSpotToCentre(double centreX, double centreY)
{
    double beamX = parentImage->getBeamX();
    double beamY = parentImage->getBeamY();
    
    vec beamCentre = new_vector(beamX, beamY, 0);
    vec circleCentre = new_vector(centreX, centreY, 0);
    vec circleCentreToBeam = vector_between_vectors(circleCentre, beamCentre);
    
    return angleInPlaneOfDetector(centreX, centreY, circleCentreToBeam);
}

double Spot::angleInPlaneOfDetector(double centreX, double centreY, vec upBeam)
{
    if (centreX == 0 && centreY == 0)
    {
        centreX = parentImage->getBeamX();
        centreY = parentImage->getBeamY();
    }
    
    vec centre = new_vector(centreX, centreY, 0);
    vec spotVec = new_vector(getX(), getY(), 0);
    vec spotVecFromCentre = vector_between_vectors(centre, spotVec);
    /*
    if (centreX == 679.5 && centreY == 843.75)
    {
        std::ostringstream logged;
        logged << "spotXY:\t" << getX() << "\t" << getY() << std::endl;
        logged << "beam:\t" << parentImage->getBeamX() << "\t" << parentImage->getBeamY() << std::endl;
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

Coord Spot::getXY()
{
    Coord shift = Panel::shiftForSpot(this);
    
    return std::make_pair(x + shift.first, y + shift.second);
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
    double beamX = parentImage->getBeamX() * parentImage->getMmPerPixel();
    double beamY = parentImage->getBeamY() * parentImage->getMmPerPixel();
    
    double wavelength = parentImage->getWavelength();
    
    double height = parentImage->getYDim();
    
    double mmX = getX() * parentImage->getMmPerPixel();
 //   double mmY = (height - getY()) * parentImage->getMmPerPixel();
    double mmY = getY() * parentImage->getMmPerPixel();
    
    double detector_distance = parentImage->getDetectorDistance();
    
    vec crystalVec = new_vector(beamX, beamY, 0 - detector_distance);
    vec spotVec = new_vector(mmX, mmY, 0);
    vec reciprocalCrystalVec = new_vector(0, 0, 0 - 1 / wavelength);
    
    vec crystalToSpot = vector_between_vectors(crystalVec, spotVec);
    scale_vector_to_distance(&crystalToSpot, 1 / wavelength);
    add_vector_to_vector(&reciprocalCrystalVec, crystalToSpot);
    
    reciprocalCrystalVec.k *= -1;
    
    return reciprocalCrystalVec;
}


