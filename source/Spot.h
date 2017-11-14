/*
 * Spot.h
 *
 *  Created on: 30 Dec 2014
 *      Author: helenginn
 */

#ifndef SPOT_H_
#define SPOT_H_

#include <vector>
#include "Image.h"
#include "parameters.h"
#include "Panel.h"
#include "Vector.h"
#include "LoggableObject.h"

class Spot : LoggableObject
{
private:
        vector<vector<double> > probe;
        ImageWeakPtr parentImage;
    double angleDetectorPlane;
    bool setAngle;
    bool checked;
    int successfulCommonLines;
    double correctedX; double correctedY;
    double x; double y;
    bool rejected;
    int height;
    int length;
    int background;
    static double maxResolution;
    static double minIntensity;
    static double minCorrelation;

public:
        Spot(ImagePtr image);
        virtual ~Spot();

    static void spotsAndVectorsToResolution(double lowRes, double highRes, std::vector<SpotPtr> spots, std::vector<SpotVectorPtr> spotVectors, std::vector<SpotPtr> *lowResSpots, std::vector<SpotVectorPtr> *lowResSpotVectors);
    double weight();
        double maximumLift(ImagePtr image, int x, int y, bool ignoreCovers);
        double maximumLift(ImagePtr image, int x, int y);
        void makeProbe(int height, int background, int size);
        void setXY(double x, double y);
    void setXYFromEstimatedVector(vec hkl);
        double scatteringAngle(ImagePtr image = ImagePtr());
        bool isAcceptable(ImagePtr image);
        static void sortSpots(vector<Spot *> *spots);
        static bool spotComparison(Spot *a, Spot *b);
    double angleFromSpotToCentre(double centreX, double centreY);
    double angleInPlaneOfDetector(double centreX = 0, double centreY = 0, vec upBeam = new_vector(0, 1, 0));
    double resolution();
    bool isOnSameLineAsSpot(SpotPtr otherSpot, double toleranceDegrees);
    static void writeDatFromSpots(std::string filename, std::vector<SpotPtr> spots);
    bool isSameAs(SpotPtr spot2);
    double closeToSecondSpot(SpotPtr spot2, double squareMinDistance);

    Coord getXY();
    double getX(bool update = false);
    double getY(bool update = false);
    Coord getRawXY();
    vec estimatedVector();
    void setUpdate();
    bool focusOnNearbySpot(double maxShift, double trialX, double trialY, int round = 0);

    std::string spotLine();

    void setRejected(bool isRejected = true)
    {
        rejected = isRejected;
    }

    bool isRejected()
    {
        return rejected;
    }

    int successfulLineCount()
    {
        return successfulCommonLines;
    }

    void addSuccessfulCommonLineCount()
    {
        successfulCommonLines++;
    }

    void deleteSuccessfulCommonLineCount()
    {
        successfulCommonLines--;
    }

    bool isChecked()
    {
        return checked;
    }

    void setChecked(bool newCheck = true)
    {
        checked = newCheck;
    }

    ImagePtr getParentImage()
        {
                return parentImage.lock();
        }

        void setParentImage(ImagePtr parentImage)
        {
                this->parentImage = parentImage;
        }


};

#endif /* SPOT_H_ */
