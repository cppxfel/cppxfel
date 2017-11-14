/*
 * Holder.h
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#ifndef HOLDER_H_
#define HOLDER_H_

#include "csymlib.h"

class Reflection;

#include "parameters.h"
#include "headers/csymlib.h"
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/uctbx.h>

class Miller;

using cctbx::sgtbx::space_group;
using cctbx::sgtbx::space_group_symbols;
using cctbx::sgtbx::space_group_type;
using cctbx::sgtbx::reciprocal_space::asu;

class Reflection
{
private:
        vector<MillerPtr> millers;
        double resolution;
    static space_group spaceGroup;
    static unsigned char spgNum;
    static space_group_type spgType;
    static asu asymmetricUnit;
    static bool hasSetup;
    static std::mutex setupMutex;
    static bool setupUnitCell;
    static std::vector<MatrixPtr> flipMatrices;
    unsigned char activeAmbiguity;
    vector<long unsigned int> reflectionIds;
    static cctbx::uctbx::unit_cell unitCell;
    MutexPtr millerMutex;
public:
    Reflection(float *unitCell = NULL, CSym::CCP4SPG *group = NULL);
    void setUnitCell(float *unitCell);
        virtual ~Reflection();

    MillerPtr acceptedMiller(int num);
        MillerPtr miller(int i);
    void printDescription();
        void addMiller(MillerPtr miller);
    void addMillerCarefully(MillerPtr miller);
        int millerCount();
        Reflection *copy(bool copyMillers = false);

        static int indexForReflection(int h, int k, int l, CSym::CCP4SPG *lowspgroup, bool inverted);
    static int reflectionIdForCoordinates(int h, int k, int l);

    int checkOverlaps();
    int checkSpotOverlaps(std::vector<SpotPtr> *spots, bool actuallyDelete = true);
    void reflectionDescription();
        void calculateResolution(MtzManager *mtz);
        void clearMillers();
    void removeMiller(int index);

        double meanIntensityWithExclusion(std::string *filename, int start = 0, int end = 0);
        double meanIntensity(bool withCutoff = true, int start = 0, int end = 0);
        double meanSigma();
        double meanSigma(bool friedel);
        double meanPartiality(bool withCutoff = true);
        double mergedIntensity(WeightType weighting);
    double meanWeight(bool cutoff = true);
    double rMergeContribution(double *numerator, double *denominator);

        bool betweenResolutions(double lowAngstroms, double highAngstroms);
        bool anyAccepted();
        int acceptedCount();
        int rejectCount();

        void merge(WeightType weighting, double *intensity, double *sigma, bool calculateRejections);
        double standardDeviation(WeightType weighting);

    void detailedDescription();
    double observedPartiality(MtzManager *reference, Miller *miller);

    double mergeSigma();

    int reflectionIdForMiller(cctbx::miller::index<> cctbxMiller);
    void generateReflectionIds();

    asu *getAsymmetricUnit()
    {
        return &asymmetricUnit;
    }

    static space_group *getSpaceGroup()
    {
        return &spaceGroup;
    }

    void incrementAmbiguity();
    static int ambiguityCount();
    static MatrixPtr matrixForAmbiguity(int i);

    static void setSpaceGroup(int spaceGroupNum);
    void setUnitCellDouble(double *theUnitCell);

    void resetFlip();
    void setFlip(int i);
    void setFlipAsActiveAmbiguity();

    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }

    void setActiveAmbiguity(int newAmbiguity)
    {
        activeAmbiguity = newAmbiguity;
    }

        const vector<MillerPtr>& getMillers() const
        {
                return millers;
        }

        void setMillers(const vector<MillerPtr>& millers)
        {
                this->millers = millers;
        }

    long unsigned int getReflId()
    {
   /*     if (activeAmbiguity > ambiguityCount())
        {
            std::cout << "Active ambiguity error" << std::endl;
        }
        */
        return reflectionIds[activeAmbiguity];
    }

        double getResolution() const
        {
                return resolution;
        }

        void setResolution(double resolution)
        {
                this->resolution = resolution;
        }

    static MatrixPtr getFlipMatrix(int i);
};

#endif /* HOLDER_H_ */
