//
//  DIALSImage.h
//  cppxfel
//
//  Created by Helen Ginn on 11/03/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__DIALSImage__
#define __cppxfel__DIALSImage__

#include <stdio.h>
#include "Image.h"
#include <dxtbx/model/detector.h>
#include <scitbx/vec3.h>
#include "parameters.h"

typedef std::shared_ptr<dxtbx::model::Detector> DxtbxDetectorPtr;
using scitbx::vec3;
class Matrix;
class InputFileParser;

namespace cppxfel
{
    class DIALSImage : public Image
    {
    private:
        // personal parser object to store parameters for refinement
        static InputFileParser *personalParser;
        
        // storage of matrices leading to solutions
        std::vector<MatrixPtr> solutions;
        
        // local copy of pointer to dxtbx-style detector.
        DxtbxDetectorPtr dxtbxDetector;
        
        // convert from one format into another and also change origin
        vec3<double> cppxfelVecToScitbxVec3(vec hkl);
        vec scitbxVec3ToCppxfelVec(vec3<double> hkl3);
        
        // converter which uses dxtbxDetector object to convert coordinates.
        virtual std::pair<double, double> reciprocalCoordinatesToPixels(vec hkl);
        
        // overrides for indexing for suitability without checking solutions.
        virtual IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr);
        virtual bool checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast = false);
    public:
        // call first, before doing anything else
        static void init();
        
        // call to set space group (integer between 1 and 210) and unit cell dimensions (e.g. 160 40 125 90 90 90)
        static void setCrystalParameters(int spaceGroupNum, double params[6]);
        
        // call to set maximum reciprocal distance in inverse Angstroms
        // (e.g. 0.1 would be up to 10 Å away from beam centre, or equivalent).
        static void setMaxReciprocalDistance(double newMaxReciprocalDistance);
        
        // call to set the estimation of the domain size in inverse Angstroms
        // e.g. 0.0001 would correspond to a domain size of 1 micrometre.
        static void setInverseDomainSize(double inverseDomainSize);
        
        // Initialise an image object:
        // you MUST call all parameter settings (function outlines are above this line)
        // before attempting to initialise an image object.
        // not yet sure if the wavelength/distance values are required - requires testing
        // tidy up this comment when you know the details!
        DIALSImage(DxtbxDetectorPtr detector, double wavelength, double distance = 0)
        : Image("", wavelength, distance)
        {
            dxtbxDetector = detector;
        }

        // pass in array of crystal-as-origin rays to be turned into spots.
        // rays should start from crystal as (0, 0, 0) and pass through the reciprocal space coordinate
        // (note: ask Helen to change this if there is a more appropriate definition...
        // I went with what I thought was the DIALS model)
        void addSpots(std::vector<vec3<double> > rays);

        // at this point, call the business function
        void findIndexingSolutions() { Image::findIndexingSolutions(); };
        
        // how many solutions did I find?
        // must only be called after findIndexingSolutions() if you don't want disappointment
        int solutionCount()
        {
            return (int)solutions.size();
        }
        
        // get a rotation matrix (i.e. U of the UB matrix)
        // must only be called after findIndexingSolutions() if you don't want disappointment
        void getSolution(int i, double *matrixParams[9]);
        
        // what kind of image am I? (cppxfel or dials) - don't worry about this
        virtual ImageClass getClass()
        {
            return ImageClassDIALS;
        }
    };
};


#endif /* defined(__cppxfel__DIALSImage__) */
