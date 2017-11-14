//
//  IndexManager.h
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexManager__
#define __cppxfel__IndexManager__

#include "cmtzlib.h"
#include "csymlib.h"
#include <stdio.h>
#include <vector>
#include <map>
#include "Image.h"
#include "Vector.h"
#include "parameters.h"
#include "Logger.h"
#include "Holder.h"
#include "LoggableObject.h"
#include <mutex>

typedef std::map<int, std::pair<int, int> > PowderHistogram;

class IndexManager : LoggableObject
{
protected:
    UnitCellLatticePtr lattice;
    std::vector<ImagePtr> images;
    std::vector<ImagePtr> mergeImages;
    std::vector<double> unitCell;
    std::vector<MatrixPtr> symOperators;
    MatrixPtr unitCellOnly;
    MatrixPtr unitCellMatrix;
    MatrixPtr unitCellMatrixInverse;
    Reflection *newReflection;
    CSym::CCP4SPG *spaceGroup;
    int spaceGroupNum;
    std::vector<MtzPtr> mtzs;
    double minimumTrustDistance;
    double minimumTrustAngle;
    double solutionAngleSpread;
    double lastTime;
    ImagePtr getNextImage();
    int nextImage;
    std::mutex indexMutex;
    bool modifyParameters();

    void updateAllSpots();
    static double metrologyTarget(void *object);
    bool matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2);
    int indexOneImage(ImagePtr image, std::vector<MtzPtr> *mtzSubset);
    double maxMillerIndexTrial;
    double maxDistance;
    double smallestDistance;
    double minReciprocalDistance;
    PowderHistogram generatePowderHistogram();
    std::vector<VectorDistance> vectorDistances;
    std::vector<IOMRefinerPtr> consolidateOrientations(ImagePtr image1, ImagePtr image2, int *oneHand, int *otherHand, int *both);
public:
    ImagePtr getImage(int i)
    {
        return images[i];
    }

    std::vector<MtzPtr> getMtzs()
    {
        return mtzs;
    }

    void setMergeImages(std::vector<ImagePtr> otherImages)
    {
        mergeImages = otherImages;
    }

    void combineLists();
    void indexingParameterAnalysis();
    void refineMetrology();
    static void indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset);
    void index();
    void indexFromScratch();
    void powderPattern();
    IndexManager(std::vector<ImagePtr>images);
};


#endif /* defined(__cppxfel__IndexManager__) */
