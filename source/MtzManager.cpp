#include "MtzManager.h"

#include <string>
#include <iostream>
#include "headers/cmtzlib.h"
#include "csymlib.h"
#include "ccp4_spg.h"
#include "ccp4_general.h"
#include "ccp4_parser.h"
#include "Vector.h"
#include "FileReader.h"
#include <cerrno>
#include <fstream>
#include <iostream>
#include <boost/variant.hpp>
#include "StatisticsManager.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "definitions.h"
#include "FileParser.h"


using namespace CMtz;

using namespace CSym;

MtzManager *MtzManager::referenceManager;
//ScoreType MtzManager::scoreType = ScoreTypeCorrelation;
double MtzManager::highRes = 3.5;
double MtzManager::lowRes = 0;

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>
using cctbx::sgtbx::space_group;
using cctbx::sgtbx::space_group_symbols;
using cctbx::sgtbx::space_group_type;
using cctbx::sgtbx::reciprocal_space::asu;


std::string MtzManager::describeScoreType()
{
    switch (scoreType)
    {
        case ScoreTypeCorrelation:
            return std::string("correl");
        case ScoreTypePartialityCorrelation:
            return std::string("partCC");
        case ScoreTypePartialityLeastSquares:
            return std::string("partLQ");
        case ScoreTypeMinimizeRSplit:
            return std::string("rfactor");
        case ScoreTypeMinimizeRSplitLog:
            return std::string("logR");
        case ScoreTypeCorrelationLog:
            return std::string("logCC");
        case ScoreTypeStandardDeviation:
            return std::string("stdev");
        case ScoreTypeMinimizeRMeas:
            return std::string("rmeas");
        default:
            return std::string("unknown");
    }
}

void MtzManager::makeSuperGaussianLookupTable(double exponent)
{
    if (optimisingExponent)
        return;
    
    if (exponent == lastExponent)
        return;
    
    superGaussianScale = pow((M_PI / 2), (2 / exponent - 1));
    
    superGaussianTable.clear();
    
    const double min = 0;
    const double max = MAX_SUPER_GAUSSIAN;
    const double step = SUPER_GAUSSIAN_STEP;
    int count = 0;
    
    for (double x = min; x < max; x += step)
    {
        double value = super_gaussian(x, 0, superGaussianScale, exponent);
        
        //      std::cout << x << "\t" << value << std::endl;
        
        superGaussianTable.push_back(value);
        
        count++;
    }
    
    lastExponent = exponent;
}

std::string MtzManager::filenameRoot()
{
    vector<std::string> components = FileReader::split(filename, '.');
    
    std::string root = "";
    
    for (int i = 0; i < components.size() - 1; i++)
    {
        root += components[i] + ".";
    }
    
    return root.substr(0, root.size() - 1);
}

double MtzManager::extreme_index(MTZ *mtz, int max)
{
    double extreme = 0;
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    if (max == 0)
    {
        extreme = col_h->min;
        extreme = (extreme > col_k->min) ? col_k->min : extreme;
        extreme = (extreme > col_l->min) ? col_l->min : extreme;
    }
    else
    {
        extreme = col_h->max;
        extreme = (extreme < col_k->max) ? col_k->max : extreme;
        extreme = (extreme < col_l->max) ? col_l->max : extreme;
    }
    
    return extreme;
}

void MtzManager::findMultiplier(MTZ *mtz, int *multiplier, int *offset)
{
    int hkl_min = extreme_index(mtz, 0);
    int hkl_max = extreme_index(mtz, 1);
    
    *offset = (hkl_min < 0) ? -hkl_min : 0;
    
    *multiplier = hkl_max - ((hkl_min < 0) ? hkl_min : 0);
    
    *multiplier = MULTIPLIER;
    *offset = OFFSET;
}

void MtzManager::hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k,
                                     int *l, int *multiplier, int *offset)
{
    if (multiplier == 0)
        findMultiplier(mtz, multiplier, offset);
    
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    *h = adata[col_h->source - 1];
    *k = adata[col_k->source - 1];
    *l = adata[col_l->source - 1];
}

int MtzManager::index_for_reflection(int h, int k, int l, bool inverted)
{
    int _h = h;
    int _k = k;
    int _l = l;
    
    if (inverted)
    {
        if (low_group->spg_num == 197)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 146)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 4 || low_group->spg_num == 21)
        {
            _h = l;
            _k = k;
            _l = h;
        }
    }
    
    ccp4spg_put_in_asu(low_group, _h, _k, _l, &h, &k, &l);
    
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

bool MtzManager::reflection_comparison(Reflection *i, Reflection *j)
{
    return (i->getReflId() < j->getReflId());
}

void MtzManager::insertionSortReflections(void)
{
    std::sort(reflections.begin(), reflections.end(), reflection_comparison);
}

void MtzManager::sortLastReflection(void)
{
    int i, j;
    Reflection *tmp;
    
    i = (int)reflections.size() - 1;
    j = i;
    tmp = reflections[i];
    while (j > 0 && tmp->getReflId() < reflections[j - 1]->getReflId())
    {
        reflections[j] = reflections[j - 1];
        j--;
    }
    reflections[j] = tmp;
    
    return;
    
    
}

void MtzManager::loadParametersMap()
{
    optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH",
                                              !OPTIMISED_WAVELENGTH);
    optimisingBandwidth = FileParser::getKey("OPTIMISING_BANDWIDTH",
                                             !OPTIMISED_BANDWIDTH);
    optimisingMosaicity = FileParser::getKey("OPTIMISING_MOSAICITY",
                                             !OPTIMISED_MOSAICITY);
    optimisingOrientation = FileParser::getKey("OPTIMISING_ORIENTATION",
                                               !OPTIMISED_ROT);
    optimisingRlpSize = FileParser::getKey("OPTIMISING_RLP_SIZE",
                                           !OPTIMISED_SPOT_SIZE);
    optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT",
                                            !OPTIMISED_EXPONENT);
    
    resetDefaultParameters();
    
    stepSizeWavelength = FileParser::getKey("STEP_SIZE_WAVELENGTH",
                                            MEAN_STEP);
    stepSizeBandwidth = FileParser::getKey("STEP_SIZE_BANDWIDTH",
                                           BANDWIDTH_STEP);
    stepSizeMosaicity = FileParser::getKey("STEP_SIZE_MOSAICITY",
                                           MOS_STEP);
    stepSizeRlpSize = FileParser::getKey("STEP_SIZE_RLP_SIZE", SPOT_STEP);
    stepSizeOrientation = FileParser::getKey("STEP_SIZE_ORIENTATION",
                                             ROT_STEP);
    stepSizeOrientABC = FileParser::getKey("STEP_SIZE_ORIENTATION_ABC",
                                             ROT_STEP * 3);
    stepSizeExponent = FileParser::getKey("STEP_SIZE_EXPONENT",
                                          EXPONENT_STEP);
    
    toleranceWavelength = FileParser::getKey("TOLERANCE_WAVELENGTH",
                                             MEAN_TOLERANCE);
    toleranceBandwidth = FileParser::getKey("TOLERANCE_BANDWIDTH",
                                            BANDWIDTH_TOLERANCE);
    toleranceMosaicity = FileParser::getKey("TOLERANCE_MOSAICITY",
                                            MOS_TOLERANCE);
    toleranceRlpSize = FileParser::getKey("TOLERANCE_RLP_SIZE",
                                          SPOT_SIZE_TOLERANCE);
    toleranceOrientation = FileParser::getKey("TOLERANCE_ORIENTATION",
                                              ROT_TOLERANCE);
    toleranceExponent = FileParser::getKey("TOLERANCE_EXPONENT",
                                           EXPO_TOLERANCE);
    
    int defaultScoreInt = FileParser::getKey("DEFAULT_TARGET_FUNCTION",
                                             (int) DEFAULT_SCORE_TYPE);
    defaultScoreType = (ScoreType) defaultScoreInt;
    
    usePartialityFunction = FileParser::getKey("USE_PARTIALITY_FUNCTION",
                                               FURTHER_OPTIMISATION);
    
    maxResolutionAll = FileParser::getKey("MAX_RESOLUTION_ALL",
                                          MAX_OPTIMISATION_RESOLUTION);
    maxResolutionRlpSize = FileParser::getKey("MAX_RESOLUTION_RLP_SIZE",
                                              MAX_SPOT_SIZE_OPT_RESOLUTION);
}

MtzManager::MtzManager(void)
{
    lastReference = NULL;
    failedCount = 0;
    filename = "";
    reflections.resize(0);
    low_group = NULL;
    bandwidth = INITIAL_BANDWIDTH;
    hRot = 0;
    kRot = 0;
    aRot = 0;
    bRot = 0;
    cRot = 0;
    mosaicity = INITIAL_MOSAICITY;
    spotSize = INITIAL_SPOT_SIZE;
    wavelength = 0;
    refCorrelation = 0;
    inverse = false;
    flipped = false;
    exponent = INITIAL_EXPONENT;
    finalised = false;
    scoreType = ScoreTypeCorrelation;
    trust = TrustLevelBad;
    maxResolutionAll = MAX_OPTIMISATION_RESOLUTION;
    maxResolutionRlpSize = MAX_SPOT_SIZE_OPT_RESOLUTION;
    freePass = false;
    defaultScoreType = DEFAULT_SCORE_TYPE;
    usePartialityFunction = FURTHER_OPTIMISATION;
    rejected = false;
    scale = 1;
    externalScale = -1;
    lastExponent = 0;
    previousReference = NULL;
    previousAmbiguity = -1;
    activeAmbiguity = 0;
    bFactor = 0;
    setInitialValues = false;
    refPartCorrel = 0;
    penaltyWeight = 0.0;
    penaltyResolution = 2.5;
    int rotMode = FileParser::getKey("ROTATION_MODE", 0);
    rotationMode = (RotationMode)rotMode;
    
    optimisingWavelength = !OPTIMISED_WAVELENGTH;
    optimisingBandwidth = !OPTIMISED_BANDWIDTH;
    optimisingMosaicity = !OPTIMISED_MOSAICITY;
    optimisingOrientation = !OPTIMISED_ROT;
    optimisingRlpSize = !OPTIMISED_SPOT_SIZE;
    optimisingExponent = !OPTIMISED_EXPONENT;
    
    stepSizeWavelength = MEAN_STEP;
    stepSizeBandwidth = BANDWIDTH_STEP;
    stepSizeMosaicity = MOS_STEP;
    stepSizeRlpSize = SPOT_STEP;
    stepSizeOrientation = ROT_STEP;
    stepSizeExponent = EXPONENT_STEP;
    
    toleranceWavelength = MEAN_TOLERANCE;
    toleranceBandwidth = BANDWIDTH_TOLERANCE;
    toleranceMosaicity = MOS_TOLERANCE;
    toleranceRlpSize = SPOT_SIZE_TOLERANCE;
    toleranceOrientation = ROT_TOLERANCE;
    toleranceExponent = EXPO_TOLERANCE;
    
    usingFixedWavelength = USE_FIXED_WAVELENGTH;
    
    matrix = MatrixPtr();
}

MtzPtr MtzManager::copy()
{
    MtzPtr newManager = MtzPtr(new MtzManager());
    newManager->filename = filename;
    double lowNum = low_group->spg_num;
    newManager->setSpaceGroup(lowNum);
    newManager->bandwidth = bandwidth;
    newManager->hRot = hRot;
    newManager->kRot = kRot;
    newManager->mosaicity = mosaicity;
    newManager->spotSize = spotSize;
    newManager->wavelength = wavelength;
    newManager->refCorrelation = refCorrelation;
    newManager->flipped = flipped;
    newManager->exponent = exponent;
    newManager->matrix = matrix;
    
    for (int i = 0; i < 3; i++)
    {
        newManager->cellDim[i] = cellDim[i];
        newManager->cellAngles[i] = cellAngles[i];
        
    }
    
    for (int i = 0; i < reflections.size(); i++)
    {
        Reflection *newReflection = reflections[i]->copy();
        newManager->reflections.push_back(newReflection);
    }
    
    return newManager;
}


void MtzManager::setUnitCell(double a, double b, double c, double alpha, double beta,
                             double gamma)
{
    cellDim[0] = a;
    cellDim[1] = b;
    cellDim[2] = c;
    cellAngles[0] = alpha;
    cellAngles[1] = beta;
    cellAngles[2] = gamma;
}

void MtzManager::setUnitCell(vector<double> unitCell)
{
    cellDim[0] = unitCell[0];
    cellDim[1] = unitCell[1];
    cellDim[2] = unitCell[2];
    cellAngles[0] = unitCell[3];
    cellAngles[1] = unitCell[4];
    cellAngles[2] = unitCell[5];
}

void MtzManager::clearReflections()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        delete reflections[i];
    }
    
    reflections.clear();
    vector<Reflection *>().swap(reflections);
}

void MtzManager::removeReflection(int i)
{
    delete reflections[i];
    
    reflections.erase(reflections.begin() + i);
}

void MtzManager::addReflection(Reflection *reflection)
{
    Reflection *newReflection = NULL;
    int insertionPoint = findReflectionWithId(reflection->getReflId(), &newReflection, true);
    
    reflections.insert(reflections.begin() + insertionPoint, reflection);
    
 //   reflections.push_back(reflection);
 //   this->sortLastReflection();
}

void MtzManager::addReflections(vector<Reflection *>reflections)
{
    for (int i = 0; i < reflections.size(); i++)
    {
        addReflection(reflections[i]);
    }
}


int MtzManager::symmetryRelatedReflectionCount()
{
    int count = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (reflection(i)->millerCount() > 1)
            count++;
    }
    
    return count;
}

int MtzManager::refinedParameterCount()
{
    int count = 0;
    
    count += optimisingBandwidth;
    count += optimisingExponent;
    count += optimisingMosaicity;
    count += optimisingOrientation * 2;
    count += optimisingRlpSize;
    count += optimisingWavelength;
    
    return count;
}

void MtzManager::setMatrix(double *components)
{
    matrix = MatrixPtr(new Matrix(components));
}

void MtzManager::setMatrix(MatrixPtr newMat)
{
    matrix = newMat;
}

void MtzManager::setDefaultMatrix()
{
    double a = cellDim[0];
    double b = cellDim[1];
    double c = cellDim[2];
    double alpha = cellAngles[0];
    double beta = cellAngles[1];
    double gamma = cellAngles[2];
    
    MatrixPtr mat = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);
    matrix = mat->copy();
}

void MtzManager::getUnitCell(double *a, double *b, double *c, double *alpha,
                             double *beta, double *gamma)
{
    *a = cellDim[0];
    *b = cellDim[1];
    *c = cellDim[2];
    *alpha = cellAngles[0];
    *beta = cellAngles[1];
    *gamma = cellAngles[2];
}

void MtzManager::copySymmetryInformationFromManager(MtzPtr toCopy)
{
    double a, b, c, alpha, beta, gamma;
    toCopy->getUnitCell(&a, &b, &c, &alpha, &beta, &gamma);
    this->setUnitCell(a, b, c, alpha, beta, gamma);
    
    std::cout << "Setting unit cell to " << a << ", " << b << ", " << c << ", " << alpha << ", " << beta << ", " << gamma <<std::endl;
    
    int spgnum = toCopy->getLowGroup()->spg_ccp4_num;
    CCP4SPG *spg = ccp4spg_load_by_ccp4_num(spgnum);
    this->setLowGroup(spg);
}

void MtzManager::loadReflections(int partials)
{
    PartialityModel model =
    partials ? PartialityModelScaled : PartialityModelSimple;
    
    loadReflections(model);
}

void MtzManager::loadReflections(PartialityModel model, bool special)
{
    if (filename.length() == 0)
    {
        std::cerr
        << "Cannot load reflections as no filename has been specified."
        << std::endl;
        return;
    }
    
    MTZ *mtz = MtzGet(filename.c_str(), 0);
    
    int fromMtzNum = MtzSpacegroupNumber(mtz);
    
    int spgnum = FileParser::getKey("SPACE_GROUP", fromMtzNum);
    
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
    
    if (low_group == NULL)
    {
        setRejected(true);
        return;
    }
    
    char *hallSymbol = ccp4spg_symbol_Hall(low_group);
    
    space_group spaceGroup = space_group(hallSymbol);
    space_group_type spgType = space_group_type(spaceGroup);
    asu asymmetricUnit = asu(spgType);

    if (mtz == NULL)
        return;
    
    float *refldata = (float *) CCP4::ccp4_utils_malloc(
                                                        (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    memset(refldata, '\0',
           (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    
    float *adata = (float *) CCP4::ccp4_utils_malloc(
                                                     (mtz->ncol_read) * sizeof(float));
    
    MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);
    
    MTZCOL *col_f = MtzColLookup(mtz, "I");
    
    if (col_f == NULL)
    {
        col_f = MtzColLookup(mtz, "IMEAN");
    }
    
    if (col_f == NULL)
    {
        col_f = MtzColLookup(mtz, "FC");
        
        if (col_f == NULL)
        {
            col_f = MtzColLookup(mtz, "FP");
            
            if (col_f == NULL)
                col_f = MtzColLookup(mtz, "F");
            
            if (col_f != NULL)
            {
                std::cout << "Warning: using observed amplitude instead of intensity" << std::endl;
            }
            else
            {
                std::cout << "Warning: could not find any amplitude columns" << std::endl;
            }
        }
        else
        {
            std::cout << "Warning: using calculated amplitude instead of intensity" << std::endl;
        }
    }
    
    if (col_f == NULL)
    {
        std::cout << "Warning: intensity column not labelled I or IMEAN or amplitude" << std::endl;
        exit(1);
    }
    
    MTZCOL *col_sigf = MtzColLookup(mtz, "SIGI");
    
    if (col_sigf == NULL)
    {
        col_sigf = MtzColLookup(mtz, "SIGIMEAN");
    }
    
    if (col_sigf == NULL)
    {
        col_sigf = MtzColLookup(mtz, "SIGF");
        
        if (col_sigf != NULL)
        {
            std::cout << "Warning: using SIGF column" << std::endl;
        }
    }
    
    if (col_sigf == NULL)
    {
        std::cout << "Warning: sigma column not labelled SIGI or SIGIMEAN" << std::endl;
    }
    
    MTZCOL *col_wave = MtzColLookup(mtz, "WAVE");
    MTZCOL *col_partials = MtzColLookup(mtz, "PART");
    MTZCOL *col_phase = MtzColLookup(mtz, "PHIC");
    MTZCOL *col_shiftx = MtzColLookup(mtz, "SHIFTX");
    MTZCOL *col_shifty = MtzColLookup(mtz, "SHIFTY");
    
    int multiplier = MULTIPLIER;
    
    MTZXTAL **xtals = MtzXtals(mtz);
    float cell[6];
    double coefhkl[6];
    ccp4_lrcell(xtals[0], cell);
    MtzHklcoeffs(cell, coefhkl);
    
    setUnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    
    
    
    int num = 0;
    
    for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
    {
        memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));
        
        int h;
        int k;
        int l;
        int offset;
        
        hkls_for_reflection(mtz, adata, &h, &k, &l, &multiplier, &offset);
        
        int reflection_reflection_index = index_for_reflection(h, k, l, false);
        
        float intensity = adata[col_f->source - 1];
        
        float sigma = 1;
        if (col_sigf != NULL)
            sigma = adata[col_sigf->source - 1];
        float wavelength = 1;
        float partiality = 1;
        float phase = 0;
        float shiftX = 0;
        float shiftY = 0;
        
        if (col_wave != NULL && col_partials != NULL)
        {
            wavelength = adata[col_wave->source - 1];
            partiality = adata[col_partials->source - 1];
        }
        
        if (col_shiftx != NULL && col_shifty != NULL)
        {
            shiftX = adata[col_shiftx->source - 1];
            shiftY = adata[col_shifty->source - 1];
        }
        
        if (col_phase != NULL)
        {
            phase = adata[col_phase->source - 1];
        }
        
        if (col_partials == NULL && intensity != intensity)
            continue;
        
        MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
        miller->setData(intensity, sigma, partiality, wavelength);
        miller->setCountingSigma(sigma);
        miller->setFilename(filename);
        miller->setPartialityModel(model);
        miller->setPhase(phase);
        miller->setShift(std::make_pair(shiftX, shiftY));
        miller->matrix = this->matrix;
        
        Reflection *prevReflection;
        
        this->findReflectionWithId(reflection_reflection_index, &prevReflection);
        
        if (prevReflection != NULL)
        {
            /** Exclude unobserved reflections by checking for nan */
            if (adata[col_f->source - 1] == adata[col_f->source - 1])
            {
                prevReflection->addMiller(miller); // TODO
                miller->setParent(prevReflection);
            }
            
            // reflection is a repeat so set flag.
        }
        
        if (prevReflection == NULL)
        {
            Reflection *newReflection = new Reflection();
            reflections.push_back(newReflection);
            newReflection->setUnitCell(cell);
            newReflection->setSpaceGroup(low_group, spgType, asymmetricUnit);
            
            reflections[num]->addMiller(miller);
            miller->setParent(reflections[num]);
            reflections[num]->calculateResolution(this);
            
            this->sortLastReflection();
            
            num++;
        }
    }
    
    free(refldata);
    free(adata);
    
    //	insertionSortReflections();
    
    std::ostringstream log;
    
    log << "Loaded " << mtz->nref_filein << " reflections (" << accepted()
    << " accepted)." << std::endl;
    
    Logger::mainLogger->addStream(&log);
    
    MtzFree(mtz);
}

void MtzManager::getRefReflections(vector<Reflection *> *refPointer,
                               vector<Reflection *> *matchPointer)
{
    if (lastReference != referenceManager)
    {
        refReflections.clear();
        matchReflections.clear();
        
        for (int i = 0; i < reflectionCount(); i++)
        {
            int reflId = reflection(i)->getReflId();
            
            Reflection *refReflection = NULL;
            referenceManager->findReflectionWithId(reflId, &refReflection);
            
            if (refReflection != NULL)
            {
                matchReflections.push_back(reflection(i));
                refReflections.push_back(refReflection);
            }
        }
        
        lastReference = referenceManager;
    }
    
    refPointer->reserve(refReflections.size());
    refPointer->insert(refPointer->begin(), refReflections.begin(),
                       refReflections.end());
    
    matchPointer->reserve(matchReflections.size());
    matchPointer->insert(matchPointer->begin(), matchReflections.begin(),
                         matchReflections.end());
    
}

void MtzManager::setReference(MtzManager *reference)
{
    if (reference != NULL)
        Logger::mainLogger->addString("Setting reference to " + reference->getFilename());
    MtzManager::referenceManager = reference;
}

void MtzManager::setFilename(std::string name)
{
    filename = name;
}

std::string MtzManager::getFilename(void)
{
    return filename;
}

void MtzManager::setSpaceGroup(int spgnum)
{
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
}

int MtzManager::findReflectionWithId(long unsigned int refl_id, Reflection **reflection, bool insertionPoint)
{
    if (reflectionCount() == 0)
    {
        *reflection = NULL;
        return 0;
    }
    
    int lower = 0;
    int higher = reflectionCount() - 1;
    int new_bound = (higher + lower) / 2;
    
    if ((refl_id < this->reflection(lower)->getReflId())
        || (refl_id > this->reflection(higher)->getReflId()))
    {
        if (insertionPoint)
        {
            if (refl_id < this->reflection(lower)->getReflId())
                return 0;
            
            if (refl_id > this->reflection(higher)->getReflId())
                return reflectionCount();
        }
        
        *reflection = NULL;
        return -1;
    }
    
    while (this->reflection(new_bound)->getReflId() != refl_id)
    {
        if (new_bound == higher || new_bound == lower)
        {
            if (this->reflection(higher)->getReflId() == refl_id)
            {
                (*reflection) = this->reflection(higher);
                return -1;
            }
            
            if (this->reflection(lower)->getReflId() == refl_id)
            {
                (*reflection) = this->reflection(lower);
                return -1;
            }
            
            if (insertionPoint)
            {
                int lowest = 0;
                
                int start = lower - 2 >= 0 ? lower - 2 : 0;
                
                for (int i = start; i < higher + 2 && i < reflectionCount(); i++)
                {
                    if (this->reflection(i)->getReflId() < refl_id)
                        lowest = i;
                }
                
                return lowest + 1;
            }
            
            *reflection = NULL;
            return -1;
        }
        
        if (this->reflection(new_bound)->getReflId() > refl_id)
        {
            higher = new_bound;
        }
        else if (this->reflection(new_bound)->getReflId() < refl_id)
        {
            lower = new_bound;
        }
        
        new_bound = (higher + lower) / 2;
    }
    
    (*reflection) = this->reflection(new_bound);
    
    return -1;
}

void MtzManager::findCommonReflections(MtzManager *other,
                                       vector<Reflection *> &reflectionVector1, vector<Reflection *> &reflectionVector2,
                                       int *num, bool force)
{
    if (other == previousReference && previousAmbiguity == activeAmbiguity && !force && false)
    {
        reflectionVector1.reserve(matchReflections.size());
        reflectionVector2.reserve(refReflections.size());
        
        reflectionVector1.insert(reflectionVector1.begin(), matchReflections.begin(), matchReflections.end());
        reflectionVector2.insert(reflectionVector2.begin(), refReflections.begin(), refReflections.end());
        
        if (num != NULL)
            *num = (int)reflectionVector1.size();
        
        return;
    }

    matchReflections.clear();
    refReflections.clear();
    
    previousReference = other;
    previousAmbiguity = activeAmbiguity;
    
    for (int i = 0; i < reflectionCount(); i++)
    {        
        int refl_id = reflection(i)->getReflId();
        
        Reflection *otherReflection = NULL;
        
        other->findReflectionWithId(refl_id, &otherReflection);
        
        if (otherReflection != NULL && otherReflection->millerCount() > 0)
        {
            reflectionVector1.push_back(reflection(i));
            matchReflections.push_back(reflection(i));
            reflectionVector2.push_back(otherReflection);
            refReflections.push_back(otherReflection);
        }
    }
    
    if (num != NULL)
    {
        *num = (int)reflectionVector1.size();
    }
}

void MtzManager::applyScaleFactorsForBins(int binCount)
{
    double highRes = FileParser::getKey("HIGH_RESOLUTION", 0.0);
    double lowRes = FileParser::getKey("LOW_RESOLUTION", 0.0);
    
    if (highRes == 0)
    {
        highRes = 1 / this->maxResolution();
    }
    
    vector<double> bins;
    StatisticsManager::generateResolutionBins(lowRes, highRes, binCount, &bins);
    
    for (int shell = 0; shell < bins.size() - 1; shell++)
    {
        double low = bins[shell];
        double high = bins[shell + 1];
        
        vector<Reflection *> refReflections, imgReflections;
        
        this->findCommonReflections(referenceManager, imgReflections, refReflections,
                                    NULL);
        
        double weights = 0;
        double refMean = 0;
        double imgMean = 0;
        int count = 0;
        
        for (int i = 0; i < imgReflections.size(); i++)
        {
            if (imgReflections[i]->betweenResolutions(low, high))
            {
                double refIntensity = refReflections[i]->meanIntensity();
                double imgIntensity = imgReflections[i]->meanIntensity();
                
                if (refIntensity != refIntensity || imgIntensity != imgIntensity)
                    continue;
                
                weights++;
                
                refMean += refIntensity;
                imgMean += imgIntensity;
                count++;
                
            }
        }
        
        refMean /= weights;
        imgMean /= weights;
        
        std::cout << refMean << ", " << imgMean << ", " << weights << ", " << low << ", " << high << std::endl;
        
        double ratio = refMean / imgMean;
        
        applyScaleFactor(ratio, low, high);
    }
}

void MtzManager::bFactorAndScale(double *scale, double *bFactor, double exponent,
                                 vector<std::pair<double, double> > *dataPoints)
{
    clearScaleFactor();
    MtzManager *reference = MtzManager::getReferenceManager();
    
    double grad = this->gradientAgainstManager(reference);
    this->applyScaleFactor(grad);
    
    vector<Reflection *> refReflections, imageReflections;
    this->findCommonReflections(reference, imageReflections, refReflections, NULL);
    
    vector<double> smoothBins;
    StatisticsManager::generateResolutionBins(50, 1.8, 24, &smoothBins);
    
    if (dataPoints == NULL)
    {
        vector<std::pair<double, double> > data = vector<std::pair<double, double> >();
        dataPoints = &data;
    }
    
    vector<boost::tuple<double, double, double> > pointsToFit;
    
    for (int shell = 0; shell < smoothBins.size() - 3; shell++)
    {
        double low = smoothBins[shell];
        double high = smoothBins[shell + 3];
        
        vector<Reflection *> refReflections, imgReflections;
        
        this->findCommonReflections(referenceManager, imgReflections, refReflections);
        
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
        
        double resolution = StatisticsManager::midPointBetweenResolutions(
                                                                          smoothBins[shell], smoothBins[shell + 3]);
        
        double res_squared = pow(1 / resolution, 2);
        
        double four_d_squared = 4 * pow(resolution, 2);
        
        double right_exp = pow(1 / four_d_squared, exponent);
        
        double four_d_to_exp = pow(2, exponent) * right_exp;
        
        double intensityRatio = (refMean / imgMean);
        double logIntensityRatio = log(intensityRatio);
        
        if (logIntensityRatio != logIntensityRatio)
            continue;
        
        double weight = 1;
        weight /= res_squared;
        
        boost::tuple<double, double, double> point = boost::make_tuple(four_d_to_exp,
                                                              logIntensityRatio, weight);
        
        std::pair<double, double> showPoint = std::make_pair(res_squared,
                                                        1 / intensityRatio);
        
        dataPoints->push_back(showPoint);
        pointsToFit.push_back(point);
    }
    
    double gradient = 0;
    double intercept = 0;
    
    regression_line(pointsToFit, intercept, gradient);
    
    double k = exp(intercept);
    
    double b = pow(fabs(gradient), 1 / (double) exponent);
    
    if (gradient < 0)
        b = -b;
    
    if (b != b)
        b = 0;
    
    *scale = k;
    *bFactor = b;
}

double MtzManager::minimizeRFactor(MtzManager *otherManager)
{
    return minimizeGradient(otherManager, false);
}

double MtzManager::minimizeGradient(MtzManager *otherManager, bool leastSquares)
{
    double resolution = 0.001;
    double step = 0.1;
    double gradient = gradientAgainstManager(otherManager);
    
    vector<Reflection *> reflections1;
    vector<Reflection *> reflections2;
    vector<double> ints1;
    vector<double> ints2;
    vector<double> weights;
    int num = 0;
    
    MtzManager::findCommonReflections(otherManager, reflections1, reflections2, &num);
    
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < reflections1[i]->millerCount(); j++)
        {
            if (!reflections1[i]->miller(j)->accepted())
                continue;
            
            double mean1 = reflections1[i]->miller(j)->intensity();
            double mean2 = reflections2[i]->meanIntensity();
            double weight = reflections2[i]->meanPartiality();
            
            if (mean1 != mean1 || mean2 != mean2 || weight != weight)
                continue;
            
            //		std::cout << mean1 << "\t" << mean2 << std::endl;
            
            ints1.push_back(mean1);
            ints2.push_back(mean2);
            weights.push_back(weight);
        }
    }
    
    double best = r_factor_between_vectors(&ints1, &ints2, &weights, gradient);
    ;
    
    while (step > resolution)
    {
        double gradUp = gradient + step;
        double gradDown = gradient - step;
        
        double upScore = r_factor_between_vectors(&ints1, &ints2, &weights,
                                                  gradUp);
        double downScore = r_factor_between_vectors(&ints1, &ints2, &weights,
                                                    gradDown);
        
        //	std::cout << upScore << "\t" << best << "\t" << downScore << std::endl;
        
        if (upScore < best)
        {
            gradient = gradUp;
            best = upScore;
        }
        else if (downScore < best)
        {
            gradient = gradDown;
            best = downScore;
        }
        else
        {
            step /= 2;
        }
    }
    
    return gradient;
}

double MtzManager::gradientAgainstManager(MtzManager *otherManager,
                                          bool withCutoff, double lowRes, double highRes)
{
    vector<Reflection *> reflections1;
    vector<Reflection *> reflections2;
    
    int count = 0;
    int num = 0;
    
    this->findCommonReflections(otherManager, reflections1, reflections2, &num);
    
    if (num <= 1)
        return 1;
    
    double x_squared = 0;
    double x_y = 0;
    
    for (int i = 0; i < num; i++)
    {
        if (reflections1[i]->acceptedCount() == 0 && withCutoff)
            continue;
        
        if (!reflections1[i]->betweenResolutions(lowRes, highRes))
            continue;
        
        double int1 = reflections1[i]->meanIntensity(withCutoff);
        double int2 = reflections2[i]->meanIntensity();
        double weight = reflections1[i]->meanPartiality();
        
        if ((int1 != int1) || (int2 != int2) || (weight != weight))
            continue;
        
        x_squared += int1 * int1 * weight;
        x_y += int1 * int2 * weight;
        
        count++;
    }
    
    double grad = x_y / x_squared;
    
    if (grad < 0)
        grad = -1;
    
    return grad;
}

void MtzManager::clearScaleFactor()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflections[i]->miller(j)->setScale(1);
            reflections[i]->miller(j)->setBFactor(0);
        }
    }
}

void MtzManager::makeScalesPermanent()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            reflection(i)->miller(j)->makeScalesPermanent();
        }
    }
}

void MtzManager::applyBFactor(double bFactor)
{
 //   if (bFactor < 0)
 //       bFactor = 0 - bFactor;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflection(i)->miller(j)->setBFactor(bFactor);
        }
    }
}

void MtzManager::applyScaleFactor(double scaleFactor,
                                  double lowRes, double highRes, bool absolute)
{
    if (scaleFactor == scaleFactor)
    {
        if (absolute)
            scale = scaleFactor;
        else
            scale *= scaleFactor;
    }
    
    logged << "Applying scale factor " << scaleFactor << " now ";
    
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            if (reflection(i)->betweenResolutions(lowRes, highRes))
            {
                if (absolute)
                    reflection(i)->miller(j)->setScale(scale);
                else
                    reflection(i)->miller(j)->applyScaleFactor(scaleFactor);
                
                if (i == 0 && j == 0)
                    logged << reflection(i)->miller(j)->getScale() << std::endl;
            }
        }
    }
    
    sendLog(LogLevelDebug);

}

double MtzManager::averageIntensity(void)
{
    double total_intensity = 0;
    double total = 0;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        double intensity = reflections[i]->meanIntensity();
        if (intensity == intensity)
        {
            double weight = reflections[i]->meanPartiality();
            total_intensity += intensity * weight;
            total += weight;
        }
    }
    
    total_intensity /= total;
    return total_intensity;
}

void MtzManager::applyPolarisation(void)
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflections[i]->miller(j)->applyPolarisation(1.459);
        }
    }
}

void MtzManager::writeToFile(std::string newFilename, bool announce, bool shifts, bool includeAmbiguity)
{
    int columns = 7;
    if (includeAmbiguity)
        flipToActiveAmbiguity();
    
    if (shifts) columns += 2;
    
    float cell[6], wavelength;
    float *fdata = new float[columns];
    
    /* variables for symmetry */
    CCP4SPG *mtzspg = low_group;
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[9];
    
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
    cell[0] = cellDim[0];
    cell[1] = cellDim[1];
    cell[2] = cellDim[2];
    cell[3] = cellAngles[0];
    cell[4] = cellAngles[1];
    cell[5] = cellAngles[2];
    wavelength = this->getWavelength();
    
    mtzout = MtzMalloc(0, 0);
    ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
    mtzout->refs_in_memory = 0;
    mtzout->fileout = MtzOpenForWrite(newFilename.c_str());
    
    // then add symm headers...
    for (int i = 0; i < mtzspg->nsymop; ++i)
        CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
    strncpy(ltypex, mtzspg->symbol_old, 1);
    ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
                mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);
    
    // then add xtals, datasets, cols
    xtal = MtzAddXtal(mtzout, "XFEL crystal", "XFEL project", cell);
    set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
    colout[0] = MtzAddColumn(mtzout, set, "H", "H");
    colout[1] = MtzAddColumn(mtzout, set, "K", "H");
    colout[2] = MtzAddColumn(mtzout, set, "L", "H");
    colout[3] = MtzAddColumn(mtzout, set, "I", "J");
    colout[4] = MtzAddColumn(mtzout, set, "SIGI", "Q");
    colout[5] = MtzAddColumn(mtzout, set, "PART", "R");
    colout[6] = MtzAddColumn(mtzout, set, "WAVE", "R");
    
    if (shifts)
    {
        colout[7] = MtzAddColumn(mtzout, set, "SHIFTX", "R");
        colout[8] = MtzAddColumn(mtzout, set, "SHIFTY", "R");
    }
    
    int num = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            if (!reflection(i)->miller(j))
                std::cout << "!miller(j) in mtz manager" << std::endl;
            
            if (reflections[i]->miller(j)->isRejected())
                continue;
            
            double intensity = reflections[i]->miller(j)->getRawIntensity();
            double sigma = reflections[i]->miller(j)->getSigma();
            double partiality = reflections[i]->miller(j)->getPartiality();
            double bFactor = reflections[i]->miller(j)->getBFactorScale();
            
            if (intensity != intensity)
            {
                continue;
            }
            
            num++;
            
            int h = reflections[i]->miller(j)->getH();
            int k = reflections[i]->miller(j)->getK();
            int l = reflections[i]->miller(j)->getL();
            int _h, _k, _l;
            ccp4spg_put_in_asu(low_group, h, k, l, &_h, &_k, &_l);
            
            // set fdata
            fdata[0] = h;
            fdata[1] = k;
            fdata[2] = l;
            fdata[3] = intensity / bFactor;
            fdata[4] = sigma;
            fdata[5] = partiality;
            fdata[6] = reflections[i]->miller(j)->getWavelength();
            
            if (shifts)
            {
                fdata[7] = reflections[i]->miller(j)->getShift().first;
                fdata[8] = reflections[i]->miller(j)->getShift().second;
            }
            
            ccp4_lwrefl(mtzout, fdata, colout, columns, num);
        }
    }
    
    MtzPut(mtzout, " ");
    MtzFree(mtzout);
    
    LogLevel shouldAnnounce = announce ? LogLevelNormal : LogLevelDebug;
    
    std::ostringstream logged;
    logged << "Written to file " << newFilename << std::endl;
    Logger::mainLogger->addStream(&logged, shouldAnnounce);
    
    if (includeAmbiguity)
        resetFlip();
    
    delete [] fdata;
}


MtzManager::~MtzManager(void)
{
    for (int i = 0; i < reflections.size(); i++)
    {
        delete reflections[i];
    }
    
    reflections.clear();
    vector<Reflection *>().swap(reflections);
    
    if (low_group != NULL)
        ccp4spg_free(&low_group);
    
}

void MtzManager::description(void)
{
    logged << "Filename: " << filename << std::endl;
    logged << "Number of reflections: " << reflectionCount() << std::endl;
    logged << "Number of accepted Millers: " << accepted() << std::endl;
    logged << "Average intensity: " << this->averageIntensity() << std::endl;
    
    sendLog();
}

void MtzManager::setSigmaToUnity()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflections[i]->miller(j)->setSigma(1);
            
            // TAKE NEXT LINE OUT
        }
    }
}

void MtzManager::setPartialityToUnity()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflections[i]->miller(j)->setPartiality(1);
            
            // TAKE NEXT LINE OUT
        }
    }
}

int MtzManager::accepted(void)
{
    int acceptedCount = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            if (reflection(i)->miller(j)->accepted())
                acceptedCount++;
        }
    }
    
    return acceptedCount;
}


void MtzManager::writeToDat()
{
    std::string name = filename;
    int lastindex = (int)name.find_last_of(".");
    std::string rootName = name.substr(0, lastindex);
    std::string datName = rootName + ".dat";
    
    std::ofstream datOutput;
    datOutput.open(datName);
    
    datOutput << this->matrix->description() << std::endl;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        Reflection *aReflection = reflection(i);
        for (int j = 0; j < aReflection->millerCount(); j++)
        {
            MillerPtr miller = aReflection->miller(j);
            double shiftX = miller->getShift().first;
            double shiftY = miller->getShift().second;
            double lastX = miller->getLastX();
            double lastY = miller->getLastY();
            
            double combinedX = shiftX + lastX;
            double combinedY = shiftY + lastY;
            
            datOutput << miller->getH() << "\t" << miller->getK() << "\t" << miller->getL();
            
            datOutput << "\t" << miller->getWavelength();
            datOutput << "\t" << miller->getRawIntensity();
            datOutput << "\t" << miller->getCountingSigma();
            datOutput << "\t" <<  combinedX << "\t" << combinedY << "\t"
            << miller->getPartiality() << "\t" << miller->getBFactorScale() << "\t" << miller->getResolution() << std::endl;
        }
    }
    
    logged << "Wrote " << datName << " file" << std::endl;
    
    datOutput.close();
}

void MtzManager::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

int MtzManager::rejectOverlaps()
{
    int count = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        count += reflection(i)->checkOverlaps();
    }
    
    return count;
}

int MtzManager::removeStrongSpots(std::vector<SpotPtr> *spots)
{
    int before = (int)spots->size();
    
    int count = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        Reflection *ref = reflection(i);
        
        count += ref->checkSpotOverlaps(spots);
    }
    
    int after = (int)spots->size();
    
    logged << "Before spot removal: " << before << " - after: " << after << std::endl;
    sendLog(LogLevelDetailed);
    
    return count;
}

int MtzManager::ambiguityCount()
{
    int count = reflection(0)->ambiguityCount();
    
    return count;
}

void MtzManager::incrementActiveAmbiguity()
{
    int count = ambiguityCount();
    
    activeAmbiguity++;
    activeAmbiguity = (activeAmbiguity % count);
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setActiveAmbiguity(activeAmbiguity);
    }
}

void MtzManager::flipToActiveAmbiguity()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setFlip(activeAmbiguity);
    }
    
    this->insertionSortReflections();
}


void MtzManager::resetFlip()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->resetFlip();
    }
    
    this->insertionSortReflections();
}

void MtzManager::setActiveAmbiguity(int newAmbiguity)
{
    activeAmbiguity = newAmbiguity;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setActiveAmbiguity(newAmbiguity);
    }
}

double MtzManager::maxResolution()
{
    double maxResolution = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (reflection(i)->getResolution() > maxResolution)
        {
            maxResolution = reflection(i)->getResolution();
        }
    }
    
    return maxResolution;
}

void MtzManager::cutToResolutionWithSigma(double acceptableSigma)
{
    double maxRes = maxResolution();
    vector<double> resolutions;
    StatisticsManager::generateResolutionBins(0, maxRes, 10, &resolutions);

    double cutoffResolution = maxResolution();
    
    for (int i = (int)resolutions.size() - 1; i > 0 && cutoffResolution == 0; i--)
    {
        double isigiSum = 0;
        double weights = 0;
        
        for (int i = 0; i < reflectionCount(); i++)
        {
            double highRes = resolutions[i];
            double lowRes = resolutions[i - 1];
            
            if (!reflection(i)->betweenResolutions(lowRes, highRes))
                continue;

            for (int j = 0; j < reflection(i)->millerCount(); j++)
            {
                MillerPtr miller = reflection(i)->miller(j);
                
                double weight = 1;
                double isigi = miller->getRawIntensity() / miller->getCountingSigma();
                
                isigiSum += isigi;
                weights += weight;
            }
        }
        
        double isigi = isigiSum / weights;
        
        if (isigi > acceptableSigma)
        {
            cutoffResolution = highRes;
        }
    }
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (reflection(i)->getResolution() > cutoffResolution)
        {
            this->removeReflection(i);
            i--;
        }
    }
    
    logged << filename << ": rejected reflections over " << 1 / cutoffResolution << " Å." << std::endl;
    sendLog();
}

void MtzManager::setAdditionalWeight(double weight)
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setAdditionalWeight(weight);
    }
}

void MtzManager::denormaliseMillers()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            reflection(i)->miller(j)->denormalise();
        }
    }
}

bool MtzManager::checkUnitCell(double trueA, double trueB, double trueC, double tolerance)
{
    double *cellDims = new double[3];
    matrix->unitCellLengths(&cellDims);
    
    if (cellDims[0] < trueA - tolerance || cellDims[0] > trueA + tolerance)
        return false;
    if (cellDims[1] < trueB - tolerance || cellDims[1] > trueB + tolerance)
        return false;
    if (cellDims[2] < trueC - tolerance || cellDims[2] > trueC + tolerance)
        return false;
    
    return true;
}

void MtzManager::incrementFailedCount()
{
    failedCount++;
}

void MtzManager::resetFailedCount()
{
    failedCount = 0;
}

#define PARAM_HROT 0
#define PARAM_KROT 1
#define PARAM_MOS 2
#define PARAM_SPOT_SIZE 3
#define PARAM_WAVELENGTH 4
#define PARAM_BANDWIDTH 5
#define PARAM_EXPONENT 6
#define PARAM_B_FACTOR 7
#define PARAM_SCALE_FACTOR 8
#define PARAM_UNIT_CELL_A 9
#define PARAM_UNIT_CELL_B 10
#define PARAM_UNIT_CELL_C 11
#define PARAM_NUM 12

std::string MtzManager::getParamLine()
{
    std::ostringstream line;
    line << "params ";
    
    line << mosaicity << " " << spotSize << " " << wavelength << " ";
    line << bFactor << " " << scale << " ";
    line << cellDim[0] << " " << cellDim[1] << " " << cellDim[2];
    
    return line.str();
}

void MtzManager::setParamLine(std::string line)
{
    std::vector<std::string> components = FileReader::split(line, ' ');
    
    if (components.size() < 8)
    {
        return;
    }
    
    mosaicity = atof(components[1].c_str());
    spotSize = atof(components[2].c_str());
    wavelength = atof(components[3].c_str());
    bFactor = atof(components[4].c_str());
    scale = atof(components[5].c_str());
    cellDim[0] = atof(components[6].c_str());
    cellDim[1] = atof(components[7].c_str());
    cellDim[2] = atof(components[8].c_str());
    
    finalised = true;
    setInitialValues = true;
}