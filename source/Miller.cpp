/*
 * Miller.cpp
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#include "Panel.h"
#include "Miller.h"

#include "definitions.h"
#include "Vector.h"
#include <cmath>
#include "parameters.h"
#include <cmath>
#include "Shoebox.h"
#include <memory>
#include <cctbx/miller/asu.h>
#include <cctbx/miller.h>
#include "FileParser.h"



using cctbx::sgtbx::reciprocal_space::asu;

asu Miller::p1_asu = asu();
bool Miller::initialised_p1 = false;
space_group Miller::p1_spg = space_group();

double Miller::integrate_sphere_gaussian(double p, double q)
{
    double mean = 0.5;
    double sigma = 0.25;
    double exponent = 1.5;
    
    double partiality_add = super_gaussian(q, mean, sigma, exponent);
    double partiality_remove = super_gaussian(p, mean, sigma, exponent);
    
    double proportion = (partiality_add + partiality_remove) / 2;
    proportion *= (q - p);
    
    
    return proportion;
}

double Miller::integrate_sphere_uniform(double p, double q)
{
    double partiality_add = 4 * q * (1 - q);
    double partiality_remove = 4 * p * (1 - p);
    
    double proportion = (partiality_add + partiality_remove) / 2;
    proportion *= (q - p);
    
    return proportion;
}

double Miller::integrate_sphere(double p, double q, double radius, double sphere_volume, double circle_surface_area)
{
    if (rlpModel == RlpModelGaussian)
        return integrate_sphere_gaussian(p, q);
    if (rlpModel == RlpModelUniform)
        return integrate_sphere_uniform(p, q);
    
    else return integrate_sphere_uniform(p, q);
}

double Miller::superGaussian(double bandwidth, double mean,
                            double sigma, double exponent)
{
    if (mtzParent == NULL || mtzParent->getLastExponent() == 0)
    {
        return super_gaussian(bandwidth, mean, sigma, exponent);
    }
    else
    {
        if (bandwidth != bandwidth || mean != mean)
            return 0;
        
        double standardisedX = fabs((bandwidth - mean) / sigma);
        
        if (standardisedX > MAX_SUPER_GAUSSIAN)
            return 0;
        
        const double step = SUPER_GAUSSIAN_STEP;
        
        int lookupInt = standardisedX / step;
        return mtzParent->superGaussianTable[lookupInt];
    }
}

double Miller::integrate_beam_slice(double pBandwidth, double qBandwidth, double mean,
                                   double sigma, double exponent)
{
    double pValue = superGaussian(pBandwidth, mean, sigma, exponent);
    double qValue = superGaussian(qBandwidth, mean, sigma, exponent);
    
    double pX = (pBandwidth - mean) / sigma;
    double qX = (qBandwidth - mean) / sigma;
    
    double width = fabs(qX - pX);
    
    double area = (pValue + qValue) / 2 * width;
    
    return area;
}

double integrate_beam(double pBandwidth, double qBandwidth, double mean,
                     double sigma)
{
    double normal_add = cdf(pBandwidth, mean, sigma);
    double normal_remove = cdf(qBandwidth, mean, sigma);
    
    double normal_slice = (normal_add - normal_remove);
    
    return normal_slice;
}

double Miller::sliced_integral(double low_wavelength, double high_wavelength,
                              double spot_size_radius, double maxP, double maxQ, double mean, double sigma,
                              double exponent)
{
    if (resolution() < 1 / trickyRes) slices = maxSlices;
    
    double bandwidth_span = high_wavelength - low_wavelength;
    
    double fraction_total = 1;
    double currentP = 0;
    double currentQ = 1 / slices;
    
    double total_integral = 0;
    
    double total_normal = 0;
    double total_sphere = 0;
    
    double sphere_volume = (4 * M_PI / 3) * pow(spot_size_radius, 3);
    double circle_surface_area = pow(spot_size_radius, 2) * M_PI;
    
    for (int i = 0; i < slices; i++)
    {
        double pBandwidth = high_wavelength - bandwidth_span * currentP;
        double qBandwidth = high_wavelength - bandwidth_span * currentQ;
        
        double normalSlice = integrate_beam_slice(pBandwidth, qBandwidth, mean,
                                                 sigma, exponent);
        
        double sphereSlice = 0;
        
        if (normalSlice > 0)
        {
            sphereSlice = integrate_sphere(currentP, currentQ,
                                              spot_size_radius, sphere_volume, circle_surface_area);
        }
        
        total_integral += normalSlice * sphereSlice;
        total_normal += normalSlice;
        total_sphere += sphereSlice;
        
        currentP = currentQ;
        currentQ += fraction_total / slices;
    }
    
    return total_integral;
}

Miller::Miller(MtzManager *parent, int _h, int _k, int _l)
{
    h = _h;
    k = _k;
    l = _l;
    lastX = 0;
    lastY = 0;
    rawIntensity = 0;
    gainScale = 1;
    sigma = 1;
    wavelength = 0;
    lastNormPartiality = 0;
    partiality = -1;
    model = PartialityModelNone;
    filename = "";
    countingSigma = 0;
    latestHRot = 0;
    latestKRot = 0;
    normalised = FileParser::getKey("NORMALISE_PARTIALITIES", true);
    correctingPolarisation = FileParser::getKey("POLARISATION_CORRECTION", false);
    polarisationFactor = FileParser::getKey("POLARISATION_FACTOR", 0.0);
    polarisationCorrection = 0;
    lastWavelength = 0;
    rejectedReasons["merge"] = false;
    scale = 1;
    bFactor = 0;
    bFactorScale = 0;
    resol = 0;
    shift = std::make_pair(0, 0);
    shoebox = ShoeboxPtr();
    selfPtr = MillerPtr();
    fakeFriedel = -1;
    rejected = false;
    calculatedRejected = true;
    denormaliseFactor = 1;
    excluded = false;
    lastRlpSize = 0;
    lastMosaicity = 0;
    lastVolume = 0;
    lastSurfaceArea = 0;
    slices = FileParser::getKey("PARTIALITY_SLICES", 8);
    trickyRes = FileParser::getKey("CAREFUL_RESOLUTION", 8.0);
    maxSlices = FileParser::getKey("MAX_SLICES", 100);
    
    partialCutoff = FileParser::getKey("PARTIALITY_CUTOFF",
                                       PARTIAL_CUTOFF);
    
    int rlpInt = FileParser::getKey("RLP_MODEL", 0);
    rlpModel = (RlpModel)rlpInt;
    
    mtzParent = parent;
    parentReflection = NULL;
    matrix = MatrixPtr();
    flipMatrix = MatrixPtr(new Matrix());
}

vec Miller::hklVector(bool shouldFlip)
{
    vec newVec = new_vector(h, k, l);
    
    if (shouldFlip)
    {
        flipMatrix->multiplyVector(&newVec);
    }
    
    return newVec;
}

int Miller::getH()
{
    return hklVector().h;
}

int Miller::getK()
{
    return hklVector().k;

}

int Miller::getL()
{
    return hklVector().l;
}

void Miller::setFlipMatrix(MatrixPtr flipMat)
{
    flipMatrix = flipMat;
}

void Miller::setParent(Reflection *reflection)
{
    parentReflection = reflection;
}

void Miller::setPartialityModel(PartialityModel _model)
{
    model = _model;
}

void Miller::setData(double _intensity, double _sigma, double _partiality,
                     double _wavelength)
{
    rawIntensity = _intensity;
    sigma = _sigma;
    partiality = _partiality;
    wavelength = _wavelength;
}

void Miller::printHkl(void)
{
    std::cout << "h k l " << h << " " << k << " " << l << std::endl;
}

bool Miller::accepted(void)
{
    if (this->model == PartialityModelNone)
        return true;

    if (this->partiality < partialCutoff)
    {
        return false;
    }
    
    if (this->isRejected())
        return false;
    
    if (this->rawIntensity != this->rawIntensity)
        return false;
    
    double sigma = this->getSigma();
    
    if (sigma == 0)
    {
        return false;
    }
    
    return true;
}

bool Miller::isRejected()
{
    if (calculatedRejected)
        return rejected;
    
    rejected = false;
    
    for (std::map<std::string, bool>::iterator it = rejectedReasons.begin();
         it != rejectedReasons.end(); ++it)
    {
        if (rejectedReasons[it->first])
            rejected = true;
    }
    
    calculatedRejected = true;
    
    return rejected;
}

double Miller::getBFactorScale()
{
    if (bFactorScale != 0)
        return bFactorScale;
    
    double resn = getResolution();
    
    double sinThetaOverLambda = resn / 2;// 1 / (2 / resn);
    
    double factor = exp(- bFactor * pow(sinThetaOverLambda, 2));
    
    bFactorScale = factor;
    
    return factor;
}

double Miller::scaleForScaleAndBFactor(double scaleFactor, double bFactor,
                                      double resn, double exponent_exponent)
{
    double four_d_squared = 4 * pow(resn, 2);
    
    double right_exp = pow(1 / four_d_squared, exponent_exponent);
    
    double four_d_to_exp = pow(2, exponent_exponent) * right_exp;
    
    double exponent = pow(bFactor, exponent_exponent) * four_d_to_exp;
    
    double scale = scaleFactor * exp(pow(exponent, exponent_exponent));
    
    return scale;
}

double Miller::intensity(bool withCutoff)
{
    if ((this->accepted() && withCutoff) || !withCutoff)
    {
        double modifier = scale * gainScale;
        modifier /= getBFactorScale();
        modifier /= correctingPolarisation ? getPolarisationCorrection() : 1;
        
        if (model == PartialityModelScaled)
        {
            modifier /= partiality * denormaliseFactor;
        }
        
        if (partiality == 0)
            modifier = 0;
        
        return modifier * rawIntensity;
    }
    else
    {
        std::cout << "Warning, allowing an unaccepted intensity!" << std::endl;
        return 0;
    }
}

void Miller::makeScalesPermanent()
{
    double scale = scaleForScaleAndBFactor(this->scale, this->bFactor,
                                          1 / this->resol);
    
    this->rawIntensity *= scale;
    this->sigma *= scale;
    
    this->scale = 1;
    bFactor = 0;
}

double Miller::getSigma(void)
{
    // bigger sigma - larger error
    
    double correction = gainScale * scale;
    /*
    double bFactorScale = getBFactorScale();
    correction /= bFactorScale;
    
    if (bFactorScale == 0) correction = 0;*/
    
    return sigma * correction;
}

double Miller::getPartiality()
{
    return partiality;
}

double Miller::getWavelength()
{
    return wavelength;
}

vec Miller::getTransformedHKL(MatrixPtr matrix)
{
    vec hkl = new_vector(h, k, l);
    
    matrix->multiplyVector(&hkl);
    
    return hkl;
}

vec Miller::getTransformedHKL(double hRot, double kRot)
{
    MatrixPtr newMatrix = MatrixPtr();
    rotateMatrixHKL(hRot, kRot, 0, matrix, &newMatrix);
    
    vec hkl = getTransformedHKL(newMatrix);
    return hkl;
}

double Miller::getWavelength(double hRot, double kRot)
{
    vec hkl = getTransformedHKL(hRot, kRot);
    
    wavelength = getEwaldSphereNoMatrix(hkl);
    
    return wavelength;
}

double Miller::getWavelength(MatrixPtr transformedMatrix)
{
    vec hkl = getTransformedHKL(transformedMatrix);
    wavelength = getEwaldSphereNoMatrix(hkl);
    
    return wavelength;
}

double Miller::getEwaldWeight(double hRot, double kRot, bool isH)
{
    double weight = 0;
    vec hkl = getTransformedHKL(hRot, kRot);
    weight = getEwaldWeightForAxis(hkl, isH);
    
    return weight;
}

double Miller::getWeight(bool cutoff, WeightType weighting)
{
    if (!this->accepted() && cutoff)
    {
        return 0;
    }
    
    double weight = 1;
    double partiality = this->getPartiality();
    
    if (weighting == WeightTypePartiality
        || weighting == WeightTypePartialityCorrelation
        || weighting == WeightTypePartialitySigma)
    {
        weight = partiality;
    }
    
    if (weighting == WeightTypePartialitySigma)
    {
        double sigma = this->getSigma();
        
        if (sigma > 0)
            weight /= sigma;
    }
    
    if (weighting == WeightTypeISigI)
    {
        double intensity = this->getRawIntensity();
        double sigma = this->getSigma();
        
        weight = intensity / sigma;
    }

    return weight;
}

void Miller::rotateMatrixABC(double aRot, double bRot, double cRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix)
{
    (*newMatrix) = oldMatrix->copy();
    
    (*newMatrix)->rotateABC(oldMatrix, aRot, bRot, cRot);
}

void Miller::rotateMatrixHKL(double hRot, double kRot, double lRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix)
{
    (*newMatrix) = oldMatrix->copy();
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    double lRad = lRot * M_PI / 180;
    
    (*newMatrix)->rotate(hRad, kRad, lRad);
}

double Miller::expectedRadius(double spotSize, double mosaicity, vec *hkl)
{
    vec usedHKL;
    
    if (hkl == NULL)
    {
        MatrixPtr newMatrix = MatrixPtr();
        rotateMatrixHKL(latestHRot, latestKRot, 0, matrix, &newMatrix);
        
        usedHKL = new_vector(h, k, l);
        hkl = &usedHKL;
        
        newMatrix->multiplyVector(hkl);
    }
    
    
    spotSize = fabs(spotSize);
    double radMos = fabs(mosaicity) * M_PI / 180;
    
    double distanceFromOrigin = length_of_vector(*hkl);
    
    double spotSizeIncrease = fabs(radMos * distanceFromOrigin);
    
    double radius = (spotSize + spotSizeIncrease);
    
    return radius;
}

double Miller::partialityForHKL(vec hkl, double mosaicity,
                               double spotSize, double wavelength, double bandwidth, double exponent)
{
    bandwidth = fabs(bandwidth);
    
    
    double radius = expectedRadius(spotSize, mosaicity, &hkl);
    
    vec mean_wavelength_ewald_centre = new_vector(0, 0, -1 / wavelength);
    
    vec centred_spot_position = vector_between_vectors(
                                                       mean_wavelength_ewald_centre, hkl);
    
    vec move_inwards_position = copy_vector(centred_spot_position);
    scale_vector_to_distance(&move_inwards_position,
                             length_of_vector(move_inwards_position) - radius);
    
    vec move_outwards_position = copy_vector(centred_spot_position);
    scale_vector_to_distance(&move_outwards_position,
                             length_of_vector(move_outwards_position) + radius);
    
    vec central_wavelength_outwards_position = copy_vector(
                                                           centred_spot_position);
    scale_vector_to_distance(&central_wavelength_outwards_position,
                             1 / wavelength + radius);
    
    vec central_wavelength_inwards_position = copy_vector(
                                                          centred_spot_position);
    
    scale_vector_to_distance(&central_wavelength_inwards_position,
                             1 / wavelength - radius);
    
    add_vector_to_vector(&move_inwards_position, mean_wavelength_ewald_centre);
    add_vector_to_vector(&move_outwards_position, mean_wavelength_ewald_centre);
    
    add_vector_to_vector(&central_wavelength_inwards_position,
                         mean_wavelength_ewald_centre);
    add_vector_to_vector(&central_wavelength_outwards_position,
                         mean_wavelength_ewald_centre);
    
    double inwards_bandwidth = getEwaldSphereNoMatrix(move_inwards_position);
    double outwards_bandwidth = getEwaldSphereNoMatrix(move_outwards_position);
    
    double stdev = wavelength * bandwidth / 2;
    
    double limit = 4 * stdev;
    double inwardsLimit = wavelength + limit;
    double outwardsLimit = wavelength - limit;
    
    if ((outwards_bandwidth > inwardsLimit || inwards_bandwidth < outwardsLimit) && resolution() > 1 / trickyRes)
    {
        return 0;
    }
    
    double thisPartiality = sliced_integral(outwards_bandwidth,
                                                        inwards_bandwidth, radius, 0, 1, wavelength, stdev, exponent);
    
    return thisPartiality;
    
}

double Miller::calculateDefaultNorm()
{
    double defaultSpotSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
    double wavelength = 1.75;
    
    double d = getResolution();
    
    double newH = 0;
    double newK = sqrt((4 * pow(d, 2) - pow(d, 4) * pow(wavelength, 2)) / 4);
    double newL = 0 - pow(d, 2) * wavelength / 2;
    
    vec newHKL = new_vector(newH, newK, newL);
    
    double normPartiality = partialityForHKL(newHKL, 0,
                                             defaultSpotSize, wavelength, INITIAL_BANDWIDTH, INITIAL_EXPONENT);
    
    return normPartiality;
}

double Miller::calculateNormPartiality(double mosaicity,
                                      double spotSize, double wavelength, double bandwidth, double exponent)
{
    vec hkl = new_vector(h, k, l);
    
    matrix->multiplyVector(&hkl);
    
    double d = length_of_vector(hkl);
    
    double newH = 0;
    double newK = sqrt((4 * pow(d, 2) - pow(d, 4) * pow(wavelength, 2)) / 4);
    double newL = 0 - pow(d, 2) * wavelength / 2;
    
    vec newHKL = new_vector(newH, newK, newL);
    
    double normPartiality = partialityForHKL(newHKL, mosaicity,
                                            spotSize, wavelength, bandwidth, exponent);
    
    return normPartiality;
}

double Miller::resolution()
{
    if (resol == 0)
    {
        vec newVec = getTransformedHKL(0, 0);
        
        resol = length_of_vector(newVec);
    }
    
    return resol;
}

void Miller::recalculatePartiality(MatrixPtr rotatedMatrix, double mosaicity,
                                   double spotSize, double wavelength, double bandwidth, double exponent)
{
    if (model == PartialityModelFixed)
    {
        return;
    }
    
    vec hkl = new_vector(h, k, l);
    
    rotatedMatrix->multiplyVector(&hkl);
    
    this->wavelength = getEwaldSphereNoMatrix(hkl);
    
    double tempPartiality = partialityForHKL(hkl, mosaicity,
                                            spotSize, wavelength, bandwidth, exponent);
    
    double normPartiality = 1;
    bool recalculatingNorm = FileParser::getKey("RECALCULATE_NORM", true);
    
    if (normalised && tempPartiality > 0)
    {
        normPartiality = lastNormPartiality;
        if (recalculatingNorm || lastNormPartiality == 0)
        {
            normPartiality = calculateNormPartiality(mosaicity, spotSize, wavelength, bandwidth, exponent);
            lastNormPartiality = normPartiality;
        }
    }
    
    lastWavelength = wavelength;
    lastBandwidth = bandwidth;
    lastRlpSize = spotSize;
    lastMosaicity = mosaicity;
    lastRadius = expectedRadius(spotSize, mosaicity, &hkl);

    partiality = tempPartiality / normPartiality;
    
    if ((!std::isfinite(partiality)) || (partiality != partiality))
    {
        partiality = 0;
    }
}

double Miller::twoTheta(bool horizontal)
{
    double usedWavelength = wavelength;
    if (lastWavelength != 0)
        usedWavelength = lastWavelength;
    
    double rlpCoordinates[2];
    double beamCoordinates[2];
    
    vec hkl = new_vector(h, k, l);
    matrix->multiplyVector(&hkl);
    
    
    vec crystalCentre = new_vector(0, 0, -1 / usedWavelength);
    vec crystalToRlp = vector_between_vectors(crystalCentre, hkl);
    
    rlpCoordinates[0] = (horizontal ? crystalToRlp.h : crystalToRlp.k);
    rlpCoordinates[1] = crystalToRlp.l;
    
    beamCoordinates[0] = 0;
    beamCoordinates[1] = 1 / usedWavelength;
    
  //  double cosTwoTheta = dotProduct / (rlpMagnitude * beamMagnitude);
    
    double hOrK = (horizontal ? hkl.h : hkl.k);
    double effectiveResolution = 1 / sqrt(hOrK * hOrK + hkl.l * hkl.l);
    double sinTheta = usedWavelength / (2 * effectiveResolution);
    double theta = asin(sinTheta);
    
    return 2 * theta;
}

void Miller::setHorizontalPolarisationFactor(double newFactor)
{
    polarisationCorrection = 0;
    
    // not finished
}

double Miller::getPolarisationCorrection()
{
    if (!matrix)
    {
        return 1;
    }
    
    if (polarisationCorrection == 0)
    {
        // P = (1 + cos2 2θM cos2 2θ) / (1 + cos2 2θM)
        
        // calculate polarisation factor

        double horizontalFactor = polarisationFactor;
        double verticalFactor = 1 - polarisationFactor;
        
        double horizontalTwoTheta = twoTheta(true);
        double verticalTwoTheta = twoTheta(false);
        
        double cosSquaredHorizontal = pow(cos(horizontalTwoTheta), 2);
        double cosSquaredVertical = pow(cos(verticalTwoTheta), 2);
        
        double numerator1 = 1 * verticalFactor + cosSquaredHorizontal * horizontalFactor;
        double denominator1 = verticalFactor + horizontalFactor;
        
        double numerator2 = 1 * horizontalFactor + cosSquaredVertical * verticalFactor;
        double denominator2 = horizontalFactor + verticalFactor;
        
        double numerator = (numerator1 + numerator2);
        double denominator = (denominator1 + denominator2);
        
        polarisationCorrection = numerator / denominator;
    }
    
    return polarisationCorrection;
}

MillerPtr Miller::copy(void)
{
    MillerPtr newMiller = MillerPtr(new Miller(mtzParent));
    
    newMiller->h = h;
    newMiller->k = k;
    newMiller->l = l;
    newMiller->rawIntensity = rawIntensity;
    newMiller->model = model;
    newMiller->sigma = sigma;
    newMiller->partiality = partiality;
    newMiller->wavelength = wavelength;
    newMiller->matrix = matrix;
    newMiller->filename = std::string(filename);
    newMiller->countingSigma = countingSigma;
    newMiller->lastX = lastX;
    newMiller->lastY = lastY;
    newMiller->shift = shift;
    
    return newMiller;
}

void Miller::flip(void)
{
    int tmp = l;
    l = k;
    k = tmp;
}

void Miller::applyScaleFactor(double scaleFactor)
{
    setScale(scale * scaleFactor);
}

void Miller::applyPolarisation(double wavelength)
{
    double components[] =
    { -0.0046462837103, 0.00749971430969, -0.00347981101231, 0.00531577207698,
        0.00576709258158, 0.00533160033255, 0.00633224815006,
        0.000661573996986, -0.00702906145495 };
    
    Matrix *mat = new Matrix(components);
    
    vec hkl = new_vector(h, k, l);
    mat->multiplyVector(&hkl);
    double resolution = length_of_vector(hkl);
    
    double sintheta = wavelength / (2 / resolution);
    double theta = asin(sintheta);
    
    double l = sin(0.5 * theta);
    
    if (l == l && l != 0 && std::isfinite(l))
        rawIntensity *= l + 1;
}

bool Miller::positiveFriedel(bool *positive, int *_isym)
{
    int h = getH();
    int k = getK();
    int l = getL();
    
    cctbx::miller::index<> newMiller = cctbx::miller::index<>(h, k, l);
    space_group *spaceGroup = parentReflection->getSpaceGroup();
    asu *asymmetricUnit = parentReflection->getAsymmetricUnit();
    
    cctbx::miller::asym_index asymmetricMiller = cctbx::miller::asym_index(*spaceGroup, *asymmetricUnit, newMiller);
    
    int isym = asymmetricMiller.isym();
    
    //  std::cout << isym << std::endl;
    
    *positive = ((isym > 0) == 1);
    
    return isym != 0;
    /*
    bool fake = FileParser::getKey("FAKE_ANOMALOUS", false);
    
    if (fake && fakeFriedel == -1)
    {
        fakeFriedel = rand() % 2;
    }
    
    if (fake)
    {
        return fakeFriedel;
    }
    
    int positives = 0;
    
    if (h > 0)
        positives++;
    if (k > 0)
        positives++;
    if (l > 0)
        positives++;
    
    return (positives == 1 || positives == 3);*/
}

double Miller::scatteringAngle(ImagePtr image)
{
    double beamX = getImage()->getBeamX();
    double beamY = getImage()->getBeamY();
    
    double distance_from_centre = sqrt(
                                      pow(lastX - beamX, 2) + pow(lastY - beamY, 2));
    double distance_pixels = distance_from_centre * getImage()->getMmPerPixel();
    double detector_distance = getImage()->getDetectorDistance();
    
    double sinTwoTheta = sin(distance_pixels / detector_distance);
    
    return sinTwoTheta;
}

void flattenAngle(double *radians)
{
    while (*radians < 0)
    {
        *radians += M_PI;
    }
    
    while (*radians > 2 * M_PI)
    {
        *radians -= M_PI;
    }
}

void Miller::calculatePosition(double distance, double wavelength, double beamX, double beamY, double mmPerPixel, MatrixPtr transformedMatrix, double *x, double *y)
{
    vec hkl = new_vector(h, k, l);
    transformedMatrix->multiplyVector(&hkl);
    
    double x_mm = (hkl.k * distance / (1 / wavelength + hkl.l));
    double y_mm = (hkl.h * distance / (1 / wavelength + hkl.l));
    
    double x_coord = beamX - x_mm / mmPerPixel;
    double y_coord = beamY - y_mm / mmPerPixel;
    
    lastX = x_coord;
    lastY = y_coord;
    
    *x = lastX;
    *y = lastY;
}

void Miller::positionOnDetector(MatrixPtr transformedMatrix, int *x,
                                int *y)
{
    double x_coord = 0;
    double y_coord = 0;
    
    double distance = getImage()->getDetectorDistance();
    double wavelength = getImage()->getWavelength();
    bool even = shoebox->isEven();

    int beamX = getImage()->getBeamX();
    int beamY = getImage()->getBeamY();
    double mmPerPixel = getImage()->getMmPerPixel();
    
    calculatePosition(distance, wavelength, beamX, beamY, mmPerPixel,
                      transformedMatrix, &x_coord, &y_coord);
    
    int intLastX = int(x_coord);
    int intLastY = int(y_coord);
    
    if (even)
    {
        intLastX = round(x_coord);
        intLastY = round(y_coord);
    }
    
    if (!Panel::shouldUsePanelInfo())
    {
        int search = indexer->getSearchSize();
        getImage()->focusOnAverageMax(&intLastX, &intLastY, search, 1, even);
        
        shift = std::make_pair(intLastX + 0.5 - x_coord, intLastY + 0.5 - y_coord);
        
        *x = intLastX;
        *y = intLastY;
    }
    else
    {
        int search = indexer->getSearchSize();
        Coord bestShift = Panel::shiftForMiller(this);
        
        if (bestShift.first != FLT_MAX)
        {
            double shiftedX = lastX + bestShift.first;
            double shiftedY = lastY + bestShift.second;
            
            int xInt = shiftedX;
            int yInt = shiftedY;
            
            getImage()->focusOnAverageMax(&xInt, &yInt, search, 1, even);
            
            shift = std::make_pair(xInt + 0.5 - x_coord, yInt + 0.5 - y_coord);
            
            *x = xInt;
            *y = yInt;
        }
        else
        {
            int search = indexer->getSearchSize();
            
            getImage()->focusOnAverageMax(&intLastX, &intLastY, search, 1, even);
            
            *x = intLastX;
            *y = intLastY;
        }
    }
}

void Miller::makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize)
{
    double radius = expectedRadius(rlpSize, mosaicity, NULL);
    
    shoebox->complexShoebox(wavelength, bandwidth, radius);
}

void Miller::integrateIntensity(MatrixPtr transformedMatrix)
{
    if (!getImage())
        throw 1;
    
    std::ostringstream logged;
    
    if (!shoebox)
    {
        MillerPtr strongSelf = selfPtr.lock();
        
        shoebox = ShoeboxPtr(new Shoebox(strongSelf));
        
        int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                                  SHOEBOX_FOREGROUND_PADDING);
        int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                               SHOEBOX_NEITHER_PADDING);
        int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                                  SHOEBOX_BACKGROUND_PADDING);
        bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);
        
        logged << "Shoebox created from values " << foregroundLength << ", " << neitherLength << ", " << backgroundLength << std::endl;
        
        shoebox->simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);
    }
    
    int x = 0;
    int y = 0;
    
    positionOnDetector(transformedMatrix, &x, &y);
    
    rawIntensity = getImage()->intensityAt(x, y, shoebox, &countingSigma, 0);

    
    if (rawIntensity > 1000 && false)
    {
        logged << "Raw intensity " << rawIntensity << ", counting sigma " << countingSigma << " for position " << x << ", " << y << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelDebug);
        getImage()->printBox(x, y, 5);
        shoebox->printShoebox();
    }
}

void Miller::incrementOverlapMask(double hRot, double kRot)
{
    int x = lastX;
    int y = lastY;
    
    getImage()->incrementOverlapMask(x, y, shoebox);
}


bool Miller::isOverlappedWithSpots(std::vector<SpotPtr> *spots)
{
    double x = lastX;
    double y = lastY;
    int count = 0;
    int tolerance = FileParser::getKey("METROLOGY_SEARCH_SIZE", 1) + 2;
    
    for (int i = 0; i < spots->size(); i++)
    {
        int x2 = (*spots)[i]->getX();
        int y2 = (*spots)[i]->getY();
        
        double xDiff = fabs(x2 - x);
        double yDiff = fabs(y2 - y);
        
        if (xDiff < tolerance && yDiff < tolerance)
        {
            spots->erase(spots->begin() + i);
            i--;
            count++;
        }
    }
    
    return (count > 0);
}

bool Miller::isOverlapped()
{
    int x = lastX;
    int y = lastY;
    
    unsigned char max = getImage()->maximumOverlapMask(x, y, shoebox);
    
    return (max >= 2);
}

void Miller::setRejected(std::string reason, bool rejection)
{
    rejectedReasons[reason] = rejection;
    
    if (!rejection)
    {
        rejectedReasons.erase(reason);
    }
    
    calculatedRejected = false;
}

bool Miller::isRejected(std::string reason)
{
    if (rejectedReasons.count(reason) == 0)
        return false;
    
    return rejectedReasons[reason];
}

double Miller::averageRawIntensity(vector<MillerPtr> millers)
{
    double allIntensities = 0;
    int num = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
        allIntensities += miller->getRawIntensity();
        num++;
    }
    
    allIntensities /= num;
    
    return allIntensities;
}

double Miller::observedPartiality(double reference)
{
    return rawIntensity / reference;
}

double Miller::observedPartiality(MtzManager *reference)
{
    return parentReflection->observedPartiality(reference, this);
}

Miller::~Miller()
{
//    if (shoebox) std::cout << "Shoebox still exists, deallocating Miller" << std::endl;
}

void Miller::denormalise()
{
    denormaliseFactor = calculateDefaultNorm();
}

// Vector stuff

