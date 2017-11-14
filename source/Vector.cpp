/*
 * Vector.c
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#include "Vector.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include "gaussianfit.h"
#include <tuple>
#include "Matrix.h"
#include "Logger.h"
#include <fstream>

vec reverseVector(vec vec1)
{
    vec vec2 = copy_vector(vec1);
    vec2.h = -vec2.h;
    vec2.k = -vec2.k;
    vec2.l = -vec2.l;

    return vec2;
}

bool vectors_are_equal(vec vec1, vec vec2)
{
    return (vec1.h == vec2.h && vec1.k == vec2.k && vec1.l == vec2.l);
}

double distance_between_vectors(vec vec1, vec vec2)
{
    take_vector_away_from_vector(vec1, &vec2);

    return length_of_vector(vec2);
}

bool within_vicinity(vec vec1, vec vec2, double maxD)
{
    if (fabs(vec1.h - vec2.h) > maxD)
        return false;

    if (fabs(vec1.k - vec2.k) > maxD)
        return false;

    if (fabs(vec1.l - vec2.l) > maxD)
        return false;

    take_vector_away_from_vector(vec1, &vec2);

    if (length_of_vector(vec2) > maxD)
        return false;

    return true;
}

void setFloatingPointErrorZerosToZero(vec *vec1, double limit)
{
    if (fabs(vec1->h) < limit)
        vec1->h = 0;

    if (fabs(vec1->k) < limit)
        vec1->k = 0;

    if (fabs(vec1->l) < limit)
        vec1->l = 0;
}

vec vector_between_vectors(vec vec1, vec vec2)
{
        vec vec;
        vec.h = vec2.h - vec1.h;
        vec.k = vec2.k - vec1.k;
        vec.l = vec2.l - vec1.l;

        return vec;
}

void take_vector_away_from_vector(vec vec1, vec *vec2)
{
        (*vec2).h -= (vec1).h;
        (*vec2).k -= (vec1).k;
        (*vec2).l -= (vec1).l;
}

void add_vector_to_vector(vec *vec1, vec vec2)
{
        (*vec1).h += vec2.h;
        (*vec1).k += vec2.k;
        (*vec1).l += vec2.l;
}

vec perpendicular_for_vectors(vec vec1, vec vec2)
{
    vec perp = cross_product_for_vectors(vec1, vec2);
    scale_vector_to_distance(&perp, 1);

    return perp;
}

MatrixPtr rotation_between_vectors(vec vec1, vec vec2)
{
    MatrixPtr matrix = MatrixPtr(new Matrix());
   // MatrixPtr matrix;

    // Find closest angle between vectors (straightest sweep between the two)
    double cosine = cosineBetweenVectors(vec1, vec2);

    vec crossVector = cross_product_for_vectors(vec1, vec2);
    scale_vector_to_distance(&crossVector, 1);

    double angle = acos(cosine);

    matrix->rotateRoundUnitVector(crossVector, angle);

    return matrix;
}

MatrixPtr closest_rotation_matrix(vec vec1, vec vec2, vec chosenCrossProduct, double *resultantAngle)
{
    bool close = false;

    // we want to minimise the angle between the vectors rotating round chosen axis. This is the starting value
    double lastCosAngle = fabs(angleBetweenVectors(vec1, vec2));
    MatrixPtr mat = MatrixPtr(new Matrix());

    // we step by this amount on each iteration.
    double step = 0.5 * M_PI / 180;

    // in case we start going in the wrong direction, we can switch direction once
    bool switchedOnce = false;

    // for very fine angles when we are closer to the solution
    bool divided = false;

    int cycles = 0;
    /*
    double biggestCosAngle = 0;
    double bestRotation = 0;

    std::ostringstream logged;

    for (double i = 0; i < 2 * M_PI; i += step)
    {
        mat->rotateRoundUnitVector(chosenCrossProduct, i);
        vec vec1Copy = copy_vector(vec1);
        mat->multiplyVector(&vec1Copy);
        double cosTheta = cosineBetweenVectors(vec1Copy, vec2);

        if (cosTheta > biggestCosAngle)
        {
            bestRotation = i;
            biggestCosAngle = cosTheta;
        }

        mat->setIdentity();
    }

    MatrixPtr mat2 = MatrixPtr(new Matrix());
    mat2->rotateRoundUnitVector(chosenCrossProduct, bestRotation);
    logged << mat2->description() << std::endl;
    return mat2;

    mat->setIdentity();
    */
    double totalStep = 0;

    while (!close && cycles < 20000)
    {
        mat->rotateRoundUnitVector(chosenCrossProduct, step);
        vec vec1Copy = copy_vector(vec1);
        mat->multiplyVector(&vec1Copy);
        double cosTheta = cosineBetweenVectors(vec1Copy, vec2);
        //double angleDiff = fabs(angleBetweenVectors(vec1Copy, vec2)); // checked

        if (cosTheta < lastCosAngle)
        {
            if (switchedOnce)
                close = true;

            step = -step;
            switchedOnce = true;
        }

        if (cosTheta > 0.95 && !divided)
        {
            step /= 10;
            divided = true;
        }

        totalStep += step;
        lastCosAngle = cosTheta;
        cycles++;
    }

    *resultantAngle = acos(lastCosAngle);

    return mat;

}

MatrixPtr rotation_between_vectors_custom_cross(vec vec1, vec vec2, vec chosenCrossProduct)
{
    scale_vector_to_distance(&vec1, 1);
    scale_vector_to_distance(&vec2, 1);
    scale_vector_to_distance(&chosenCrossProduct, 1);

    MatrixPtr vx = MatrixPtr(new Matrix());
    MatrixPtr vxSquared = MatrixPtr(new Matrix());
    double c = cosineBetweenVectors(vec1, vec2);


    vec v = cross_product_for_vectors(vec1, vec2);// chosenCrossProduct;

    vx = MatrixPtr(new Matrix());
    vx->components[0] = 0;
    vx->components[1] = v.l;
    vx->components[2] = -v.k;
    vx->components[4] = -v.l;
    vx->components[5] = 0;
    vx->components[6] = v.h;
    vx->components[8] = v.k;
    vx->components[9] = -v.h;
    vx->components[10] = 0;

    vxSquared = vx->copy();
    vxSquared->multiply(*vx);

    vec vxvec1 = copy_vector(vec1);
    vx->multiplyVector(&vxvec1);

    vec vxSquaredVec1 = copy_vector(vec1);
    vxSquared->multiplyVector(&vxSquaredVec1);

    vec vec3 = vec2;
    take_vector_away_from_vector(vxvec1, &vec3);
    take_vector_away_from_vector(vec1, &vec3);


    double scale = (1 - c);
    vxSquared->multiply(scale);

    MatrixPtr rotation = MatrixPtr(new Matrix());
    rotation->add(vx);
    rotation->add(vxSquared);

    return rotation;
}

vec cross_product_for_vectors(vec vec1, vec vec2)
{
        double new_h = vec1.k * vec2.l - vec1.l * vec2.k;
        double new_k = vec1.l * vec2.h - vec1.h * vec2.l;
        double new_l = vec1.h * vec2.k - vec1.k * vec2.h;

        return new_vector(new_h, new_k, new_l);
}

double dot_product_for_vectors(vec vec1, vec vec2)
{
    return vec1.h * vec2.h + vec1.k * vec2.k + vec1.l * vec2.l;
}

vec copy_vector(vec old_vec)
{
        vec vec;
        vec.h = old_vec.h;
        vec.k = old_vec.k;
        vec.l = old_vec.l;

        return vec;
}

double cosineBetweenVectors(vec vec1, vec vec2)
{
    double dotProduct = vec1.h * vec2.h + vec1.k * vec2.k + vec1.l * vec2.l;
    double vec1_length = length_of_vector(vec1);
    double vec2_length = length_of_vector(vec2);

    double cosTheta = dotProduct / (vec1_length * vec2_length);

    return cosTheta;
}

double angleBetweenVectors(vec vec1, vec vec2)
{
    double cosTheta = cosineBetweenVectors(vec1, vec2);

    if (cosTheta > 1)
        cosTheta = 1;

    if (cosTheta < -1)
        cosTheta = -1;

        double angle = acos(cosTheta);

    if (angle != angle && (cosTheta < 1.0001))
        angle = 0;

        return angle;
}

double length_of_vector_squared(vec vec)
{
    return vec.h * vec.h + vec.k * vec.k + vec.l * vec.l;
}

double length_of_vector(vec vec)
{
        return pow(vec.h * vec.h + vec.k * vec.k + vec.l * vec.l, 0.5);
}

vec new_vector(double h, double k, double l)
{
        vec vec;
        vec.h = h;
        vec.k = k;
        vec.l = l;

        return vec;
}

void scale_vector_to_distance(vec *vector, double new_distance)
{
        double distance = sqrt(
                        pow(vector->h, 2) + pow(vector->k, 2) + pow(vector->l, 2));

        double scale = new_distance / distance;

        (*vector).h *= scale;
        (*vector).k *= scale;
        (*vector).l *= scale;
}

double getEwaldSphereNoMatrix(vec index)
{
    if (index.l == 0)
        return 0;

    double ewald_radius = index.h * index.h + index.k * index.k
                        + index.l * index.l;

        ewald_radius /= (0 - 2 * index.l);
        double ewald_wavelength = 1 / ewald_radius;

        return ewald_wavelength;
}

double getEwaldWeightForAxis(vec index, bool isH)
{
    vec axisVec = isH ? new_vector(1, 0, 0) : new_vector(0, 1, 0);
    vec indexVec = new_vector(index.h, index.k, 0);

    double cos = cosineBetweenVectors(axisVec, indexVec);

    return 1 - fabs(cos);
}

double cartesian_to_distance(double x, double y)
{
    double distance = sqrt(pow(x, 2) + pow(y, 2));

    return distance;
}

double cartesian_to_angle(double x, double y)
{
    double angle = atan(y / x);

    if ((x < 0 && y > 0) ||
        (x < 0 && y < 0))
        angle += M_PI;

    if (x > 0 && y < 0)
        angle += M_PI * 2;

    return angle;
}

double cdf(double x, double mean, double sigma)
{
        double y = (x - mean) / sigma;

        double cumulative = _cdf(y);

        return cumulative;
}

double _cdf(double x)
{
        // constants
        double a1 = 0.254829592;
        double a2 = -0.284496736;
        double a3 = 1.421413741;
        double a4 = -1.453152027;
        double a5 = 1.061405429;
        double p = 0.3275911;

        // Save the sign of x
        int sign = 1;
        if (x < 0)
                sign = -1;
        x = fabs(x) / sqrt(2.0);

        // A&S formula 7.1.26
        double t = 1.0 / (1.0 + p * x);
        double y = 1.0
                        - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
                                        * exp(-x * x);

        return 0.5 * (1.0 + sign * y);
}

double normal_distribution(double x, double mean, double sigma)
{
        double power = 0 - pow((x - mean), 2) / (2 * sigma * sigma);
        double exp = pow(M_E, power);

        double denominator = sigma * sqrt(2 * M_PI);

        return exp / denominator;
}

double super_gaussian(double x, double mean, double sigma_0, double exponent)
{
        double correction_sigma = 2 / exponent - 1;
        double sigma = pow(M_PI / 2, correction_sigma) * sigma_0;

        double power = 0 - pow(fabs((x - mean)), exponent) / (2 * pow(sigma, exponent));

    double exp = pow(M_E, power);

        double denominator = sigma * sqrt(2 * M_PI);

        return exp / denominator;
}

double correlation_between_vectors(vector<double> *vec1,
                vector<double> *vec2)
{
        return correlation_between_vectors(vec1, vec2, NULL);
}

double correlation_between_vectors(vector<double> *vec1,
                vector<double> *vec2, vector<double> *weights)
{
        return correlation_between_vectors(vec1, vec2, weights, -1);
}

double correlation_through_origin(vector<double> *vec1,
                vector<double> *vec2, vector<double> *weights)
{
        double mean_y = weighted_mean(vec2, weights);
        double grad = gradient_between_vectors(vec1, vec2);

        double residuals_squared = 0;
        double denominator = 0;

        for (int i = 0; i < vec1->size(); i++)
        {
                double mean1 = (*vec1)[i];
                double mean2 = (*vec2)[i];

                if (mean1 != mean1 || mean2 != mean2)
                        continue;

                residuals_squared += pow(mean2 - grad * mean1, 2);
                double addition = pow(mean2 - mean_y, 2);

                denominator += addition;
        }

        double R_squared = 1 - residuals_squared / denominator;

        if (R_squared < 0)
                R_squared = 0;
        double R = sqrt(R_squared);
        if (grad < 0)
                R = 0;
        if (R != R)
                R = -1;

        return R;
}

double sum(vector<double> values)
{
    double sum = 0;
    for (int i = 0; i < values.size(); i++)
    {
        sum += values[i];
    }

    return sum;
}

double correlation_between_vectors(vector<double> *vec1,
                vector<double> *vec2, vector<double> *weights, int exclude)
{
        double sum_x = 0;
        double sum_y = 0;
        double num = 0;
        bool needWeights = (weights == NULL);

        for (int i = 0; i < vec1->size(); i++)
        {
                if (i == exclude)
                {
                        continue;
                }

                if (weights != NULL && (*weights)[i] == 0)
                        continue;

                double addition = (*vec1)[i];
                if (weights != NULL)
                        addition *= (*weights)[i];
                sum_x += addition;

                addition = (*vec2)[i];
                if (weights != NULL)
                        addition *= (*weights)[i];
                sum_y += addition;

                double weight = (weights == NULL) ? 1 : (*weights)[i];
                num += weight;
        }
        double mean_x = sum_x / num;
        double mean_y = sum_y / num;

        if (mean_x != mean_x || mean_y != mean_y)
                return 0;

        double sum_x_y_minus_mean_x_y = 0;
        double sum_x_minus_mean_x_sq = 0;
        double sum_y_minus_mean_y_sq = 0;

        for (int i = 0; i < vec1->size(); i++)
        {
                if (i == exclude)
                        continue;

                if (weights != NULL && (*weights)[i] == 0)
                        continue;

                double addition = ((*vec1)[i] - mean_x) * ((*vec2)[i] - mean_y);
                if (weights != NULL)
                        addition *= (*weights)[i];
                sum_x_y_minus_mean_x_y += addition;

                addition = std::pow((*vec1)[i] - mean_x, 2);
                if (weights != NULL)
                        addition *= (*weights)[i];
                sum_x_minus_mean_x_sq += addition;

                addition = std::pow((*vec2)[i] - mean_y, 2);
                if (weights != NULL)
                        addition *= (*weights)[i];

                sum_y_minus_mean_y_sq += addition;

        }

        sum_x_y_minus_mean_x_y /= num;
        sum_x_minus_mean_x_sq /= num;
        sum_y_minus_mean_y_sq /= num;

        double r = sum_x_y_minus_mean_x_y
                        / (sqrt(sum_x_minus_mean_x_sq) * sqrt(sum_y_minus_mean_y_sq));

        if (needWeights)
                delete weights;

        return r;
}

double gradient_between_vectors(vector<double> *vec1,
                vector<double> *vec2)
{

        double sum_x_y = 0;
        double sum_x_squared = 0;

        for (int i = 0; i < vec1->size(); i++)
        {
                sum_x_y += (*vec1)[i] * (*vec2)[i];
                sum_x_squared += std::pow((*vec1)[i], 2);
        }

        double grad = sum_x_y / sum_x_squared;

        return grad;
}

double minimize_gradient_between_vectors(vector<double> *vec1,
                vector<double> *vec2)
{
        double firstSlope = gradient_between_vectors(vec1, vec2);

        double bestSlope = firstSlope;
        double bestScore = least_squares_between_vectors(vec1, vec2, bestSlope);

        int count = 0;
        bool optimised = false;
        double step = 0.05;
        double finished = 0.001;

        while (!optimised && count < 10)
        {
                double higherSlope = bestSlope + step;
                double lowerSlope = bestSlope - step;

                double higherScore = least_squares_between_vectors(vec1, vec2, higherSlope);
                double lowerScore = least_squares_between_vectors(vec1, vec2, lowerSlope);

                if (bestScore < higherScore && bestScore < lowerScore)
                {
                        step /= 2;
                        continue;
                }

                double newScore = (higherScore < lowerScore) ? higherScore : lowerScore;
                double newSlope = (higherScore < lowerScore) ? higherSlope : lowerSlope;

                bestSlope = newSlope;
                bestScore = newScore;

                count++;

                if (step < finished)
                        break;
        }

        return bestSlope;
}


double least_squares_between_vectors(vector<double> *vec1,
                vector<double> *vec2, double slope)
{
        if (slope == 0)
        {
        slope = gradient_between_vectors(vec1, vec2);
        }

        if (vec1->size() == 0)
                return FLT_MAX;

        double total = 0;
        double num = 0;

        for (int i = 0; i < vec1->size(); i++)
        {
                double a = (*vec1)[i];
                double b = (*vec2)[i] * slope;

                total += pow(a - b, 2);
                num++;
        }

        return total / num;
}

double r_factor_between_vectors(vector<double> *vec1,
                vector<double> *vec2, vector<double> *weights, double scale)
{
        double sum_numerator = 0;
        double sum_denominator = 0;
        double absolute = 0;
        double allWeights = 0;

        for (int i = 0; i < vec1->size(); i++)
        {
                double weight = (*weights)[i];
                double int1 = (*vec1)[i] * scale;
                double int2 = (*vec2)[i];

                if (weight != weight || int1 != int1 || int2 != int2)
                        continue;

                allWeights += weight;

                double sqrtD = 1;

                double absAddition = fabs(int1 - int2) * sqrtD;
                absAddition /= (int1 + int2) * sqrtD / 2;

                absolute += fabs(absAddition) * weight;

                sum_numerator += fabs(int1 - int2) * weight;
                sum_denominator += (int1 + int2) * weight / 2;
        }

        double r_split = sum_numerator / (sum_denominator * sqrt(2));

        return r_split;
}

double weighted_mean(vector<double> *means, vector<double> *weights)
{
        double sum = 0;
        double weight_sum = 0;

        for (int i = 0; i < means->size(); i++)
        {
                double weight = (weights == NULL ? 1 : (*weights)[i]);

                sum += (*means)[i] * weight;
                weight_sum += weight;
        }

        return sum / weight_sum;
}

bool higher(double mean1, double mean2)
{
    return mean1 > mean2;
}

double median(vector<double> *means)
{
    std::sort(means->begin(), means->end(), higher);
    int mid = (int)means->size() / 2;
    double midPoint = 0;

    if (means->size() % 2 == 0)
    {
        midPoint = ((*means)[mid] + (*means)[mid + 1]) / 2;
    }
    else
        midPoint = (*means)[mid];

    return midPoint;
}

void histogram_gaussian(vector<double> *means, vector<int> *freq,
                double &mean, double &stdev)
{
        double sum = 0;
        double weight_sum = 0;

        for (int i = 0; i < means->size(); i++)
        {
        double weight = (*freq)[i];
                sum += (*means)[i] * weight;
                weight_sum += weight;
        }

        mean = sum / weight_sum;

        double squaredSum = 0;
        double weightSqSum = 0;

        for (int i = 0; i < means->size(); i++)
        {
        double weight = (*freq)[i];
                squaredSum += pow(mean - (*means)[i], 2) * weight;
                weightSqSum += weight;
        }

        stdev = sqrt(squaredSum / weightSqSum);
}

double least_squares_gaussian_fit(vector<double> *means,
                vector<int> *freq)
{
        double mean = 0;
        double stdev = 0;
        double score = 0;

        gaussian_fit(*means, *freq, (int)means->size(), &mean, &stdev, &score, false);

        return score;
}

double standard_deviation(vector<double> *values, vector<double> *weights)
{
    double mean = weighted_mean(values, weights);

    return standard_deviation(values, weights, mean);
}

double standard_deviation(vector<double> *values, vector<double> *weights, double mean)
{
        double squaredSum = 0;
        double weightSqSum = 0;

        for (int i = 0; i < values->size(); i++)
        {
                double value = (*values)[i];

                if (value != value || value == FLT_MAX)
                        continue;

                squaredSum += pow(mean - value, 2);

        double weight = (weights == NULL) ? 1 : (*weights)[i];

                weightSqSum += weight;
        }

        double stdev = sqrt(squaredSum / weightSqSum);

        return stdev;
}

void regression_line(vector<boost::tuple<double, double, double> > values, double &intercept, double &gradient)
{
        double sigma_x = 0;
        double sigma_y = 0;
        double sigma_x_y = 0;
        double sigma_x_2 = 0;
        double weight_n = 0;

        for (int i=0; i < values.size(); i++)
        {
                double x = boost::get<0>(values[i]);
                double y = boost::get<1>(values[i]);
                double weight = boost::get<2>(values[i]);

                sigma_x += x * weight;
                sigma_y += y * weight;
                sigma_x_y += x * y * weight;
                sigma_x_2 += x * x * weight;
                weight_n += weight;
        }

        double mean_x = sigma_x / weight_n;
        double mean_y = sigma_y / weight_n;

        double sxy = sigma_x_y - sigma_x * sigma_y / weight_n;
        double sxx = sigma_x_2 - pow(sigma_x, 2) / weight_n;

        gradient = sxy / sxx;

        intercept = mean_y - gradient * mean_x;
}

double minimizeParam(double &step, double &param, double (*score)(void *object),
                     void *object)
{
    return minimizeParameter(step, &param, score, object);
}

double minimizeParameter(double &step, double *param, double (*score)(void *object),
                                    void *object)
{
    double param_trials[3];
    double param_scores[3];

    int j = 0;
    int param_min_num = 1;

    double bestParam = *param;

    for (double i = bestParam - step; j < 3; i += step)
    {
        *param = i;

        double aScore = (*score)(object);

        if (aScore != aScore) aScore = FLT_MAX;

        param_scores[j] = aScore;
        param_trials[j] = i;
        j++;
    }

    double param_min_score = param_scores[1];

    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }

    *param = param_trials[param_min_num];
    (*score)(object);

    if (param_min_num == 1)
        step /= 2;

    return param_scores[param_min_num];
}

std::map<double, int> histogram(std::vector<double> values, double step)
{
    std::map<int, int> tempHistogram;

    for (int i = 0; i < values.size(); i++)
    {
        if (values[i] != values[i])
            continue;

        int category = values[i] / step;

        if (values[i] == 0)
            category = 0;

        tempHistogram[category]++;
    }

    std::map<double, int> realHistogram;

    if (values.size() == 0)
        return realHistogram;

   // Logger::mainLogger->addString("Histogram has " + std::to_string((int)tempHistogram.size()) + " categories.");

    int minCategory = INT_MAX;
    int maxCategory = -INT_MAX;

    for (std::map<int, int>::iterator it = tempHistogram.begin(); it != tempHistogram.end(); it++)
    {
        int category = it->first;

        if (category < minCategory)
            minCategory = category;
        if (category > maxCategory)
            maxCategory = category;
    }

    for (int i = minCategory; i <= maxCategory; i++)
    {
        double realStep = i * step;
        int frequency = tempHistogram[i];

        realHistogram[realStep] = frequency;
    }

    return realHistogram;
}

void histogramCSV(std::string filename, std::map<double, int> map1, std::map<double, int> map2)
{
    std::ofstream stream;
    stream.open(filename.c_str());

    for (std::map<double, int>::iterator it = map1.begin(); it != map1.end(); it++)
    {
        double category = it->first;
        int first = it->second;
        int second = 0;

        if (map2.count(category))
        {
            second = map2[category];
        }

        stream << category << "," << first << "," << second << std::endl;
    }

    stream.close();
}
