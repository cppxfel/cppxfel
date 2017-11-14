#include <string>
#include <iostream>
#include "MtzManager.h"
#include "ScalingManager.h"
#include "ReflectionManager.h"
#include <new>
#include <vector>
#include <cmath>
#include <algorithm>
#include "misc.h"
#include "parameters.h"

#include "Holder.h"
#include "Miller.h"


int primary_call = 0;
int primary_call_2 = 0;
int primary_call_3 = 0;

void ReflectionManager::GForImage(ImageReflection *image,
                vector<Scale_factor> *Gs, Scale_factor **G, int *num)
{
        primary_call_2++;
        double l = image->l;

    if (num != NULL)
        *num = l;

        *G = &(*Gs)[l];
}

void ReflectionManager::reflection_for_image(double l, Reflection **reflection)
{
        primary_call_3++;
        int lower = 0;
        int higher = (int)reflections.size() - 1;
        int new_bound = (higher + lower) / 2;

        int i = 0;

        if (reflections.size() == 2)
        {
                if (reflections[0].l == l)
                        *reflection = reflections[0].reflection;

                else if (reflections[1].l == l)
                        *reflection = reflections[1].reflection;

                else
                        *reflection = NULL;

                return;
        }

        if ((l < reflections[lower].l) || (l > reflections[higher].l))
        {
                *reflection = NULL;
                return;
        }

        while (reflections[new_bound].l != l)
        {
                i++;
                if (higher == lower + 1)
                {
                        if (reflections[lower].l == l)
                        {
                                *reflection = reflections[lower].reflection;
                                return;
                        }
                        else if (reflections[higher].l == l)
                        {
                                *reflection = reflections[higher].reflection;
                                return;
                        }
                        else
                        {
                                *reflection = NULL;
                                return;
                        }
                }

                if (reflections[new_bound].l > l)
                {
                        higher = new_bound;
                }
                else if (reflections[new_bound].l < l)
                {
                        lower = new_bound;
                }

                new_bound = (higher + lower) / 2;
        }

        (*reflection) = reflections[new_bound].reflection;
}

void ReflectionManager::alpha_beta_hj(vector<Scale_factor> Gs, int j,
                double *alpha_hj, double *beta_hj)
{
        primary_call++;
        double ref_intensity = reflections[j].reflection->meanIntensity();

        double ref_sigma = 1;

        Scale_factor *G_ref;
        GForImage(&reflections[j], &Gs, &G_ref);

        double G1 = Gs[0].G;

        *beta_hj = ref_intensity / ref_sigma;
        *alpha_hj = (*G_ref).G / (G1 * ref_sigma);
}

double ReflectionManager::E_alpha(vector<Scale_factor> Gs, int l,
                double *sum_alpha_beta, double *sum_alpha_alpha)
{
        double beta_hl, alpha_hl = 0;
        alpha_beta_hj(Gs, l, &alpha_hl, &beta_hl);

        for (int m = 0; m < reflections.size(); m++)
        {
                double beta_hm, alpha_hm;
                alpha_beta_hj(Gs, m, &alpha_hm, &beta_hm);

                (*sum_alpha_beta) += alpha_hm * beta_hm;
                (*sum_alpha_alpha) += alpha_hm * alpha_hm;
        }

        double E_alpha = beta_hl
                        - (alpha_hl * (*sum_alpha_beta)) / (*sum_alpha_alpha);

        return E_alpha;
}

double ReflectionManager::l_sum_E_alpha_2(vector<Scale_factor> Gs)
{
        double sum_over_images = 0;

        /* when a reflection is not present in an image the E_alpha
         term is equal to 0 due to the beta terms being 0
         (intensity = 0) */

        for (int l = 0; l < reflections.size(); l++)
        {
                double sum_alpha_m_squared = 0;
                double sum_alpha_m_beta_m = 0;

                double e_alpha = E_alpha(Gs, l, &sum_alpha_m_beta_m,
                                &sum_alpha_m_squared);
                double E_alpha_squared = e_alpha * e_alpha;
                sum_over_images += E_alpha_squared;
        }

        return sum_over_images;
}

double ReflectionManager::contributionToGradient(vector<Scale_factor> Gs, int l)
{
        if (reflections.size() == 1)
                return 0;

        if (l == 0)
        {
                double sum_over_images = 0;

                for (int j = 0; j < reflections.size(); j++)
                {
                        double sum_alpha_m_squared = 0;
                        double sum_alpha_m_beta_m = 0;

                        double e_alpha = E_alpha(Gs, j, &sum_alpha_m_beta_m,
                                        &sum_alpha_m_squared);

                        if (e_alpha == 0)
                                continue;

                        double beta_hj, alpha_hj;
                        alpha_beta_hj(Gs, j, &alpha_hj, &beta_hj);

                        double m2_first_sum = 0;
                        double m2_second_sum = 0;

                        double partiality = reflections[j].reflection->meanPartiality();
                        double sigma_j = 1 / sqrt(partiality);
                        sigma_j = 1;

                        Reflection *lth_ref;
                        reflection_for_image(l, &lth_ref);

                        double epsilon = (reflections[j].reflection != lth_ref);
                        Scale_factor *G_j;
                        GForImage(&reflections[j], &Gs, &G_j);

                        double forgotten_term = 0;

                        if (sigma_j != 0 && sum_alpha_m_squared != 0)
                        {
                                forgotten_term = epsilon * (*G_j).G
                                                / (pow(Gs[0].G, 2) * sigma_j);
                                forgotten_term *= sum_alpha_m_beta_m / sum_alpha_m_squared;
                        }

                        for (int m = 1; m < reflections.size(); m++)
                        {
                                double beta_hm, alpha_hm;
                                alpha_beta_hj(Gs, m, &alpha_hm, &beta_hm);

                                Scale_factor *G_ref;
                                GForImage(&reflections[m], &Gs, &G_ref);
                                double Gm = (*G_ref).G;
                                double partiality = reflections[m].reflection->meanPartiality();
                                double sigma_m = 1 / sqrt(partiality);
                                sigma_m = 1;

                                m2_first_sum += beta_hm * Gm / sigma_m;

                                m2_second_sum += pow(Gm, 2) / pow(sigma_m, 2);
                        }

                        double first_term = e_alpha;
                        double second_term = alpha_hj / pow(sum_alpha_m_squared, 2);

                        if (alpha_hj == 0)
                                second_term = 0;

                        double left_side_third = sum_alpha_m_squared / pow(Gs[0].G, 2)
                                        * m2_first_sum;

                        double right_side_third = 2 * sum_alpha_m_beta_m / pow(Gs[0].G, 3)
                                        * m2_second_sum;

                        double third_term = left_side_third - right_side_third;

                        double addition = 2 * first_term
                                        * (forgotten_term + second_term * third_term);

                        if (first_term != first_term || second_term != second_term
                                        || third_term != third_term)
                                addition = 0;

                        if (addition != addition)
                                addition = 0;

                        sum_over_images += addition;
                }

                return sum_over_images;
        }
        else
        {
                Reflection *lth_ref;
                reflection_for_image(l, &lth_ref);

                if (lth_ref == NULL)
                {
                        return 0;
                }

                double sum_over_images = 0;

                for (int j = 0; j < reflections.size(); j++)
                {
                        double sum_alpha_m_squared = 0;
                        double sum_alpha_m_beta_m = 0;

                        double e_alpha = E_alpha(Gs, j, &sum_alpha_m_beta_m,
                                        &sum_alpha_m_squared);

                        if (e_alpha == 0)
                                continue;

                        double beta_hj, alpha_hj;
                        alpha_beta_hj(Gs, j, &alpha_hj, &beta_hj);

                        double sigma_il = 0;
                        double intensity_il = 0;

                        int delta_lj = (reflections[j].reflection == lth_ref);

                        double partiality = lth_ref->meanPartiality();
                        sigma_il = 1 / sqrt(partiality);
                        sigma_il = 1;

                        intensity_il = lth_ref->meanIntensity();

                        double alpha_il = Gs[l].G / (Gs[0].G * sigma_il);
                        double beta_il = intensity_il / sigma_il;

                        double first_l_term = 0;
                        double second_l_term = 0;
                        double third_l_term = 0;

                        if (sigma_il != 0)
                        {
                                first_l_term = delta_lj / (Gs[0].G * sigma_il);
                                second_l_term = beta_il / (Gs[0].G * sigma_il);
                                third_l_term = 2 * alpha_il / (Gs[0].G * sigma_il);
                        }

                        double second_term_left = first_l_term * sum_alpha_m_beta_m
                                        / sum_alpha_m_squared;
                        if (sum_alpha_m_squared == 0)
                                second_term_left = 0;

                        double second_term_right_left = alpha_hj
                                        / pow(sum_alpha_m_squared, 2);

                        if (sum_alpha_m_squared == 0)
                                second_term_right_left = 0;

                        double second_term_right_right = sum_alpha_m_squared * second_l_term
                                        - sum_alpha_m_beta_m * third_l_term;

                        double second_term = second_term_left
                                        + second_term_right_left * second_term_right_right;
                        double addition = -2 * e_alpha * second_term;

                        if (addition != addition)
                        {
                                addition = 0;
                        }

                        sum_over_images += addition;
                }

                return sum_over_images;
        }
}

bool reflection_comparison(ImageReflection ref1, ImageReflection ref2)
{
        return (ref1.l < ref2.l);
}

void ReflectionManager::sortReflections(void)
{
        std::sort(reflections.begin(), reflections.end(), reflection_comparison);

        /*      double mean_intensity = 0;

         for (int i=0; i < reflections.size(); i++)
         {
         mean_intensity += (*reflections[i].reflection).mean_intensity;
         }

         mean_intensity /= reflections.size();

         for (int i=0; i < reflections.size(); i++)
         {
         if ((*reflections[i].reflection).mean_intensity < mean_intensity)
         {
         reflections.erase(reflections.begin()+i);
         i--;
         }
         }*/
}

void ReflectionManager::addReflectionForImage(Reflection *reflection, MtzManager *manager,
                double l)
{
        double intensity = reflection->meanIntensity();

        if (intensity != intensity)
                return;

        int num = (int)reflections.size();

        reflections.resize(num + 1);

        reflections[num].reflection = reflection;
        reflections[num].manager = manager;
        reflections[num].l = l;
}

double ReflectionManager::intensity(vector<Scale_factor> *Gs, double *sigma)
{
    double intensity = 0;
    double weights = 0;
    int num = 0;

    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i].reflection->millerCount(); j++)
        {
            int l = 0;

            Scale_factor *G = NULL;
            GForImage(&reflections[i], Gs, &G, &l);

            MillerPtr miller = reflections[i].reflection->miller(j);

            if (!miller->accepted())
                continue;

            double singleIntensity = miller->intensity();
            double weight = miller->getWeight();

            if (weight != weight || weight == 0)
                continue;

            double scale = fabs(G->G);

            if (scale < 0.1)
                continue;

            weight *= scale;
            singleIntensity *= scale;

            intensity += singleIntensity * weight;
            weights += weight;
            num++;
        }
    }

    intensity /= weights;
    *sigma = weights / num;

    return intensity;
}

Reflection *ReflectionManager::mergedReflection(vector<Scale_factor> *Gs, bool half, bool all)
{
    Reflection *newReflection = reflections[0].reflection->copy(false);

    newReflection->setFlipAsActiveAmbiguity();

    int _h = newReflection->miller(0)->getH();
    int _k = newReflection->miller(0)->getK();
    int _l = newReflection->miller(0)->getL();

    int h, k, l;

    ccp4spg_put_in_asu(group, _h, _k, _l, &h, &k, &l);

    MillerPtr miller = MillerPtr(new Miller(NULL, h, k, l));

    miller->setParent(newReflection);

    double newSigma = 0;
    double newIntensity = 0;
    double int1, int2;
    splitIntensities(Gs, &int1, &int2);

    if (!all)
    {
        newIntensity = half ? int1 : int2;
    }
    else
    {
        newIntensity = intensity(Gs, &newSigma);
    }

    if ((newIntensity < int1 && newIntensity < int2) || (newIntensity > int1 && newIntensity > int2))
        std::cout << "Warning intensities!" << std::endl;

    miller->setData(newIntensity, 1, 1, 0);

    newReflection->clearMillers();
    newReflection->addMiller(miller);

    return newReflection;
}

void ReflectionManager::splitIntensities(vector<Scale_factor> *Gs,
                                 double *int1, double *int2)
{
    double intensity1 = 0;
    double intensity2 = 0;
    double weights1 = 0;
    double weights2 = 0;

    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i].reflection->millerCount(); j++)
        {
            int l = 0;

            Scale_factor *G = NULL;
            GForImage(&reflections[i], Gs, &G, &l);

            double &intensity = (l % 2 == 0) ? intensity2 : intensity1;
            double &weights = (l % 2 == 0) ? weights2 : weights1;

            MillerPtr miller = reflections[i].reflection->miller(j);

            if (!miller->accepted())
                continue;

            double singleIntensity = miller->intensity();
            double weight = miller->getWeight();

            if (weight != weight || weight == 0)
                continue;

            double scale = fabs(G->G);

            if (scale < 0.1)
                continue;

            weight *= scale;
            singleIntensity *= scale;

            intensity += singleIntensity * weight;
            weights += weight;
        }
    }

    intensity1 /= weights1;
    intensity2 /= weights2;

    *int1 = intensity1;
    *int2 = intensity2;
}

double ReflectionManager::rMerge(vector<Scale_factor> *Gs,
                double *numContribution, double *denomContribution, double resolution)
{
        if (reflections.size() < 2)
                return 0;

    double mean = 0;
    double count = 0;
    int num = 0;

    if (intensities.size() == 0)
    {
        for (int i = 0; i < reflections.size(); i++)
        {
            if (!reflections[i].reflection->betweenResolutions(0, resolution))
                continue;

            for (int j = 0; j < reflections[i].reflection->millerCount(); j++)
            {
                Scale_factor *G = NULL;
                GForImage(&reflections[i], Gs, &G);
                MillerPtr miller = reflections[i].reflection->miller(j);

                if (!miller->accepted())
                    continue;

                double intensity = miller->intensity();
                double *scale = &(G->G);

                intensities.push_back(intensity);
                scales.push_back(scale);
            }
        }
    }

    for (int i = 0; i < intensities.size(); i++)
    {
        double intensity = intensities[i];
        double scale = fabs(*(scales[i]));

        if (scale < 0.1)
            continue;

        intensity *= scale;

        mean += intensity;
        count++;
        num++;
    }


    if (num <= 1)
                return 0;

        mean /= intensities.size();

        if (mean != mean || std::isinf(mean))
                return 0;

        double numerator = 0;
        double denominator = 0;


    for (int i = 0; i < intensities.size(); i++)
    {
        double intensity = intensities[i];
        double scale = *scales[i];

        intensity *= scale;

        double numAddition = fabs((mean - intensity));
        double denomAddition = fabs(intensity);

        numerator += numAddition;
        denominator += denomAddition;
        }

        double correction_term = sqrt(num / (num - 1));

        if (numerator != numerator || denominator != denominator)
                return 0;

        *numContribution += correction_term * numerator;
        *denomContribution += denominator;

        return correction_term * numerator / denominator;
}

ReflectionManager::ReflectionManager(void)
{
        refl_id = 0;
        prev_num = -1;
        reflections.resize(0);
}

ReflectionManager::~ReflectionManager(void)
{
//      std::cout << "deallocating RefManager" << std::endl;

//      print_trace();
}
