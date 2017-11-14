#include <string>
#include <iostream>
#include "MtzManager.h"
#include "ScalingManager.h"
#include <new>
#include <limits>
#include <cmath>
#include "FileParser.h"
#include "definitions.h"
#include "Holder.h"
#include "Miller.h"

//#include <boost/thread/thread.hpp>

#define MULT_G 1
#define MULT_WAVELENGTH 1




double ScalingManager::evaluate_psi()
{
    double total = 0;

    for (int h = 0; h < (*refs).size(); h++)
    {
        double addition = ((*refs)[h]).l_sum_E_alpha_2(Gs);
        total += addition;
    }

    return total;
}

void ScalingManager::evaluate_gradient_at_l(int l)
{
    double total = 0;

    int j = 0;

    for (int i = 0; i < (*refs).size(); i++)
    {
        double addition = ((*refs)[i]).contributionToGradient(Gs, l);
        total += addition;

        if (addition != 0)
            j++;
    }

    //  std::cout << "G for (l=" << l << ") = " << total << " using "
    //                  << j << " reflections." << std::endl;

    /*if ((*Gs[l].reflections).size() > 0) {
     std::cout << "Reflections ";
     for (int i = 0; i < (*Gs[l].reflections).size(); i++) {
     // std::cout << (*Gs[l].reflections)[i] << ", ";
     }

     std::cout << std::endl;
     }
     */
    gradients[l] = total;
}


double ScalingManager::set_Gs(scitbx::af::shared<double> new_Gs)
{
    int count = 0;

    for (int i = 0; i < Gs.size(); i++)
    {
        Gs[i].G = new_Gs[i];
        count++;
    }

    return Gs.size();
}

int ScalingManager::orderFinalReflection(void)
{
    int j = (int)(*refs).size() - 1;

    while (j >= 1 && ((*refs)[j]).refl_id < ((*refs)[j - 1]).refl_id)
    {
        ReflectionManager tmp = (*refs)[j - 1];
        (*refs)[j - 1] = (*refs)[j];
        (*refs)[j] = tmp;
        j--;
    }

    return j;
}

void ScalingManager::findReflectionInVector(int refl_id,
                                            ReflectionManager **ref, int *l)
{
    if ((*refs).size() == 0)
    {
        *ref = NULL;
        return;
    }

    int lower = 0;
    int higher = (int)refs->size() - 1;
    int new_bound = (higher + lower) / 2;

    if (((*refs)[higher]).refl_id == refl_id)
    {
        *ref = &(*refs)[higher];
        return;
    }

    if ((refl_id < ((*refs)[lower]).refl_id)
        || (refl_id > ((*refs)[higher]).refl_id))
    {
        *ref = NULL;
        return;
    }

    while (((*refs)[new_bound]).refl_id != refl_id)
    {
        if (higher == lower + 1)
        {
            *ref = NULL;
            return;
        }

        if (((*refs)[new_bound]).refl_id > refl_id)
        {
            higher = new_bound;
        }
        else if (((*refs)[new_bound]).refl_id < refl_id)
        {
            lower = new_bound;
        }

        new_bound = (higher + lower) / 2;
    }

    (*ref) = &(*refs)[new_bound];
    *l = new_bound;
}

void ScalingManager::processReflections(void)
{
    if (mtz_num == 0)
    {
        std::cout << "No images for processing." << std::endl;
        return;
    }

    refs = new vector<ReflectionManager>;

    for (int i = 0; i < mtz_num; i++)
    {
        Gs[i].reflections = new vector<int>;
        Gs[i].multipliers.assign(PARAM_NUM, 1);
    }

    int unique_refl_num = 0;

    for (int i = 0; i < mtzs.size(); i++) // image
    {
        for (int j = 0; j < mtzs[i]->reflectionCount(); j++) // reflection in image
        {
            if (!mtzs[i]->reflection(j)->anyAccepted())
                continue;

            int refl_id = mtzs[i]->reflection(j)->getReflId();

            if (mtzs[i]->reflection(j)->getResolution() > 1 / 2.0)
                        continue;

            int l = 0;
            ReflectionManager *matching_ref;
            findReflectionInVector(refl_id, &matching_ref, &l);

            if (matching_ref != NULL)
            {
                matching_ref->addReflectionForImage(mtzs[i]->reflection(j), &*mtzs[i],
                                                i);
            }
            else
            {
                ReflectionManager *newRef = new ReflectionManager();
                newRef->setGroup(group);
                refs->push_back(*newRef);

                delete newRef;

                ((*refs)[unique_refl_num]).refl_id = refl_id;

                ((*refs)[unique_refl_num]).addReflectionForImage(
                                                             mtzs[i]->reflection(j), &*mtzs[i], i);

                orderFinalReflection();
                unique_refl_num++;
            }
        }
    }

    double totalRefl = 0;
    int count = 0;

    for (int i = 0; i < (*refs).size(); i++)
    {
        if (((*refs)[i]).reflections.size() > 1)
        {
            ((*refs)[i]).sortReflections();

            totalRefl += ((*refs)[i]).reflections.size();
            count++;

            for (int j = 0; j < ((*refs)[i]).reflections.size(); j++)
            {
                int l = ((*refs)[i]).reflections[j].l;
                (*Gs[l].reflections).push_back(i);
                Gs[l].refManagers.push_back(&((*refs)[i]));
            }
        }
    }

    totalRefl /= count;

    std::cout << "Average reflection count: " << totalRefl << std::endl;
}

ScalingManager::ScalingManager(vector<MtzPtr> mtzs)
{
    this->mtzs = mtzs;
    Gs.resize(mtzs.size());
    gradients.resize(mtzs.size());

    for (int i = 0; i < mtzs.size(); i++)
    {
        Gs[i].G = 1;
        Gs[i].image = &*mtzs[i];
        gradients[i] = 0;
        Gs[i].params.resize(PARAM_NUM);
        Gs[i].multipliers.resize(PARAM_NUM);
    }

    if (mtzs.size() == 0)
        return;

    group = mtzs[0]->getLowGroup();

    mtz_num = (int)mtzs.size();
    refs = NULL;
}

ScalingManager::ScalingManager(char **filenames, int filenum, int partiality)
{
    mtzs.resize(filenum);

    Gs.resize(filenum);
    gradients.resize(filenum);

    for (int i = 0; i < filenum; i++)
    {
        Gs[i].G = 1;
        Gs[i].image = (&*mtzs[i]);
        gradients[i] = 0;
    }

    for (int i = 0; i < filenum; i++)
    {
        std::cout << "Loading file " << filenames[i] << std::endl;
        std::string filename_string = filenames[i];

        mtzs[i]->setFilename(filename_string);
        mtzs[i]->loadReflections(partiality);
    }

    mtz_num = filenum;
    refs = NULL;
}

double ScalingManager::residualsForImage(MtzManager *reference,
                                         MtzManager *image, double scale)
{
    double total_residual = 0;

    for (int j = 0; j < image->reflectionCount(); j++)
    {
        int refl_id = (int)image->reflection(j)->getReflId();

        Reflection *mainImage;
        reference->findReflectionWithId(refl_id, &mainImage);

        if (mainImage == NULL)
            continue;

        double residual = image->reflection(j)->meanIntensity() * scale
        - mainImage->meanIntensity();

        residual *= residual;

        residual = sqrt(residual);

        total_residual += residual;
    }

    return total_residual;
}

void ScalingManager::referenceScaling(MtzManager *reference)
{
    for (int i = 0; i < mtz_num; i++)
    {
        double scale = 1;
        double step = 0.5;

        double bestResidual = std::numeric_limits<double>::max();
        double bestScale = scale;

        while (step > 0.001)
        {
            for (double j = scale - step; j <= scale + step; j += step)
            {
                double residual = residualsForImage(reference, &*mtzs[i], j);

                if (residual < bestResidual)
                {
                    bestResidual = residual;
                    bestScale = j;
                }

            }

            if (bestScale == scale)
                step /= 2;

            scale = bestScale;

            printf("Residual = %.2f, scale = %.4f\n", bestResidual, bestScale);
        }

        std::cout << mtzs[i]->getFilename() << " " << 1 / scale << std::endl;
    }
}

double generalGradient(double minX, double maxX, double minY, double maxY)
{
    if (maxY != maxY || minY != minY)
        return 0;

    return (maxY - minY) / (maxX - minX);
}

double ScalingManager::multiplierForParam(int l, int paramNum, double currentR,
                                          double resolution)
{
    double param = Gs[l].params[paramNum];

    double minParam = param + resolution;
    double maxParam = param - resolution;

    MtzManager *image = Gs[l].image;

    Gs[l].params[paramNum] = minParam;
    image->refreshPartialities(&Gs[l].params[0]);

    double rMergeBehind = rMerge();

    Gs[l].params[paramNum] = maxParam;
    image->refreshPartialities(&Gs[l].params[0]);

    double rMergeInFront = rMerge();

    double derivative1 = generalGradient(minParam, param, rMergeBehind,
                                         currentR);
    double derivative2 = generalGradient(param, maxParam, currentR,
                                         rMergeInFront);

    double secondDerivative = generalGradient(minParam, param, derivative1,
                                              derivative2);

    Gs[l].params[paramNum] = param;
    image->refreshPartialities(&Gs[l].params[0]);

    if (secondDerivative == 0)
        std::cout << minParam << "\t" << maxParam << std::endl;

    //  std::cout << "Second derivative: " << rMergeBehind << "\t" << rMergeInFront
    //                  << "\t" << secondDerivative << std::endl;

    if (secondDerivative == 0)
        secondDerivative = 1;

    return fabs(secondDerivative);
}

double ScalingManager::gradientForParam(int l, int paramNum, double currentR)
{
    double param = Gs[l].params[paramNum];

    double maxParam = param + 0.00001;

    MtzManager *image = Gs[l].image;

    Gs[l].params[paramNum] = maxParam;
    image->refreshPartialities(&Gs[l].params[0]);

    double rMergeInFront = rMerge();

    double firstDerivative = generalGradient(param, maxParam, currentR,
                                             rMergeInFront);

    Gs[l].params[paramNum] = param;
    image->refreshPartialities(&Gs[l].params[0]);

    return firstDerivative;
}

double ScalingManager::multiplierForL(int l, double currentR)
{
    double *G = &Gs[l].G;

    double G1 = *G - 0.0001;
    double G2 = *G + 0.0001;
    double G0 = *G;

    *G = G1;
    double rMerge1 = rMerge();

    *G = G2;
    double rMerge2 = rMerge();

    *G = G0;

    double deriv1 = generalGradient(G1, G0, rMerge1, currentR);
    double deriv2 = generalGradient(G0, G2, currentR, rMerge2);

    double final = generalGradient(G1, G0, deriv1, deriv2);

    if (final == 0)
        final = 1;

    return fabs(final);
}

double ScalingManager::gradientForL(int l, double currentR)
{
    double *G = &Gs[l].G;

    double G1 = *G - 0.005;
    double G2 = *G + 0.005;
    double G0 = *G;

    *G = G1;
    double rMerge1 = rMerge(l);

    *G = G2;
    double rMerge2 = rMerge(l);

    double grad = rMerge2 - rMerge1;
    grad /= G2 - G1;

    *G = G0;

    return grad;
}

void ScalingManager::mergedReflections(MtzManager *templateMtz, vector<Reflection *> &reflections, bool half, bool all)
{
    for (int i = 0; i < refs->size(); i++)
    {
        Reflection *reflection = (*refs)[i].mergedReflection(&Gs, half, all);
        reflection->calculateResolution(templateMtz);

        reflections.push_back(reflection);
    }
}

void ScalingManager::rMergeThreaded(ScalingManager *self, int offset, double *numerator, double *denominator, int imageNumber)
{
    double resolution = 2.0;

    int maxThreads = FileParser::getMaxThreads();

    if (imageNumber == -1)
    {
        for (int i = offset; i < self->refs->size(); i += maxThreads)
        {
            (*self->refs)[i].rMerge(&self->Gs, numerator, denominator, resolution);
        }

        return;
    }
}

double ScalingManager::rMerge(int imageNumber)
{
    double resolution = 2.0;

    double numerator = 0;
    double denominator = 0;

    int maxThreads = FileParser::getMaxThreads();

    if (imageNumber >= 0)
    {
        Scale_factor *g = &Gs[imageNumber];

        for (int i = 0; i < g->refManagers.size(); i++)
        {
           g->refManagers[i]->rMerge(&Gs, &numerator, &denominator, resolution);
        }

        double rMerge = numerator / denominator;
        return rMerge;
    }

    vector<double> numerators, denominators;
    numerators.resize(maxThreads);
    denominators.resize(maxThreads);

    boost::thread_group threads;

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(rMergeThreaded, this, i, &numerators[i], &denominators[i], imageNumber);
        threads.add_thread(thr);
    }

    threads.join_all();

    for (int i = 0; i < maxThreads; i++)
    {
        numerator += numerators[i];
        denominator += denominators[i];
    }

    if (denominator == 0)
        return 1000;

    return numerator / denominator;
}

double ScalingManager::rSplit(void)
{
    double numerator = 0;
    double denominator = 0;

    for (int i = 0; i < refs->size(); i++)
    {
        double int1 = 0;
        double int2 = 0;

        (*refs)[i].splitIntensities(&Gs, &int1, &int2);

        if (int1 != int1 || int2 != int2)
            continue;

        double mean = (int1 + int2) / 2;

        numerator += fabs(int1 - int2);
        denominator += mean;
    }

    double rsplit = numerator / (denominator * sqrt(2));

    return rsplit;
}

double ScalingManager::rMergeWrapper(void *object)
{
    return static_cast<ScalingManager *>(object)->rMerge();
}

void ScalingManager::gridSearch(void)
{
    //  outputSimpleScaling(NULL);

    std::cout << "Beginning f(x): " << rMerge() << std::endl;

    for (int num = 0; num < 2; num++)
    {
        for (int i = 0; i < mtzs.size(); i++)
        {
            /// mtzs[i]->gridSearch(rMergeWrapper, this);
        }

        std::cout << "Ending cycle R-merge: " << rMerge() << std::endl;
    }
}
ScalingManager::~ScalingManager(void)
{
 //   std::cout << "Deallocating scaling manager" << std::endl;

    delete refs;

    //  for (int i=0; i < (*refs).size(); i++)
    //          delete &(*refs)[i];
}
