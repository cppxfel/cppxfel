//
//  AmbiguityBreaker.cpp
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "AmbiguityBreaker.h"
#include <scitbx/lbfgsb.h>
#include "StatisticsManager.h"
#include "FileParser.h"
#include <vector>

double AmbiguityBreaker::dotProduct(int imageNumI, int imageNumJ)
{
    int arrayStartI = imageNumI * ambiguityCount;
    int arrayStartJ = imageNumJ * ambiguityCount;
    
    double fx = 0;
    
    for (int i = 0; i < ambiguityCount; i++)
    {
        double a = x[arrayStartI + i];
        double b = x[arrayStartJ + i];
        
        double addition = a * b;
        
        fx += addition;
    }
    
    return fx;
}

double AmbiguityBreaker::evaluation()
{
    double fx = 0;
    
    int imageCount = (int)mtzs.size();
    
    for (int i = 0; i < imageCount - 1; i++)
    {
        for (int j = i + 1; j < imageCount; j++)
        {
            double correl = gridCorrelation(i, j);

            if (correl == -1)
                continue;
            
            double dot = dotProduct(i, j);
            
            fx += pow(correl - dot, 2);
        }
    }
    
    return fx;
}

double AmbiguityBreaker::gradientForImage(int imageNum, int axis)
{
    double g = 0;
    
    for (int j = 0; j < mtzs.size(); j++)
    {
        if (j == imageNum)
            continue;
        
        double correl = gridCorrelation(imageNum, j);
        
        if (correl == -1)
            continue;

        int arrayStartJ = j * ambiguityCount;
        
        double dot = dotProduct(imageNum, j);
        
        double jAxis = x[arrayStartJ + axis];
        
        g += (-2) * jAxis * (correl - dot);
    }
    
    return g;
}

double AmbiguityBreaker::gridCorrelation(int imageNumI, int imageNumJ)
{
    return statsManager->gridCorrelation(imageNumI, imageNumJ);
}

void AmbiguityBreaker::makeCorrelationGrid()
{
    statsManager = new StatisticsManager();
    statsManager->setMtzs(mtzs);
    
    statsManager->generate_cc_grid();
}

// Call constructor and then run()

AmbiguityBreaker::AmbiguityBreaker(vector<MtzPtr> newMtzs)
{
    setMtzs(newMtzs);
}

void AmbiguityBreaker::setMtzs(vector<MtzPtr> newMtzs)
{
    mtzs = newMtzs;
    
    if (!mtzs.empty())
    {
        ambiguityCount = mtzs[0]->ambiguityCount();
    }
    
    
    bool shouldApplyUnrefinedPartiality = FileParser::getKey("APPLY_UNREFINED_PARTIALITY", false);
    
    if (shouldApplyUnrefinedPartiality)
    {
        for (int i = 0; i < mtzs.size(); i++)
        {
            mtzs[i]->applyUnrefinedPartiality();
        }
    }
    
    
    if (ambiguityCount > 1)
        makeCorrelationGrid();
}


void AmbiguityBreaker::run()
{
    if (ambiguityCount > 1)
    {
        breakAmbiguity();
        printResults();
        split();
    }
    merge();
}

void AmbiguityBreaker::split()
{
    if (ambiguityCount == 2)
    {
        for (int i = 0; i < mtzs.size(); i++)
        {
            int startI = i * ambiguityCount;
            
            if (x[startI] > x[startI + 1])
            {
                mtzs[i]->setActiveAmbiguity(1);
            }
            else
            {
                mtzs[i]->setActiveAmbiguity(0);
            }
        }
    }
    else if (ambiguityCount > 2)
    {
        Logger::mainLogger->addString("WARNING! cppxfel has not been configured to split a set of images with more than two possible indexing solutions. The resulting structure will be twinned. Instead, feel free to provide a reference MTZ under the \"INITIAL_MTZ\" keyword which will avoid the requirement to break the indexing ambiguity. (Sorry)");
    }
}

void AmbiguityBreaker::merge()
{
    for (int i = 0; i < mtzs.size(); i++)
    {
    //    mtzs[i]->applyUnrefinedPartiality();
    }
    
    MtzGrouper *idxGrouper = new MtzGrouper();
    idxGrouper->setWeighting(WeightTypeAverage);
    idxGrouper->setExcludeWorst(false);
    idxGrouper->setMtzManagers(mtzs);
    idxGrouper->merge(&merged);
    delete idxGrouper;
}

void AmbiguityBreaker::printResults()
{
    std::ostringstream logged;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        int start = i * ambiguityCount;
        
        std::string filename = mtzs[i]->getFilename();
        logged << filename << "\t";
        
        for (int j = 0; j < ambiguityCount; j++)
        {
            logged << x[start + j] << "\t";
        }
        
        logged << std::endl;
    }
    
    Logger::mainLogger->addStream(&logged);
}

void AmbiguityBreaker::breakAmbiguity()
{
    int n = ambiguityCount * (int)mtzs.size();

    std::ostringstream logged;
    
    scitbx::af::shared<int> nbd(n);
    scitbx::af::shared<double> lowerLims(n);
    scitbx::af::shared<double> upperLims(n);
    
    x = scitbx::af::shared<double>(n);
    
    logged << "Number of variables: " << x.size() << std::endl;
    
    for (int i = 0; i < n; i++)
        nbd[i] = 0;
    
    double factr = 1.e+7;
    double pgtol = 0;
    int iprint = 0;
    
    scitbx::lbfgsb::minimizer<double> *minimizer =
    new scitbx::lbfgsb::minimizer<double>(n, 5, lowerLims, upperLims,
                                          nbd, false, factr, pgtol, iprint);
    
    scitbx::lbfgs::traditional_convergence_test<double, int> is_converged(n);
    
    scitbx::af::ref<double> xb(x.begin(), n);
    
    for (int i = 0; i < n; i++)
        x[i] = 1;
    
    for (int i = 0; i < x.size(); i++)
    {
        double random = (double)rand() / RAND_MAX;
        x[i] = random;
    }
    
    Logger::mainLogger->addStream(&logged);
    
    scitbx::af::shared<double> g(n);
    scitbx::af::ref<double> gb(g.begin(), n);
    
    for (int i = 0; i < n; i++)
        g[i] = 0;
    
    try
    {
        for (int num = 0;; num++)
        {
            float fx = evaluation(); // calculate fx
            
            // gradient values (set to 0 initially, then populate arrays)
            
            for (int i = 0; i < n; i++)
            {
                g[i] = 0;
            }
            
            for (int i = 0; i < mtzs.size(); i++)
            {
                for (int j = 0; j < ambiguityCount; j++)
                {
                    g[i * ambiguityCount + j] = gradientForImage(i, j);
                }
            }
            
            std::cout << "fx: " << fx << std::endl;
            
            minimizer->process(xb, fx, gb);
            
            if (minimizer->is_terminated())
                break;
            
            if (minimizer->n_iteration() > 50)
                break;
        }
    } catch (scitbx::lbfgs::error &e)
    {
        std::cout << e.what() << std::endl;
    } catch (void *e)
    {
        std::cout << "Unknown error!" << std::endl;
    }
}