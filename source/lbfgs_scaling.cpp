#include "lbfgs_scaling.h"
#include <vector>
#include <string>
#include <iostream>
#include "stdio.h"
#include "definitions.h"
#include <scitbx/lbfgsb.h>

ScalingManager *scaling;

Lbfgs_Scaling::Lbfgs_Scaling(char **filenames, int filenum)
{
        scaling = new ScalingManager(filenames, filenum, 1);

        (*scaling).processReflections();
}

Lbfgs_Scaling::Lbfgs_Scaling(vector<MtzPtr> mtzs)
{
        scaling = new ScalingManager(mtzs);

        scaling->processReflections();
}

void Lbfgs_Scaling::run(void)
{
        int n = scaling->mtz_num;

    std::ostringstream logged;

        // set up lbfgs algorithm

        double fx = scaling->rMerge();
        std::cout << "R meas after gradient scaling: " << fx << std::endl;

    scitbx::af::shared<int> nbd(n);
    scitbx::af::shared<double> lowerLims(n);
    scitbx::af::shared<double> upperLims(n);

    scitbx::af::shared<double> x = scitbx::af::shared<double>(n);

    logged << "Number of variables: " << x.size() << std::endl;

    for (int i = 0; i < n; i++)
    {
        nbd[i] = 0;
        lowerLims[i] = 0;
        upperLims[i] = 0;
    }

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


    scitbx::af::shared<double> g(n);
    scitbx::af::ref<double> gb(g.begin(), n);

    for (int i = 0; i < n; i++)
        g[i] = 0;

    try
    {
        for (int num = 0;; num++)
        {
            float fx = scaling->rMerge(); // calculate fx

            // gradient values (set to 0 initially, then populate arrays)

            for (int i = 0; i < n; i++)
            {
                g[i] = scaling->gradientForL(i, fx);
            }

            std::cout << "fx: " << fx << std::endl;

            minimizer->process(xb, fx, gb);

            scaling->set_Gs(x);

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


    // after refinement


    for (int i = 0; i < n; i++)
        {
        std::string filename = scaling->mtzs[i]->getFilename();

        //      std::cout << filename << " " << scaling->Gs[i].G << std::endl;

                if (scaling->Gs[i].G > 0.5)
                        scaling->mtzs[i]->applyScaleFactor(scaling->Gs[i].G);
                else
                        scaling->mtzs[i]->setRefCorrelation(0);

        x[i] = 1;
        }

    fx = scaling->rMerge();
    scaling->set_Gs(x);

    std::cout << "R meas after correcting scales: " << fx << std::endl;

    double rsplit = scaling->rSplit();
    std::cout << "Estimated R split: " << rsplit << std::endl;
}

Lbfgs_Scaling::~Lbfgs_Scaling(void)
{
        delete scaling;
}
