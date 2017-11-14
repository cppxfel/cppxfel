#ifndef statistics
#define statistics


#include <string>
#include <iostream>
#include "MtzManager.h"

struct Partial
{
    MillerPtr miller;
        double partiality;
        double percentage;
        double resolution;
        double wavelength;
};

class StatisticsManager
{
private:

public:
        StatisticsManager(void);
        ~StatisticsManager(void);

    vector<MtzPtr> mtzs;
        void loadFiles(char **filenames, int filenum, int partiality);
    void generate_cc_grid();
    void ccGridThreaded(int offset, int calculationsPerThread, std::map<int, int> *histogram, int histogramCount, double slice, int *num_cc, int *num_inv_cc);
    static void ccGridThreadedWrapper(StatisticsManager *object, int offset, int calculationsPerThread, std::map<int, int> *histogram, int histogramCount, double slice, int *num_cc, int *num_inv_cc);
    double gridCorrelation(int imageNumI, int imageNumJ);

        void printGradientsAgainstRef(MtzManager *reference);

        void setMtzs(vector<MtzPtr> mtzs);

        void twoImagePartialityStats(int num);
    void partialityStats(int num, double threshold = -100, int h = 0, int k = 0, int l = 0);
        static void twoImagePartialityStatsWritten(vector<Partial> *partials,
                        MtzManager **image, MtzManager **test_image);
        void write_refls(int num);

        static double midPointBetweenResolutions(double minD, double maxD);
        static void generateResolutionBins(double minD, double maxD, int binCount,
                        vector<double> *bins);
        static void convertResolutions(double lowAngstroms, double highAngstroms,
                        double *lowReciprocal, double *highReciprocal);

        double cc_through_origin(int num1, int num2, int silent, int inverted,
                        int *hits);
        double cc_through_origin(int num1, int num2, int silent, int inverted,
                        int *hits, double lowResolution, double highResolution, bool log);

        static double cc_pearson(MtzManager *shot1, MtzManager *shot2, int silent,
            int *hits, double *multiplicity, double lowResolution,
                        double highResolution, bool log = false, bool freeOnly = false);
        double cc_pearson(int num1, int num2, int silent, int *hits,
                        double *multiplicity, double lowResolution = 0, double highResolution =
                                        0, bool log = false, bool freeOnly = false);

        static double r_factor(RFactorType rFactor, MtzManager *shot1, int *hits,
                        double *multiplicity, double lowResolution, double highResolution, bool freeOnly = false);

        static double r_split(MtzManager *shot1, MtzManager *shot2, int silent,
                        int *hits, double *multiplicity, double lowResolution,
                        double highResolution, bool log, bool freeOnly = false);

        double **cc_array;
        double **inv_cc_array;
        int **hits;
        int **inv_hits;
        vector<vector<double> > gradient_array;

        int mtz_num;
};

#endif // statistics
