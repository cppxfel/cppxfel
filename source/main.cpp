//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <string>

#include "InputFileParser.h"
#include "MtzManager.h"
#include "StatisticsManager.h"
#include "misc.h"
#include "GraphDrawer.h"
#include "Wiki.h"
#include "Logger.h"
#include <fstream>

void finishJobNotification(int argc, char *argv[], int minutes)
{
    const char *jobNotificationFileStr = std::getenv("JOB_NOTIFICATION_FILE");
    
    if (jobNotificationFileStr == NULL)
    {
        return;
    }
    
    std::ostringstream command;
    command << "cppxfel.run ";
    for (int i = 1; i < argc; i++)
    {
        command << argv[i] << " ";
    }
    
    std::ostringstream notification;
    notification << "osascript -e 'display notification \"" << command.str() << "\" with title \"Job finished\" subtitle \"" << minutes << " minutes to complete\" sound name \"Glass\"'" << std::endl;
    
    std::ofstream jobNotificationFile;
    jobNotificationFile.open(jobNotificationFileStr, std::ofstream::out | std::ofstream::app);
    jobNotificationFile << notification.str();
    jobNotificationFile.close();
    
    std::cout << "Job notification posted." << std::endl;
}

void new_main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    new_main(argc, argv);
}

void new_main(int argc, char *argv[])
{
    time_t startcputime;
    time(&startcputime);
    
	if (argc == 1)
	{
        std::cout << "Welcome to cppxfel version 1.1!" << std::endl;
        std::cout << "Please refer to & cite paper in Journal of Applied Crystallography (unpublished)" << std::endl << std::endl;
        std::cout << "Command order for regular structure solution:" << std::endl;
        std::cout << "\tcppxfel.run_dials shot*.pickle" << std::endl;
        std::cout << "\tcppxfel.input_gen" << std::endl;
        std::cout << "\tcppxfel.run -i integrate.txt" << std::endl;
        std::cout << "\tcppxfel.run -i refine.txt" << std::endl;
        std::cout << "\tcppxfel.run -i merge.txt" << std::endl << std::endl;;
        std::cout << "Other functions for assessing data quality:" << std::endl << std::endl;
        
        std::cout << "Correlation between two MTZ files:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -cc firstFile.mtz secondFile.mtz [ambiguity] [lowRes] [highRes] [bins]" << std::endl << std::endl;
        std::cout << "ambiguity: 0, 1, 2 or 3 - will compare different indexing solutions where the Bravais lattice symmetry is higher than that of the point group for certain space groups. Default 0" << std::endl;
        std::cout << "lowRes and highRes: set to resolution in Angstroms to bound the results, or set to 0 to take lowest/highest resolution data. Default 0, 0" << std::endl;
        std::cout << "bins: number of bins to report correlation statistics. Default 20." << std::endl << std::endl;
        
        std::cout << "Partiality CSV files:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -partiality reference.mtz ref-img-shot-number.mtz [highRes]" << std::endl << std::endl;
        std::cout << "highRes: 0, 1, 2 or 3 - highest resolution reflection to report results on. Default 1.4" << std::endl;
        std::cout << "This outputs partiality_[m].csv where m is bin number, which can be imported into other graphing softwares such as R." << std::endl << std::endl;;

        
        std::cout << "Merging statistics:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -rpim unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl;
        std::cout << "\tcppxfel.run -rmeas unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl;
        std::cout << "\tcppxfel.run -rmerge unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl << std::endl;
        std::cout << "lowRes and highRes: set to resolution in Angstroms to bound the results, or set to 0 to take lowest/highest resolution data. Default 0, 0" << std::endl;
        std::cout << "bins: number of bins to report correlation statistics. Default 20." << std::endl << std::endl;
        
        exit(1);
	}

	if (strcmp(argv[1], "-wiki") == 0)
	{
		if (argc <= 2)
		{
			std::cout << "arguments: -wiki <logfile>" << std::endl;
			exit(1);
		}

        Wiki wiki = Wiki(std::string(argv[2]));
		wiki.process();

		exit(1);
	}
    
	std::cout << "Welcome to cppxfel!" << std::endl;
    
	if (strcmp(argv[1], "-i") == 0)
	{
		if (argc < 3)
		{
			std::cout << "arguments: -i <input_script>" << std::endl;
			exit(1);
		}

        std::vector<std::string> extras;
        
        for (int i = 3; i < argc; i++)
        {
            extras.push_back(std::string(argv[i]));
        }
        
        InputFileParser *parser = new InputFileParser(std::string(argv[2]), extras);
		parser->parse(false);
        
        delete parser;
	}
    else
    {
        Logger::mainLogger = LoggerPtr(new Logger());
        boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);
    }
    
    if (strcmp(argv[1], "-b") == 0)
    {
        float bFactor = atof(argv[2]);
        
        MtzManager *mtz1 = new MtzManager();
        mtz1->setFilename(std::string(argv[3]));
        mtz1->loadReflections(0);
        
        mtz1->applyBFactor(bFactor);
        
        mtz1->writeToFile("b-" + std::string(argv[3]));

    }
    
	if (strcmp(argv[1], "-rmerge") == 0 || strcmp(argv[1], "-rpim") == 0
			|| strcmp(argv[1], "-rmeas") == 0)
	{
		if (argc < 3)
		{
			std::cout << "arguments: -r{merge} <file1> [lowRes] [highRes] [bins]."
					<< std::endl;
			exit(1);
		}

		RFactorType rFactor = RFactorTypeMerge;

		if (strcmp(argv[1], "-rpim") == 0)
		{
			rFactor = RFactorTypePim;
		}
		else if (strcmp(argv[1], "-rmeas") == 0)
		{
			rFactor = RFactorTypeMeas;
		}

		double highRes = 0;
		double lowRes = 0;
		int bins = 20;

		if (argc > 5)
		{
			lowRes = atof(argv[4]);
			highRes = atof(argv[5]);
		}
		if (argc > 6)
		{
			bins = atoi(argv[6]);
		}

        MtzManager *mtz = new MtzManager();
		mtz->setFilename(std::string(argv[2]));
		mtz->loadReflections(1);

		mtz->rFactorWithManager(rFactor, false, false, lowRes, highRes, bins);
	}

	if (strcmp(argv[1], "-cc") == 0 || strcmp(argv[1], "-rsplit") == 0)
	{
		if (argc < 4)
		{
			std::cout
					<< "arguments: -cc <file1> <file2> (0 / 1) [lowRes] [highRes] [bins]."
					<< std::endl;
			exit(1);
		}

		bool rsplit = (strcmp(argv[1], "-rsplit") == 0);

		int inverted = 0;
		double highRes = 0;
		double lowRes = 0;
		int bins = 20;

		if (argc > 4)
		{
			inverted = atoi(argv[4]);
		}
		if (argc > 6)
		{
			lowRes = atof(argv[5]);
			highRes = atof(argv[6]);
		}
		if (argc > 7)
		{
			bins = atoi(argv[7]);
		}

		MtzManager *mtz1 = new MtzManager();
		mtz1->setFilename(std::string(argv[2]));
		mtz1->loadReflections(1);

		MtzManager *mtz2 = new MtzManager();
		mtz2->setFilename(std::string(argv[3]));
		mtz2->loadReflections(PartialityModelScaled, true);

		if (inverted)
            mtz1->setActiveAmbiguity(1);

		if (rsplit)
		{
			mtz1->rSplitWithManager(mtz2, 1, 0, lowRes, highRes, bins);
		}
		else
		{
			mtz1->correlationWithManager(mtz2, 1, 0, lowRes, highRes, bins);

		}

		delete mtz1;
		delete mtz2;

		exit(1);
	}

	if (strcmp(argv[1], "-gradscaling") == 0)
	{
		if (argc <= 3)
		{
			std::cout << "arguments: -gradscaling <ref> <file2> ... <filen>."
					<< std::endl;
			exit(1);
		}

		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *image = new MtzManager();
			image->setFilename(argv[i]);
			image->loadReflections(1);

			double gradientOld = image->gradientAgainstManager(reference);
			double gradientBefore = image->minimizeRFactor(reference);

			std::cout << image->getFilename() << "\t" << gradientOld << "\t"
					<< gradientBefore << std::endl;
		}
	}

	if (strcmp(argv[1], "-inv") == 0)
	{
		if (argc <= 2)
		{
			std::cout << "arguments: -inv <file>." << std::endl;
			exit(1);
		}
        
        MtzManager *mtz = new MtzManager();
		mtz->setFilename(std::string(argv[2]));
		mtz->loadReflections(1);
        mtz->setActiveAmbiguity(1);

		mtz->writeToFile(std::string("inv-") + argv[2], false, false, true);
	}

	if (strcmp(argv[1], "-stats") == 0)
	{
		if (argc < 3)
		{
			std::cout << "arguments: -stats <filein> <threshold>." << std::endl;
			exit(1);
		}

		if (argc >= 3)
		{
			StatisticsManager stats;
			stats.loadFiles(&argv[2], 1, 0);
            
            double threshold = -100;
            int h = 0; int k = 0; int l = 0;
            
            if (argc >= 4)
            {
                threshold = atof(argv[3]);
            }
            
            if (argc >= 7)
            {
                h = atoi(argv[4]);
                k = atoi(argv[5]);
                l = atoi(argv[6]);
            }

            GraphDrawer drawer = GraphDrawer(&*stats.mtzs[0]);
            drawer.plotPartialityStats();
		}
	}
    
    if (strcmp(argv[1], "-intensities") == 0)
    {
        if (argc <= 5)
        {
            std::cout << "arguments: -intensities h k l <file1> {<file2> ...}." << std::endl;
            exit(1);
        }
        
        int h = atoi(argv[2]);
        int k = atoi(argv[3]);
        int l = atoi(argv[4]);
        
        if (h == 0 && k == 0 && l == 0)
        {
            h = (double)(rand() / RAND_MAX) * 40;
            k = (double)rand() / RAND_MAX * (40 - h) + h;
            l = (double)rand() / RAND_MAX * (40 - k) + k;
        }
        
        std::cout << "Plotting intensities for (" << h << ", " << k << ", " << l << ")" << std::endl;
        
        
        std::vector<MtzPtr> mtzs;
        
        for (int i = 5; i < argc; i++)
        {
            MtzPtr mtz = MtzPtr(new MtzManager());
            mtz->setFilename(argv[i]);
            mtz->loadReflections(1);
            
            mtzs.push_back(mtz);
        }
        
        GraphDrawer drawer = GraphDrawer(&*mtzs[0]);
        drawer.plotReflectionFromMtzs(mtzs, h, k, l);
    }

#ifdef MAC
    if (strcmp(argv[1], "-partiality") == 0)
	{
		if (argc <= 3)
		{
			std::cout << "arguments: -partiality <ref> <filein> {<maxres>}." << std::endl;
			exit(1);
		}

		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);
		MtzManager::setReference(reference);
        
        double maxRes = 1.6;
        
        if (argc > 4)
        {
            maxRes = atof(argv[4]);
        }
        
        
        MtzManager *mtz = new MtzManager();
        mtz->setFilename(argv[3]);
        std::cout << "Partiality plot for " << argv[3] << std::endl;
        mtz->loadReflections(PartialityModelNone, true);
        
        vector<Reflection *>refReflections, imageReflections;
        
        GraphDrawer graph = GraphDrawer(mtz);
        
        graph.partialityPlot("partiality", GraphMap(), maxRes);
        
        delete mtz;

		delete reference;
	}
    
    if (strcmp(argv[1], "-scale") == 0)
    {
        if (argc <= 3)
        {
            std::cout << "arguments: -scale <ref> <file1>." << std::endl;
            exit(1);
        }
        
        MtzManager *reference = new MtzManager();
        reference->setFilename(argv[2]);
        reference->loadReflections(1);
        MtzManager::setReference(reference);
        
        MtzManager *mtz = new MtzManager();
        mtz->setFilename(argv[3]);
        mtz->loadReflections(1);
        
        mtz->applyScaleFactorsForBins(50);
        mtz->writeToFile("scaled-" + std::string(argv[3]));
    }

	if (strcmp(argv[1], "-bfactor") == 0)
	{
		if (argc <= 2)
		{
			std::cout << "arguments: -bfactor <ref> <file1> {<file2> ...}." << std::endl;
			exit(1);
		}

		vector<MtzManager *> managers = vector<MtzManager *>();

		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);
		MtzManager::setReference(reference);
        FileParser::setKey("REFINE_B_FACTOR", true);
        
		for (int i = 3; i < argc; i++)
		{
			MtzManager *mtz = new MtzManager();
			mtz->setFilename(argv[i]);
			mtz->loadReflections(1);
            managers.push_back(mtz);
		}

		GraphDrawer graph = GraphDrawer(reference);

		graph.resolutionStatsPlot(managers, "intensity_bins", GraphMap(), true,
				false);
		graph.resolutionStatsPlot(managers, "intensity_bins_2", GraphMap(),
				true, true);
		graph.resolutionStatsPlot(managers);
        graph.bFactorPlot(managers);
        
		for (int i = 0; i < managers.size(); i++)
			delete managers[i];

	}

	if (strcmp(argv[1], "-ccplot") == 0)
	{
		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);
		MtzManager::setReference(reference);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *image = new MtzManager();
			image->setFilename(argv[i]);
			image->loadReflections(1);

			GraphDrawer graph = GraphDrawer(image);
			graph.correlationPlot("correl");

			delete image;
		}

		delete reference;
	}
#endif
    
    time_t endcputime;
    time(&endcputime);
    
    clock_t difference = endcputime - startcputime;
    double seconds = difference;
    
    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;
    
    std::ostringstream logged;
    logged << "N: Total time: " << minutes << " minutes, "
    << finalSeconds << " seconds (" << seconds << " seconds)." << std::endl;

	logged << "Done" << std::endl;
    Logger::mainLogger->addStream(&logged);
    
    if (strcmp(argv[1], "-i") == 0)
        finishJobNotification(argc, argv, minutes);
    
    sleep(2);
}
