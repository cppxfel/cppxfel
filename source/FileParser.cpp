/*
 * FileParser.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "FileParser.h"

#include "FileReader.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <locale>
#include <stdio.h>

#define FILE_PARSER_CPP_

#include "MtzRefiner.h"

ParametersMap FileParser::parameters;
std::ostringstream FileParser::log;
int FileParser::threadsFound = 0;
char FileParser::splitCharMajor = ' ';
char FileParser::splitCharMinor = ' ';

int FileParser::getMaxThreads()
{
    if (threadsFound != 0)
        return threadsFound;

    std::ostringstream logged;

    char *nslots;
    nslots = getenv("NSLOTS");
    int maxThreads = 0;

    if (nslots != NULL)
    {
        maxThreads = atoi(nslots);
        logged << "Using environment variable NSLOTS: " << maxThreads << " threads." << std::endl;
    }

    if (maxThreads == 0)
    {
        maxThreads = getKey("MAX_THREADS", MAX_THREADS);

        if (hasKey("MAX_THREADS"))
        {
            logged << "Getting number of threads from user input: " << maxThreads << " threads." << std::endl;
        }
        else
        {
            logged << "Using default number of threads: " << maxThreads << " threads." << std::endl;
        }
    }

    Logger::mainLogger->addStream(&logged);

    threadsFound = maxThreads;

    return maxThreads;
}

bool FileParser::hasKey(std::string key)
{
        return (parameters.count(key) > 0);
}


void FileParser::simpleFloat(ParametersMap *map, std::string command,
                std::string rest)
{
        double theFloat = atof(rest.c_str());

        log << "Setting double " << command << " to " << theFloat << std::endl;

        (*map)[command] = theFloat;
}

void FileParser::simpleBool(ParametersMap *map, std::string command,
                std::string rest)
{
        bool on = (rest == "ON" || rest == "on" ? 1 : 0);

        log << "Setting bool " << command << " to " << on << std::endl;

        (*map)[command] = on;
}

void FileParser::simpleString(ParametersMap *map, std::string command,
                std::string rest)
{
        (*map)[command] = rest;

        log << "Setting string " << command << " to " << rest << std::endl;

}

void FileParser::simpleInt(ParametersMap *map, std::string command,
                std::string rest)
{
        int theInt = atoi(rest.c_str());

        log << "Setting int " << command << " to " << theInt << std::endl;

        (*map)[command] = theInt;
}

void FileParser::doubleVector(ParametersMap *map, std::string command,
                              std::string rest)
{
        vector<std::string> components = FileReader::split(rest, splitCharMinor);
        vector<double> doubleVector;

        log << "Setting " << command << " to ";

        for (int i = 0; i < components.size(); i++)
        {
                double theFloat = atof(components[i].c_str());
                log << theFloat << " ";
                doubleVector.push_back(theFloat);
        }

        log << std::endl;

        (*map)[command] = doubleVector;
}

void FileParser::intVector(ParametersMap *map, std::string command,
                std::string rest)
{
        vector<std::string> components = FileReader::split(rest, splitCharMinor);
        vector<int> intVector;

        log << "Setting " << command << " to ";

        for (int i = 0; i < components.size(); i++)
        {
                int theInt = atoi(components[i].c_str());
                log << theInt << " ";
                intVector.push_back(theInt);
        }

        log << std::endl;

        (*map)[command] = intVector;
}

void FileParser::generateFunctionList()
{
        parserMap = ParserMap();

    parserMap["VERBOSITY_LEVEL"] = simpleInt;
    parserMap["MAX_THREADS"] = simpleInt;

        // Refinement parameters
        parserMap["REMOVE_WEDGE"] = simpleFloat;

    parserMap["MINIMUM_CYCLES"] = simpleInt;
    parserMap["MAXIMUM_CYCLES"] = simpleInt;
    parserMap["STOP_REFINEMENT"] = simpleBool;

    parserMap["SET_SIGMA_TO_UNITY"] = simpleBool;
    parserMap["APPLY_UNREFINED_PARTIALITY"] = simpleBool;
    parserMap["BINARY_PARTIALITY"] = simpleBool;
    parserMap["MINIMIZATION_METHOD"] = simpleInt;
    parserMap["NELDER_MEAD_CYCLES"] = simpleInt;
    parserMap["MEDIAN_WAVELENGTH"] = simpleBool;
    parserMap["WAVELENGTH_RANGE"] = doubleVector;
    parserMap["LANDSCAPE_DIVISIONS"] = simpleInt;
    parserMap["EXCLUSION_BY_CC_HALF"] = simpleBool;
    parserMap["ACCEPTABLE_UNIT_CELL_TOLERANCE"] = simpleFloat;
    parserMap["ALLOW_TRUST"] = simpleBool;
    parserMap["EXCLUDE_OWN_REFLECTIONS"] = simpleBool;
    parserMap["PARTIALITY_CUTOFF"] = simpleFloat;
        parserMap["DEFAULT_TARGET_FUNCTION"] = simpleInt;
    parserMap["TARGET_FUNCTIONS"] = intVector;
    parserMap["USE_PARTIALITY_FUNCTION"] = simpleBool;
    parserMap["RLP_MODEL"] = simpleInt;
        parserMap["CORRELATION_THRESHOLD"] = simpleFloat;
    parserMap["PARTIALITY_CORRELATION_THRESHOLD"] = simpleFloat;
        parserMap["MAX_RESOLUTION_ALL"] = simpleFloat;
        parserMap["MAX_RESOLUTION_RLP_SIZE"] = simpleFloat;
    parserMap["MIN_REFINED_RESOLUTION"] = simpleFloat;
    parserMap["INITIAL_CORRELATION_THRESHOLD"] = simpleFloat;
        parserMap["THRESHOLD_SWAP"] = simpleInt;
        parserMap["OUTLIER_REJECTION_SIGMA"] = simpleFloat;
        parserMap["OUTLIER_REJECTION"] = simpleBool;
        parserMap["CORRELATION_REJECTION"] = simpleBool;
        parserMap["PARTIALITY_REJECTION"] = simpleBool;
    parserMap["POLARISATION_CORRECTION"] = simpleBool;
    parserMap["POLARISATION_FACTOR"] = simpleFloat;
    parserMap["REFINEMENT_INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["TRUST_INDEXING_SOLUTION"] = simpleBool;
    parserMap["REFINE_B_FACTOR"] = simpleBool;
    parserMap["INITIAL_GRID_SEARCH"] = simpleBool;
    parserMap["MASS_SCALING"] = simpleBool;
    parserMap["R_FACTOR_THRESHOLD"] = simpleFloat;
    parserMap["REINITIALISE_WAVELENGTH"] = simpleBool;
    parserMap["PENALTY_WEIGHT"] = simpleFloat;
    parserMap["PENALTY_RESOLUTION"] = simpleFloat;
    parserMap["PARTIALITY_SLICES"] = simpleInt;
    parserMap["MAX_SLICES"] = simpleInt;
    parserMap["CAREFUL_RESOLUTION"] = simpleFloat;
    parserMap["SMOOTH_FUNCTION"] = simpleBool;
    parserMap["NORMALISE_PARTIALITIES"] = simpleBool;
    parserMap["REPLACE_REFERENCE"] = simpleBool;
    parserMap["REFINE_ENERGY_SPECTRUM"] = simpleBool;
    parserMap["FREE_MILLER_LIST"] = simpleString;
    parserMap["FREE_MILLER_PROPORTION"] = simpleFloat;

        parserMap["INITIAL_WAVELENGTH"] = simpleFloat;
        parserMap["INITIAL_BANDWIDTH"] = simpleFloat;
        parserMap["INITIAL_MOSAICITY"] = simpleFloat;
        parserMap["INITIAL_EXPONENT"] = simpleFloat;
        parserMap["INITIAL_RLP_SIZE"] = simpleFloat;

        parserMap["STEP_SIZE_WAVELENGTH"] = simpleFloat;
        parserMap["STEP_SIZE_BANDWIDTH"] = simpleFloat;
        parserMap["STEP_SIZE_MOSAICITY"] = simpleFloat;
        parserMap["STEP_SIZE_EXPONENT"] = simpleFloat;
        parserMap["STEP_SIZE_ORIENTATION"] = simpleFloat;
    parserMap["STEP_SIZE_ORIENTATION_ABC"] = simpleFloat;
    parserMap["STEP_SIZE_RLP_SIZE"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_A"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_B"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_C"] = simpleFloat;

        parserMap["TOLERANCE_WAVELENGTH"] = simpleFloat;
        parserMap["TOLERANCE_BANDWIDTH"] = simpleFloat;
        parserMap["TOLERANCE_MOSAICITY"] = simpleFloat;
        parserMap["TOLERANCE_EXPONENT"] = simpleFloat;
        parserMap["TOLERANCE_ORIENTATION"] = simpleFloat;
        parserMap["TOLERANCE_RLP_SIZE"] = simpleFloat;

        parserMap["OPTIMISING_WAVELENGTH"] = simpleBool;
        parserMap["OPTIMISING_BANDWIDTH"] = simpleBool;
        parserMap["OPTIMISING_MOSAICITY"] = simpleBool;
        parserMap["OPTIMISING_EXPONENT"] = simpleBool;
        parserMap["OPTIMISING_ORIENTATION"] = simpleBool;
        parserMap["OPTIMISING_RLP_SIZE"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_A"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_B"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_C"] = simpleBool;

        parserMap["ORIENTATION_MATRIX_LIST"] = simpleString;
    parserMap["SECOND_MATRIX_LIST"] = simpleString;
    parserMap["MATRIX_LIST_VERSION"] = simpleFloat;
        parserMap["INITIAL_MTZ"] = simpleString;
    parserMap["IMAGE_LIMIT"] = simpleInt;
    parserMap["IMAGE_SKIP"] = simpleInt;
    parserMap["NEW_MATRIX_LIST"] = simpleString;

    parserMap["RECALCULATE_SIGMA"] = simpleBool;
    parserMap["MERGE_ANOMALOUS"] = simpleBool;
    parserMap["FAKE_ANOMALOUS"] = simpleBool;
    parserMap["SCALING_STRATEGY"] = simpleInt;
        parserMap["MINIMUM_REFLECTION_CUTOFF"] = simpleInt;
    parserMap["APPLY_INFLATION"] = simpleBool;
    parserMap["MINIMUM_MULTIPLICITY"] = simpleInt;
    parserMap["THREADED_MERGE"] = simpleBool;
    parserMap["FAST_MERGE"] = simpleBool;
    parserMap["INDEX_STARTING_DATA"] = doubleVector;

        // Indexing parameters

    parserMap["DETECTOR_GAIN"] = simpleFloat;
    parserMap["BITS_PER_PIXEL"] = simpleInt;
        parserMap["SPACE_GROUP"] = simpleInt;
        parserMap["INTEGRATION_WAVELENGTH"] = simpleFloat;
        parserMap["DETECTOR_DISTANCE"] = simpleFloat;
    parserMap["BEAM_CENTRE"] = doubleVector;
    parserMap["MM_PER_PIXEL"] = simpleFloat;
    parserMap["DETECTOR_SIZE"] = doubleVector;
        parserMap["OVER_PRED_BANDWIDTH"] = simpleFloat;
        parserMap["OVER_PRED_RLP_SIZE"] = simpleFloat;
    parserMap["REFINE_ORIENTATIONS"] = simpleBool;
    parserMap["REFINE_DISTANCES"] = simpleBool;
    parserMap["INDEXING_ORIENTATION_TOLERANCE"] = simpleFloat;
    parserMap["LOW_INTENSITY_PENALTY"] = simpleBool;
        parserMap["INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["ABSOLUTE_INTENSITY"] = simpleBool;
        parserMap["METROLOGY_SEARCH_SIZE"] = simpleInt;
    parserMap["METROLOGY_MOVE_THRESHOLD"] = simpleFloat;
    parserMap["FOCUS_ON_PEAK_SIZE"] = simpleInt;
        parserMap["SHOEBOX_FOREGROUND_PADDING"] = simpleInt;
        parserMap["SHOEBOX_NEITHER_PADDING"] = simpleInt;
        parserMap["SHOEBOX_BACKGROUND_PADDING"] = simpleInt;
    parserMap["SHOEBOX_MAKE_EVEN"] = simpleBool;
    parserMap["COMPLEX_SHOEBOX"] = simpleBool;
    parserMap["MIN_INTEGRATED_RESOLUTION"] = simpleFloat;
    parserMap["MAX_INTEGRATED_RESOLUTION"] = simpleFloat;
        parserMap["UNIT_CELL"] = doubleVector;
    parserMap["FIX_UNIT_CELL"] = simpleBool;
    parserMap["ADD_MASK"] = intVector;
    parserMap["INITIAL_ORIENTATION_STEP"] = simpleFloat;
    parserMap["SHOEBOX_BANDWIDTH_MULTIPLIER"] = simpleFloat;
    parserMap["PIXEL_LEAK"] = simpleFloat;
    parserMap["ORIENTATION_SCORE"] = simpleInt;
    parserMap["ORIENTATION_CORRECTION"] = doubleVector;
    parserMap["IMAGE_MASKED_VALUE"] = simpleInt;
    parserMap["IMAGE_IGNORE_UNDER_VALUE"] = simpleInt;
    parserMap["SPHERE_THICKNESS"] = simpleFloat;
    parserMap["SIGMA_RESOLUTION_CUTOFF"] = simpleFloat;
    parserMap["PIXEL_COUNT_CUTOFF"] = simpleInt;
    parserMap["EXPECTED_SPOTS"] = simpleInt;
    parserMap["INDEXING_SLICE_ANGLE"] = simpleFloat;
    parserMap["REFINE_UNIT_CELL_A"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_B"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_C"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_ALPHA"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_BETA"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_GAMMA"] = simpleBool;
    parserMap["STEP_UNIT_CELL_A"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_B"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_C"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_ALPHA"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_BETA"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_GAMMA"] = simpleFloat;
    parserMap["FROM_DIALS"] = simpleBool;
    parserMap["DO_NOT_REJECT_REFLECTIONS"] = simpleBool;
    parserMap["REFINE_IN_PLANE_OF_DETECTOR"] = simpleBool;
    parserMap["ROTATION_MODE"] = simpleInt;
    parserMap["FIT_BACKGROUND_AS_PLANE"] = simpleBool;
    parserMap["SKIP_BAD_PIXELS"] = simpleBool;
    parserMap["ROUGH_CALCULATION"] = simpleBool;

    parserMap["LEARNING_TO_INDEX"] = simpleBool;
    parserMap["MINIMUM_SPOTS_EXPLAINED"] = simpleInt;
    parserMap["MINIMUM_TRUST_ANGLE"] = simpleFloat;
    parserMap["MINIMUM_TRUST_DISTANCE"] = simpleFloat;
    parserMap["SOLUTION_ANGLE_SPREAD"] = simpleFloat;
    parserMap["REJECT_CLOSE_SPOTS"] = simpleBool;
    parserMap["THOROUGH_SOLUTION_SEARCHING"] = simpleBool;
    parserMap["MAX_SEARCH_NUMBER_MATCHES"] = simpleInt;
    parserMap["MAX_SEARCH_NUMBER_SOLUTIONS"] = simpleInt;
    parserMap["MAX_MILLER_INDEX_TRIAL"] = simpleInt;
    parserMap["ACCEPT_ALL_SOLUTIONS"] = simpleBool;
    parserMap["INDEXING_MIN_RESOLUTION"] = simpleFloat;
    parserMap["SPOTS_PER_LATTICE"] = simpleInt;
    parserMap["RECIPROCAL_TOLERANCE"] = simpleFloat;
    parserMap["GOOD_SOLUTION_ST_DEV"] = simpleFloat;
    parserMap["BAD_SOLUTION_ST_DEV"] = simpleFloat;
    parserMap["GOOD_SOLUTION_SUM_RATIO"] = simpleFloat;
    parserMap["GOOD_SOLUTION_HIGHEST_PEAK"] = simpleInt;
    parserMap["SOLUTION_ATTEMPTS"] = simpleInt;
    parserMap["ONE_INDEXING_CYCLE_ONLY"] = simpleBool;
    parserMap["NEW_INDEXING_METHOD"] = simpleBool;
    parserMap["MAX_RECIPROCAL_DISTANCE"] = simpleFloat;
    parserMap["ALWAYS_FILTER_SPOTS"] = simpleBool;
    parserMap["MINIMUM_NEIGHBOURS"] = simpleInt;
    parserMap["MINIMUM_SOLUTION_NETWORK_COUNT"] = simpleInt;
    parserMap["SCRAMBLE_SPOT_VECTORS"] = simpleBool;
    parserMap["NETWORK_TRIAL_LIMIT"] = simpleInt;
    parserMap["INDEXING_TIME_LIMIT"] = simpleInt;
    parserMap["MAX_LATTICES_PER_IMAGE"] = simpleInt;
    parserMap["CHECKING_COMMON_SPOTS"] = simpleBool;
    parserMap["REJECT_IF_SPOT_COUNT"] = simpleInt;
    parserMap["POWDER_PATTERN_STEP"] = simpleFloat;
    parserMap["POWDER_PATTERN_STEP_ANGLE"] = simpleFloat;
    parserMap["BAD_SOLUTION_HIGHEST_PEAK"] = simpleInt;

    parserMap["IMAGE_MIN_SPOT_INTENSITY"] = simpleFloat;
    parserMap["IMAGE_MIN_CORRELATION"] = simpleFloat;
    parserMap["IMAGE_PIXEL_JUMP"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_HEIGHT"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_BACKGROUND"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_PADDING"] = simpleInt;
    parserMap["PROBE_DISTANCES"] = doubleVector;
    parserMap["RECIPROCAL_UNIT_CELL"] = doubleVector;

    parserMap["IGNORE_MISSING_IMAGES"] = simpleBool;

    parserMap["PIXEL_TOLERANCE"] = simpleFloat;
    parserMap["MINIMUM_CIRCLE_SPOTS"] = simpleInt;
    parserMap["COMMON_CIRCLE_THRESHOLD"] = simpleFloat;
    parserMap["MAX_UNIT_CELL"] = simpleFloat;
    parserMap["COMMON_CIRCLE_ANGLE_RANGE"] = doubleVector;
    parserMap["RANDOM_SEED"] = simpleInt;

        parserMap["PANEL_LIST"] = simpleString;
    parserMap["SKIP_LINES"] = simpleInt;
}

ParserFunction FileParser::splitLine(std::string line, std::string &command,
                std::string &rest)
{
    int space_index = (int)line.find_first_of(splitCharMajor);

    if (space_index == std::string::npos)
    {
        command = "NULL";
        return NULL;
    }

        command = line.substr(0, space_index);

        std::ostringstream stream;

        std::locale theLocale;
        for (std::string::size_type j = 0; j < command.length(); ++j)
                stream << std::toupper(command[j], theLocale);

        std::string upperCommand = stream.str();

        rest = line.substr(space_index + 1, std::string::npos);

    if (parserMap.count(upperCommand) == 0 && upperCommand != "PANEL" && upperCommand != "MASK")
    {
        std::cout << "Error: do not understand command " << upperCommand << std::endl;
        exit(1);
    }

    ParserFunction function = parserMap[upperCommand];
        command = upperCommand;

        return function;
}

bool FileParser::checkSpaces(std::string line)
{
    int space_index = (int)line.find_first_of(splitCharMajor);

        if (space_index == std::string::npos)
        {
                log << "Warning: " << line << " has no assignment" << std::endl;
                return false;
        }

        return true;
}

FileParser::FileParser(void)
{

}

FileParser::FileParser(std::string name, std::vector<std::string> someExtras)
{
        std::cout << "Initialising parser" << std::endl;

        this->filename = name;
        generateFunctionList();

    extras = someExtras;
}

FileParser::~FileParser()
{
        // TODO Auto-generated destructor stub
}
