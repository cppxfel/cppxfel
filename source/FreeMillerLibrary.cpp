//
//  FreeMillerLibrary.cpp
//  cppxfel
//
//  Created by Helen Ginn on 17/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "FreeMillerLibrary.h"
#include "FileReader.h"
#include "FileParser.h"
#include "csymlib.h"
#include "parameters.h"
#include "UnitCellLattice.h"
#include <sstream>
#include "CSV.h"
#include "Miller.h"

FreeMillerLibrary FreeMillerLibrary::library;

FreeMillerLibrary::FreeMillerLibrary(std::string filename)
{
    std::string contents = FileReader::get_file_contents(filename.c_str());
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    for (int i = 1; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
        std::vector<std::string> components = FileReader::split(line, ',');
        
        if (components.size() >= 3)
        {
            int h = atoi(components[0].c_str());
            int k = atoi(components[1].c_str());
            int l = atoi(components[2].c_str());
            
            if (h == 0 && k == 0 && l == 0)
                continue;
            
            vec3<int> hkl = scitbx::vec3<int>(h, k, l);
            
            freeIndices.push_back(hkl);
        }
    }
}

FreeMillerLibrary::FreeMillerLibrary(std::string filename, double maxResolution)
{
    double freeProportion = FileParser::getKey("FREE_MILLER_FRACTION", 0.05);
    int spaceGroup = FileParser::getKey("SPACE_GROUP", 0);
    std::vector<double> unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    if (spaceGroup == 0)
    {
        std::cout << "Space group is unassigned! Please assign your space group using the command SPACE_GROUP in order to generate a free Miller set during refinement." << std::endl;
        exit(1);
    }
    
    if (unitCell.size() == 0)
    {
        std::cout << "Please specify the UNIT_CELL in order to generate a free Miller set during refinement." << std::endl;
        exit(1);
    }
    
    UnitCellLatticePtr lattice = UnitCellLatticePtr(new UnitCellLattice(unitCell[0], unitCell[1], unitCell[2],
                                                     unitCell[3], unitCell[4], unitCell[5], spaceGroup, maxResolution));
    
    for (int i = 0; i < lattice->standardVectorCount(); i++)
    {
        vec3<int> oneVec = lattice->intVector(i);
        
        double randomNumber = rand() / (double)RAND_MAX;
        
        if (randomNumber < freeProportion)
        {
            freeIndices.push_back(oneVec);
        }
    }
    
    if (filename != "")
    {
        CSV csv = CSV(3, "h", "k", "l");
        
        for (int i = 0; i < freeIndexCount(); i++)
        {
            csv.addEntry(0, (double)freeIndices[i][0], (double)freeIndices[i][1], (double)freeIndices[i][2]);
        }
        
        csv.writeToFile(filename);
    }
}

void FreeMillerLibrary::setup()
{
    double resolution = FileParser::getKey("MAX_RESOLUTION_ALL", 1.0);
    std::string filename = FileParser::getKey("FREE_MILLER_LIST", std::string("freeMillers.txt"));
    
    if (filename == "")
    {
        library = FreeMillerLibrary(filename, resolution);
        return;
    }
    
    std::string contents = "";
    
    try
    {
        contents = FileReader::get_file_contents(filename.c_str());
    }
    catch (int)
    {
        std::ostringstream logged;
        logged << "Could not find " << filename << " so will make a new free Miller list." << std::endl;
        Logger::mainLogger->addStream(&logged);
    }
    
    if (contents.length() > 0)
    {
        std::ostringstream logged;
        logged << "Loading pre-determined free set from " << filename << "." << std::endl;
        Logger::mainLogger->addStream(&logged);
        
        library = FreeMillerLibrary(filename);
    }
    else
    {
        library = FreeMillerLibrary(filename, resolution);
        return;
    }
}

void FreeMillerLibrary::printSummary()
{
    std::ostringstream logged;
    logged << "Total number of free Miller indices: " << freeIndexCount() << std::endl;
    Logger::mainLogger->addStream(&logged);
}

bool FreeMillerLibrary::isCertainMillerFree(Miller *miller)
{
    if (freeIndexCount() == 0)
        return false;
    
    int h = miller->getH();
    int k = miller->getK();
    int l = miller->getL();
    
    for (int i = 0; i < freeIndexCount(); i++)
    {
        if (freeIndices[i][0] == h && freeIndices[i][1] == k && freeIndices[i][2] == l)
        {
            return true;
        }
    }
    
    return false;
}