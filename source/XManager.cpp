//
//  XManager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 03/02/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#include "XManager.h"
#include "FileReader.h"
#include "Miller.h"
#include "parameters.h"
#include "FileParser.h"



XManager::XManager()
{
    lineSplitters.push_back(0);
    lineSplitters.push_back(4);
    lineSplitters.push_back(8);
    lineSplitters.push_back(12);
    lineSplitters.push_back(14);
    lineSplitters.push_back(22);
    lineSplitters.push_back(30);
    lineSplitters.push_back(37);
    lineSplitters.push_back(43);
    lineSplitters.push_back(49);
    lineSplitters.push_back(56);
    lineSplitters.push_back(63);
    lineSplitters.push_back(69);
    lineSplitters.push_back(77);
}

XManager::~XManager()
{
    
}

void XManager::setFilenames(vector<std::string> newFiles)
{
    filenames = newFiles;
    filename = newFiles[0];
}

void XManager::loadReflections(PartialityModel model)
{
    int fileCount = (int)filenames.size();
    double singleRefPartiality = 1 / (double)fileCount;
    int skipLines = FileParser::getKey("SKIP_LINES", 0);
    
    logged << "First filename: " << filename;
    logged << " containing X file components: " << fileCount << std::endl;

    sendLog();
    
    // set unit cell and space group
    
    int spgNum = FileParser::getKey("SPACE_GROUP", 0);
    
    if (spgNum == 0)
    {
        std::cout << "Please enter space group under keyword SPACE_GROUP using designated space group number" << std::endl;
        exit(1);
    }
    
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    if (unitCell.size() == 0)
    {
        std::cout << "Please set unit cell under keyword UNIT_CELL followed by six parameters." << std::endl;
        exit(1);
    }
    
    setUnitCell(unitCell);
    
    setDefaultMatrix();
    
    CCP4SPG *spg = ccp4spg_load_by_ccp4_num(spgNum);
    this->setLowGroup(spg);
    
    for (int i = 0; i < fileCount; i++)
    {
        std::string filename = filenames[i];
        std::string contents = FileReader::get_file_contents(filename.c_str());
        vector<std::string> lines = FileReader::split(contents, '\n');
        
        for (int j = skipLines; j < lines.size(); j++)
        {
            vector<std::string> components;
            
            int success = FileReader::splitAtIndices(lines[j], lineSplitters, components);
            
            if (success == 0 && reflections.size() == 0)
            {
                std::cout << "Cannot parse line " << j << " of file, perhaps you need to skip some lines. Include SKIP_LINES keyword." << std::endl;
                exit(1);
            }
            else if (success == 0)
            {
                continue;
            }
            
            int h = atoi(components[0].c_str());
            int k = atoi(components[1].c_str());
            int l = atoi(components[2].c_str());
            
            double intensity = atof(components[4].c_str());
            double sigi = atof(components[5].c_str());
        //    double correction = atof(components[11].c_str());
            
      //      if (sigi < 0)
       //         continue;
            
            int reflid = this->index_for_reflection(h, k, l, false);
            
            Reflection *currentReflection = NULL;
            
            this->findReflectionWithId(reflid, &currentReflection);
            
            if (currentReflection == NULL)
            {
                currentReflection = new Reflection();
                
                MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
                miller->setData(0, 0, 0, 0);
                miller->setParent(currentReflection);
                miller->setMatrix(matrix);
                miller->setPartialityModel(PartialityModelFixed);
                miller->setFilename(filename);
                currentReflection->addMiller(miller);
                
                currentReflection->calculateResolution(this);
                
                reflections.push_back(currentReflection);
                sortLastReflection();
            }
            
            MillerPtr miller = currentReflection->miller(0);
            
            double rawInt = miller->getRawIntensity();
            double sigma = miller->getSigma();
            
            double partiality = miller->getPartiality();
            
            miller->setRawIntensity(rawInt + intensity);
            
            miller->setSigma(sigma + sigi);
            miller->setPartiality(partiality + singleRefPartiality);
        }
    }
}