//
//  FileReader.cpp
//  GameDriver
//
//  Created by Helen Ginn on 21/05/2014.
//  Copyright (c) 2014 Helen Ginn. All rights reserved.
//

#include "FileReader.h"

#include <fstream>
#include <sstream>
#include <cerrno>
#include <sys/stat.h>

#include "Logger.h"

vector<std::string> &FileReader::split(const std::string &s, char delim, vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<std::string> FileReader::split(const std::string s, const std::string &delim)
{
    vector<std::string> elems;
    std::string rest = s;
    
    int count = 0;
    bool finished = false;
    
    while (!finished)
    {
        count++;
        size_t index = rest.substr(1, rest.length() - 1).find(delim);
        
        if (index == std::string::npos)
        {
            index = rest.length() - 1;
            finished = true;
        }
        
        std::string cutout = rest.substr(1, index + 1);
        elems.push_back(cutout);
        
        rest = rest.substr(index + 1, s.length() - index - 1);
        
        if (index == 0)
            break;
    }
    
    return elems;
}

vector<std::string> FileReader::split(const std::string &s, char delim) {
    vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int FileReader::splitAtIndices(const std::string &s, vector<int> &positions, vector<std::string> &elems) {
    
    for (int i = 0; i < positions.size() - 1; i++)
    {
        int start = positions[i];
        int length = positions[i + 1] - start;
        
        if (s.size() < positions[i + 1])
            return 0;
        
        std::string segment = s.substr(start, length);
        elems.push_back(segment);
    }
    
    return 1;
}

bool FileReader::exists(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

std::string FileReader::get_file_contents(const char *filename)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned long)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return(contents);
    }
    
    sleep(1);
    
    std::string errString = "Could not get file contents for file " + std::string(filename);
    Logger::mainLogger->addString(errString);
    Logger::mainLogger->addString(strerror(errno));
    
    throw(errno);
}
