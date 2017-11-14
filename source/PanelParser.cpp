/*
 * PanelParser.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#include "PanelParser.h"
#include "parameters.h"
#include "FileReader.h"
#include <boost/variant.hpp>
#include "Panel.h"
#include "Logger.h"

PanelParser::PanelParser(std::string filename) :
                FileParser(filename)
{

}

PanelParser::~PanelParser()
{
    panels.clear();
}

void PanelParser::addPanel(std::string rest, PanelTag tag)
{
    std::string keyword = (tag == PanelTagNormal ? "PANEL" : "MASK");

        ParametersMap singleParameter;

    this->doubleVector(&singleParameter, keyword, rest);

        vector<double> panelCoords = boost::get<vector<double> >(singleParameter[keyword]);

        PanelPtr ptr = PanelPtr(new Panel(panelCoords));
    ptr->setTag(tag);
        Panel::setupPanel(ptr);
        panels.push_back(ptr);
}

void PanelParser::parse(bool fromPython)
{
    if (panels.size())
    {
        Logger::mainLogger->addString("Already loaded panels.");
        return;
    }

        std::string fileContents = FileReader::get_file_contents(filename.c_str());
        vector<std::string> fileLines = FileReader::split(fileContents, '\n');

    splitCharMajor = ' ';
    splitCharMinor = ' ';

        for (int i = 0; i < fileLines.size(); i++)
        {
                std::string line = fileLines[i];

                if (line.length() == 0)
                        continue;

                if (line.substr(0, 1) == "#")
                        continue;

                int space_index = (int)line.find_first_of(" ");

                if (space_index == std::string::npos)
                {
                        std::cout << "Warning: " << line << " has no assignment, ignored"
                                        << std::endl;
                        continue;
                }

                if (!checkSpaces(line))
                        continue;

                std::string command;
                std::string rest;

                splitLine(line, command, rest);

                if (command == "PANEL")
                {
                        addPanel(rest);
                }

        if (command == "MASK")
        {
            addPanel(rest, PanelTagBad);
        }
        }

    Logger::mainLogger->addStream(&log);
    log.str("");
}
