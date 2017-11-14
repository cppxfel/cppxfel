//
//  CSV.cpp
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "CSV.h"
#include <fstream>
#include <float.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string.h>

void CSV::histogram(std::map<double, int> histogram)
{
    for (std::map<double, int>::iterator it = histogram.begin(); it != histogram.end(); it++)
    {
        addEntry(0, it->first, (double)it->second);
    }
}

void CSV::addEntry(int dummy, ...)
{
    va_list arguments;
    va_start(arguments, dummy);
    Entry newEntry;

    for (int i = 0; i < headers.size(); i++)
    {
        double value = va_arg(arguments, double);
        newEntry.push_back(value);
    }

    entries.push_back(newEntry);
}

void CSV::writeToFile(std::string filename)
{
    std::ofstream csv;
    csv.open(filename.c_str());

    for (int i = 0; i < headers.size(); i++)
    {
        csv << headers[i] << ",";
    }

    csv << std::endl;

    for (int i = 0; i < entries.size(); i++)
    {
        Entry anEntry = entries[i];

        for (int j = 0; j < anEntry.size(); j++)
        {
            csv << anEntry[j] << ",";
        }

        csv << std::endl;
    }

    csv.close();
}

double CSV::valueForEntry(std::string header, int entry)
{
    for (int i = 0; i < headerCount(); i++)
    {
        if (headers[i] == header)
        {
            return entries[entry][i];
        }
    }

    return 0;
}

std::string CSV::mapToAscii(Plot plot)
{
    std::ostringstream stream;

    for (Plot::iterator it = plot.begin(); it != plot.end(); it++)
    {
        Row row = it->second;

        stream << "N: ";

        for (Row::iterator it2 = row.begin(); it2 != row.end(); it2++)
        {
            char charray[2];
            charray[0] = it2->second;
            charray[1] = '\0';

            stream << charray;
        }

        stream << std::endl;
    }

    return stream.str();
}

void CSV::minMaxCol(int col, double *min, double *max)
{
    *min = FLT_MAX;
    *max = -FLT_MAX;

    for (int i = 0; i < entries.size(); i++)
    {
        if (entries[i][col] > *max)
            *max = entries[i][col];

        if (entries[i][col] < *min)
            *min = entries[i][col];
    }
}

void CSV::writeStringToPlot(std::string text, Plot *plot, int y, int x)
{
    for (int i = 0; i < text.length(); i++)
    {
        (*plot)[y][x + i] = text[i];
    }
}

std::string CSV::plotColumns(int col1, int col2)
{
    char *cols = std::getenv("COLUMNS");
    char *linesstr = getenv("LINES");

 //   std::cout << "Plotting " << entries.size() << " points" << std::endl;

    int columns = 100;

    if (cols && strlen(cols))
    {
        columns = atoi(cols);
    }

    int lines = 50;

    if (cols && strlen(linesstr))
    {
        lines = atoi(linesstr);
    }

    if (columns < 80) columns = 80;
    if (lines < 50) lines = 50;

    const int leftMargin = 8;
    const int bottomMargin = 4;

    int graphCols = columns - leftMargin;
    int graphLines = lines - bottomMargin;

    Plot plot;

    for (int i = 0; i < lines; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            plot[i][j] = PlotBlank;
        }
    }

    double xMin, xMax, yMin, yMax;
    minMaxCol(col1, &xMin, &xMax);
    minMaxCol(col2, &yMin, &yMax);

  //  double xStep = (xMax - xMin) / (double)graphCols;
  //  double yStep = (yMax - yMin) / (double)graphLines;

    for (int i = 0; i < graphLines; i++)
    {
        plot[i][leftMargin] = PlotVerticalLine;
    }

    for (int i = 0; i < graphLines; i += ((double)graphLines / 10.))
    {
        plot[i][leftMargin - 1] = PlotHorizontalTickMark;
        std::ostringstream valueStream;
        double value = ((double)i / (double)graphLines) * (yMin - yMax) + yMax;

        valueStream << std::fixed << std::setprecision(2) << value;

        writeStringToPlot(valueStream.str(), &plot, i, 0);
    }

    for (int i = 0; i < graphCols; i++)
    {
        plot[graphLines][leftMargin + i] = PlotHorizontalLine;
    }

    writeStringToPlot(headers[col1], &plot, graphLines + bottomMargin - 1, graphCols / 2);

    for (int i = 0; i < graphCols; i += ((double)graphCols / 10.))
    {
        plot[graphLines + 1][leftMargin + i] = PlotVerticalTickMark;

        std::ostringstream valueStream;
        double value = ((double)i / (double)graphCols) * (xMax - xMin) + xMin;

        valueStream << std::fixed << std::setprecision(2) << value;

        writeStringToPlot(valueStream.str(), &plot, graphLines + 2, leftMargin + i - 2);
    }

    for (int i = 0; i < entries.size(); i++)
    {
        double xValue = entries[i][col1];
        double yValue = entries[i][col2];

        double xProportion = (xValue - xMin) / (xMax - xMin);
        double yProportion = (yValue - yMin) / (yMax - yMin);

        int xCharsIn = xProportion * (double)(graphCols - 1);
        int yCharsIn = yProportion * (double)graphLines;

        int xFullCharsIn = xCharsIn + leftMargin + 2;
        int yFullCharsIn = (graphLines - yCharsIn);

        if (yFullCharsIn == graphLines)
            yFullCharsIn -= 1;

        char currentChar = plot[yFullCharsIn][xFullCharsIn];

        if (currentChar == ' ')
        {
            plot[yFullCharsIn][xFullCharsIn] = '.';
        }
        else if (currentChar == '.')
        {
            plot[yFullCharsIn][xFullCharsIn] = '*';
        }
        else if (currentChar == '*')
        {
            plot[yFullCharsIn][xFullCharsIn] = '#';
        }
    }

    std::ostringstream ascii;
    ascii << std::endl << headers[col2] << std::endl << std::endl;

    ascii << mapToAscii(plot);

    std::ostringstream logged;
    logged << ascii.str() << std::endl;
    Logger::mainLogger->addStream(&logged);

    return ascii.str();
}
