//
//  CSV.h
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__CSV__
#define __cppxfel__CSV__

#include <stdio.h>
#include <vector>
#include <string>
#include "parameters.h"
#include <cstdarg>
#include "LoggableObject.h"

typedef enum
{
    PlotVerticalLine = '|',
    PlotHorizontalLine = '_',
    PlotHorizontalTickMark = '-',
    PlotVerticalTickMark = '\'',
    PlotBlank = ' ',

} PlotChar;

typedef std::vector<double> Entry;
typedef std::map<int, char> Row;
typedef std::map<int, Row > Plot;

class CSV
{
private:
    std::vector<std::string> headers;
    std::vector<Entry> entries;
    void minMaxCol(int col, double *min, double *max);
    std::string mapToAscii(Plot plot);
    void writeStringToPlot(std::string text, Plot *plot, int x, int y);
public:
    CSV(int count, ...)
    {
        va_list arguments;
        va_start(arguments, count);

        for (int i = 0; i < count; i++)
        {
            std::string header = std::string(va_arg(arguments, char *));
            headers.push_back(header);
        }

        va_end(arguments);
    }

    void addEntry(int dummy, ...);
    void writeToFile(std::string filename);
    double valueForEntry(std::string header, int entry);
    void histogram(std::map<double, int> histogram);
    std::string plotColumns(int col1, int col2);

    int entryCount()
    {
        return (int)entries.size();
    }

    int headerCount()
    {
        return (int)headers.size();
    }
};

#endif /* defined(__cppxfel__CSV__) */
