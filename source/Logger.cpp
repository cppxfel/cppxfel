//
//  Logger.cpp
//  cppxfel
//
//  Created by Helen Ginn on 18/01/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#include "Logger.h"

LoggerPtr Logger::mainLogger;
bool Logger::ready;

Logger::Logger()
{
    ready = false;
    printedLogLevel = LogLevelNormal;
}

Logger::~Logger()
{

}

bool Logger::isReady()
{
    return mainLogger->ready;
}

void Logger::addString(std::string theString, LogLevel level)
{
    std::ostringstream stream;
    stream << theString << std::endl;
    addStream(&stream, level);
}

bool Logger::tryLock(std::mutex &lock, int maxTries)
{
    int tryCount = 0;
    bool locked = false;

    while (!locked)
    {
        try
        {
            locked = lock.try_lock();
        }
        catch (std::system_error err)
        {
            locked = false;
        }

        tryCount ++;
    }

    return locked;
}

void Logger::addStream(std::ostringstream *stream, LogLevel level)
{
    if (level > printedLogLevel)
        return;

    bool success = tryLock(mtx);

    if (!success)
        return;

    boost::thread::id thread_id = boost::this_thread::get_id();

    StreamPtr ptr = StreamPtr(new std::ostringstream());
    *ptr << stream->str();

    if (stringsToOutput.count(thread_id) == 0)
    {
        stringsToOutput[thread_id] = vector<LogAndLevel>();
        stringsToOutput[thread_id].reserve(100);
    }

    LogAndLevel logAndLevel = std::make_pair(ptr, level);

    stringsToOutput[thread_id].push_back(logAndLevel);

    ready = true;

    printBlock.notify_all();

    mtx.unlock();
}

void Logger::awaitPrinting()
{
    std::unique_lock<std::mutex> lck(mtx);

    while (true)
    {
        printBlock.wait(lck, isReady);

        std::unique_lock<std::mutex> writeLock(writing, std::defer_lock);

        while (!writeLock.try_lock()) {}

        int count = 0;
        int num = 0;

        for (StringMap::iterator it = stringsToOutput.begin();
             it != stringsToOutput.end(); ++it)
        {
            for (int i = 0; i < stringsToOutput[it->first].size(); i++)
            {
                if (stringsToOutput[it->first][i].second > printedLogLevel)
                    continue;

                std::cout << stringsToOutput[it->first][i].first->str() << std::flush;
                count++;

                stringsToOutput[it->first][i].first->str("");
                stringsToOutput[it->first][i].first->clear();
            }

            num++;

            stringsToOutput[it->first].clear();
            vector<LogAndLevel>().swap(stringsToOutput[it->first]);
        }

        ready = false;

        writing.unlock();
    }
}

void Logger::changePriorityLevel(LogLevel newLevel)
{
    printedLogLevel = newLevel;
}

void Logger::awaitPrintingWrapper(LoggerPtr logger)
{
    logger->awaitPrinting();
}
