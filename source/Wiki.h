/*
 * Wiki.h
 *
 *  Created on: 17 Nov 2014
 *      Author: helenginn
 */

#ifndef WIKI_H_
#define WIKI_H_

#include <string>

class Wiki
{
private:
        std::string contents;
        std::string filename;
public:
        Wiki();
        Wiki(std::string filename);
        void process();
        virtual ~Wiki();
};

#endif /* WIKI_H_ */
