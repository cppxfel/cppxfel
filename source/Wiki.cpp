/*
 * Wiki.cpp
 *
 *  Created on: 17 Nov 2014
 *      Author: helenginn
 */

#include "Wiki.h"
#include "FileReader.h"
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "definitions.h"

std::string directory()
{
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
	{
		return std::string(cwd);
	}
	else
	{
		perror("getcwd() error");
		return "Unknown";
	}
}

Wiki::Wiki()
{
	// TODO Auto-generated constructor stub
	contents = "";
	filename = "";
}

Wiki::Wiki(std::string filename)
{
	this->filename = filename;
	contents = FileReader::get_file_contents(filename.c_str());
}

void Wiki::process()
{
	vector<std::string> lines = FileReader::split(contents, '\n');

	bool firstTabSeparated = false;

	/* TAGS */

	vector<std::string> tags;

#ifdef MAC
	tags.push_back("Apple");
#else
	tags.push_back("Linux");
#endif

	vector<std::string> directories = FileReader::split(filename, '/');

	int num = (int)directories.size() - 2;

	tags.push_back(directories[num]);

	std::cout << "{{tag>";

	for (int i=0; i < tags.size(); i++)
	{
		std::cout << tags[i] << " ";
	}

	std::cout << "}}" << std::endl;

	/* HEAD TITLE */

	struct tm* clock;				// create a time structure

	struct stat attrib;			// create a file attribute structure

	stat(filename.c_str(), &attrib);		// get the attributes of afile.txt

	clock = gmtime(&(attrib.st_mtime));	// Get the last modified time and put it into the time structure

	std::cout << "====== " << directories[directories.size() - 1] << ", ";

	std::cout << clock->tm_mday << " ";

	switch (clock->tm_mon)
	{
	case 0:
		std::cout << "January";
		break;
	case 1:
		std::cout << "February";
		break;
	case 2:
		std::cout << "March";
		break;
	case 3:
		std::cout << "April";
		break;
	case 4:
		std::cout << "May";
		break;
	case 5:
		std::cout << "June";
		break;
	case 6:
		std::cout << "July";
		break;
	case 7:
		std::cout << "August";
		break;
	case 8:
		std::cout << "September";
		break;
	case 9:
		std::cout << "October";
		break;
	case 10:
		std::cout << "November";
		break;
	case 11:
		std::cout << "December";
		break;
	}

	std::cout << " " << clock->tm_year + 1900 << " ======" << std::endl;

	/* DATA */

	for (int i = 0; i < lines.size(); i++)
	{
		std::string *line = &lines[i];
		if (line->substr(0, 3) == "N: ")
		{
			std::string cutLine = line->substr(3, std::string::npos);

			vector<std::string> tabSeparated = FileReader::split(cutLine,
					'\t');

			if (tabSeparated.size() >= 2)
			{
				std::string separator = "|";

				if (firstTabSeparated == false)
				{
					// make headings
					separator = "^";
					firstTabSeparated = true;
				}

				for (int i = 0; i < tabSeparated.size(); i++)
				{
					std::cout << separator << " " << tabSeparated[i] << " ";
				}
				std::cout << " " << separator << std::endl;
			}
			else
			{
				firstTabSeparated = false;
				std::cout << cutLine << std::endl;
				std::cout << std::endl;
			}
		}
	}
}

Wiki::~Wiki()
{
	// TODO Auto-generated destructor stub
}

