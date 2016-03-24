#include <string>
#include <sstream>
#include "misc.h"
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>
#include <math.h>

std::string f_to_str(double val)
{
	std::ostringstream ss;
	ss << val;
	std::string temp = ss.str();

	return temp;
}

std::string i_to_str(int val)
{
	std::ostringstream ss;
	ss << val;
	std::string temp = ss.str();
	
	return temp;
}

bool replace(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string getBaseFilename(std::string filename)
{
    std::string fName = getFilename(filename);
    size_t pos = fName.rfind(".");
    if(pos == std::string::npos)  //No extension.
        return fName;
    
    if(pos == 0)    //. is at the front. Not an extension.
        return fName;
    
    return fName.substr(0, pos);
}

std::string getFilename(std::string filename)
{
    size_t pos = filename.rfind("/");
    if(pos == std::string::npos)  //No path.
        return filename;
    
    return filename.substr(pos + 1, filename.length());
}

std::string getPath(std::string filename)
{
    size_t pos = filename.rfind("/");
    if(pos == std::string::npos)  //No path.
        return filename;
    
    return filename.substr(0, pos + 1);
}

unsigned long factorial(unsigned long n)
{
    unsigned long ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

unsigned int choose(unsigned long n, unsigned long choose)
{
    if (n > 18) n = 18;
    
    unsigned long nFac = factorial(n);
    unsigned long cFac = factorial(choose);
    unsigned long ncFac = factorial(n - choose);
    
    unsigned int value = (unsigned int)(nFac / (cFac * ncFac));
    
    return value;
}

double proportion(int n)
{
    double prop = n / pow(2, n);
    
    return prop;
}