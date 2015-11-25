#include <string>
#include <stdio.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>

std::string i_to_str(int val);
std::string f_to_str(double val);

bool replace(std::string& str, const std::string& from, const std::string& to);

std::string getBaseFilename(std::string filename);
std::string getFilename(std::string filename);
std::string getPath(std::string filename);

unsigned long factorial(unsigned long n);
unsigned int choose(unsigned long n, unsigned long choose);
double proportion(int n);