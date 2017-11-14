//#ifndef lbfgs_scaling
//#define lbfgs_scaling

#include <vector>
#include <string>
#include <iostream>
#include "ScalingManager.h"
#include "parameters.h"

class Lbfgs_Scaling
{

public:
        Lbfgs_Scaling(char **filenames, int filenum);
        ~Lbfgs_Scaling(void);
        void run(void);
        Lbfgs_Scaling(vector<MtzPtr>mtzs);
};

//#endif
