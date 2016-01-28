//
//  Shoebox.cpp
//  cppxfel
//
//  Created by Helen Ginn on 26/01/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#include "Vector.h"
#include "Shoebox.h"
#include "Miller.h"
#include "Image.h"
#include "FileParser.h"
#include "Logger.h"

Shoebox::Shoebox(MillerPtr parent)
{
    miller = parent;
    
    if (parent)
    {
        ImagePtr image = parent->getImage();
        mmPerPixel = image->getMmPerPixel();
        detectorDistance = image->getDetectorDistance();
    }
    
    pixelLeak = FileParser::getKey("PIXEL_LEAK", PIXEL_LEAK);
    bandwidthMultiplier = FileParser::getKey("SHOEBOX_BANDWIDTH_MULTIPLIER", SHOEBOX_BANDWIDTH_MULTIPLIER);
}

void Shoebox::centre(int *centreX, int *centreY)
{
    int slowSide = 0; int fastSide = 0;
    sideLengths(&slowSide, &fastSide);
    
    *centreX = (slowSide % 2 == 0) ? slowSide / 2 : (slowSide - 1) / 2;
    *centreY = (fastSide % 2 == 0) ? fastSide / 2 : (fastSide - 1) / 2;
}

void Shoebox::sideLengths(int *slowSide, int *fastSide)
{
    *slowSide = (int)shoebox.size();
    if (slowSide == 0)
    {
        fastSide = 0;
        return;
    }
    
    *fastSide = (int)shoebox[0].size();
}

void Shoebox::simpleShoebox(int foregroundLength, int neitherLength, int backgroundLength, bool shoeboxEven)
{
    int centre = backgroundLength;
    int adjustment1 = shoeboxEven ? 0 : 1;
    int adjustment2 = shoeboxEven ? -1 : 0;
    even = shoeboxEven;
    
    sendLog();
    
    clearShoebox();
    
    for (int i = 0; i < backgroundLength * 2 + adjustment1; i++)
    {
        shoebox.push_back(vector<double>());
        
        for (int j = 0; j < backgroundLength * 2 + adjustment1; j++)
        {
            shoebox[i].push_back(-1);
        }
    }
    
    for (int i = centre - neitherLength; i <= centre + neitherLength + adjustment2; i++)
    {
        for (int j = centre - neitherLength; j <= centre + neitherLength + adjustment2; j++)
        {
            shoebox[i][j] = 0;
        }
    }
    
    for (int i = centre - foregroundLength; i <= centre + foregroundLength + adjustment2; i++)
    {
        for (int j = centre - foregroundLength; j <= centre + foregroundLength  + adjustment2; j++)
        {
            shoebox[i][j] = 1;
        }
    }
}

void shoeboxCoordinates(double largeDim, double resolution, int i, int j, double *x, double *y)
{
    double midCoord = (largeDim) / 2;
    
    double extremeCoord = largeDim / resolution;
    extremeCoord /= 2;
    
    double iMinus = (double)i - midCoord;
    double jMinus = (double)j - midCoord;
    
    *x = iMinus / midCoord;
    *y = jMinus / midCoord;
    
    *x *= extremeCoord;
    *y *= extremeCoord;
}

void Shoebox::printShoebox()
{
    std::ostringstream logged;
    
    logged << "Shoebox for Miller" << std::endl;
    
    for (int i = 0; i < shoebox.size(); i++)
    {
        for (int j = 0; j < shoebox[i].size(); j++)
        {
            logged << shoebox[i][j] << "\t";
        }
        
        logged << std::endl;
    }
    
    logged << std::endl;
    
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

void Shoebox::putBoxOnBackground(Box &smallBox, Box &newShoebox)
{
    int slowSize = (int)smallBox.size();
    int fastSize = (int)smallBox[0].size();
    
    int newSlow = slowSize + 2;
    int newFast = fastSize + 2;
    
    newShoebox.resize(newSlow);
    for (int i = 0; i < newSlow; i++)
    {
        newShoebox[i].resize(newFast);
        for (int j = 0; j < newFast; j++)
        {
            newShoebox[i][j] = -1;
            
            if (i > 0 && i < newSlow - 1 && j > 0 && j < newFast - 1)
            {
                int x = i - 1;
                int y = j - 1;
                
                newShoebox[i][j] = smallBox[x][y];
            }
        }
    }
}

void Shoebox::putBoxOnPaddedBackground(Box &smallBox, Box &newShoebox)
{
    int slowSize = (int)smallBox.size();
    int fastSize = (int)smallBox[0].size();
    
    int newSlow = slowSize + 4;
    int newFast = fastSize + 4;
    
    newShoebox.resize(newSlow);
    for (int i = 0; i < newSlow; i++)
    {
        newShoebox[i].resize(newFast);
        for (int j = 0; j < newFast; j++)
        {
            newShoebox[i][j] = 0;
            
            if (i == 0 || i == newSlow - 1 || j == 0 || j == newFast - 1)
                newShoebox[i][j] = -1;
            
            if (i > 1 && i < newSlow - 2 && j > 1 && j < newFast - 2)
            {
                int x = i - 2;
                int y = j - 2;
                
                newShoebox[i][j] = smallBox[x][y];
            }
        }
    }
}

void Shoebox::chopBox(Box &smallBox)
{
    for (int i = 0; i < smallBox.size(); i++)
    {
        bool allZero = true;
        
        for (int j = 0; j < smallBox[i].size(); j++)
        {
            if (smallBox[i][j] != 0)
                allZero = false;
        }
        
        if (allZero)
        {
            smallBox.erase(smallBox.begin() + i);
            i--;
        }
    }
    
    if (smallBox.size() == 0)
        return;
    
    for (int i = 0; i < smallBox[0].size(); i++)
    {
        bool allZero = true;
        
        for (int j = 0; j < smallBox.size(); j++)
        {
            if (smallBox[j][i] != 0)
                allZero = false;
        }
        
        if (allZero)
        {
            for (int j = 0; j < smallBox.size(); j++)
            {
                smallBox[j].erase(smallBox[j].begin() + i);
            }
            i--;
        }
    }
}

void Shoebox::centreShoebox(Box &smallBox)
{
    if (smallBox.size() == 0)
        return;
    
    int slowSize = (int)smallBox.size();
    int fastSize = (int)smallBox[0].size();
    
    if (slowSize % 2 == 0)
    {
        smallBox.push_back(vector<double>());
        
        smallBox[slowSize].resize(fastSize);
        
        for (int i = 0; i < fastSize; i++)
        {
            smallBox[slowSize][i] = 0;
        }
        
        for (int i = 0; i < fastSize; i++)
        {
            for (int j = slowSize; j > 0; j--)
            {
                smallBox[j][i] = 0.5 * smallBox[j][i] + 0.5 * smallBox[j - 1][i];
            }
            
            smallBox[0][i] *= 0.5;
        }
        
        slowSize = (int)smallBox.size();
        fastSize = (int)smallBox[0].size();
    }
    
    if (fastSize % 2 == 0)
    {
        for (int i = 0; i < slowSize; i++)
        {
            smallBox[i].push_back(0.0);
        }
        
        for (int i = 0; i < slowSize; i++)
        {
            for (int j = fastSize; j > 0; j--)
            {
                smallBox[i][j] = 0.5 * smallBox[i][j] + 0.5 * smallBox[i][j - 1];
            }
            
            smallBox[i][0] *= 0.5;
        }
    }

}

void Shoebox::compressBigShoebox(double width, double height, double angle, Box &smallBox)
{
    double maxForegroundDimension = (width > height) ? width : height;
    int maxForegroundDimensionInt = maxForegroundDimension + 1;
    
    if (maxForegroundDimensionInt % 2 == 0)
        maxForegroundDimensionInt++;
    
    double resolution = 10;
    
    double bigMaxForeDimFloat = resolution * maxForegroundDimensionInt;
    
    int bigMaxForeDim = int(bigMaxForeDimFloat);
    
    int bufferedMax = bigMaxForeDim + resolution;
    
    Box bigBox;
    bigBox.resize(bufferedMax + 1);
    for (int i = 0; i < bufferedMax + 1; i++)
        bigBox[i].resize(bufferedMax + 1);
    
    Matrix shoeboxMatrix = Matrix();
    
    shoeboxMatrix.scale(width / 2, height / 2, 1);
    shoeboxMatrix.rotate2D(angle);
    
    Matrix inverseMatrix = shoeboxMatrix.inverse2DMatrix();
    
    for (int i = 0; i < bufferedMax + 1; i++)
    {
        for (int j = 0; j < bufferedMax + 1; j++)
        {
            double x; double y;
            shoeboxCoordinates(bigMaxForeDim, resolution, i, j, &x, &y);
            
            vec coordVec = new_vector(x, y, 1);
            inverseMatrix.multiplyVector(&coordVec);
            
            //      double origDistance = sqrt(x * x + y * y);
            double distance = sqrt(pow(coordVec.h, 2) + pow(coordVec.k, 2));
            
            bigBox[i][j] = distance < 1;
        }
    }
    
    int smallMaxForeDim = int(maxForegroundDimensionInt);
    
    int bufferedSmall = smallMaxForeDim;
    
    smallBox.resize(bufferedSmall);
    for (int i = 0; i < bufferedSmall; i++)
        smallBox[i].resize(bufferedSmall);
    
    
    for (int i = 0; i < bufferedSmall; i++)
    {
        for (int j = 0; j < bufferedSmall; j++)
        {
            double minI = i * resolution;
            double maxI = (i + 1) * resolution;
            
            double minJ = j * resolution;
            double maxJ = (j + 1) * resolution;
            
            double total = 0;
            double count = 0;
            
            for (int m = minI; m < maxI; m++)
            {
                for (int n = minJ; n < maxJ; n++)
                {
                    count++;
                    total += bigBox[m][n];
                }
            }
            
            total /= count;
            
            smallBox[i][j] = total;
        }
    }
}

void Shoebox::complexShoebox(double wavelength, double bandwidth, double radius)
{
    MillerPtr tempMiller = getMiller();
    
    if (!tempMiller)
    {
        return;
    }
    
    MatrixPtr mat;
    mat = tempMiller->getMatrix();
    
    vec hkl = new_vector(tempMiller->getH(), tempMiller->getK(), tempMiller->getL());
    mat->multiplyVector(&hkl);
    
    double hk = sqrt(hkl.h * hkl.h + hkl.k * hkl.k);
    
    double maxHK = hk + radius;
    double minHK = hk - radius;
    
    double minWave = wavelength - bandwidth * bandwidthMultiplier;
    double maxWave = wavelength + bandwidth * bandwidthMultiplier;
    
    double max_mm = (maxHK * detectorDistance / (1 / maxWave + hkl.l));
    double min_mm = (minHK * detectorDistance / (1 / minWave + hkl.l));
    
    double width_mm = (radius * detectorDistance / (1 / wavelength));
    
    double mmDiff = max_mm - min_mm;
    
    double widthLeak = 1 * pixelLeak;

    double pixelHeight = mmDiff / mmPerPixel + pixelLeak;
    double pixelWidth = 2 * width_mm / mmPerPixel + widthLeak;
    
    double angle = cartesian_to_angle(hkl.h, hkl.k);
    angle += M_PI / 2;
    
    Box smallBox;
    Box newShoebox;
    compressBigShoebox(pixelWidth, pixelHeight, angle, smallBox);
    chopBox(smallBox);
    centreShoebox(smallBox);
    putBoxOnPaddedBackground(smallBox, newShoebox);

    shoebox = newShoebox;
    printShoebox();
}

void Shoebox::clearShoebox()
{
    for (int i = 0; i < shoebox.size(); i++)
    {
        shoebox[i].clear();
        vector<double>().swap(shoebox[i]);
    }
    
    shoebox.clear();
    vector<vector<double> >().swap(shoebox);
}

Shoebox::~Shoebox()
{
//    std::cout << "Deallocating shoebox" << std::endl;
}
