//
//  Matrix.h
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#ifndef __RaddoseViewer__Matrix__
#define __RaddoseViewer__Matrix__

#include <iostream>
#include "Vector.h"
#include "parameters.h"
#include <cctbx/miller.h>
#include <scitbx/mat3.h>
#include "csymlib.h"

class Matrix
{
private:
    class Proxy
    {
        Matrix &a;
        int idx;
    public:
        Proxy(Matrix &a, int idx) : a(a), idx(idx) {}
        double operator= (double x) { a.components[idx] = x; return a.components[idx]; }
        double operator* (double x) { return a.components[idx] * x; }
        double operator* (Proxy x) { return a.components[idx] * x.a.components[idx]; }
        double operator*= (double x) { a.components[idx] *= x; return a.components[idx]; }
        double operator*= (Proxy x) { a.components[idx] *= x.a.components[idx]; return a.components[idx]; }
        
    };
    
    MatrixPtr unitCell;
    MatrixPtr rotation;
    double eulerA;
    double eulerB;
    double eulerC;
    
public:
    double components[16];
    
    bool isIdentity();
    void multiply(double scalar);
    void add(MatrixPtr secondMatrix);
    void subtract(MatrixPtr secondMatrix);
    Matrix(void);
    Matrix(double *components);
    Matrix(scitbx::mat3<double> unitcell, scitbx::mat3<double> rotation);
    MatrixPtr copy(void);
    void printDescription(bool detailed = false);
    std::string description(bool detailed = false, bool submatrix = false);
    Matrix inverse2DMatrix();
    MatrixPtr inverse3DMatrix();
    Matrix transpose();
    cctbx::miller::index<double> multiplyIndex(cctbx::miller::index<> *index);
    static void symmetryOperatorsForSpaceGroup(std::vector<MatrixPtr> *matrices, CSym::CCP4SPG *spaceGroup);
    
    void translate(double x, double y, double z);
    void rotateHK(double hRot, double kRot);
    void rotateABC(MatrixPtr oldMatrix, double aRot, double bRot, double cRot);
    void rotate(double alpha, double beta, double gamma);
    void rotateRoundUnitVector(vec unitVector, double radians);
    void rotateRoundUnitVector(double *unitVector, double radians);
    void multiply(Matrix &b);
    void multiplyVector(vec *vector);
    void preMultiply(Matrix &b);
    void scale(double scale);
    void scale(double a, double b, double c);
    void identity(void);
    void rotateModelAxes(double alpha, double beta, double gamma);
    void newMultiplyVector(double *vector[]);
    static MatrixPtr matrixFromUnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    void orientationMatrixUnitCell(double *a, double *b, double *c);
    void changeOrientationMatrixDimensions(double newA, double newB, double newC, double alpha, double beta, double gamma);
    void scaleUnitCellAxes(double aScale, double bScale, double cScale);
    void setComplexMatrix(MatrixPtr unitCell, MatrixPtr rotation);
    
    void rotate2D(double angle);
    void translation(double **vector);
    
    double getEwaldSphere(vec *vector);
    double getEwaldSphereNoMatrix(vec index);
    
    void eulerAngles(double *theta, double *phi, double *psi);
    double similarityToRotationMatrix(MatrixPtr mat2, double tolerance);
    void unitCellLengths(double **lengths);
    scitbx::mat3<double> cctbxMatrix(MatrixPtr theMatrix = MatrixPtr());
    void threeDimComponents(double **componentArray);
    void assignFromCctbxMatrix(scitbx::mat3<double> newMat);
    void assignFromCctbxMatrix(Matrix *changeMat, scitbx::mat3<double> newMat);
    double *array();
    void print(void);
    void recalculateOrientationMatrix();
    std::string summary();
    
    bool isComplex()
    {
        if (unitCell)
            return true;
        
        return false;
    }
    
    double determinant();
    Matrix operator*=(Matrix &b);
    Matrix operator*(Matrix &b);
    double &operator[](int index) {return components[index]; };
};

#endif /* defined(__RaddoseViewer__Matrix__) */
