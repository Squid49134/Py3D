// C LIBRARY TO EXTEND FIELD TRACER
// William Ransom, University of Delaware 2017

// gcc -o TracerFunctions.so -shared -fPIC TracerFunctions.c

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// MACRO DEFINITIONS FOR CONSOLIDATED SYNTAX

// 3D
#define Bx(i,j,k)  Bx[(i) + (j)*Xsize + (k)*Xsize*Ysize]
#define By(i,j,k)  By[(i) + (j)*Xsize + (k)*Xsize*Ysize]
#define Bz(i,j,k)  Bz[(i) + (j)*Xsize + (k)*Xsize*Ysize]

// 2D
#define Bx2(i,j)  Bx[(i) + (j)*Xsize]
#define By2(i,j)  By[(i) + (j)*Xsize]


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// SLOPE CALCULATING METHOD FOR E

// Works identically to CalcSlopes3 accept it does not normalize the B
// components
void CalcSlopes3E(double x, double y, double z, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, double * B){

    double Wx, Wy, Wz;
    int i, j, k, i1, j1, k1;
    
    if (x >= 0 && y >= 0 && z >= 0){
        Wx = fmod(x, 1);
        Wy = fmod(y, 1);
        Wz = fmod(z, 1);
        i = (int) x;
        j = (int) y;
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
    }
    else if (y >= 0 && z >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
    }
    else if (x >= 0 && z >= 0){
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wx = fmod(x, 1);
        i = (int) x;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }
    else if (y >= 0 && x >= 0){
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wx = fmod(x, 1);
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else if (z >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }
    else if (y >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else if (x >= 0){
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wx = fmod(x, 1);
        i = (int) x;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else{
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }

    B[0] = (1-Wx)*(1-Wy)*(1-Wz)*Bx(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bx(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bx(i1,j,k) + Wx*Wy*(1-Wz)*Bx(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bx(i,j,k1) + (1-Wx)*Wy*Wz*Bx(i,j1,k1) + Wx*(1-Wy)*Wz*Bx(i1,j,k1) + Wx*Wy*Wz*Bx(i1,j1,k1);
    B[1] = (1-Wx)*(1-Wy)*(1-Wz)*By(i,j,k) + (1-Wx)*Wy*(1-Wz)*By(i,j1,k) + Wx*(1-Wy)*(1-Wz)*By(i1,j,k) + Wx*Wy*(1-Wz)*By(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*By(i,j,k1) + (1-Wx)*Wy*Wz*By(i,j1,k1) + Wx*(1-Wy)*Wz*By(i1,j,k1) + Wx*Wy*Wz*By(i1,j1,k1);
    B[2] = (1-Wx)*(1-Wy)*(1-Wz)*Bz(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bz(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bz(i1,j,k) + Wx*Wy*(1-Wz)*Bz(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bz(i,j,k1) + (1-Wx)*Wy*Wz*Bz(i,j1,k1) + Wx*(1-Wy)*Wz*Bz(i1,j,k1) + Wx*Wy*Wz*Bz(i1,j1,k1);

}


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// STANDARD 3D TRACING METHODS

// Method to determine Bx/|B|, By/|B| and Bz/|B| at given point [x, y, z] in
// simulation of size Xsize, Ysize, Zsize, returning the results in vector B
// for the 3D RK4 process
void CalcSlopes3(double x, double y, double z, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, double * B){ 
    // Variables for interpolation of B at points in between grid spaces
    double Wx, Wy, Wz;
    int i, j, k, i1, j1, k1;
    
    // If x, y and z are greater than 0, one method may be used to interpolate
    // data, if x, y, or z is less than 0, a different method must be used
    // This distinction is necessary since the simulation has periodic
    // boundaries, and therefore points may be negative
    if (x >= 0 && y >= 0 && z >= 0){
        // Wx represents the distance between x and the x value of the closest
        // lower left data grid space, such that if x = 12.7, Wx = .7
        // Wy and Wz follow the same pattern
        // if x, y, and z are positive this value can be found using modulus
        Wx = fmod(x, 1);
        Wy = fmod(y, 1);
        Wz = fmod(z, 1);
        // i is the x value of the closest lower left data grid space, such 
        // that if x = 12.7, i = 12, j and k follow the same pattern
        i = (int) x;
        j = (int) y;
        k = (int) z;
        // i1, j1 and k1 are i + 1, j + 1 and k + 1
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        // It is important to notice that due to the periodic boundary
        // conditions of the data, it would normally be possible for x, y, and
        // z to range from  -.5 to (Size - .5), however since RK4 checks the
        // slopes at 4 points per step, and the RK4 method only checks that
        // x, y, and z are within the allowed range once per step, it is
        // possible for x, y, and z to be outside the range (-.5, Size - .5)
        // Therefore if x, y and z are greater than 0, it is possible for i, j,
        //and k to be to be >= Size
        // if i, j or k = Size - 1, i, j, or k +1 must = 0
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        // if i, j or k > Size - 1, i, j, k, i1, j1, and k1 will all move to
        // the other side of the simulation, so Size is subtracted from their
        // value
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
    }
    // if x is < 0 conditions are different, y and z follow similar patterns
    else if (y >= 0 && z >= 0){
        // finding Wx for negative x such that if x = -4.7 Wx = .3
        Wx = 1 - ((int) x - x);   
        // And i = -5     
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        // if i is = -1, i + 1 = 0 and i = Xsize - 1
        if (i == -1){
            i = Xsize - 1;      
        }
        // if i < -1, both i and i + 1 move to the other side of the simulation
        // so Xsize is added to their value
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
    }
    else if (x >= 0 && z >= 0){
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wx = fmod(x, 1);
        i = (int) x;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }
    else if (y >= 0 && x >= 0){
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wx = fmod(x, 1);
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else if (z >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (k == (Zsize - 1)){
            k1 = 0;
        }
        if (k > (Zsize - 1)){
            k = k - Zsize;            
            k1 = k1 - Zsize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }
    else if (y >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else if (x >= 0){
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        Wx = fmod(x, 1);
        i = (int) x;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    else{
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wz = 1 - ((int) z - z);
        k = (int) z;
        i1 = i + 1;
        j1 = j + 1;
        k1 = k + 1;
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
        if (k == -1){
            k = Zsize - 1;      
        }
        if (k < -1){
            k = k + Zsize;
            k1 = k1 + Zsize;        
        }
    }
    
    // With Wx, Wy, Wz, i, j, k, i1, j1, and k1 determined the interpolated
    // values of components of B at point [x, y, z] can be found from their
    // values at the nearest 8 data points
    double B_Wx = (1-Wx)*(1-Wy)*(1-Wz)*Bx(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bx(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bx(i1,j,k) + Wx*Wy*(1-Wz)*Bx(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bx(i,j,k1) + (1-Wx)*Wy*Wz*Bx(i,j1,k1) + Wx*(1-Wy)*Wz*Bx(i1,j,k1) + Wx*Wy*Wz*Bx(i1,j1,k1);
    double B_Wy = (1-Wx)*(1-Wy)*(1-Wz)*By(i,j,k) + (1-Wx)*Wy*(1-Wz)*By(i,j1,k) + Wx*(1-Wy)*(1-Wz)*By(i1,j,k) + Wx*Wy*(1-Wz)*By(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*By(i,j,k1) + (1-Wx)*Wy*Wz*By(i,j1,k1) + Wx*(1-Wy)*Wz*By(i1,j,k1) + Wx*Wy*Wz*By(i1,j1,k1);
    double B_Wz = (1-Wx)*(1-Wy)*(1-Wz)*Bz(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bz(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bz(i1,j,k) + Wx*Wy*(1-Wz)*Bz(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bz(i,j,k1) + (1-Wx)*Wy*Wz*Bz(i,j1,k1) + Wx*(1-Wy)*Wz*Bz(i1,j,k1) + Wx*Wy*Wz*Bz(i1,j1,k1);
    // Finding |B| at the interpolated point
    double B_Wm = sqrt(pow(B_Wx, 2) + pow(B_Wy, 2) + pow(B_Wz, 2));
    
    // normalization of components
    B[0] = B_Wx/B_Wm;
    B[1] = B_Wy/B_Wm;
    B[2] = B_Wz/B_Wm;

}


// Fourth Order Runge-Kutta procedure for 3D
int RK4_3D(double * Line_X, double * Line_Y, double * Line_Z, double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float ds, int steps, double * Ex, double * Ey, double * Ez, double * EI){
    // Value to hold interpolated E along trace path    
    double interp = 0;
    // Determining whether or not to interpolate E
    int Eint = 0;
    if ((Ex[0] == 0) && (Ex[1] == 0) && (Ex[2] == 0)){
        Eint = 0;
    }
    else{
        Eint = 1;
    }

    // Starting point
    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;

    // Slope checking variables
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];

    // Loop to step forward from initial point
    int step;
    for (step = 0; step < steps; step++){
        
        // Ensuring X, Y, and Z are within the allowed range of 
        // (-.5, Size - .5)
        if (X < -.5){
            X = X + Xsize;
        }
        if (Y < -.5){
            Y = Y + Ysize;
        }
        if (Z < -.5){
            Z = Z + Ysize;
        }
        if (X > (Xsize - .5)){
            X = X - Xsize;
        }
        if (Y > (Ysize - .5)){
            Y = Y - Ysize;
        }
        if (Z > (Zsize - .5)){
            Z = Z - Zsize;
        }
        
        // E interpolation
        if (Eint == 1){
            CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
            CalcSlopes3E(X, Y, Z, Ex, Ey, Ez, Xsize, Ysize, Zsize, Slopes2);
            // interp total = 
            // Sum(((Ex * Bx/|B|) + (Ey * By/|B|) + (Ez * Bz/|B|)) * ds)
            interp = interp + (((Slopes[0]*Slopes2[0]) + (Slopes[1]*Slopes2[1]) + (Slopes[2]*Slopes2[2]))*ds);
            EI[step] = interp;
        }

        // Adding data points
        Line_X[step] = X;
        Line_Y[step] = Y;
        Line_Z[step] = Z;
        
        // RK4

        //CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        
        //CalcSlopes3(X + (ds/2)*K1x , Y + (ds/2)*K1y, Z + (ds/2)*K1z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        CalcSlopes3(X + (ds/2)*Slopes[0] , Y + (ds/2)*Slopes[1], Z + (ds/2)*Slopes[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        
        //CalcSlopes3(X + (ds/2)*K2x , Y + (ds/2)*K2y, Z + (ds/2)*K2z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        CalcSlopes3(X + (ds/2)*Slopes2[0] , Y + (ds/2)*Slopes2[1], Z + (ds/2)*Slopes2[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        
        //CalcSlopes3(X + ds*K3x , Y + ds*K3y, Z + ds*K3z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);
        CalcSlopes3(X + ds*Slopes3[0] , Y + ds*Slopes3[1], Z + ds*Slopes3[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);

        //X = X + (ds/6)*(K1x + 2*K2x + 2*K3x + K4x);
        //Y = Y + (ds/6)*(K1y + 2*K2y + 2*K3y + K4y);
        //Z = Z + (ds/6)*(K1z + 2*K2z + 2*K3z + K4z);
        X = X + (ds/6)*(Slopes[0] + 2*Slopes2[0] + 2*Slopes3[0] + Slopes4[0]);
        Y = Y + (ds/6)*(Slopes[1] + 2*Slopes2[1] + 2*Slopes3[1] + Slopes4[1]);
        Z = Z + (ds/6)*(Slopes[2] + 2*Slopes2[2] + 2*Slopes3[2] + Slopes4[2]);

    }
    // Returns E interpolation total
    return interp;
}


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// PUNCTURE METHODS

// Puncture point locating method, PunctAxis is the line trace data of the
// axis to be punctured along
// For instance, if the user wants to plot punctures of the plane Z = 100,
// PunctAxis would contain the Z line points, Val would = 100, OtherAxis1 and 
// 2 will contain the X and Y line points, Points1 and 2 hold the X and Y 
// values of the puncture points
int Punct(double * PunctAxis, double Val, float ds, int Steps, double * OtherAxis1, double * OtherAxis2, double * Points1, double * Points2){
    // Variable to count total puncture points    
    int Points = 0;
    int i;
    // Loop to find puncture points
    for (i = 0; i < Steps; i++){
        // If a point in the line data is within  ds/2 of the value of the
        // plane to be punctured, record that points locationin the puncture
        // plane
        if (((Val - (ds/2)) < PunctAxis[i]) && (PunctAxis[i] < (Val + (ds/2)))){
            Points1[Points] = OtherAxis1[i];
            Points2[Points] = OtherAxis2[i];
            // And increment the Points total
            Points++;
        }
    }
    // Returns the total number of points found
    return Points;
}


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// SEPARATOR MAPPING METHODS

// Tracing method used to find separator points
// Identical to the RK4_3D method, accept this one will cut off the trace
// if reconnection is detected, and no E interpolation for speed
int RK4_SepUP(double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float ds, int steps){

    // Variable to check for reconnection
    double SlopeCheck[3];

    // Variable to determine if the inflow region above the separator has
    // positive or negative X component
    CalcSlopes3(Xinit, Yinit + 20, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, SlopeCheck);
    double InitXComp = SlopeCheck[0];

    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;
    
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];

    int step;
    for (step = 0; step < steps; step++){

        if (X < -.5){
            X = X + Xsize;
        }
        if (Y < -.5){
            Y = Y + Ysize;
        }
        if (Z < -.5){
            Z = Z + Ysize;
        }
        if (X > (Xsize - .5)){
            X = X - Xsize;
        }
        if (Y > (Ysize - .5)){
            Y = Y - Ysize;
        }
        if (Z > (Zsize - .5)){
            Z = Z - Zsize;
        }
        
        // The reconnection check
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, SlopeCheck);
        // If at any point the X component of B at point [X, Y, Z] is the
        // opposite of InitXComp, the line is either reconnecting, or has gone
        // past the reconnection zone into the other inflow region, return 0
        if (InitXComp < 0){
            if (SlopeCheck[0] > 0){
               return 0;
            }
        }
        else{
            if (SlopeCheck[0] < 0){
               return 0;
            }
        }

        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        CalcSlopes3(X + (ds/2)*Slopes[0] , Y + (ds/2)*Slopes[1], Z + (ds/2)*Slopes[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        CalcSlopes3(X + (ds/2)*Slopes2[0] , Y + (ds/2)*Slopes2[1], Z + (ds/2)*Slopes2[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        CalcSlopes3(X + ds*Slopes3[0] , Y + ds*Slopes3[1], Z + ds*Slopes3[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);

        X = X + (ds/6)*(Slopes[0] + 2*Slopes2[0] + 2*Slopes3[0] + Slopes4[0]);
        Y = Y + (ds/6)*(Slopes[1] + 2*Slopes2[1] + 2*Slopes3[1] + Slopes4[1]);
        Z = Z + (ds/6)*(Slopes[2] + 2*Slopes2[2] + 2*Slopes3[2] + Slopes4[2]);

    }
    // If the loop complete without the X component of B ever being the
    // opposite of InitXComp, the line is in the inflow region, return 1
    return 1;

}


// Identical to the RK4_SepUP accept the reconnection check is reversed
int RK4_SepLOW(double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float ds, int steps){

    double SlopeCheck[3];
    
    // Variable to determine if the inflow region below the separator has
    // positive or negative X component
    CalcSlopes3(Xinit, Yinit - 20, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, SlopeCheck);
    double InitXComp = SlopeCheck[0];

    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;
    
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];

    int step;
    for (step = 0; step < steps; step++){

        if (X < -.5){
            X = X + Xsize;
        }
        if (Y < -.5){
            Y = Y + Ysize;
        }
        if (Z < -.5){
            Z = Z + Ysize;
        }
        if (X > (Xsize - .5)){
            X = X - Xsize;
        }
        if (Y > (Ysize - .5)){
            Y = Y - Ysize;
        }
        if (Z > (Zsize - .5)){
            Z = Z - Zsize;
        }
     
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, SlopeCheck);
        if (InitXComp < 0){
            if (SlopeCheck[0] > 0){
               return 0;
            }
        }
        else{
            if (SlopeCheck[0] < 0){
               return 0;
            }
        }

        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        CalcSlopes3(X + (ds/2)*Slopes[0] , Y + (ds/2)*Slopes[1], Z + (ds/2)*Slopes[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        CalcSlopes3(X + (ds/2)*Slopes2[0] , Y + (ds/2)*Slopes2[1], Z + (ds/2)*Slopes2[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        CalcSlopes3(X + ds*Slopes3[0] , Y + ds*Slopes3[1], Z + ds*Slopes3[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);

        X = X + (ds/6)*(Slopes[0] + 2*Slopes2[0] + 2*Slopes3[0] + Slopes4[0]);
        Y = Y + (ds/6)*(Slopes[1] + 2*Slopes2[1] + 2*Slopes3[1] + Slopes4[1]);
        Z = Z + (ds/6)*(Slopes[2] + 2*Slopes2[2] + 2*Slopes3[2] + Slopes4[2]);

    }
    
    return 1;

}


// Method used to find points on the separator surface by tracing points
// above and below and testing if they reconnect
double SepPoints(int N, double * Start, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float ds, int Steps, double * SeparatorX, double * SeparatorY, double * SeparatorZ, int UPorLOW){
    // Increment specifying the distance between the two points to be tested
    float inc = 4;
    // Variable to hold the initial Y value of the starting point
    double Y0;
    Y0 = Start[1];
    // Variables to specify whether tested lines 1 or 2 are in the inflow
    // region, inflow = 1 means it is inflow, inflow = 0 means reconnecting
    int inflow;
    int inflow2;

    // Testing if the program should be seaching for a separator above or
    // below the reconnection zone, above = 1, below = 0
    if (UPorLOW == 1){
        // Loop to repeatedly test and adjust the lower point until it is in
        // the reconnection zone
        while (1 == 1){
            // Testing the lower point
            inflow = RK4_SepUP(Start[0], Start[1], Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, ds, (int) Steps);
            // if it is reconnecting, break, move on to the upper point
            if (inflow == 0){
                break;
            }
            // Fails if the reconnection zone cannnot be found within 20 grid
            // spaces of the starting Y value
            if (Start[1] < (Y0 - 20)){
                printf("Failed most likely because N was too big or Yinit was not close enough");
                return;    
            }
            // If inflow = 1  you know the lower point is in the inflow region,
            // so shrink the increment, and test 1 grid space lower until the
            // lower point reconnects
            inc = 1;
            Start[1] = Start[1] - 1;
                  
        }
        // Once the lower point is in the reconnection zone, test the upper
        // point until the upper point is in the inflow region, the lower point
        // is in the reconnection zone, and the increment between them is <= .5
        while (1 == 1){
            // Testing the upper point
            inflow2 = RK4_SepUP(Start[0], Start[1] + inc, Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, ds, (int) Steps);
            // If it is in the inflow region, either reduce the size of the
            // increment, or break if the increment is <= .5, this means
            // the separator sheet is between the upper and lower point
            if (inflow2 == 1){
                if (inc <= .5){
                    // If the lower point is reconnecting, the upper point
                    //  is inflow, and the inc < = .5, pick a Y value half
                    // way between the two points, ensuring this Y value is
                    // within .25 grid spaces of the separator sheets actual
                    // location
                    Start[1] = Start[1] + (inc/2);
                    break;
                }
                inc = inc/2;        
            }
            // if the upper point is brought down too low, and reconnects,
            // move the lower point up to the upper points position
            else if (inflow2 == 0){
                Start[1] = Start[1] + inc;
            }
        }
        // Return the Y value of the separator sheet at this point
        return Start[1];
    }
    // If looking for a separator sheet below the reconnection zone, do
    // the same as above upside down
    else{
        while (1 == 1){
            inflow = RK4_SepLOW(Start[0], Start[1], Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, ds, (int) Steps);
            if (inflow == 0){
                break;
            }
            if (Start[1] > (Y0 + 20)){
                printf("Failed most likely because N was too big or Yinit was not close enough");
                return;    
            }
            inc = 1;
            Start[1] = Start[1] + 1;
                  
        }
        while (1 == 1){
            inflow2 = RK4_SepLOW(Start[0], Start[1] - inc, Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, ds, (int) Steps);
            if (inflow2 == 1){
                if (inc <= .5){
                    Start[1] = Start[1] - (inc/2);
                    break;
                }
                inc = inc/2;        
            }
            else if (inflow2 == 0){
                Start[1] = Start[1] - inc;
            }
        }
        return Start[1];
    }
}


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// STANDARD 2D TRACING METHODS

// 2D interpolation / slope calculating method, works identically to
// CalcSlopes3 accept only interpolates nearest 4 points, not 8 since 2D
void CalcSlopes2(double x, double y, float * Bx, float * By, int Xsize, int Ysize, double * B){ 

    double Wx, Wy;
    int i, j, i1, j1;
    
    double Bx_ij, Bx_i1j, Bx_ij1, Bx_i1j1, By_ij, By_i1j, By_ij1, By_i1j1;
    
    if ((x >= 0) && (y >= 0)){
        Wx = fmod(x, 1);
        Wy = fmod(y, 1);
        i = (int) x;
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
    }
    else if (y >= 0){
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        if (j == (Ysize - 1)){
            j1 = 0;
        }
        if (j > (Ysize - 1)){
            j = j - Ysize;            
            j1 = j1 - Ysize;
        }
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
    }
    else if (x >= 0){
        Wy = 1 - ((int) y - y);
        j = (int) y;
        Wx = fmod(x, 1);
        i = (int) x;
        i1 = i + 1;
        j1 = j + 1;
        if (i == (Xsize - 1)){
            i1 = 0;
        }
        if (i > (Xsize - 1)){
            i = i - Xsize;            
            i1 = i1 - Xsize;
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }
    else{
        Wx = 1 - ((int) x - x);        
        i = (int) x;
        Wy = 1 - ((int) y - y);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        if (i == -1){
            i = Xsize - 1;      
        }
        if (i < -1){
            i = i + Xsize;
            i1 = i1 + Xsize;        
        }
        if (j == -1){
            j = Ysize - 1;      
        }
        if (j < -1){
            j = j + Ysize;
            j1 = j1 + Ysize;        
        }
    }

    Bx_ij = Bx2(i,j);
    Bx_i1j = Bx2(i1,j);
    Bx_ij1 = Bx2(i,j1);
    Bx_i1j1 = Bx2(i1,j1);
    
    By_ij = By2(i,j);
    By_i1j = By2(i1,j);
    By_ij1 = By2(i,j1);
    By_i1j1 = By2(i1,j1);

    double B_Wx = (1-Wx)*(1-Wy)*Bx_ij + (1-Wx)*Wy*Bx_ij1 + Wx*(1-Wy)*Bx_i1j + Wx*Wy*Bx_i1j1;
    double B_Wy = (1-Wx)*(1-Wy)*By_ij + (1-Wx)*Wy*By_ij1 + Wx*(1-Wy)*By_i1j + Wx*Wy*By_i1j1;
    double B_Wm = sqrt(pow(B_Wx, 2) + pow(B_Wy, 2));

    B[0] = B_Wx/B_Wm;
    B[1] = B_Wy/B_Wm;

}


// Fourth Order Runge-Kutta procedure for 2D, works similarly to 3D accept
// no Z component, and the trace does not run for the total number of steps
// instead it is cut off when the line bites its own tail, so this method
// returns the number of steps required to complete the line rather than 
// the interpolation total
int RK4_2D(double* Line_X, double * Line_Y, double Xinit, double Yinit, float * Bx, float * By, int Xsize, int Ysize, float ds, int steps){

    double K1x, K2x, K3x, K4x, K1y, K2y, K3y, K4y = 0;

    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];
  
    double X;
    double Y;
    
    X = Xinit;
    Y = Yinit;
    
    int step;    

    // Variables to determine when to cut off the trace
    int killin10000 = 0;
    int kill = 0;
    // Variable to trim line data arrays to proper length
    int ArrayLength;

    for (step = 0; step < steps; step++){
        if (X < -.5){
            X = X + Xsize;
        }
        if (Y < -.5){
            Y = Y + Ysize;
        }
        if (X > (Xsize - .5)){
            X = X - Xsize;
        }
        if (Y > (Ysize - .5)){
            Y = Y - Ysize;
        }
     
        Line_X[step] = X;
        Line_Y[step] = Y;

        // Lines with low curvature will usually bite their own tail within ds
        // of their starting point
        if ((step > 5) && (fabs(X - Line_X[0]) < ds) && (fabs(Y - Line_Y[0]) < ds)){
            // Setting the length to trim the line data to
            ArrayLength = step + 1;
            // returning this value 
            return ArrayLength;
            break;
        }
        // Lines of greater curvature will not bite their own tail exactly.
        // Instead they will pass close to their starting point,
        // there is however no way to distinguish between a line which has
        // missed its starting point by a small distance, and one which will
        // perfectly bite its own tail but simply will not do so for a few more
        // steps.  To solve this when any line passes within a small distance
        // of its starting location, kill is set to 1, such that after each
        // step killin10000 will be incremented by 1, if killin10000 reaches 
        // 10000 before the line gets any closer to its starting point, the
        // trace will be cut off, and the data arrays will be trimmed to the
        // number of steps taken before the killin10000 stage began (when the 
        // line was the closest it will get to its starting position).  
        // However, if the line passes closer to its starting position before 
        // killin10000 reaches 10000, the trace will cut off there.
        else if ((step > 100) && (fabs(X - Line_X[0]) < (25*ds)) && (fabs(Y - Line_Y[0]) < (25*ds))){
            ArrayLength = step + 1;
            killin10000 = 0;
            kill = 1;
        }
        // Even the highest curvature lines should pass within 250*ds of their
        // starting position, found by trial and error, not certain
        else if ((step > 1000) && (fabs(X - Line_X[0]) < (250*ds)) && (fabs(Y - Line_Y[0]) < (250*ds))){
            ArrayLength = step + 1;
            killin10000 = 0;
            kill = 1;
            
        }
        // if the line does not pass closer to its starting position, return
        // the value of Array length when the killin10000 process began
        if ((kill == 1) && (killin10000 == 10000)){
            return ArrayLength;
            break;
        }

        //RK4

        CalcSlopes2(X, Y, Bx, By, Xsize, Ysize, Slopes);
        K1x = Slopes[0];
        K1y = Slopes[1];
        
        CalcSlopes2(X + (ds/2)*K1x, Y + (ds/2)*K1y, Bx, By, Xsize, Ysize, Slopes2);
        K2x = Slopes2[0];
        K2y = Slopes2[1];
        
        CalcSlopes2(X + (ds/2)*K2x, Y + (ds/2)*K2y, Bx, By, Xsize, Ysize, Slopes3);
        K3x = Slopes3[0];
        K3y = Slopes3[1];
        
        CalcSlopes2(X + ds*K3x, Y + ds*K3y, Bx, By, Xsize, Ysize, Slopes4);
        K4x = Slopes4[0];
        K4y = Slopes4[1];
        
        X = X + (ds/6)*(K1x + 2*K2x + 2*K3x + K4x);
        Y = Y + (ds/6)*(K1y + 2*K2y + 2*K3y + K4y);

        // incrementing killin10000
        killin10000++;
    }
    // if the max number of steps is reached, the line has not passed within
    // 250*ds of its starting position and therefore is considered incomplete
    printf("Warning, line not completed");
    return steps;
}


