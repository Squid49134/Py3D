/* C functions for the test particle module */
/* Incase you forget how to compile the shared lib file
 * gcc -o TracerFunctions.so -shared -fPIC TracerFunctions.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Bx(i,j,k)  Bx[(i) + (j)*Xsize + (k)*Xsize*Ysize]
#define By(i,j,k)  By[(i) + (j)*Xsize + (k)*Xsize*Ysize]
#define Bz(i,j,k)  Bz[(i) + (j)*Xsize + (k)*Xsize*Ysize]

#define Bx2(i,j)  Bx[(i) + (j)*Xsize]
#define By2(i,j)  By[(i) + (j)*Xsize]


// 3D slope calculator for RK4 procedure
//double * CalcSlopes3(double x, double y, double z, double * Bx, double * By,
void CalcSlopes3(double x, double y, double z, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, double * B){
    // first finds closest LOWER LEFT data grid point (i,j,k) to field line 
    // point (x,y,z) and distance between them such that Wx = x - i, Wy = y - j, 
    // Wz = z - k for interpolation of field between grid points  

    double Wx, Wy, Wz;
    int i, j, k, i1, j1, k1;
    
    if (x >= 0 && y >= 0 && z >= 0){
        // if both x, y and z are greater than 0, modulus is used to find point (i,j,k)
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
        // if x, y or z is below 0, modulus wont work properly
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

    // finding average Bx, By, Bz and Bm at field line point (x,y,z) from 8 nearest 
    // neighboring data points (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1), 
    // (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (x,y,z)
    double B_Wx = (1-Wx)*(1-Wy)*(1-Wz)*Bx(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bx(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bx(i1,j,k) + Wx*Wy*(1-Wz)*Bx(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bx(i,j,k1) + (1-Wx)*Wy*Wz*Bx(i,j1,k1) + Wx*(1-Wy)*Wz*Bx(i1,j,k1) + Wx*Wy*Wz*Bx(i1,j1,k1);
    double B_Wy = (1-Wx)*(1-Wy)*(1-Wz)*By(i,j,k) + (1-Wx)*Wy*(1-Wz)*By(i,j1,k) + Wx*(1-Wy)*(1-Wz)*By(i1,j,k) + Wx*Wy*(1-Wz)*By(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*By(i,j,k1) + (1-Wx)*Wy*Wz*By(i,j1,k1) + Wx*(1-Wy)*Wz*By(i1,j,k1) + Wx*Wy*Wz*By(i1,j1,k1);
    double B_Wz = (1-Wx)*(1-Wy)*(1-Wz)*Bz(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bz(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bz(i1,j,k) + Wx*Wy*(1-Wz)*Bz(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bz(i,j,k1) + (1-Wx)*Wy*Wz*Bz(i,j1,k1) + Wx*(1-Wy)*Wz*Bz(i1,j,k1) + Wx*Wy*Wz*Bz(i1,j1,k1);
    double B_Wm = sqrt(pow(B_Wx, 2) + pow(B_Wy, 2) + pow(B_Wz, 2));

    // returns slopes of field and magnitude at point (x, y, z)
    B[0] = B_Wx/B_Wm;
    B[1] = B_Wy/B_Wm;
    B[2] = B_Wz/B_Wm;

}

void CalcSlopes3E(double x, double y, double z, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, double * B){
    // first finds closest LOWER LEFT data grid point (i,j,k) to field line 
    // point (x,y,z) and distance between them such that Wx = x - i, Wy = y - j, 
    // Wz = z - k for interpolation of field between grid points  

    double Wx, Wy, Wz;
    int i, j, k, i1, j1, k1;
    
    if (x >= 0 && y >= 0 && z >= 0){
        // if both x, y and z are greater than 0, modulus is used to find point (i,j,k)
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
        // if x, y or z is below 0, modulus wont work properly
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

    // finding average Bx, By, Bz and Bm at field line point (x,y,z) from 8 nearest 
    // neighboring data points (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1), 
    // (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (x,y,z)
    B[0] = (1-Wx)*(1-Wy)*(1-Wz)*Bx(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bx(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bx(i1,j,k) + Wx*Wy*(1-Wz)*Bx(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bx(i,j,k1) + (1-Wx)*Wy*Wz*Bx(i,j1,k1) + Wx*(1-Wy)*Wz*Bx(i1,j,k1) + Wx*Wy*Wz*Bx(i1,j1,k1);
    B[1] = (1-Wx)*(1-Wy)*(1-Wz)*By(i,j,k) + (1-Wx)*Wy*(1-Wz)*By(i,j1,k) + Wx*(1-Wy)*(1-Wz)*By(i1,j,k) + Wx*Wy*(1-Wz)*By(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*By(i,j,k1) + (1-Wx)*Wy*Wz*By(i,j1,k1) + Wx*(1-Wy)*Wz*By(i1,j,k1) + Wx*Wy*Wz*By(i1,j1,k1);
    B[2] = (1-Wx)*(1-Wy)*(1-Wz)*Bz(i,j,k) + (1-Wx)*Wy*(1-Wz)*Bz(i,j1,k) + Wx*(1-Wy)*(1-Wz)*Bz(i1,j,k) + Wx*Wy*(1-Wz)*Bz(i1,j1,k) + (1-Wx)*(1-Wy)*Wz*Bz(i,j,k1) + (1-Wx)*Wy*Wz*Bz(i,j1,k1) + Wx*(1-Wy)*Wz*Bz(i1,j,k1) + Wx*Wy*Wz*Bz(i1,j1,k1);

}

void CalcSlopes2(double x, double y, float * Bx, float * By, int Xsize, int Ysize, double * B){ 

    double Wx, Wy;
    int i, j, i1, j1;
    
    double Bx_ij, Bx_i1j, Bx_ij1, Bx_i1j1, By_ij, By_i1j, By_ij1, By_i1j1;
    
    if ((x >= 0) && (y >= 0)){
        // if both x, y and z are greater than 0, modulus is used to find point (i,j,k)
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
        // if x, y or z is below 0, modulus wont work properly
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

    // identifying Bx, By and Bz at closest 4 data grid points
    Bx_ij = Bx2(i,j);
    Bx_i1j = Bx2(i1,j);
    Bx_ij1 = Bx2(i,j1);
    Bx_i1j1 = Bx2(i1,j1);
    
    By_ij = By2(i,j);
    By_i1j = By2(i1,j);
    By_ij1 = By2(i,j1);
    By_i1j1 = By2(i1,j1);

    // finding average Bx, By, Bz and Bm at field line point (x,y,z) from 8 nearest 
    // neighboring data points (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1), 
    // (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (x,y,z)
    double B_Wx = (1-Wx)*(1-Wy)*Bx_ij + (1-Wx)*Wy*Bx_ij1 + Wx*(1-Wy)*Bx_i1j + Wx*Wy*Bx_i1j1;
    double B_Wy = (1-Wx)*(1-Wy)*By_ij + (1-Wx)*Wy*By_ij1 + Wx*(1-Wy)*By_i1j + Wx*Wy*By_i1j1;
    double B_Wm = sqrt(pow(B_Wx, 2) + pow(B_Wy, 2));

    // returns slopes of field and magnitude at point (x, y, z)
    B[0] = B_Wx/B_Wm;
    B[1] = B_Wy/B_Wm;

}

int FieldLine2D(double* Line_X, double * Line_Y, double Xinit, double Yinit, float * Bx, float * By, int Xsize, int Ysize, float dx, int steps){
    
    // RK4 variables
    double K1x, K2x, K3x, K4x, K1y, K2y, K3y, K4y = 0;

    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];
  
    double X;
    double Y;
    // initial point for a field line
    X = Xinit;
    Y = Yinit;
    
    // loop to step forward line from initial point
    int step;    
    int killin10000 = 0;
    int kill = 0;
    int ArrayLength;
    for (step = 0; step < steps; step++){  
        // for periodic boundaries checking if X, or Y has moved outside 
        // range (-.5 to size-.5) and if so switches to other side of data grid
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
        
        // adding points (X,Y) to field line
        Line_X[step] = X;
        Line_Y[step] = Y;
        
        // breaks out of loop when X and Y are within dx of their values at 
        // the second step (first step not chosen since it has higher error)
        if ((step > 5) && (fabs(X - Line_X[0]) < dx) && (fabs(Y - Line_Y[0]) < dx)){
            ArrayLength = step + 1;
            return ArrayLength;
            break;
        }
        else if ((step > 100) && (fabs(X - Line_X[0]) < (25*dx)) && (fabs(Y - Line_Y[0]) < (25*dx))){
            ArrayLength = step + 1;
            killin10000 = 0;
            kill = 1;
        }
        else if ((step > 1000) && (fabs(X - Line_X[0]) < (250*dx)) && (fabs(Y - Line_Y[0]) < (250*dx))){
            ArrayLength = step + 1;
            killin10000 = 0;
            kill = 1;
            
        }
    
        if ((kill == 1) && (killin10000 == 10000)){
            return ArrayLength;
            break;
        }
        
        // RK4
        CalcSlopes2(X, Y, Bx, By, Xsize, Ysize, Slopes);
        K1x = Slopes[0];
        K1y = Slopes[1];
        
        CalcSlopes2(X + (dx/2)*K1x, Y + (dx/2)*K1y, Bx, By, Xsize, Ysize, Slopes2);
        K2x = Slopes2[0];
        K2y = Slopes2[1];
        
        CalcSlopes2(X + (dx/2)*K2x, Y + (dx/2)*K2y, Bx, By, Xsize, Ysize, Slopes3);
        K3x = Slopes3[0];
        K3y = Slopes3[1];
        
        CalcSlopes2(X + dx*K3x, Y + dx*K3y, Bx, By, Xsize, Ysize, Slopes4);
        K4x = Slopes4[0];
        K4y = Slopes4[1];
        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x);
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y);

        killin10000++;
    }
    printf("Warning, line not completed");
    return steps;
}

// Core 3D tracing method (RK4)
int FieldLine3D(double * Line_X, double * Line_Y, double * Line_Z, double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float dx, int steps, double * Ex, double * Ey, double * Ez, double * EI){
    double interp = 0;
    int Eint = 0;
    if ((Ex[0] == 0) && (Ex[1] == 0) && (Ex[2] == 0)){
        Eint = 0;
    }
    else{
        Eint = 1;
    }

    // initial starting points for field line
    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;
    
    // RK4 variables
    double K1x, K2x, K3x, K4x, K1y, K2y, K3y, K4y, K1z, K2z, K3z, K4z = 0;
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];
    
    // loop to step forward line from initial point
    int step;
    for (step = 0; step < steps; step++){
        
        // for periodic boundaries checking if X, Y, or Z has moved outside 
        // range (-.5 to size-.5) and if so switches to other side of data grid
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
        
        if (Eint == 1){
            CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
            CalcSlopes3E(X, Y, Z, Ex, Ey, Ez, Xsize, Ysize, Zsize, Slopes2);
            interp = interp + (((Slopes[0]*Slopes2[0]) + (Slopes[1]*Slopes2[1]) + (Slopes[2]*Slopes2[2]))*dx);
            EI[step] = interp;
        }

        // adding points (X,Y) to field line
        Line_X[step] = X;
        Line_Y[step] = Y;
        Line_Z[step] = Z;
        
        // RK4, slightly different implementation from 2D for speed
        //Slopes = CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize);
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        K1x = Slopes[0];
        K1y = Slopes[1];
        K1z = Slopes[2];
        
        CalcSlopes3(X + (dx/2)*K1x , Y + (dx/2)*K1y, Z + (dx/2)*K1z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        K2x = Slopes2[0];
        K2y = Slopes2[1];
        K2z = Slopes2[2];
        
        CalcSlopes3(X + (dx/2)*K2x , Y + (dx/2)*K2y, Z + (dx/2)*K2z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        K3x = Slopes3[0];
        K3y = Slopes3[1];
        K3z = Slopes3[2];
        
        CalcSlopes3(X + dx*K3x , Y + dx*K3y, Z + dx*K3z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);
        K4x = Slopes4[0];
        K4y = Slopes4[1];
        K4z = Slopes4[2];

        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x);
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y);
        Z = Z + (dx/6)*(K1z + 2*K2z + 2*K3z + K4z);

    }
    return interp;

}

int Punct(double * PunctAxis, double Val, float ds, int Steps, double * OtherAxis1, double * OtherAxis2, double * Points1, double * Points2){
    int i;
    int Points = 0;
    for (i = 0; i < Steps; i++){
        if (((Val - (ds/2)) < PunctAxis[i]) && (PunctAxis[i] < (Val + (ds/2)))){
            Points1[Points] = OtherAxis1[i];
            Points2[Points] = OtherAxis2[i];
            Points++;
        }
    }
    return Points;
}

int FieldLineSepUp(double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float dx, int steps){

    double SlopeCheck[3];

    // initial starting points for field line
    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;
    
    // RK4 variables
    double K1x, K2x, K3x, K4x, K1y, K2y, K3y, K4y, K1z, K2z, K3z, K4z = 0;
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];

    // loop to step forward line from initial point
    int step;
    for (step = 0; step < steps; step++){
        
        // for periodic boundaries checking if X, Y, or Z has moved outside 
        // range (-.5 to size-.5) and if so switches to other side of data grid
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
        if (SlopeCheck[0] < 0){
           return 0;
        }

        // RK4, slightly different implementation from 2D for speed
        //Slopes = CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize);
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        CalcSlopes3(X + (dx/2)*Slopes[0] , Y + (dx/2)*Slopes[1], Z + (dx/2)*Slopes[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        CalcSlopes3(X + (dx/2)*Slopes2[0] , Y + (dx/2)*Slopes2[1], Z + (dx/2)*Slopes2[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        CalcSlopes3(X + dx*Slopes3[0] , Y + dx*Slopes3[1], Z + dx*Slopes3[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);

        X = X + (dx/6)*(Slopes[0] + 2*Slopes2[0] + 2*Slopes3[0] + Slopes4[0]);
        Y = Y + (dx/6)*(Slopes[1] + 2*Slopes2[1] + 2*Slopes3[1] + Slopes4[1]);
        Z = Z + (dx/6)*(Slopes[2] + 2*Slopes2[2] + 2*Slopes3[2] + Slopes4[2]);

    }
    
    return 1;

}

int FieldLineSepLow(double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float dx, int steps){

    double SlopeCheck[3];
    // initial starting points for field line
    double X = Xinit;
    double Y = Yinit;
    double Z = Zinit;
    
    // RK4 variables
    double K1x, K2x, K3x, K4x, K1y, K2y, K3y, K4y, K1z, K2z, K3z, K4z = 0;
    double Slopes[3];
    double Slopes2[3];
    double Slopes3[3];
    double Slopes4[3];

    // loop to step forward line from initial point
    int step;
    for (step = 0; step < steps; step++){
        
        // for periodic boundaries checking if X, Y, or Z has moved outside 
        // range (-.5 to size-.5) and if so switches to other side of data grid
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
        if (SlopeCheck[0] > 0){
           return 0;
        }

        // RK4, slightly different implementation from 2D for speed
        //Slopes = CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize);
        CalcSlopes3(X, Y, Z, Bx, By, Bz, Xsize, Ysize, Zsize, Slopes);
        CalcSlopes3(X + (dx/2)*Slopes[0] , Y + (dx/2)*Slopes[1], Z + (dx/2)*Slopes[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes2);
        CalcSlopes3(X + (dx/2)*Slopes2[0] , Y + (dx/2)*Slopes2[1], Z + (dx/2)*Slopes2[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes3);
        CalcSlopes3(X + dx*Slopes3[0] , Y + dx*Slopes3[1], Z + dx*Slopes3[2], Bx, By, Bz, Xsize, Ysize, Zsize, Slopes4);

        X = X + (dx/6)*(Slopes[0] + 2*Slopes2[0] + 2*Slopes3[0] + Slopes4[0]);
        Y = Y + (dx/6)*(Slopes[1] + 2*Slopes2[1] + 2*Slopes3[1] + Slopes4[1]);
        Z = Z + (dx/6)*(Slopes[2] + 2*Slopes2[2] + 2*Slopes3[2] + Slopes4[2]);

    }
    
    return 1;

}

double SepPoints(int N, double * Start, double * Bx, double * By, double * Bz, int Xsize, int Ysize, int Zsize, float dx, int Steps, double * SeparatorX, double * SeparatorY, double * SeparatorZ, int UPorLOW){
    int i;
    int j;
    float inc = 4;
    double Y0;
    Y0 = Start[1];
    int inflow;
    int inflow2;

    if (UPorLOW == 1){
        while (1 == 1){
            inflow = FieldLineSepUp(Start[0], Start[1], Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, dx, (int) Steps);
            if (inflow == 0){
                break;
            }
            if (Start[1] < (Y0 - 20)){
                printf("Failed most likely because N was too big or Yinit was not close enough");
                return;    
            }
            inc = 1;
            Start[1] = Start[1] - 1;
                  
        }
        while (1 == 1){
            inflow2 = FieldLineSepUp(Start[0], Start[1] + inc, Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, dx, (int) Steps);
            if (inflow2 == 1){
                if (inc <= .5){
                    Start[1] = Start[1] + (inc/2);
                    break;
                }
                inc = inc/2;        
            }
            else if (inflow2 == 0){
                Start[1] = Start[1] + inc;
            }
        }
        return Start[1];
    }
    else{
        while (1 == 1){
            inflow = FieldLineSepLow(Start[0], Start[1], Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, dx, (int) Steps);
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
            inflow2 = FieldLineSepLow(Start[0], Start[1] - inc, Start[2], Bx, By, Bz, Xsize, Ysize, Zsize, dx, (int) Steps);
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


//----------------------------------------------------------------------------//

//int main(){
//    printf("\n");
//    printf("\n");
//    printf("Cleared");
//    printf("\n");
//    printf("\n");
//    return 0;
//}
