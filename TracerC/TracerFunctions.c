/* C functions for the test particle module */
/* Incase you forget how to compile the shared lib file
 * gcc -o TracerFunctions.so -shared -fPIC TracerFunctions.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Bx(i,j,k)  Bx[(i) + (j)*SizeX + (k)*SizeX*SizeY]
#define By(i,j,k)  By[(i) + (j)*SizeX + (k)*SizeX*SizeY]
#define Bz(i,j,k)  Bz[(i) + (j)*SizeX + (k)*SizeX*SizeY]


// 3D slope calculator for RK4 procedure
//double * CalcSlopes3(double x, double y, double z, double * Bx, double * By,
void CalcSlopes3(double x, double y, double z, double * Bx, double * By, double * Bz, int SizeX, int SizeY, int SizeZ, double * B){
    // first finds closest LOWER LEFT data grid point (i,j,k) to field line 
    // point (x,y,z) and distance between them such that Wx = x - i, Wy = y - j, 
    // Wz = z - k for interpolation of field between grid points
    double Wx, Wy, Wz;
    int i, j, k, i1, j1, k1;
    
    double Bx_ijk, Bx_i1jk, Bx_ij1k, Bx_i1j1k, Bx_ijk1, Bx_i1jk1, Bx_ij1k1, Bx_i1j1k1, By_ijk, By_i1jk, By_ij1k, By_i1j1k, By_ijk1, By_i1jk1, By_ij1k1, By_i1j1k1, Bz_ijk, Bz_i1jk, Bz_ij1k, Bz_i1j1k, Bz_ijk1, Bz_i1jk1, Bz_ij1k1, Bz_i1j1k1;
    
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
        if (i == (SizeX - 1)){
            i1 = 0;
        }
        if (j == (SizeY - 1)){
            j1 = 0;
        }
        if (k == (SizeZ - 1)){
            k1 = 0;
        }
    }
    else if (y >= 0 && z >= 0){
        i = (SizeX - 1);
        Wx = 1 + x;
        Wy = fmod(y, 1);
        j = (int) y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = 0;
        j1 = j + 1;
        k1 = k + 1;
        if (j == (SizeY - 1)){
            j1 = 0;
        }
        if (k == (SizeZ - 1)){
            k1 = 0;
        }
    }
    else if (x >= 0 && z >= 0){
        j = (SizeY - 1);
        Wy = 1 + y;
        Wx = fmod(x, 1);
        i = (int) x;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = i + 1;
        j1 = 0;
        k1 = k + 1;
        if (i == (SizeX - 1)){
            i1 = 0;
        }
        if (k == (SizeZ - 1)){
            k1 = 0;
        }
    }
    else if (y >= 0 && x >= 0){
        k = (SizeZ - 1);
        Wz = 1 + z;
        Wx = fmod(x, 1);
        i = (int) x;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = i + 1;
        j1 = j + 1;
        k1 = 0;
        if (i == (SizeX - 1)){
            i1 = 0;
        }
        if (j == (SizeY - 1)){
            j1 = 0;
        }
    }
    else if (z >= 0){
        i = (SizeX - 1);
        Wx = 1 + x;
        j = (SizeY - 1);
        Wy = 1 + y;
        Wz = fmod(z, 1);
        k = (int) z;
        i1 = 0;
        j1 = 0;
        k1 = k + 1;
        if (k == (SizeZ - 1)){
            k1 = 0;
        }
    }
    else if (y >= 0){
        i = (SizeX - 1);
        Wx = 1 + x;
        k = (SizeZ - 1);
        Wz = 1 + z;
        Wy = fmod(y, 1);
        j = (int) y;
        i1 = 0;
        j1 = j + 1;
        k1 = 0;
        if (j == (SizeY - 1)){
            j1 = 0;
        }
    }
    else if (x >= 0){
        k = (SizeZ - 1);
        Wz = 1 + z;
        j = (SizeY - 1);
        Wy = 1 + y;
        Wx = fmod(x, 1);
        i = (int) x;
        i1 = i + 1;
        j1 = 0;
        k1 = 0;
        if (i == (SizeX - 1)){
            i1 = 0;
        }
    }
    else{
        // if x, y or z is below 0, modulus wont work properly
        i = (SizeX - 1);
        Wx = 1 + x;
        j = (SizeY - 1);
        Wy = 1 + y;
        k = (SizeZ - 1);
        Wz = 1 + z;
        i1 = 0;
        j1 = 0;
        k1 = 0;
    }
        
    // identifying Bx, By and Bz at closest 4 data grid points
    Bx_ijk = Bx(i,j,k);
    Bx_i1jk = Bx(i1,j,k);
    Bx_ij1k = Bx(i,j1,k);
    Bx_i1j1k = Bx(i1,j1,k);
    Bx_ijk1 = Bx(i,j,k1);
    Bx_i1jk1 = Bx(i1,j,k1);
    Bx_ij1k1 = Bx(i,j1,k1);
    Bx_i1j1k1 = Bx(i1,j1,k1);
    
    By_ijk = By(i,j,k);
    By_i1jk = By(i1,j,k);
    By_ij1k = By(i,j1,k);
    By_i1j1k = By(i1,j1,k);
    By_ijk1 = By(i,j,k1);
    By_i1jk1 = By(i1,j,k1);
    By_ij1k1 = By(i,j1,k1);
    By_i1j1k1 = By(i1,j1,k1);
    
    Bz_ijk = Bz(i,j,k);
    Bz_i1jk = Bz(i1,j,k);
    Bz_ij1k = Bz(i,j1,k);
    Bz_i1j1k = Bz(i1,j1,k);
    Bz_ijk1 = Bz(i,j,k1);
    Bz_i1jk1 = Bz(i1,j,k1);
    Bz_ij1k1 = Bz(i,j1,k1);
    Bz_i1j1k1 = Bz(i1,j1,k1);

    // finding average Bx, By, Bz and Bm at field line point (x,y,z) from 8 nearest 
    // neighboring data points (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1), 
    // (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (x,y,z)
    double B_Wx = (1-Wx)*(1-Wy)*(1-Wz)*Bx_ijk + (1-Wx)*Wy*(1-Wz)*Bx_ij1k + Wx*(1-Wy)*(1-Wz)*Bx_i1jk + Wx*Wy*(1-Wz)*Bx_i1j1k + (1-Wx)*(1-Wy)*Wz*Bx_ijk1 + (1-Wx)*Wy*Wz*Bx_ij1k1 + Wx*(1-Wy)*Wz*Bx_i1jk1 + Wx*Wy*Wz*Bx_i1j1k1;
    double B_Wy = (1-Wx)*(1-Wy)*(1-Wz)*By_ijk + (1-Wx)*Wy*(1-Wz)*By_ij1k + Wx*(1-Wy)*(1-Wz)*By_i1jk + Wx*Wy*(1-Wz)*By_i1j1k + (1-Wx)*(1-Wy)*Wz*By_ijk1 + (1-Wx)*Wy*Wz*By_ij1k1 + Wx*(1-Wy)*Wz*By_i1jk1 + Wx*Wy*Wz*By_i1j1k1;
    double B_Wz = (1-Wx)*(1-Wy)*(1-Wz)*Bz_ijk + (1-Wx)*Wy*(1-Wz)*Bz_ij1k + Wx*(1-Wy)*(1-Wz)*Bz_i1jk + Wx*Wy*(1-Wz)*Bz_i1j1k + (1-Wx)*(1-Wy)*Wz*Bz_ijk1 + (1-Wx)*Wy*Wz*Bz_ij1k1 + Wx*(1-Wy)*Wz*Bz_i1jk1 + Wx*Wy*Wz*Bz_i1j1k1;
    double B_Wm = sqrt(pow(B_Wx, 2) + pow(B_Wy, 2) + pow(B_Wz, 2));
    
    // returns slopes of field and magnitude at point (x, y, z)
    B[0] = B_Wx/B_Wm;
    B[1] = B_Wy/B_Wm;
    B[2] = B_Wz/B_Wm;
}


// Core 3D tracing method (RK4)
void FieldLine3D(double * Line_X, double * Line_Y, double * Line_Z, double Xinit, double Yinit, double Zinit, double * Bx, double * By, double * Bz, int SizeX, int SizeY, int SizeZ, float dx, int steps){
    
    printf("\n tracing...");

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
            X = X + SizeX;
        }
        if (Y < -.5){
            Y = Y + SizeY;
        }
        if (Z < -.5){
            Z = Z + SizeY;
        }
        if (X > (SizeX - .5)){
            X = X - SizeX;
        }
        if (Y > (SizeY - .5)){
            Y = Y - SizeY;
        }
        if (Z > (SizeZ - .5)){
            Z = Z - SizeZ;
        }
        
        // adding points (X,Y) to field line
        Line_X[step] = X;
        Line_Y[step] = Y;
        Line_Z[step] = Z;
        
        // RK4, slightly different implementation from 2D for speed
        //Slopes = CalcSlopes3(X, Y, Z, Bx, By, Bz, SizeX, SizeY, SizeZ);
        CalcSlopes3(X, Y, Z, Bx, By, Bz, SizeX, SizeY, SizeZ, Slopes);
        K1x = Slopes[0];
        K1y = Slopes[1];
        K1z = Slopes[2];
        
        CalcSlopes3(X + (dx/2)*K1x , Y + (dx/2)*K1y, Z + (dx/2)*K1z, Bx, By, Bz, SizeX, SizeY, SizeZ, Slopes2);
        K2x = Slopes2[0];
        K2y = Slopes2[1];
        K2z = Slopes2[2];
        
        CalcSlopes3(X + (dx/2)*K2x , Y + (dx/2)*K2y, Z + (dx/2)*K2z, Bx, By, Bz, SizeX, SizeY, SizeZ, Slopes3);
        K3x = Slopes3[0];
        K3y = Slopes3[1];
        K3z = Slopes3[2];
        
        CalcSlopes3(X + dx*K3x , Y + dx*K3y, Z + dx*K3z, Bx, By, Bz, SizeX, SizeY, SizeZ, Slopes4);
        K4x = Slopes4[0];
        K4y = Slopes4[1];
        K4z = Slopes4[2];

        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x);
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y);
        Z = Z + (dx/6)*(K1z + 2*K2z + 2*K3z + K4z);

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
