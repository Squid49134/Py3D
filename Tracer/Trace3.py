
# 3D Field Line Tracing Program
# takes about 74 seconds per Million Steps

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# loading data
print('loading data...')
Bx =  np.load('/scratch-fast/asym030/bx.npy')
By =  np.load('/scratch-fast/asym030/by.npy')
Bz =  np.load('/scratch-fast/asym030/bz.npy')
xx = np.linspace(0,51.2, num = 2048)
yy = np.linspace(0,25.6, num = 1024)
zz = np.linspace(0,25.6, num = 1024)

def CalcSlopes(x, y, z, Bx, By, Bz):
# distance between field line point (x,y,z) and lower left data grid 
    #point (i,j,k), such that Wx = x - i, Wy = y - j, Wz = z - k
    if x > 0 and y > 0 and z > 0:
        Wx = x % 1
        Wy = y % 1
        Wz = z % 1
        # corresponding closest lower left data grid point (i,j,k) to field line
        # point (x,y,z)
        i = int(x)
        j = int(y)
        k = int(z)
        i1 = i + 1
        j1 = j + 1
        k1 = k + 1
        if i == 2047:
            i1 = 0
        if j == 1023:
            j1 = 0
        if k == 1023:
            k1 = 0
    elif y > 0 and z > 0:
        i = 2047
        Wx = 1 + x
        Wy = y % 1
        j = int(y)
        Wz = z % 1
        k = int(z)
        i1 = 0
        j1 = j + 1
        k1 = k + 1
        if j == 1023:
            j1 = 0
        if k == 1023:
            k1 = 0
    elif x > 0 and z > 0:
        j = 1023
        Wy = 1 + y
        Wx = x % 1
        i = int(x)
        Wz = z % 1
        k = int(z)
        i1 = i + 1
        j1 = 0
        k1 = k + 1
        if i == 2047:
            i1 = 0
        if k == 1023:
            k1 = 0
    elif y > 0 and x > 0:
        k = 1023
        Wz = 1 + z
        Wx = x % 1
        i = int(x)
        Wy = y % 1
        j = int(y)
        i1 = i + 1
        j1 = j + 1
        k1 = 0
        if i == 2047:
            i1 = 0
        if j == 1023:
            j1 = 0
    elif z > 0:
        i = 2047
        Wx = 1 + x
        j = 1023
        Wy = 1 + y
        Wz = z % 1
        k = int(z)
        i1 = 0
        j1 = 0
        k1 = k + 1
        if k == 1023:
            k1 = 0
    elif y > 0:
        i = 2047
        Wx = 1 + x
        k = 1023
        Wz = 1 + z
        Wy = y % 1
        j = int(y)
        i1 = 0
        j1 = j + 1
        k1 = 0
        if j == 1023:
            j1 = 0
    elif x > 0:
        k = 1023
        Wz = 1 + z
        j = 1023
        Wy = 1 + y
        Wx = x % 1
        i = int(x)
        i1 = i + 1
        j1 = 0
        k1 = 0
        if i == 2047:
            i1 = 0
    else:
        i = 2047
        Wx = 1 + x
        j = 1023
        Wy = 1 + y
        k = 1023
        Wz = 1 + z
        i1 = 0
        j1 = 0
        k1 = 0
        
        # identifying Bx, By and Bz at closest 4 data grid points
    Bx_ijk = Bx[i,j,k]
    Bx_i1jk = Bx[i1,j,k]
    Bx_ij1k = Bx[i,j1,k]
    Bx_i1j1k = Bx[i1,j1,k]
    Bx_ijk1 = Bx[i,j,k1]
    Bx_i1jk1 = Bx[i1,j,k1]
    Bx_ij1k1 = Bx[i,j1,k1]
    Bx_i1j1k1 = Bx[i1,j1,k1]
    
    By_ijk = By[i,j,k]
    By_i1jk = By[i1,j,k]
    By_ij1k = By[i,j1,k]
    By_i1j1k = By[i1,j1,k]
    By_ijk1 = By[i,j,k1]
    By_i1jk1 = By[i1,j,k1]
    By_ij1k1 = By[i,j1,k1]
    By_i1j1k1 = By[i1,j1,k1]
    
    Bz_ijk = Bz[i,j,k]
    Bz_i1jk = Bz[i1,j,k]
    Bz_ij1k = Bz[i,j1,k]
    Bz_i1j1k = Bz[i1,j1,k]
    Bz_ijk1 = Bz[i,j,k1]
    Bz_i1jk1 = Bz[i1,j,k1]
    Bz_ij1k1 = Bz[i,j1,k1]
    Bz_i1j1k1 = Bz[i1,j1,k1]

    # finding average Bx, By, Bz and Bm from 4 nearest neighboring data 
    # points (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), (i,j,k+1), 
    # (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (x,y,z)
    B_Wx = (1-Wx)*(1-Wy)*(1-Wz)*Bx_ijk + (1-Wx)*Wy*(1-Wz)*Bx_ij1k + Wx*(1-Wy)*(1-Wz)*Bx_i1jk + Wx*Wy*(1-Wz)*Bx_i1j1k + (1-Wx)*(1-Wy)*Wz*Bx_ijk1 + (1-Wx)*Wy*Wz*Bx_ij1k1 + Wx*(1-Wy)*Wz*Bx_i1jk1 + Wx*Wy*Wz*Bx_i1j1k1
    B_Wy = (1-Wx)*(1-Wy)*(1-Wz)*By_ijk + (1-Wx)*Wy*(1-Wz)*By_ij1k + Wx*(1-Wy)*(1-Wz)*By_i1jk + Wx*Wy*(1-Wz)*By_i1j1k + (1-Wx)*(1-Wy)*Wz*By_ijk1 + (1-Wx)*Wy*Wz*By_ij1k1 + Wx*(1-Wy)*Wz*By_i1jk1 + Wx*Wy*Wz*By_i1j1k1
    B_Wz = (1-Wx)*(1-Wy)*(1-Wz)*Bz_ijk + (1-Wx)*Wy*(1-Wz)*Bz_ij1k + Wx*(1-Wy)*(1-Wz)*Bz_i1jk + Wx*Wy*(1-Wz)*Bz_i1j1k + (1-Wx)*(1-Wy)*Wz*Bz_ijk1 + (1-Wx)*Wy*Wz*Bz_ij1k1 + Wx*(1-Wy)*Wz*Bz_i1jk1 + Wx*Wy*Wz*Bz_i1j1k1
    B_Wm = np.sqrt(B_Wx**2 + B_Wy**2 + B_Wz**2)
    
    return B_Wx/B_Wm, B_Wy/B_Wm, B_Wz/B_Wm

# line tracing method
def Line(X, Y, Z): 
    
    # maximum number of steps along line
    MaxSteps = 500000
    
    K1x = 0
    K2x = 0
    K3x = 0
    K4x = 0
    
    K1y = 0
    K2y = 0
    K3y = 0
    K4y = 0
    
    K1z = 0
    K2z = 0
    K3z = 0
    K4z = 0
    
    # arrays to hold X, Y and Z coordinates of field line points
    Line_X = np.zeros(MaxSteps)
    Line_Y = np.zeros(MaxSteps)
    Line_Z = np.zeros(MaxSteps)
    
    # differential step towards next point on field line
    # dx has little effect on accuracy of tracing, accuracy is most strongly
    # affected by the curl of the field line
    dx = .1

    # loop to step forward line from initial point
    for step in range(0,MaxSteps):
        
        # checking if X, Y, or Z has moved outside range (-.5 to size-.5) and if so
        # switches to other side of data grid
        if X < -.5:
            X = X + 2048
        if Y < -.5:
            Y = Y + 1024
        if Z < -.5:
            Z = Z + 1024
        if X > 2047.5:
            X = X - 2048
        if Y > 1023.5:
            Y = Y - 1024
        if Z > 1023.5:
            Z = Z - 1024
        
        # just a status update, sometimes it takes a long time
        if (step % 100000) == 0 and step != 0:
            print("Step = " + str(step))
        
        # adding points (X,Y) to field line divided by 20 for real space
        Line_X[step] = X/40
        Line_Y[step] = Y/40
        Line_Z[step] = Z/40
            
        # finding next point on Line
        Slopes = CalcSlopes(X, Y, Z, Bx, By, Bz)
        K1x = Slopes[0]
        K1y = Slopes[1]
        K1z = Slopes[2]
        
        Slopes2 = CalcSlopes(X + (dx/2)*K1x , Y + (dx/2)*K1y, Z + (dx/2)*K1z, Bx, By, Bz)
        K2x = Slopes2[0]
        K2y = Slopes2[1]
        K2z = Slopes2[2]
        
        Slopes3 = CalcSlopes(X + (dx/2)*K2x , Y + (dx/2)*K2y, Z + (dx/2)*K2z, Bx, By, Bz)
        K3x = Slopes3[0]
        K3y = Slopes3[1]
        K3z = Slopes3[2]
        
        Slopes4 = CalcSlopes(X + dx*K3x , Y + dx*K3y, Z + dx*K3z, Bx, By, Bz)
        K4x = Slopes4[0]
        K4y = Slopes4[1]
        K4z = Slopes4[2]

        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x)
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y)
        Z = Z + (dx/6)*(K1z + 2*K2z + 2*K3z + K4z)

    # returns line as three arrays of X, Y and Z points
    return Line_X, Line_Y, Line_Z


print('tracing...')
# the figure
fig1 = plt.figure(1)
fig1.set_size_inches(30,8, forward = True)
# making 3D plot
ax = fig1.add_subplot(131, projection = '3d')
# 500, 180, 250 is inside reconn zone
Contour = Line(1277,275,500)
print('plotting...')
ax.plot(Contour[0],Contour[1],Contour[2], linestyle = 'none', marker = '.', markersize = .1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# elev = 0 and azim = 0 looks down X axis
# elev = 90 and azim = 270 looks down Z axis
# elev = 0 and azim = 270 looks down Y axis
# BE CAREFUL AXIS LIMITS ARE BACKWARDS ie 51.2 to 0 not 0 to 51.2, this is done
# to match plots of J.transpose
# this plot looks down X axis
ax.view_init(elev = 0, azim = 0)
ax.set_zlim([25.6,0])
ax.set_ylim([25.6,0])
ax.set_xlim([51.2,0])

ax2 = fig1.add_subplot(132, projection = '3d')
ax2.plot(Contour[0],Contour[1],Contour[2], linestyle = 'none', marker = '.', markersize = .1)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
# looking down Y axis
ax2.view_init(elev = 0, azim = 270)
ax2.set_zlim([25.6,0])
ax2.set_ylim([25.6,0])
ax2.set_xlim([51.2,0])

ax3 = fig1.add_subplot(133, projection = '3d')
ax3.plot(Contour[0],Contour[1],Contour[2], linestyle = 'none', marker = '.', markersize = .1)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')
# looking down Z axis
ax3.view_init(elev = 90, azim = 270)
ax3.set_zlim([25.6,0])
ax3.set_ylim([25.6,0])
ax3.set_xlim([51.2,0])

fig1.tight_layout()

plt.show()

# TODO:
# C implementation
# figure out debye length

