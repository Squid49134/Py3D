
# 3D Field Line Tracing Program
# 500,000 steps in 30 seconds at .04
# RK4 took 37 seconds for 500000 at .02

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.ndimage.filters import gaussian_filter



# loading data
print('loading data...')
Bx =  np.load('/scratch-fast/asym030/bx.npy')
By =  np.load('/scratch-fast/asym030/by.npy')
Bz =  np.load('/scratch-fast/asym030/bz.npy')
xx = np.linspace(0,51.2, num = 2048)
yy = np.linspace(0,25.6, num = 1024)
zz = np.linspace(0,25.6, num = 1024)

# the magnetic field components, gaussian filter is applied to imporve
# tracing accuracy
#print('filtering...')
#Bx = gaussian_filter(bx, sigma = 3, mode = 'wrap')
#By = gaussian_filter(by, sigma = 3, mode = 'wrap')
#Bz = gaussian_filter(bz, sigma = 3, mode = 'wrap')

def CalcB(x, y, z):
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
    # (i+1,j,k+1), (i,j+1,k+1) and (i+1,j+1,k+1) at field line point (X,Y,Z)
    Bwx = (1-Wx)*(1-Wy)*(1-Wz)*Bx_ijk + (1-Wx)*Wy*(1-Wz)*Bx_ij1k + Wx*(1-Wy)*(1-Wz)*Bx_i1jk + Wx*Wy*(1-Wz)*Bx_i1j1k + (1-Wx)*(1-Wy)*Wz*Bx_ijk1 + (1-Wx)*Wy*Wz*Bx_ij1k1 + Wx*(1-Wy)*Wz*Bx_i1jk1 + Wx*Wy*Wz*Bx_i1j1k1
    Bwy = (1-Wx)*(1-Wy)*(1-Wz)*By_ijk + (1-Wx)*Wy*(1-Wz)*By_ij1k + Wx*(1-Wy)*(1-Wz)*By_i1jk + Wx*Wy*(1-Wz)*By_i1j1k + (1-Wx)*(1-Wy)*Wz*By_ijk1 + (1-Wx)*Wy*Wz*By_ij1k1 + Wx*(1-Wy)*Wz*By_i1jk1 + Wx*Wy*Wz*By_i1j1k1
    Bwz = (1-Wx)*(1-Wy)*(1-Wz)*Bz_ijk + (1-Wx)*Wy*(1-Wz)*Bz_ij1k + Wx*(1-Wy)*(1-Wz)*Bz_i1jk + Wx*Wy*(1-Wz)*Bz_i1j1k + (1-Wx)*(1-Wy)*Wz*Bz_ijk1 + (1-Wx)*Wy*Wz*Bz_ij1k1 + Wx*(1-Wy)*Wz*Bz_i1jk1 + Wx*Wy*Wz*Bz_i1j1k1
    Bmag = np.sqrt(Bwx**2 + Bwy**2 + Bwz**2)
    return Bmag, Bwx, Bwy, Bwz
        
def CalcRadCurv(x, y, z, S, Bx, By, Bz, Bm):
    # make a step even smaller than the differential step
    ds = S/2
    # find slope(B) at (X,Y,Z) +/- ds, this is a tangent
    B2 = CalcB(x + (ds*Bx/Bm), y + (ds*By/Bm), z + (ds*Bz/Bm))
    B1 = CalcB(x - (ds*Bx/Bm), y - (ds*By/Bm), z - (ds*Bz/Bm))  
    # find the change in the slope of the field (the tangent)
    dTds = [(B2[1] - B1[1])/(2*ds), (B2[2] - B1[2])/(2*ds), (B2[3] - B1[3])/(2*ds)]
    # find the magnitude off this change
    MagdTds = np.sqrt(dTds[0]**2 + dTds[1]**2 + dTds[2]**2)
    # radius of curvature is the reciprocal of this magnitude, must limit
    # to avoid radiuses equal to infinity
    if MagdTds > ds:
        r = (MagdTds)**(-1)
    else:
        r = 1/ds
    return r

# line tracing method
def Line(InitX, InitY, InitZ): 
    
    # maximum number of steps along line
    MaxSteps = 5000000
    
    # slope variable to hold Bx/|B| at point (X,Y,Z)
    DeltaX = 0
    # and By/|B|
    DeltaY = 0
    # and Bz/|B|
    DeltaZ = 0
    
    # arrays to hold X, Y and Z coordinates of field line points
    Line_X = np.zeros(MaxSteps)
    Line_Y = np.zeros(MaxSteps)
    Line_Z = np.zeros(MaxSteps)    
    
    # initial point for a field line MUST BE FLOATS
    X = float(InitX)
    Y = float(InitY)
    Z = float(InitZ)
     
    # s is the differential step
    s = .04
    # da is a place holder
    da = 0
    
    # loop to step forward line from initial point
    for step in range(0,MaxSteps - 1):

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
            
        Line_X[step] = X/40
        Line_Y[step] = Y/40
        Line_Z[step] = Z/40
        
        # just a status update, sometimes it takes a long time
        if (step % 100000) == 0 and step != 0:
            print("Step = " + str(step))
            
        B = CalcB(X,Y,Z)
        B_Wm = B[0]        
        B_Wx = B[1] 
        B_Wy = B[2]
        B_Wz = B[3]
                
        # finding next point on line
        if step == 0:
            # first step is linear
            DeltaX = ((s/2)*B_Wx/B_Wm)
            DeltaY = ((s/2)*B_Wy/B_Wm)
            DeltaZ = ((s/2)*B_Wz/B_Wm)
            X = X + DeltaX
            Y = Y + DeltaY
            Z = Z + DeltaZ
        else:
            # after first step, use curvature
            # first calculate radius of curvature of the field at X,Y,Z
            R = CalcRadCurv(X,Y,Z,s, B_Wx, B_Wy, B_Wz, B_Wm)
            # then take an arc length, s, of a circle with this radius and make
            # da equal to the length of a chord from the beginning of the arc 
            # to the end.  Picture slicing an arc of length s off the top of a 
            # circle, da is then equal to the linear cut made in the circle.  
            # In a perfect circle, the slope of this cut/chord will be tangent 
            # to the circle at a point at the center of the arc
            da = 2*R*np.sin((.5*s)/R)
            # finally the X,Y,Z steps are determined by taking a step of length
            # da from point N-1 (the beginning of the arc/chord) in the direction 
            # of Slope(B) at point N (the center of the arc), towards point
            # N+1 (the end of the arc/chord)
            DeltaX = (da*B_Wx/B_Wm)
            DeltaY = (da*B_Wy/B_Wm)
            DeltaZ = (da*B_Wz/B_Wm)
            X = 40*Line_X[step - 1] + DeltaX
            Y = 40*Line_Y[step - 1] + DeltaY
            Z = 40*Line_Z[step - 1] + DeltaZ

    # returns line as two arrays of X and Y points
    return Line_X, Line_Y, Line_Z


print('tracing...')
# the figure
fig1 = plt.figure(1)
fig1.set_size_inches(10,10, forward = True)
# making 3D plot
ax = fig1.add_subplot(111, projection = '3d')
Contour = Line(1000,200,500)
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
ax.set_zlim([0,25.6])
ax.set_ylim([0,25.6])
ax.set_xlim([0,51.2])

fig2 = plt.figure(2)
fig2.set_size_inches(10,10, forward = True)
ax2 = fig2.add_subplot(111, projection = '3d')
ax2.plot(Contour[0],Contour[1],Contour[2], linestyle = 'none', marker = '.', markersize = .1)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
# looking down Y axis
ax2.view_init(elev = 0, azim = 270)
ax2.set_zlim([0,25.6])
ax2.set_ylim([0,25.6])
ax2.set_xlim([0,51.2])

fig3 = plt.figure(3)
fig3.set_size_inches(10,10, forward = True)
ax3 = fig3.add_subplot(111, projection = '3d')
ax3.plot(Contour[0],Contour[1],Contour[2], linestyle = 'none', marker = '.', markersize = .1)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')
# looking down Z axis
ax3.view_init(elev = 90, azim = 270)
ax3.set_zlim([0,25.6])
ax3.set_ylim([0,25.6])
ax3.set_xlim([0,51.2])

plt.show()

# TODO:
# C implementation
# figure out debye length

