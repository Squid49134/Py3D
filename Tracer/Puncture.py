
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

# data in
# /share/gander/asym030

# loading data

print('loading data...')
Bx =  np.load('/scratch-fast/asym030/bx.npy')
By =  np.load('/scratch-fast/asym030/by.npy')
Bz =  np.load('/scratch-fast/asym030/bz.npy')
xx = np.linspace(0,51.2, num = 2048)
yy = np.linspace(0,25.6, num = 1024)
zz = np.linspace(0,25.6, num = 1024)
Bm = np.sqrt(Bx**2 + By**2 + Bz**2)

# maximum number of steps along line
MaxSteps = 100000

# differential step towards next point on field line
# dx has little effect on accuracy of tracing, accuracy is most strongly
# affected by the curl of the field line
dx = .1

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

Xinit = 480
Yinit = 240
Zinit = 500

Contour = Line(Xinit,Yinit,Zinit)

Xz = []
Yz = []

print('puncturing...')
for i in range(0,MaxSteps):
    if (Zinit - (dx/2))/40 < Contour[2][i] < (Zinit + (dx/2))/40:
        Xz.append(Contour[0][i])
        Yz.append(Contour[1][i])

Xz = np.array(Xz)
Yz = np.array(Yz)

Xy = []
Zy = []

for i in range(0,MaxSteps):
    if (Yinit - dx/2)/40 < Contour[1][i] < (Yinit + dx/2)/40:
        Xy.append(Contour[0][i])
        Zy.append(Contour[2][i])

Xy = np.array(Xy)
Zy = np.array(Zy)

Zx = []
Yx = []

for i in range(0,MaxSteps):
    if (Xinit - dx/2)/40 < Contour[0][i] < (Xinit + dx/2)/40:
        Zx.append(Contour[2][i])
        Yx.append(Contour[1][i])

Zx = np.array(Zx)
Yx = np.array(Yx)

print('plotting...')
fig1 = plt.figure(1)
fig1.set_size_inches(6, 12, forward = True)
fig1.suptitle('Puncture Plots', fontsize=20)
fig1.subplots_adjust(hspace = .4)
fig1.subplots_adjust(top = .9)

ax1 = fig1.add_subplot(311)
ax1.pcolormesh(xx, yy, Bm[:,:,Zinit].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
ax1.scatter(Xz, Yz, c = 'b', s = 3)
ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
ax1.set_aspect('equal')
ax1.set_title('Z = ' + str(Zinit/40))
ax1.locator_params(nbins = 6, axis = 'y')
ax1.locator_params(nbins = 8, axis = 'x')
ax1.set_ylim([0,25.6])
ax1.set_xlim([0,51.2])

ax2 = fig1.add_subplot(312)
ax2.pcolormesh(xx, zz, Bm[:,Yinit,:].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
ax2.scatter(Xy, Zy, c = 'b', s = 3)
ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
ax2.set_aspect('equal')
ax2.set_title('Y = ' + str(Yinit/40))
ax2.locator_params(nbins = 6, axis = 'y')
ax2.locator_params(nbins = 8, axis = 'x')
ax2.set_ylim([0,25.6])
ax2.set_xlim([0,51.2])

ax3 = fig1.add_subplot(313)
ax3.pcolormesh(yy, zz, Bm[Xinit,:,:].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
ax3.scatter(Yx, Zx, c = 'b', s = 3)
ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
ax3.set_aspect('equal')
ax3.set_title('X = ' + str(Xinit/40))
ax3.locator_params(nbins = 6, axis = 'y')
ax3.locator_params(nbins = 8, axis = 'x')
ax3.set_ylim([0,25.6])
ax3.set_xlim([0,25.6])

plt.show()

