#-----------------------------------------------------------------------------#
# 2D and 3D magnetic field line tracing program

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from types import StringType
from Py3D.sub import load_movie
from matplotlib.ticker import AutoMinorLocator




def TraceField(Xinit = None, Yinit = None, Zinit = None, Bx = None, By = None, Bz = None, ds = .1, passes = 10):
    if Xinit == None:
        while True:        
            try:
                Xinit = abs(int(raw_input('\n Enter X value of starting position: \n')))
                break
            except:
                print('invalid input try again \n')
                continue
    if Yinit == None:
        while True:        
            try:
                Yinit = abs(int(raw_input('Enter Y value of starting position: \n')))
                break
            except:
                print('invalid input try again \n')
                continue
    if Zinit == None:
        while True:        
            try:
                Zinit = abs(int(raw_input('Enter Z value of starting position: \n')))
                break
            except:
                print('invalid input try again \n')
                continue
    if Bx == None or By == None or Bz == None:
        while True:
            Mov = raw_input('\n Is data in Movie file? Y or N: \n')
            if Mov == 'Y' or Mov == 'y':
                d = load_movie()
                Bx = d['bx']
                By = d['by']
                Bz = d['bz']
                break            
            elif Mov == 'N' or Mov == 'n':
                while True:        
                    PathX = raw_input('\n Specify file path for X component field data: \n' )
                    try:
                        print('loading...')
                        Bx = np.load(PathX)
                        break
                    except:
                        print('no file found try again \n')
                        continue
                while True:    
                    PathY = raw_input('Specify file path for Y component field data: \n' )
                    try:
                        print('loading...')
                        By = np.load(PathY)
                        break
                    except:
                        print('no file found try again \n')
                        continue
                while True:   
                    PathZ = raw_input('Specify file path for Z component field data: \n' )
                    try:
                        print('loading...')
                        Bz = np.load(PathZ)
                        break
                    except:
                        print('no file found try again \n')
                        continue
                break
            else:
                print('invalid input try again \n')
                continue
            
    Dimension = 0
    if len(Bx.shape) == 3:
        Dimension = 3
        if Bx.shape[0] >= Bx.shape[1] and Bx.shape[0] >= Bx.shape[2]:
            Steps = passes * (Bx.shape[0] / ds)
        elif Bx.shape[1] >= Bx.shape[0] and Bx.shape[1] >= Bx.shape[2]:
            Steps = passes * (Bx.shape[1] / ds)
        else:
            Steps = passes * (Bx.shape[2] / ds)
    elif len(Bx.shape) == 2:
        Dimension = 2
        if Bx.shape[0] >= Bx.shape[1]:
            Steps = passes * (Bx.shape[0] / ds)
        else:
            Steps = passes * (Bx.shape[1] / ds)
            
    if Dimension == 3:
        Xsize = Bx.shape[0]
        Ysize = Bx.shape[1]
        Zsize = Bx.shape[2]
        while True:
            FieldLine3D(Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps)
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                continue            
            elif cont == 'N' or cont == 'n':
                break
            else:
                print('invalid input try again \n')
    elif Dimension == 2:
        Xsize = Bx.shape[0]
        Ysize = Bx.shape[1]
        FieldLine2D(Xinit, Yinit, Bx, By, Xsize, Ysize, ds, Steps)
        #Generate Multiline plot?
    else:
        print('invalid data file, aborting...')


def FieldLine2D(Xinit, Yinit, Bx, By, SizeX, SizeY, dx = .1, steps = 100000):

    return 0
    
def FieldLine3D(Xinit, Yinit, Zinit, Bx, By, Bz, SizeX, SizeY, SizeZ, dx = .1, steps = 100000):
    
    print('tracing...')

    X = Xinit
    Y = Yinit
    Z = Zinit
    
    steps = int(steps)
    
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
    Line_X = np.zeros(steps)
    Line_Y = np.zeros(steps)
    Line_Z = np.zeros(steps)

    # loop to step forward line from initial point
    for step in range(0,steps):
        
        # checking if X, Y, or Z has moved outside range (-.5 to size-.5) and if so
        # switches to other side of data grid
        if X < -.5:
            X = X + SizeX
        if Y < -.5:
            Y = Y + SizeY
        if Z < -.5:
            Z = Z + SizeY
        if X > (SizeX - .5):
            X = X - SizeX
        if Y > (SizeY - .5):
            Y = Y - SizeY
        if Z > (SizeZ - .5):
            Z = Z - SizeZ
        
        # just a status update, sometimes it takes a long time
        if (100*float(step)/float(steps))%10 == 0:
            print(str(100*float(step)/float(steps)) + '%')
        
        # adding points (X,Y) to field line divided by 20 for real space
        Line_X[step] = X
        Line_Y[step] = Y
        Line_Z[step] = Z
            
        # finding next point on Line
        Slopes = CalcSlopes(X, Y, Z, Bx, By, Bz, SizeX, SizeY, SizeZ)
        K1x = Slopes[0]
        K1y = Slopes[1]
        K1z = Slopes[2]
        
        Slopes2 = CalcSlopes(X + (dx/2)*K1x , Y + (dx/2)*K1y, Z + (dx/2)*K1z, Bx, By, Bz, SizeX, SizeY, SizeZ)
        K2x = Slopes2[0]
        K2y = Slopes2[1]
        K2z = Slopes2[2]
        
        Slopes3 = CalcSlopes(X + (dx/2)*K2x , Y + (dx/2)*K2y, Z + (dx/2)*K2z, Bx, By, Bz, SizeX, SizeY, SizeZ)
        K3x = Slopes3[0]
        K3y = Slopes3[1]
        K3z = Slopes3[2]
        
        Slopes4 = CalcSlopes(X + dx*K3x , Y + dx*K3y, Z + dx*K3z, Bx, By, Bz, SizeX, SizeY, SizeZ)
        K4x = Slopes4[0]
        K4y = Slopes4[1]
        K4z = Slopes4[2]

        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x)
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y)
        Z = Z + (dx/6)*(K1z + 2*K2z + 2*K3z + K4z)
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(30,8, forward = True)
    # making 3D plot
    ax = fig1.add_subplot(131, projection = '3d')
    ax.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$X$' + ' ' + '$axis$', fontsize=20)
    ax.plot(Line_X,Line_Y,Line_Z, linestyle = 'none', marker = '.', markersize = .1)
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
    ax.set_xlim([SizeX-1,0])
    ax.set_ylim([SizeY-1,0])
    ax.set_zlim([SizeZ-1,0])
    
    ax2 = fig1.add_subplot(132, projection = '3d')
    ax2.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Y$' + ' ' + '$axis$', fontsize=20)
    ax2.plot(Line_X,Line_Y,Line_Z, linestyle = 'none', marker = '.', markersize = .1)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    # looking down Y axis
    ax2.view_init(elev = 0, azim = 270)
    ax2.set_xlim([SizeX-1,0])
    ax2.set_ylim([SizeY-1,0])
    ax2.set_zlim([SizeZ-1,0])
    
    ax3 = fig1.add_subplot(133, projection = '3d')
    ax3.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Z$' + ' ' + '$axis$', fontsize=20)
    ax3.plot(Line_X,Line_Y,Line_Z, linestyle = 'none', marker = '.', markersize = .1)
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    # looking down Z axis
    ax3.view_init(elev = 90, azim = 270)
    ax3.set_xlim([SizeX-1,0])
    ax3.set_ylim([SizeY-1,0])
    ax3.set_zlim([SizeZ-1,0])
    
    fig1.tight_layout()

    while True:
        Punct = raw_input('Generate puncture plots? Y or N \n')
        if Punct == 'Y' or Punct == 'y':
            Puncture(Xinit, Yinit, Zinit, Line_X, Line_Y, Line_Z, Bx, By, Bz, SizeX, SizeY, SizeZ, steps, dx)   
            break
        elif Punct == 'N' or Punct == 'n':
            print('Puncture plots cancelled \n')
            break
        else:
            print('invalid input try again \n')
            
    print('plotting...')

    plt.show()


def CalcSlopes(x, y, z, Bx, By, Bz, SizeX, SizeY, SizeZ):
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
        if i == (SizeX - 1):
            i1 = 0
        if j == (SizeY - 1):
            j1 = 0
        if k == (SizeZ - 1):
            k1 = 0
    elif y > 0 and z > 0:
        i = (SizeX - 1)
        Wx = 1 + x
        Wy = y % 1
        j = int(y)
        Wz = z % 1
        k = int(z)
        i1 = 0
        j1 = j + 1
        k1 = k + 1
        if j == (SizeY - 1):
            j1 = 0
        if k == (SizeZ - 1):
            k1 = 0
    elif x > 0 and z > 0:
        j = (SizeY - 1)
        Wy = 1 + y
        Wx = x % 1
        i = int(x)
        Wz = z % 1
        k = int(z)
        i1 = i + 1
        j1 = 0
        k1 = k + 1
        if i == (SizeX - 1):
            i1 = 0
        if k == (SizeZ - 1):
            k1 = 0
    elif y > 0 and x > 0:
        k = (SizeZ - 1)
        Wz = 1 + z
        Wx = x % 1
        i = int(x)
        Wy = y % 1
        j = int(y)
        i1 = i + 1
        j1 = j + 1
        k1 = 0
        if i == (SizeX - 1):
            i1 = 0
        if j == (SizeY - 1):
            j1 = 0
    elif z > 0:
        i = (SizeX - 1)
        Wx = 1 + x
        j = (SizeY - 1)
        Wy = 1 + y
        Wz = z % 1
        k = int(z)
        i1 = 0
        j1 = 0
        k1 = k + 1
        if k == (SizeZ - 1):
            k1 = 0
    elif y > 0:
        i = (SizeX - 1)
        Wx = 1 + x
        k = (SizeZ - 1)
        Wz = 1 + z
        Wy = y % 1
        j = int(y)
        i1 = 0
        j1 = j + 1
        k1 = 0
        if j == (SizeY - 1):
            j1 = 0
    elif x > 0:
        k = (SizeZ - 1)
        Wz = 1 + z
        j = (SizeY - 1)
        Wy = 1 + y
        Wx = x % 1
        i = int(x)
        i1 = i + 1
        j1 = 0
        k1 = 0
        if i == (SizeX - 1):
            i1 = 0
    else:
        i = (SizeX - 1)
        Wx = 1 + x
        j = (SizeY - 1)
        Wy = 1 + y
        k = (SizeZ - 1)
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

def Puncture(Xinit, Yinit, Zinit, LineX, LineY, LineZ, Bx, By, Bz, SizeX, SizeY, SizeZ, Steps, dx):
    
    print('WARNING, COLORMAPS MAY NEED ADJUSTING')

    print('\n' + 'calculating |B|...')
    Bm = np.sqrt(Bx**2 + By**2 + Bz**2)    
    
    print('puncturing...')
    Xz = []
    Yz = []
    
    Steps = int(Steps)
    for i in range(0,Steps):
        if (100*float(i)/float(3*Steps))%10 == 0:
            print(str(100*float(i)/float(3*Steps)) + '%')
        if (Zinit - (dx/2)) < LineZ[i] < (Zinit + (dx/2)):
            Xz.append(LineX[i])
            Yz.append(LineY[i])
    
    Xz = np.array(Xz)
    Yz = np.array(Yz)
    
    Xy = []
    Zy = []
    
    for i in range(0,Steps):
        if (100*((Steps+float(i))/float(3*Steps)))%10 == 0 and i != 0:
            print(str(100*((Steps+float(i))/float(3*Steps))) + '%')
        if (Yinit - dx/2) < LineY[i] < (Yinit + dx/2):
            Xy.append(LineX[i])
            Zy.append(LineZ[i])
    
    Xy = np.array(Xy)
    Zy = np.array(Zy)
    
    Zx = []
    Yx = []
    
    for i in range(0,Steps):
        if (100*((2*float(Steps)+float(i))/float(3*Steps)))%10 == 0 and i != 0:
            print(str(100*((2*float(Steps)+float(i))/float(3*Steps))) + '%')
        if (Xinit - dx/2) < LineX[i] < (Xinit + dx/2):
            Zx.append(LineZ[i])
            Yx.append(LineY[i])
    
    Zx = np.array(Zx)
    Yx = np.array(Yx)
    
    fig2 = plt.figure(2)
    fig2.set_size_inches(6, 12, forward = True)
    fig2.suptitle('$Puncture$' + ' ' + '$Plots$' + ' ' + '$of$' + ' ' + '$\mid B\mid$', fontsize=20)
    fig2.subplots_adjust(hspace = .4)
    fig2.subplots_adjust(top = .9)
    
    ax1 = fig2.add_subplot(311)
    ax1.pcolormesh(Bm[:,:,Zinit].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
    ax1.scatter(Xz, Yz, c = 'b', s = 3)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.set_aspect('equal')
    ax1.set_title('$Z$' + ' ' + '$=$'+ ' ' + str(Zinit))
    ax1.locator_params(nbins = 6, axis = 'y')
    ax1.locator_params(nbins = 8, axis = 'x')
    ax1.set_ylim([0,SizeY-1])
    ax1.set_xlim([0,SizeX-1])
    
    ax2 = fig2.add_subplot(312)
    ax2.pcolormesh(Bm[:,Yinit,:].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
    ax2.scatter(Xy, Zy, c = 'b', s = 3)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_aspect('equal')
    ax2.set_title('$Y$' + ' ' + '$=$'+ ' ' + str(Yinit))
    ax2.locator_params(nbins = 6, axis = 'y')
    ax2.locator_params(nbins = 8, axis = 'x')
    ax2.set_ylim([0,SizeZ-1])
    ax2.set_xlim([0,SizeX-1])
    
    ax3 = fig2.add_subplot(313)
    ax3.pcolormesh(Bm[Xinit,:,:].T, vmin=1, vmax=2, cmap=plt.cm.bwr)
    ax3.scatter(Yx, Zx, c = 'b', s = 3)
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_aspect('equal')
    ax3.set_title('$X$' + ' ' + '$=$'+ ' ' + str(Xinit))
    ax3.locator_params(nbins = 6, axis = 'y')
    ax3.locator_params(nbins = 8, axis = 'x')
    ax3.set_ylim([0,SizeZ-1])
    ax3.set_xlim([0,SizeY-1])
    
    plt.show()

TraceField()
    
