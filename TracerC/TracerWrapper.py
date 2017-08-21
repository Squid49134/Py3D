#-----------------------------------------------------------------------------#
# 2D and 3D magnetic field line tracing program

import numpy as np
import matplotlib as mpt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import AutoMinorLocator
#from Py3D.sub import load_movie
import os

# now for linking C to python
from numpy.ctypeslib import ndpointer
import ctypes

# This is called here so realpath is where the .so file is
_pathlibdir   = os.path.dirname(os.path.realpath(__file__))

# linking to the C library
pathlib   = os.path.join(_pathlibdir,'TracerFunctions.so')
_lib = ctypes.cdll.LoadLibrary(pathlib)

#-----------------------------------------------------------------------------#

# TRACER
def TraceField(B, Start, ds = .1, passes = 10):
    
    try:
        if (len(B) == 2 or len(B) == 3):
            0
        else:
            print('invalid B')
            return 0
    except:
        print('invalid B')
        return 0
        
    # 3D
    if len(B) == 3:
        # checks validity of provided args
        try:
            assert(B[0].shape == B[1].shape == B[2].shape)
            assert(0 < Start[0] < B[0].shape[0] - 1)
            assert(0 < Start[1] < B[0].shape[1] - 1)
            assert(0 < Start[2] < B[0].shape[2] - 1)
        except:
            print('invalid arguments, B components must be equal sized arrays,')
            print('X, Y, Z must be within the size of B arrays')
            return 0
        
        # calculates max number of differential steps along line
        if B[0].shape[0] >= B[0].shape[1] and B[0].shape[0] >= B[0].shape[2]:
            Steps = passes * (B[0].shape[0] / ds)
        elif B[0].shape[1] >= B[0].shape[0] and B[0].shape[1] >= B[0].shape[2]:
            Steps = passes * (B[0].shape[1] / ds)
        else:
            Steps = passes * (B[0].shape[2] / ds)
        
        # sizes of data files are imoprtant to boundary conditions in 
        # the tracing routines
        Xsize = B[0].shape[0]
        Ysize = B[0].shape[1]
        Zsize = B[0].shape[2]
            
        # the 3D trace call
        RK4_3D(Start[0], Start[1], Start[2], B[0], B[1], B[2], Xsize,
        Ysize, Zsize, ds, int(Steps))
        # checks if user would like to trace a different point with differnt
        # X, Y, Z starting position
        while True:    
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                while True:        
                    try:
                        Start[0] = abs(float(raw_input('\n' + 'Enter X value of starting position: \n')))
                        Start[1] = abs(float(raw_input('Enter Y value of starting position: \n')))
                        Start[2] = abs(float(raw_input('Enter Z value of starting position: \n')))
                        ds =  abs(float(raw_input('Enter ds value: \n')))
                        passes = abs(float(raw_input('Enter number of passes: \n')))
                        break
                    except:
                        print('invalid input try again \n')
                        continue
                plt.close('all')
                TraceField(B, Start, ds, passes)
                break
            elif cont == 'N' or cont == 'n':
                print('finished')
                break
            else:
                print('invalid input try again \n')
                continue
            
    else:
        # checks validity of provided args
        try:
            assert(B[0].shape == B[1].shape)
            assert(0 < Start[0] < B[0].shape[0] - 1)
            assert(0 < Start[1] < B[0].shape[1] - 1)
        except:
            print('invalid arguments, B components must be equal sized arrays,')
            print('X, Y must be within the size of B arrays')
            return 0
        
        # calculates max number of differential steps along line
        if B[0].shape[0] >= B[0].shape[1]:
            Steps = passes * (B[0].shape[0] / ds)
        else:
            Steps = passes * (B[0].shape[1] / ds)
        
        # sizes of data files are imoprtant to boundary conditions in 
        # the tracing routines
        Xsize = B[0].shape[0]
        Ysize = B[0].shape[1]
            
        # the 2D trace call
        RK4_2D(Start[0], Start[1], B[0], B[1], Xsize, Ysize, ds, Steps)
        
        # checks if user would like to trace a different point with differnt
        # X, Y, Z starting position
        while True:    
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                while True:        
                    try:
                        Start[0] = abs(float(raw_input('\n' + 'Enter X value of starting position: \n')))
                        Start[1] = abs(float(raw_input('Enter Y value of starting position: \n')))
                        ds =  abs(float(raw_input('Enter ds value: \n')))
                        passes = abs(float(raw_input('Enter number of passes: \n')))
                        break
                    except:
                        print('invalid input try again \n')
                        continue
                plt.close('all')
                TraceField(B, Start, ds, passes)
                break
            elif cont == 'N' or cont == 'n':
                print('finished')
                break
            else:
                print('invalid input try again \n')
                continue
            
def RK4_2D(Xinit, Yinit, B1, B2, Xsize, Ysize, ds, Steps):
    func          = _lib.FieldLine2D
    func.restype  = int
    func.argtypes = [ndpointer(ctypes.c_double),  
                     ndpointer(ctypes.c_double), 
                     ctypes.c_double,
                     ctypes.c_double, 
                     ndpointer(ctypes.c_float), 
                     ndpointer(ctypes.c_float), 
                     ctypes.c_uint, 
                     ctypes.c_uint, 
                     ctypes.c_float, 
                     ctypes.c_uint]
                     
    Line_X = np.zeros(int(Steps))
    Line_Y = np.zeros(int(Steps))
    
    Bx = np.zeros(Xsize * Ysize)
    By = np.zeros(Xsize * Ysize)

    Bx = B1.reshape(Xsize*Ysize,order='F')
    By = B2.reshape(Xsize*Ysize,order='F')
    
    print('Tracing...')
    Length = func(Line_X, Line_Y, Xinit, Yinit, Bx, By, Xsize, Ysize, ds, int(Steps))
    Line_X = Line_X[:Length]
    Line_Y = Line_Y[:Length]
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(10,10, forward = True)
    ax = fig1.add_subplot(111)
    ax.set_ylim([0, Ysize-1])
    ax.set_xlim([0, Xsize-1])
    
    print('\n' + 'Calculating |B|...')
    Bm = np.sqrt(B1**2 + B2**2)

    print('\n' + 'plotting...')
    ax.pcolormesh(Bm.T)            
    ax.plot(Line_X,Line_Y, linestyle = 'none', marker = '.', markersize = .1, color='black')
    ax.set_title('$Started$' + ' ' + '$at$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(Xinit) + ', ' + '$Y$' + ' ' + '$=$' + ' ' + str(Yinit), fontsize=20)
    
    plt.show()
            
def RK4_3D(Xinit, Yinit, Zinit, B1, B2, B3, Xsize, Ysize, Zsize, ds, Steps):
    func          = _lib.FieldLine3D
    func.restype  = None
    func.argtypes = [ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ctypes.c_double, 
                     ctypes.c_double, 
                     ctypes.c_double, 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ctypes.c_uint, 
                     ctypes.c_uint, 
                     ctypes.c_uint, 
                     ctypes.c_float, 
                     ctypes.c_uint]
                     
    Line_X = np.zeros(Steps)
    Line_Y = np.zeros(Steps)
    Line_Z = np.zeros(Steps)
    
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    Bx = B1.reshape(Xsize*Ysize*Zsize,order='F')
    By = B2.reshape(Xsize*Ysize*Zsize,order='F')
    Bz = B3.reshape(Xsize*Ysize*Zsize,order='F')
    
    print('Tracing...')
    func(Line_X, Line_Y, Line_Z, Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps)
    print('Plotting!')
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(30,8, forward = True)
    # making 3D plot
    ax = fig1.add_subplot(131)
    ax.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$X$' + ' ' + '$axis$', fontsize=20) 
    ax.plot(Line_Y,Line_Z, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax.set_xlabel('Y')
    ax.set_ylabel('Z')
    ax.set_xlim([0,Ysize-1])
    ax.set_ylim([0,Zsize-1])
    
    ax2 = fig1.add_subplot(132)
    ax2.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Y$' + ' ' + '$axis$', fontsize=20) 
    ax2.plot(Line_X,Line_Z, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_xlim([0,Xsize-1])
    ax2.set_ylim([0,Zsize-1])
    
    ax3 = fig1.add_subplot(133)
    ax3.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Z$' + ' ' + '$axis$', fontsize=20)
    ax3.plot(Line_X,Line_Y, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_xlim([0,Xsize-1])
    ax3.set_ylim([0,Ysize-1])
    
    fig1.tight_layout()
    
    fig2 = plt.figure(2)
    fig2.set_size_inches(8,8, forward = True)
    # making 3D plot
    ax4 = fig2.add_subplot(111, projection = '3d')
    # 500, 180, 250 is inside reconn zone
    
    ax4.plot(Line_Z, Line_X, Line_Y, linestyle = 'none', marker = '.', markersize = .1)
    ax4.view_init(elev = 5, azim = 5)
    ax4.set_zlim([0, Ysize-1])
    ax4.set_ylim([0, Xsize-1])
    ax4.set_xlim([0, Zsize-1])
    ax4.set_xlabel('Z')
    ax4.set_ylabel('X')
    ax4.set_zlabel('Y')    
    
    plt.show()
    
    
#-----------------------------------------------------------------------------#
# SEPARATOR

def MapSeparator(B, ds = 10, passes = 2):
    N = 1
    SeparatorX = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SeparatorY = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SeparatorZ = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ)
    print('saving...')
    
    np.save('SepX5dx2', SeparatorX)
    np.save('SepY5dx2', SeparatorY)
    np.save('SepZ5dx2', SeparatorZ)
    
    print('plotting...')
    
    fig1 = plt.figure(1)
    fig1.set_size_inches(30,8, forward = True)
    # making 3D plot
    ax = fig1.add_subplot(131)
    ax.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$X$' + ' ' + '$axis$', fontsize=20) 
    ax.plot(SeparatorY,SeparatorZ, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax.set_xlabel('Y')
    ax.set_ylabel('Z')
    ax.set_xlim([0,B[0].shape[1]-1])
    ax.set_ylim([0,B[0].shape[2]-1])
    
    ax2 = fig1.add_subplot(132)
    ax2.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Y$' + ' ' + '$axis$', fontsize=20) 
    ax2.plot(SeparatorX,SeparatorZ, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_xlim([0,B[0].shape[0]-1])
    ax2.set_ylim([0,B[0].shape[2]-1])
    
    ax3 = fig1.add_subplot(133)
    ax3.set_title('$Projection$' + ' ' + '$down$' + ' ' + '$Z$' + ' ' + '$axis$', fontsize=20)
    ax3.plot(SeparatorX,SeparatorY, linestyle = 'none', marker = '.', markersize = .1, color = 'k')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_xlim([0,B[0].shape[0]-1])
    ax3.set_ylim([0,B[0].shape[1]-1])
    
    fig1.tight_layout()
    
    # the figure
    fig2 = plt.figure(2)
    fig2.set_size_inches(8,8, forward = True)
    # making 3D plot
    ax4 = fig2.add_subplot(111, projection = '3d')
    # 500, 180, 250 is inside reconn zone
    
    ax4.plot(SeparatorZ, SeparatorX, SeparatorY, linestyle = 'none', marker = '.', markersize = .1)
    ax4.view_init(elev = 5, azim = 5)
    ax4.set_zlim([0, B[0].shape[1]-1])
    ax4.set_ylim([0, B[0].shape[0]-1])
    ax4.set_xlim([0, B[0].shape[2]-1])
    ax4.set_xlabel('Z')
    ax4.set_ylabel('X')
    ax4.set_zlabel('Y')
    
    plt.show()


def SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ):
    func          = _lib.SepPoints
    func.restype  = ctypes.c_double
    func.argtypes = [ctypes.c_uint, 
                     ndpointer(ctypes.c_double),
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ctypes.c_uint, 
                     ctypes.c_uint, 
                     ctypes.c_uint, 
                     ctypes.c_float, 
                     ctypes.c_uint,
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double),]
     
    Steps = passes * (B[0].shape[0] / ds)
    Xsize = B[0].shape[0]
    Ysize = B[0].shape[1]
    Zsize = B[0].shape[2]
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    Bx = B[0].reshape(Xsize*Ysize*Zsize,order='F')
    By = B[1].reshape(Xsize*Ysize*Zsize,order='F')
    Bz = B[2].reshape(Xsize*Ysize*Zsize,order='F')
    
    Start = np.zeros(3)
    Start[0] = 0
    Start[1] = 318
    Start[2] = 0
    Yval = func(N, Start, Bx, By, Bz, Xsize, Ysize, Zsize, float(ds), int(Steps), SeparatorX, SeparatorY, SeparatorZ)
    SeparatorX[0] = 0
    SeparatorY[0] = Yval
    SeparatorZ[0] = 0

    for i in range(0, Xsize/N):
        print(i)
        for j in range(0, Zsize/N):
            if ((i == 0) and (j == 0)):
                continue
            Start[0] = (N * i)
            Start[1] = SeparatorY[(j + (Zsize/N)*(i))-1]
            Start[2] = (N * j)
            Yval = func(N, Start, Bx, By, Bz, Xsize, Ysize, Zsize, float(ds), int(Steps), SeparatorX, SeparatorY, SeparatorZ)
            SeparatorX[j + (Zsize/N)*(i)] = N * i
            SeparatorY[j + (Zsize/N)*(i)] = Yval
            SeparatorZ[j + (Zsize/N)*(i)] = N * j


def SeparatorSlice():
    print('Slice X or Z value:')
    Axis = str(raw_input())
    
    if ((Axis == 'X') or (Axis == 'x')):
        print('')
        print('Enter Path To Separator Y Value Data:')
        PathY = raw_input()
        print('Enter Path To Separator Z Value Data:')
        PathZ = raw_input()
        while True:
            print('Enter X Value For Slice:' )
            SliceX = int(raw_input())
            print('Loading')
            SepY =  np.load(PathY)
            SepZ =  np.load(PathZ)
            print('Loaded')
            
            X = SliceX
            
            First = 1024*(X)
            SepYSlice = SepY[First:First+1023]
            SepZSlice = SepZ[First:First+1023]
            
            fig = plt.figure(1)
            fig.set_size_inches(30,8, forward = True)
            # making 3D plot
            ax = fig.add_subplot(111)
            # 500, 180, 250 is inside reconn zone
            print('plotting...')
            ax.plot(SepZSlice, SepYSlice)
            
            fig.tight_layout()
            
            print('Slice another point? Y/N')
            ans = raw_input()
            if ((ans == 'Y') or (ans == 'y')):
                1
            else:
                break            
            
            plt.show()
        
    elif ((Axis == 'Z') or (Axis == 'z')):
        print('')
        print('Enter Path To Separator X Value Data:')
        PathX = raw_input()
        print('Enter Path To Separator Y Value Data:')
        PathY = raw_input()
        while True:
            print('Enter Z Value For Slice:' )
            SliceZ = int(raw_input())
            print('Loading')
            SepX =  np.load(PathX)
            SepY =  np.load(PathY)
            print('Loaded')
            
            Z = SliceZ
            
            SepYSlice = np.zeros(2048)        
            SepXSlice = np.zeros(2048)        
            
            for i in range(0, 2048):
                SepYSlice[i] = SepY[1024*i + Z]
                SepXSlice[i] = SepX[1024*i + Z]
            
            fig = plt.figure(1)
            fig.set_size_inches(30,8, forward = True)
            # making 3D plot
            ax = fig.add_subplot(111)
            # 500, 180, 250 is inside reconn zone
            print('plotting...')
            ax.plot(SepXSlice, SepYSlice)
            
            fig.tight_layout()
            
            print('Slice another point? Y/N')
            ans = raw_input()
            if ((ans == 'Y') or (ans == 'y')):
                1
            else:
                break            
            
            plt.show()
            
    else:
        print('invalid option')

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

print('loading')
B1 =  np.load('/scratch-fast/asym030/bx.npy')
B2 =  np.load('/scratch-fast/asym030/by.npy')
B3 =  np.load('/scratch-fast/asym030/bz.npy')
print('loaded')
#TraceField([B1, B2, B3], [1439, 256, 313], 2, 1.5)

#d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by'], 0)   
#Bx = d['by']
#By = d['bx']
#print('loaded')
#TraceField([Bx, By], [2000, 2000], .1, 5)

MapSeparator([B1, B2, B3], 2, 1.5)

#SeparatorSlice()
