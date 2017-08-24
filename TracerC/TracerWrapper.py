#-----------------------------------------------------------------------------#
# 2D and 3D magnetic field line tracing program

import os
import numpy as np
import matplotlib as mpt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import AutoMinorLocator
from Py3D.sub import load_movie

__all__ = ['TraceField', 'MapSeparator', 'SeparatorLoader', 'SeparatorSlice']

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
            assert(0 <= Start[0] <= B[0].shape[0] - 1)
            assert(0 <= Start[1] <= B[0].shape[1] - 1)
            assert(0 <= Start[2] <= B[0].shape[2] - 1)
            float(passes)
            float(ds)
        except:
            print('invalid arguments, B components must be equal sized arrays')
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
            
        RK4_3D(Start[0], Start[1], Start[2], B[0], B[1], B[2], Xsize,
        Ysize, Zsize, ds, int(Steps))
        
        while True:
            print('')
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                while True:        
                    try:
                        print('\n')
                        Start[0] = abs(float(raw_input('Enter X value of starting position: \n')))
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
                print('Finished.')
                break
            else:
                print('invalid input try again \n')
                continue
            
    else:
        ds2D = .1
        # checks validity of provided args
        try:
            assert(B[0].shape == B[1].shape)
            assert(0 <= Start[0] <= B[0].shape[0] - 1)
            assert(0 <= Start[1] <= B[0].shape[1] - 1)
        except:
            print('invalid arguments, B components must be equal sized arrays')
            print('X, Y must be within the size of B arrays')
            return 0
        
        # calculates max number of differential steps along line
        if B[0].shape[0] >= B[0].shape[1]:
            Steps = 10 * (B[0].shape[0] / ds)
        else:
            Steps = 10 * (B[0].shape[1] / ds)
        
        # sizes of data files are imoprtant to boundary conditions in 
        # the tracing routines
        Xsize = B[0].shape[0]
        Ysize = B[0].shape[1]
            
        # the 2D trace call
        RK4_2D(Start[0], Start[1], B[0], B[1], Xsize, Ysize, ds2D, Steps)
        
        # checks if user would like to trace a different point with differnt
        # X, Y, Z starting position
        while True:    
            print(' ')
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                while True:        
                    try:
                        Start[0] = abs(float(raw_input('\n' + 'Enter X value of starting position: \n')))
                        Start[1] = abs(float(raw_input('Enter Y value of starting position: \n')))
                        break
                    except:
                        print('invalid input try again \n')
                        continue
                plt.close('all')
                TraceField(B, Start)
                break
            elif cont == 'N' or cont == 'n':
                print('Finished.')
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
                     
    print('\n' + 'Configuring B...')                
    Line_X = np.zeros(int(Steps))
    Line_Y = np.zeros(int(Steps))
    
    Bx = np.zeros(Xsize * Ysize)
    By = np.zeros(Xsize * Ysize)

    Bx = B1.reshape(Xsize*Ysize,order='F')
    By = B2.reshape(Xsize*Ysize,order='F')
    
    print('\n' + 'Tracing...')
    Length = func(Line_X, Line_Y, Xinit, Yinit, Bx, By, Xsize, Ysize, ds, int(Steps))
    Line_X = Line_X[:Length]
    Line_Y = Line_Y[:Length]
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(10,10, forward = True)
    ax = fig1.add_subplot(111)
    ax.set_ylim([0, Ysize-1])
    ax.set_xlim([0, Xsize-1])
    ax.set_ylabel('Y', rotation = 0)
    ax.set_xlabel('X')
    ax.yaxis.set_label_position('right')
    ax.yaxis.labelpad = 10
    
    print('\n' + 'Calculating |B|...')
    Bm = np.sqrt(B1**2 + B2**2)

    print('\n' + 'plotting...')
    ax.pcolormesh(Bm.T, vmin = -.1)            
    ax.plot(Line_X,Line_Y, linestyle = 'none', marker = '.', markersize = .01, color='black')
    ax.set_title('$2D$' + ' ' + '$Line$' + ' ' + '$Trace$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$' + '\n' + '$Started$' + ' ' + '$at$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(Xinit) + ', ' + '$Y$' + ' ' + '$=$' + ' ' + str(Yinit), fontsize=20)
    
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
                     
    print('\n' + 'Configuring B...')    
    Line_X = np.zeros(Steps)
    Line_Y = np.zeros(Steps)
    Line_Z = np.zeros(Steps)
    
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    Bx = B1.reshape(Xsize*Ysize*Zsize,order='F')
    By = B2.reshape(Xsize*Ysize*Zsize,order='F')
    Bz = B3.reshape(Xsize*Ysize*Zsize,order='F')
    
    print('\n' + 'Tracing...')
    func(Line_X, Line_Y, Line_Z, Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps)
    print('Done Tracing, preparing projection plots...')    
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(5, 10, forward = True)
    plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
    fig1.suptitle('$Line$' + ' ' + '$Trace$' + ' ' + '$Projections$' + '\n' + '$Started$' + ' ' + '(' + str(Xinit) + ',' + str(Yinit) + ',' + str(Zinit) + ')', fontsize = 20)
    fig1.patch.set_facecolor('lightgrey')
    # making 3D plot
    
    ax = fig1.add_subplot(311)
    ax.set_title('$Down$' + ' ' + '$Z$' + ' ' + '$axis$', fontsize = 14)
    ax.plot(Line_X,Line_Y, linestyle = 'none', marker = '.', markersize = .01, color = 'b')
    ax.set_aspect('equal')
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set_xlabel('X')
    ax.set_ylabel('Y', rotation = 0)
    ax.set_xlim([0,Xsize-1])
    ax.set_ylim([0,Ysize-1])
    ax.yaxis.set_label_position('right')
    ax.yaxis.labelpad = 10
    
    ax2 = fig1.add_subplot(312)
    ax2.set_title('$Down$' + ' ' + '$Y$' + ' ' + '$axis$', fontsize = 14) 
    ax2.plot(Line_X,Line_Z, linestyle = 'none', marker = '.', markersize = .01, color = 'b')
    ax2.set_aspect('equal')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z', rotation = 0)
    ax2.set_xlim([0,Xsize-1])
    ax2.set_ylim([0,Zsize-1])
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.labelpad = 10
    
    ax3 = fig1.add_subplot(313)
    ax3.set_title('$Down$' + ' ' + '$X$' + ' ' + '$axis$', fontsize = 14) 
    ax3.plot(Line_Y,Line_Z, linestyle = 'none', marker = '.', markersize = .01, color = 'b')
    ax3.set_aspect('equal')
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xlabel('Y')
    ax3.set_ylabel('Z', rotation = 0)
    ax3.set_xlim([0,Ysize-1])
    ax3.set_ylim([0,Zsize-1])
    ax3.yaxis.set_label_position('right')
    ax3.yaxis.labelpad = 10
    
    while True:
        ans = raw_input('\n' + 'Create 3D and Puncture Plots? Y or N \n')
        if ((ans == 'Y') or (ans == 'y')):
            try:        
                Bm[1]
                print('\n' + '|B| already Calculated')
            except:
                print('\n' + 'Calculating |B|...')
                Bm0 = np.zeros(4)
                Bm1 = np.zeros(4)
                Bm2 = np.zeros(4)
                Bm3 = np.zeros(4)
                Bm4 = np.zeros(4)
                Bm5 = np.zeros(4)
                
                Bm0 = np.sqrt(B1[0,:,:]**2 + B2[0,:,:]**2 + B3[0,:,:]**2)
                Bm1 = np.sqrt(B1[:,0,:]**2 + B2[:,0,:]**2 + B3[:,0,:]**2)
                Bm2 = np.sqrt(B1[:,:,0]**2 + B2[:,:,0]**2 + B3[:,:,0]**2)
                Bm3 = np.sqrt(B1[int(Xinit),:,:]**2 + B2[int(Xinit),:,:]**2 + B3[int(Xinit),:,:]**2)
                Bm4 = np.sqrt(B1[:,int(Yinit),:]**2 + B2[:,int(Yinit),:]**2 + B3[:,int(Yinit),:]**2)
                Bm5 = np.sqrt(B1[:,:,int(Zinit)]**2 + B2[:,:,int(Zinit)]**2 + B3[:,:,int(Zinit)]**2)
            
            MIN = Bm2.min()
            MAX = Bm2.max()
            
            print('\n'+'Making Colormeshs...')
            Xval = np.linspace(0, Xsize - 1, Xsize)
            Yval = np.linspace(0, Ysize - 1, Ysize)
            Zval = np.linspace(0, Zsize - 1, Zsize)
        
            XZ, YZ = np.meshgrid(Xval, Yval)
            XZ = XZ.T
            YZ = YZ.T
            ZZ = np.zeros((Xsize, Ysize))
            cmp = plt.cm.bwr
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsZ = cmp(norm(Bm2))
            
            XY, ZY = np.meshgrid(Xval, Zval)
            XY = XY.T
            ZY = ZY.T
            YY = np.zeros((Xsize, Zsize))
            cmp = plt.cm.bwr
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsY = cmp(norm(Bm1))
            
            YX, ZX = np.meshgrid(Yval, Zval)
            ZX = ZX.T
            YX = YX.T
            XX = np.zeros((Zsize, Ysize))
            cmp = plt.cm.bwr
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsX = cmp(norm(Bm0))    
            
            fig2 = plt.figure(2)
            fig2.set_size_inches(10, 10, forward = True)
            fig2.patch.set_facecolor('lightgrey')
            fig2.suptitle('$Line$' + ' ' + '$Trace$' + ' ' + '$3D$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$' + '\n' + '$Started$' + ' ' + '(' + str(Xinit) + ',' + str(Yinit) + ',' + str(Zinit) + ')', fontsize = 20, y = .95)
            # making 3D plot
            ax4 = fig2.add_subplot(111, projection = '3d')
            # 500, 180, 250 is inside reconn zone
            
            ax4.plot(Line_Z, Line_X, Line_Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
            ax4.view_init(elev = 5, azim = 5)
            ax4.plot_surface(ZZ, XZ, YZ, facecolors = colorsZ, shade = False, rstride = 10, cstride = 10)
            ax4.plot_surface(ZY, XY, YY, facecolors = colorsY, shade = False, rstride = 10, cstride = 10)
            ax4.plot_surface(ZX, XX, YX, facecolors = colorsX, shade = False, rstride = 10, cstride = 10) 
            ax4.set_zlim([0, Ysize-1])
            ax4.set_ylim([0, Xsize-1])
            ax4.set_xlim([0, Zsize-1])
            ax4.set_xlabel('Z')
            ax4.set_ylabel('X')
            ax4.set_zlabel('Y')   
            
            fig2.tight_layout()
            
            print('\n' + 'Puncturing...')
            #Puncture Plot
            PunctZ = Punct(Line_Z, Zinit, ds, Steps, Line_X, Line_Y)
            Xz = PunctZ[0]
            Yz = PunctZ[1]
            
            PunctY = Punct(Line_Y, Yinit, ds, Steps, Line_X, Line_Z) 
            Xy = PunctY[0]
            Zy = PunctY[1]
            
            PunctX = Punct(Line_X, Xinit, ds, Steps, Line_Z, Line_Y)
            Zx = PunctX[0]
            Yx = PunctX[1]
            
            xx = np.linspace(0, Xsize, num = Xsize)
            yy = np.linspace(0, Ysize, num = Ysize)
            zz = np.linspace(0, Zsize, num = Zsize)
            
            fig3 = plt.figure(3)
            fig3.set_size_inches(5, 10, forward = True)
            fig3.suptitle('$Puncture$' + ' ' + '$Plots$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$', fontsize = 20)
            plt.subplots_adjust(left = .1, bottom = .1, right = .9, top = .9, wspace = None, hspace = .4)
            fig3.patch.set_facecolor('lightgrey')
            
            ax5 = fig3.add_subplot(311)
            ax5.pcolormesh(xx, yy, Bm5.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
            ax5.scatter(Xz, Yz, c = 'k', s = 1)
            ax5.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax5.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax5.set_aspect('equal')
            ax5.set_title('$Slice$' + ' ' + '$Z$' + '$=$' + str(Zinit), fontsize = 14)
            ax5.locator_params(nbins = 6, axis = 'y')
            ax5.locator_params(nbins = 8, axis = 'x')
            ax5.set_ylim([0, Ysize])
            ax5.set_xlim([0, Xsize])
            ax5.set_ylabel('Y', rotation = 0)
            ax5.set_xlabel('X')
            ax5.yaxis.set_label_position('right')
            ax5.yaxis.labelpad = 10
            
            ax6 = fig3.add_subplot(312)
            ax6.pcolormesh(xx, zz, Bm4.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
            ax6.scatter(Xy, Zy, c = 'k', s = 1)
            ax6.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax6.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax6.set_aspect('equal')
            ax6.set_title('$Slice$' + ' ' + '$Y$' + '$=$' + str(Yinit), fontsize = 14)
            ax6.locator_params(nbins = 6, axis = 'y')
            ax6.locator_params(nbins = 8, axis = 'x')
            ax6.set_ylim([0, Zsize])
            ax6.set_xlim([0, Xsize])
            ax6.set_ylabel('Z', rotation = 0)
            ax6.set_xlabel('X')
            ax6.yaxis.set_label_position('right')
            ax6.yaxis.labelpad = 10
            
            ax7 = fig3.add_subplot(313)
            ax7.pcolormesh(yy, zz, Bm3.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
            ax7.scatter(Yx, Zx, c = 'k', s = 1)
            ax7.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax7.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax7.set_aspect('equal')
            ax7.set_title('$Slice$' + ' ' + '$X$' + '$=$' + str(Xinit), fontsize = 14)
            ax7.locator_params(nbins = 6, axis = 'y')
            ax7.locator_params(nbins = 8, axis = 'x')
            ax7.set_ylim([0, Zsize])
            ax7.set_xlim([0, Ysize])
            ax7.set_ylabel('Z', rotation = 0)
            ax7.set_xlabel('Y')
            ax7.yaxis.set_label_position('right')
            ax7.yaxis.labelpad = 10
            break
        elif ((ans == 'N') or (ans == 'n')):
            break
        else:
            print('invalid input try again \n')
        
    print('\n' + 'Plotting...')
    plt.show()
   
   
def Punct(PunctAxisData, Val, ds, Steps, OtherAxisData1, OtherAxisData2):
    func          = _lib.Punct
    func.restype  = ctypes.c_int
    func.argtypes = [ndpointer(ctypes.c_double),
                     ctypes.c_double,
                     ctypes.c_float,
                     ctypes.c_uint, 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double),
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double),]
                     
    Points1 = np.zeros(Steps)
    Points2 = np.zeros(Steps)           
    Total = func(PunctAxisData, Val, ds, int(Steps), OtherAxisData1, OtherAxisData2, Points1, Points2)
    Points1 = Points1[:Total]
    Points2 = Points2[:Total]
    
    return Points1, Points2
#-----------------------------------------------------------------------------#
# SEPARATOR

def MapSeparator(B, Ystart, ds = 2, passes = 1.5):
    try:
        assert(B[0].shape == B[1].shape == B[2].shape)
        assert(float(ds))
        assert(float(passes))
        assert(float(Ystart))
    except:
        print('invalid arguments, B components must be equal sized arrays')
        return 0
    if ((ds != 2) or (passes != 1.5)):
        print('Warning, altering ds or the passes may prevent program from accurately mapping separator.')
    N = 1
    SeparatorX = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SeparatorY = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SeparatorZ = np.zeros((B[0].shape[0]/N)*(B[0].shape[2]/N))
    SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ, Ystart)
    
    print('Saving...')
    np.save('SepX5dx2', SeparatorX)
    np.save('SepY5dx2', SeparatorY)
    np.save('SepZ5dx2', SeparatorZ)


def SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ, Ystart):
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
    Start[1] = Ystart
    Start[2] = 0
    Yval = func(N, Start, Bx, By, Bz, Xsize, Ysize, Zsize, float(ds), int(Steps), SeparatorX, SeparatorY, SeparatorZ)
    SeparatorX[0] = 0
    SeparatorY[0] = Yval
    SeparatorZ[0] = 0

    for i in range(0, Xsize/N):
        print(str(i+1) + '/' + str(Xsize))
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


def SeparatorSlice(SepDataPathX, SepDataPathY, SepDataPathZ, SizeX, SizeY, SizeZ, Bx = None, By = None, Bz = None):
    try:
        assert(Bx.shape == By.shape == Bz.shape)
        assert(int(SizeX))
        assert(int(SizeY))
        assert(int(SizeZ))
    except:
        try:
            assert(Bx == None)
            assert(int(SizeX))
            assert(int(SizeY))
            assert(int(SizeZ))
        except:
            print('invalid arguments')
            return 0
            
    while True:
        print(' ')
        Axis = str(raw_input('Slice X or Z value: \n'))
        if ((Axis == 'X') or (Axis == 'x')):
            while True:
                while True:
                    try:
                        X = int(raw_input('Enter X value for slice: \n'))
                        assert(0 <= X < SizeX)
                        break;
                    except:
                        print('invalid X try again. \n')
                        continue
                
                print('\n' + 'Loading')
                SepY =  np.load(SepDataPathY)
                SepZ =  np.load(SepDataPathZ)
                print('Loaded' + '\n')
                print('Slicing...')
                
                SepYSlice = np.zeros(SizeZ)        
                SepZSlice = np.zeros(SizeZ)             
                
                First = SizeZ*(X)
                SepYSlice = SepY[First:(First + SizeZ - 1)]
                SepZSlice = SepZ[First:(First + SizeZ - 1)]
                
                print('\n' + 'Generating Figures...')
                
                fig = plt.figure(1)
                fig.set_size_inches(30,6, forward = True)
                plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                # making 3D plot
                ax = fig.add_subplot(111)
                # 500, 180, 250 is inside reconn zone
                ax.plot(SepZSlice, SepYSlice, linewidth = .5, color = 'b')
                ax.set_ylabel('Y', rotation = 0)
                ax.set_xlabel('Z')
                ax.yaxis.set_label_position('right')
                ax.yaxis.labelpad = 10
                ax.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(X), fontsize=20)
                
                try:
                    Bz[1]
                    
                    fig2 = plt.figure(2)
                    fig2.set_size_inches(30,6, forward = True)
                    plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                    # making 3D plot
                    ax2 = fig2.add_subplot(111)
                    # 500, 180, 250 is inside reconn zone
                    ax2.plot(SepZSlice, SepYSlice, linewidth = 2, color = 'g', label = 'Separator Sheet')
                    ax2.set_ylabel('Y', rotation = 0)
                    ax2.set_xlabel('Z')
                    ax2.yaxis.set_label_position('right')
                    ax2.yaxis.labelpad = 10
                    ax2.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Puncture$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(X), fontsize=20)
                    
                    Steps = 1000*Bx.shape[0]
                    Upper = RK4_3D_SepSlice(0, SepY[0] + .3, 0, Bx, By, Bz, SizeX, SizeY, SizeZ, .1, Steps)
                    Lower = RK4_3D_SepSlice(0, SepY[0] - .3, 0, Bx, By, Bz, SizeX, SizeY, SizeZ, .1, Steps)
                    UpperPunct = Punct(Upper[0], X, .1, Steps, Upper[2], Upper[1])
                    LowerPunct = Punct(Lower[0], X, .1, Steps, Lower[2], Lower[1])
                    ax2.scatter(UpperPunct[0], UpperPunct[1], s = 3, c = 'r', label = 'Line Started Above')
                    ax2.scatter(LowerPunct[0], LowerPunct[1], s = 3, c = 'b', label = 'Line Started Below')
                    ax2.legend(loc='center right', fontsize = 'large')
                                   
                except:
                    1
                
                print('\n' + 'Plotting...')
                plt.show()
                
                while True:
                    print('\n')
                    ans = raw_input('Slice another point? Y or N \n')
                    if ((ans == 'Y') or (ans == 'y')):
                        plt.close('all')
                        break
                    elif ((ans == 'N') or (ans == 'n')):
                        return
                    else:
                        print('invalid input try again... \n')
                break
            
        elif ((Axis == 'Z') or (Axis == 'z')):
            while True:
                while True:
                    try:
                        Z = int(raw_input('Enter Z value for slice: \n'))
                        assert(0 <= Z < SizeZ)
                        break
                    except:
                        print('invalid Z try again. \n')
                        continue
                
                print('\n' + 'Loading')
                SepX =  np.load(SepDataPathX)
                SepY =  np.load(SepDataPathY)
                print('Loaded' + '\n')
                print('Slicing...')
                
                SepYSlice = np.zeros(SizeX)        
                SepXSlice = np.zeros(SizeX)        
                
                for i in range(0, SizeX):
                    SepYSlice[i] = SepY[SizeZ*i + Z]
                    SepXSlice[i] = SepX[SizeZ*i + Z]
                
                print('\n' + 'Generating Figures...')            
                
                fig = plt.figure(1)
                fig.set_size_inches(30,6, forward = True)
                plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                # making 3D plot
                ax = fig.add_subplot(111)
                # 500, 180, 250 is inside reconn zone
                ax.plot(SepXSlice, SepYSlice, linewidth = .5, color = 'b')
                ax.set_ylabel('Y', rotation = 0)
                ax.set_xlabel('X')
                ax.yaxis.set_label_position('right')
                ax.yaxis.labelpad = 10
                ax.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(Z), fontsize=20)
                
                try:
                    Bz[1]
                    
                    fig2 = plt.figure(2)
                    fig2.set_size_inches(30,6, forward = True)
                    plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                    # making 3D plot
                    ax2 = fig2.add_subplot(111)
                    # 500, 180, 250 is inside reconn zone
                    ax2.plot(SepXSlice, SepYSlice, linewidth = 2, color = 'g', label = 'Separator Sheet')
                    ax2.set_ylabel('Y', rotation = 0)
                    ax2.set_xlabel('X')
                    ax2.yaxis.set_label_position('right')
                    ax2.yaxis.labelpad = 10
                    ax2.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Puncture$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(Z), fontsize=20)
                    
                    Steps = 1000*Bx.shape[0]
                    Upper = RK4_3D_SepSlice(0, SepY[0] + .3, 0, Bx, By, Bz, SizeX, SizeY, SizeZ, .1, Steps)
                    Lower = RK4_3D_SepSlice(0, SepY[0] - .3, 0, Bx, By, Bz, SizeX, SizeY, SizeZ, .1, Steps)
                    UpperPunct = Punct(Upper[2], Z, .1, Steps, Upper[0], Upper[1])
                    LowerPunct = Punct(Lower[2], Z, .1, Steps, Lower[0], Lower[1])
                    ax2.scatter(UpperPunct[0], UpperPunct[1], s = 3, c = 'r', label = 'Line Started Above')
                    ax2.scatter(LowerPunct[0], LowerPunct[1], s = 3, c = 'b', label = 'Line Started Below')
                    ax2.legend(loc='center right', fontsize = 'large')
                    
                except:
                    1
                    
                print('\n' + 'Plotting...')
                plt.show()            
                
                while True:
                    print('\n')
                    ans = raw_input('Slice another point? Y or N \n')
                    if ((ans == 'Y') or (ans == 'y')):
                        plt.close('all')
                        break
                    elif ((ans == 'N') or (ans == 'n')):
                        return
                    else:
                        print('invalid input try again... \n')
                break
                
        else:
            print('invalid option')


def RK4_3D_SepSlice(Xinit, Yinit, Zinit, B1, B2, B3, Xsize, Ysize, Zsize, ds, Steps):
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
    
    print('\n' + 'Tracing...')
    func(Line_X, Line_Y, Line_Z, Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps)
    print('Done Tracing')
    
    return Line_X, Line_Y, Line_Z
    
    
def SeparatorLoader(PathSepX, PathSepY, PathSepZ, Xsize, Ysize, Zsize, B1 = None, B2 = None, B3 = None):
    try:
        assert(B1.shape == B2.shape == B3.shape)
        assert(int(Xsize))
        assert(int(Ysize))
        assert(int(Zsize))
    except:
        try:
            assert(B1 == None)
            assert(int(Xsize))
            assert(int(Ysize))
            assert(int(Zsize))
        except:
            print('invalid arguments')
            return 0
            
    print('Loading Separator...')    
    X = np.load(PathSepX)
    Y = np.load(PathSepY)
    Z = np.load(PathSepZ)
    print('Loaded \n')
    
    print('Generating Figures...')
    fig1 = plt.figure(1)
    fig1.set_size_inches(5, 10, forward = True)
    plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
    fig1.suptitle('$Separator$' + ' ' + '$Projections$', fontsize = 20)
    fig1.patch.set_facecolor('lightgrey')
    # making 3D plot
    
    ax = fig1.add_subplot(311)
    ax.set_title('$Down$' + ' ' + '$Z$' + ' ' + '$axis$', fontsize = 14)
    ax.plot(X,Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
    ax.set_aspect('equal')
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set_xlabel('X')
    ax.set_ylabel('Y', rotation = 0)
    ax.set_xlim([0,Xsize-1])
    ax.set_ylim([0,Ysize-1])
    ax.yaxis.set_label_position('right')
    ax.yaxis.labelpad = 10
    
    ax2 = fig1.add_subplot(312)
    ax2.set_title('$Down$' + ' ' + '$Y$' + ' ' + '$axis$', fontsize = 14) 
    ax2.plot(X,Z, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
    ax2.set_aspect('equal')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z', rotation = 0)
    ax2.set_xlim([0,Xsize-1])
    ax2.set_ylim([0,Zsize-1])
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.labelpad = 10
    
    ax3 = fig1.add_subplot(313)
    ax3.set_title('$Down$' + ' ' + '$X$' + ' ' + '$axis$', fontsize = 14) 
    ax3.plot(Y,Z, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
    ax3.set_aspect('equal')
    ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax3.set_xlabel('Y')
    ax3.set_ylabel('Z', rotation = 0)
    ax3.set_xlim([0,Ysize-1])
    ax3.set_ylim([0,Zsize-1])
    ax3.yaxis.set_label_position('right')
    ax3.yaxis.labelpad = 10
    
    try:
        B3[0,0,0]
        
        print('\n'+'Calculating |B|...')
        Bm0 = np.zeros(4)
        Bm1 = np.zeros(4)
        Bm2 = np.zeros(4)
        Bm3 = np.zeros(4)
        Bm0 = np.sqrt(B1[0,:,:]**2 + B2[0,:,:]**2 + B3[0,:,:]**2)
        Bm1 = np.sqrt(B1[:,0,:]**2 + B2[:,0,:]**2 + B3[:,0,:]**2)
        Bm2 = np.sqrt(B1[:,:,0]**2 + B2[:,:,0]**2 + B3[:,:,0]**2)
        Bm3 = np.sqrt(B1[:,:,int(Zsize/2)]**2 + B2[:,:,int(Zsize/2)]**2 + B3[:,:,int(Zsize/2)]**2)

        MIN = Bm2.min()
        MAX = Bm2.max()
        
        print('\n'+'Making Colormeshs...')
        Xval = np.linspace(0, Xsize - 1, Xsize)
        Yval = np.linspace(0, Ysize - 1, Ysize)
        Zval = np.linspace(0, Zsize - 1, Zsize)
    
        XZ, YZ = np.meshgrid(Xval, Yval)
        XZ = XZ.T
        YZ = YZ.T
        ZZ = np.zeros((Xsize, Ysize))
        cmp = plt.cm.bwr
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsZ = cmp(norm(Bm2))
        
        XY, ZY = np.meshgrid(Xval, Zval)
        XY = XY.T
        ZY = ZY.T
        YY = np.zeros((Xsize, Zsize))
        cmp = plt.cm.bwr
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsY = cmp(norm(Bm1))
        
        YX, ZX = np.meshgrid(Yval, Zval)
        ZX = ZX.T
        YX = YX.T
        XX = np.zeros((Zsize, Ysize))
        cmp = plt.cm.bwr
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsX = cmp(norm(Bm0))    
        
        fig2 = plt.figure(2)
        fig2.set_size_inches(10, 10, forward = True)
        fig2.patch.set_facecolor('lightgrey')
        fig2.suptitle('$Separator$' + ' ' + '$3D$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$', fontsize = 20, y = .95)
        # making 3D plot
        ax4 = fig2.add_subplot(111, projection = '3d')
        # 500, 180, 250 is inside reconn zone
        
        ax4.plot(Z, X, Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
        ax4.view_init(elev = 10, azim = 10)
        ax4.plot_surface(ZZ, XZ, YZ, facecolors = colorsZ, shade = False, rstride = 10, cstride = 10)
        ax4.plot_surface(ZY, XY, YY, facecolors = colorsY, shade = False, rstride = 10, cstride = 10)
        ax4.plot_surface(ZX, XX, YX, facecolors = colorsX, shade = False, rstride = 10, cstride = 10) 
        ax4.set_zlim([0, Ysize-1])
        ax4.set_ylim([0, Xsize-1])
        ax4.set_xlim([0, Zsize-1])
        ax4.set_xlabel('Z')
        ax4.set_ylabel('X')
        ax4.set_zlabel('Y')   
        
        fig2.tight_layout()
        
        fig3 = plt.figure(3)
        fig3.set_size_inches(18, 9, forward = True)
        plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
        fig3.suptitle('$Separator$' + ' ' + '$Projection$', fontsize = 20)
        fig3.patch.set_facecolor('lightgrey')
        # making 3D plot

        xx = np.linspace(0, Xsize, num = Xsize)
        yy = np.linspace(0, Ysize, num = Ysize)
        
        ax5 = fig3.add_subplot(111)
        ax5.set_title('$Over$' + ' ' + '$|$' + '$B$' + '$|$' + ' ' + '$with$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(int(Zsize/2)), fontsize = 14)
        ax5.pcolormesh(xx, yy, Bm3.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)        
        ax5.plot(X,Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
        ax5.set_aspect('equal')
        ax5.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax5.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax5.set_xlabel('X')
        ax5.set_ylabel('Y', rotation = 0)
        ax5.set_xlim([0,Xsize-1])
        ax5.set_ylim([0,Ysize-1])
        ax5.yaxis.set_label_position('right')
        ax5.yaxis.labelpad = 10
    except:
        1
        
    print('\n' + 'Plotting...')
    plt.show()
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


# 2D TESTS

#d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by'], 0)   
#Bx = d['by']
#By = d['bx']
#print('Loaded')
#TraceField([Bx, By], [2000, 2000])
#TraceField([Bx, By], [4596, 4596])

# 3D TESTS
print('Loading')
BX =  np.load('/scratch-fast/asym030/bx.npy')
BY =  np.load('/scratch-fast/asym030/by.npy')
BZ =  np.load('/scratch-fast/asym030/bz.npy')
print('Loaded')
#TraceField([BX, BY, BZ], [1200, 500, 500], .01, 100)
#TraceField([BX, BY, BZ], [500, 200, 500], 20, 100)

#SeparatorSlice('SepX4dx2.npy', 'SepY4dx2.npy', 'SepZ4dx2.npy', 2048, 1024, 1024)
SeparatorSlice('SepX4dx2.npy', 'SepY4dx2.npy', 'SepZ4dx2.npy', 2048, 1024, 1024, BX, BY, BZ)

#SeparatorLoader('SepX4dx2.npy','SepY4dx2.npy','SepZ4dx2.npy', 2048, 1024, 1024)

#SeparatorLoader('SepX4dx2.npy','SepY4dx2.npy','SepZ4dx2.npy', 2048, 1024, 1024, BX, BY, BZ)

#MapSeparator([BX, BY, BZ], 318)
