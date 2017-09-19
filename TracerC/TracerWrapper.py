# FIELD TRACER

# SEE DOCUMENTATION.TXT FOR INFO ON METHOD ARGUMENTS AND FUNCTIONS

import numpy as np
import matplotlib as mpt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import AutoMinorLocator
from Py3D.sub import load_movie
import ctypes
import os
from numpy.ctypeslib import ndpointer

__all__ = ['TraceField', 'MapSeparator', 'SeparatorLoader', 'SeparatorSlice']

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# LINKING THE SHARED C LIBRARY

# The path to the current file
_pathlibdir   = os.path.dirname(os.path.realpath(__file__))

# Linking the shared object file
pathlib   = os.path.join(_pathlibdir,'TracerFunctions.so')
_lib = ctypes.cdll.LoadLibrary(pathlib)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# FIELD TRACING METHODS

# Master 2D and 3D tracing method
def TraceField(SIMds, B, Start, ds = None, passes = None, Saves = None):
    # Checking length of B vector, for 2D, 3D, or 3D with E interp
    try:
        if ((len(B) == 2) or (len(B) == 3) or (len(B) == 6)):
            0
        else:
            print('invalid B')
            return 0
    except:
        print('invalid B')
        return 0
        
    #3D
    if ((len(B) == 3) or (len(B) == 6)):
        # Checking the validitiy of provided arguments
        try:
            try:
                Saves[0] == 'String'
                Saves[1] == 'String'
                Saves[2] == 'String'
            except:
                if (Saves == None):
                    1
                else:
                    print('invalid arguments')
                    return 0
            if (len(B) == 6):
                B[3][0,0,0] * B[4][0,0,0] * B[5][0,0,0]
                assert(B[3].shape == B[4].shape == B[5].shape == B[0].shape)
            B[0][0,0,0] * B[1][0,0,0] * B[2][0,0,0]
            assert(B[0].shape == B[1].shape == B[2].shape)
            assert(0 <= Start[0]/SIMds <= B[0].shape[0] - 1)
            assert(0 <= Start[1]/SIMds <= B[0].shape[1] - 1)
            assert(0 <= Start[2]/SIMds <= B[0].shape[2] - 1)
            float(passes)
            float(ds)
        except:
            print('invalid arguments')
            return 0
            
        # Converting to grid space
        ds = ds / SIMds
        Start[0] = Start[0]/SIMds
        Start[1] = Start[1]/SIMds
        Start[2] = Start[2]/SIMds
            
        # Determining the number of steps
        if B[0].shape[0] >= B[0].shape[1] and B[0].shape[0] >= B[0].shape[2]:
            Steps = passes * (B[0].shape[0] / ds)
        elif B[0].shape[1] >= B[0].shape[0] and B[0].shape[1] >= B[0].shape[2]:
            Steps = passes * (B[0].shape[1] / ds)
        else:
            Steps = passes * (B[0].shape[2] / ds)
        
        # Setting simulation size
        Xsize = B[0].shape[0]
        Ysize = B[0].shape[1]
        Zsize = B[0].shape[2]
        
        # The 3D tracing wrapper call either with or without save paths
        try:
            Saves[0]
            FieldLine_3D(SIMds, Start[0], Start[1], Start[2], B, Xsize, Ysize, Zsize, ds, int(Steps), Saves)
        except:
            FieldLine_3D(SIMds, Start[0], Start[1], Start[2], B, Xsize, Ysize, Zsize, ds, int(Steps))
        
        # Checking if the user wants to trace another point, and if so 
        # takes new input arguments
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
                        assert(0 <= Start[0]/SIMds <= B[0].shape[0] - 1)
                        assert(0 <= Start[1]/SIMds <= B[0].shape[1] - 1)
                        assert(0 <= Start[2]/SIMds <= B[0].shape[2] - 1)
                        break
                    except:
                        print('invalid input try again \n')
                        continue
                plt.close('all')
                TraceField(SIMds, B, Start, ds, passes)
                break
            elif cont == 'N' or cont == 'n':
                print('Finished.')
                break
            else:
                print('invalid input try again \n')
                continue
    
    # 2D
    else:
        # ds is set in 2D for the purpose of checking for line completion
        # ds IS IN GRID SPACE
        ds2D = .1
        # Checking validity of arguments
        try:
            B[0][0,0] * B[1][0,0]
            assert(B[0].shape == B[1].shape)
            assert(0 <= Start[0]/SIMds <= B[0].shape[0] - 1)
            assert(0 <= Start[1]/SIMds <= B[0].shape[1] - 1)
        except:
            print('invalid arguments')
            return 0
            
        Start[0] = Start[0]/SIMds
        Start[1] = Start[1]/SIMds
        
        # The number of steps is also set in 2D for checking line completion
        if B[0].shape[0] >= B[0].shape[1]:
            Steps = 10 * (B[0].shape[0] / ds2D)
        else:
            Steps = 10 * (B[0].shape[1] / ds2D)
        
        # Setting simulation size
        Xsize = B[0].shape[0]
        Ysize = B[0].shape[1]
        
        # The 2D tracing wrapper call
        FieldLine_2D(SIMds, Start[0], Start[1], B[0], B[1], Xsize, Ysize, ds2D, int(Steps))
        
        # Checks if the user would like to trace another point, and if so takes
        # new input arguments
        while True:    
            print(' ')
            cont = raw_input('Trace another point? Y or N \n')
            if cont == 'Y' or cont == 'y':
                while True:        
                    try:
                        Start[0] = abs(float(raw_input('\n' + 'Enter X value of starting position: \n')))
                        Start[1] = abs(float(raw_input('Enter Y value of starting position: \n')))
                        assert(0 <= Start[0]/SIMds <= B[0].shape[0] - 1)
                        assert(0 <= Start[1]/SIMds <= B[0].shape[1] - 1)
                        break
                    except:
                        print('invalid input try again \n')
                        continue
                plt.close('all')
                TraceField(SIMds, B, Start)
                break
            elif cont == 'N' or cont == 'n':
                print('Finished.')
                break
            else:
                print('invalid input try again \n')
                continue
            
# 3D tracing wrapper and plotter
def FieldLine_3D(SIMds, Xinit, Yinit, Zinit, B, Xsize, Ysize, Zsize, ds, Steps, Saves = None):
    # Identifying the C tracing function
    func          = _lib.RK4_3D
    func.restype  = int
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
                     ctypes.c_uint,
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double),
                     ndpointer(ctypes.c_double)]
                     
    print('\n' + 'Configuring B...')  
    # Arrays for the traced line data points
    Line_X = np.zeros(Steps)
    Line_Y = np.zeros(Steps)
    Line_Z = np.zeros(Steps)
    
    # Arrays for restructuring B
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    # Restructuring Bx, By, and Bz into 1D arrays so they can be passed as
    # double pointers to the C function
    Bx = B[0].reshape(Xsize * Ysize * Zsize, order='F')
    By = B[1].reshape(Xsize * Ysize * Zsize, order='F')
    Bz = B[2].reshape(Xsize * Ysize * Zsize, order='F')
    
    # If E is not included Ex, Ey, and Ez are just [0,0,0] so they can be
    # distinguished form Ex, Ey, and Ez data sets when passed to the C func.
    Ex = np.zeros(3)
    Ey = np.zeros(3)
    Ez = np.zeros(3)
    # An array for the E interpolated data
    EI = np.zeros(3)
    
    # If E is included, reshape E like B
    if (len(B) == 6):
        Ex = np.zeros(Xsize * Ysize * Zsize)
        Ey = np.zeros(Xsize * Ysize * Zsize)
        Ez = np.zeros(Xsize * Ysize * Zsize)
        EI = np.zeros(Steps)

        Ex = B[3].reshape(Xsize * Ysize * Zsize, order='F')
        Ey = B[4].reshape(Xsize * Ysize * Zsize, order='F')
        Ez = B[5].reshape(Xsize * Ysize * Zsize, order='F')
    
    print('\n' + 'Tracing...')
    # Callling the C 3D tracing method which returns the total of E
    # interpolated along the trace path, or 0 if E is not passed
    interp = func(Line_X, Line_Y, Line_Z, Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps, Ex, Ey, Ez, EI)
    print('Done Tracing.')
    
    # Saving if appropriate
    try:
        Saves[0]
        print('Saving...')
        np.save(Saves[0], Line_X)
        np.save(Saves[1], Line_Y)
        np.save(Saves[2], Line_Z)
        print('Saved.')
    except:
        1
        
    # Plotting the basic line trace projections which can be done quickly 
    # since no colormeshes are involved
    print('\n' + 'Preparing Projection Plots...')    
    
    # The figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(5, 10, forward = True)
    plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
    fig1.suptitle('$Line$' + ' ' + '$Trace$' + ' ' + '$Projections$' + '\n' + '$Started$' + ' ' + '(' + str(Xinit*SIMds) + ',' + str(Yinit*SIMds) + ',' + str(Zinit*SIMds) + ')' + '\n' + '$ds$' + ' ' + '$=$' + ' ' + str(ds*SIMds) + '$,$' + ' ' + ' ' + '$Passes$' + ' ' + '$=$' + ' ' + str(float(Steps)*ds/B[0].shape[0]), fontsize = 18, y = .99)
    fig1.patch.set_facecolor('lightgrey')
    
    # The plots
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
    
    # Forcing matplot to generate the axes tick labels
    fig1.canvas.draw()
    
    # Returns an list of tick label objects
    LabelList = ax.get_xticklabels()
    # Adjusting the tick labels
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    # Resetting the tick labels
    ax.set_xticklabels(LabelList)
    
    LabelList = ax2.get_xticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax2.set_xticklabels(LabelList)
    
    LabelList = ax3.get_xticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax3.set_xticklabels(LabelList)
    
    LabelList = ax.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax.set_yticklabels(LabelList)
    
    LabelList = ax2.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax2.set_yticklabels(LabelList)
    
    LabelList = ax3.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax3.set_yticklabels(LabelList)
    
    # If E is included plot Einterp per step
    if (len(B) == 6):
        fig4 = plt.figure(4)
        Xvals = np.linspace(0, Steps - 1, Steps)
        fig4.set_size_inches(30,7, forward = True)
        fig4.patch.set_facecolor('lightgrey')
        plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
        
        ax4 = fig4.add_subplot(111)
        ax4.plot(Xvals, EI, linewidth = .5, color = 'g')
        ax4.set_xlabel('Steps')
        ax4.yaxis.set_label_position('right')
        ax4.yaxis.labelpad = 10
        ax4.set_title('$E$' + ' ' + '$Interpolated$' + ' ' + '$Per$' + ' ' + '$Step$' + ' ' + '$For$' + ' ' + str(float(Steps)*ds / B[0].shape[0]) + ' ' + '$Passes$' + '\n' + '$ds$' + ' ' + '$=$' + ' ' + str(ds*SIMds) + '$,$' + ' ' + ' ' + '$Einterp$' + ' ' + '$Total$' + ' ' + '$=$' + ' ' + str(interp), fontsize=18)

    # Checking if the user would like to take extra time creating 3D and
    # puncture plots which require developing colormeshs
    while True:
        ans = raw_input('\n' + 'Create 3D and Puncture Plots? Y or N \n')
        if ((ans == 'Y') or (ans == 'y')):
            print('\n' + 'Calculating |B|...')
            
            # Only calculates |B| for the 6 planes which are plotted to
            # avoid excessive time required to find |B| for the entire sim.
            Bm0 = np.sqrt(B[0][0,:,:]**2 + B[1][0,:,:]**2 + B[2][0,:,:]**2)
            Bm1 = np.sqrt(B[0][:,0,:]**2 + B[1][:,0,:]**2 + B[2][:,0,:]**2)
            Bm2 = np.sqrt(B[0][:,:,0]**2 + B[1][:,:,0]**2 + B[2][:,:,0]**2)
            Bm3 = np.sqrt(B[0][int(Xinit),:,:]**2 + B[1][int(Xinit),:,:]**2 + B[2][int(Xinit),:,:]**2)
            Bm4 = np.sqrt(B[0][:,int(Yinit),:]**2 + B[1][:,int(Yinit),:]**2 + B[2][:,int(Yinit),:]**2)
            Bm5 = np.sqrt(B[0][:,:,int(Zinit)]**2 + B[1][:,:,int(Zinit)]**2 + B[2][:,:,int(Zinit)]**2)
        
            # Min and Max values for normalizing colormaps
            MIN = Bm2.min()
            MAX = Bm2.max()
            
            print('\n'+'Making Colormeshs...')
            # Arrays for making meshgrids
            Xval = np.linspace(0, Xsize - 1, Xsize)
            Yval = np.linspace(0, Ysize - 1, Ysize)
            Zval = np.linspace(0, Zsize - 1, Zsize)
        
            # Meshgrids to be plotted in 3D after applying the |B| colormap
            # allows plotting 2D colormesh in 3D
            XZ, YZ = np.meshgrid(Xval, Yval)
            # Field lines actually follow the transpose of the data
            XZ = XZ.T
            YZ = YZ.T
            # This meshgrid will be on the plane Z = 0
            ZZ = np.zeros((Xsize, Ysize))
            cmp = plt.cm.bwr
            # Normalizing the colomap
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsZ = cmp(norm(Bm2))
            
            XY, ZY = np.meshgrid(Xval, Zval)
            XY = XY.T
            ZY = ZY.T
            # This meshgrid will be on the plane Y = 0
            YY = np.zeros((Xsize, Zsize))
            cmp = plt.cm.bwr
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsY = cmp(norm(Bm1))
            
            YX, ZX = np.meshgrid(Yval, Zval)
            ZX = ZX.T
            YX = YX.T
            # This meshgrid will be on the plane X = 0
            XX = np.zeros((Zsize, Ysize))
            cmp = plt.cm.bwr
            norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
            colorsX = cmp(norm(Bm0))    
            
            print('\n' + 'Puncturing...')
            # Calling the puncture plot generators which return 2 arrays
            # containing coordinates for the puncture points, this one
            # punctures the Z plane at Z = Zinit
            PunctZ = Punct(Line_Z, Zinit, ds, Steps, Line_X, Line_Y)
            Xz = PunctZ[0]
            Yz = PunctZ[1]
            
            # The Y plane puncture
            PunctY = Punct(Line_Y, Yinit, ds, Steps, Line_X, Line_Z) 
            Xy = PunctY[0]
            Zy = PunctY[1]
            
            # The X plane puncture
            PunctX = Punct(Line_X, Xinit, ds, Steps, Line_Z, Line_Y)
            Zx = PunctX[0]
            Yx = PunctX[1]
            
            # Arrays for plotting puncture plots
            xx = np.linspace(0, Xsize - 1, Xsize)
            yy = np.linspace(0, Ysize - 1, Ysize)
            zz = np.linspace(0, Zsize - 1, Zsize)
            
            print('\n' + 'Plotting...')            
            
            # The 3D plot
            fig2 = plt.figure(2)
            fig2.set_size_inches(10, 10, forward = True)
            fig2.patch.set_facecolor('lightgrey')
            fig2.suptitle('$Line$' + ' ' + '$Trace$' + ' ' + '$3D$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$' + '\n' + '$Started$' + ' ' + '(' + str(Xinit*SIMds) + ',' + str(Yinit*SIMds) + ',' + str(Zinit*SIMds) + ')' + '\n' + '$ds$' + ' ' + '$=$' + ' ' + str(ds*SIMds) + '$,$' + ' ' + ' ' + '$Passes$' + ' ' + '$=$' + ' ' + str(float(Steps)*ds/B[0].shape[0]), fontsize = 20, y = .96)
            
            ax4 = fig2.add_subplot(111, projection = '3d')
            ax4.plot(Line_Z, Line_X, Line_Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
            ax4.view_init(elev = 5, azim = 5)
            # rstride and cstride may be adjusted to change the resolution of
            # the 3D colormeshs, and the plotting speed, 10 means every 10th point,
            # 5 every 5th etc.
            ax4.plot_surface(ZZ, XZ, YZ, facecolors = colorsZ, shade = False, rstride = 10, cstride = 10)
            ax4.plot_surface(ZY, XY, YY, facecolors = colorsY, shade = False, rstride = 10, cstride = 10)
            ax4.plot_surface(ZX, XX, YX, facecolors = colorsX, shade = False, rstride = 10, cstride = 10) 
            ax4.set_zlim([0, Ysize-1])
            ax4.set_ylim([0, Xsize-1])
            ax4.set_xlim([0, Zsize-1])
            ax4.set_xlabel('Z')
            ax4.set_ylabel('X')
            ax4.set_zlabel('Y')   
            
            # Forcing matplot to generate the axes tick labels
            fig2.canvas.draw()
            
            # Returns an list of tick label objects
            LabelList = ax4.get_xticklabels()
            # Adjusting the tick labels
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            # Resetting the tick labels
            ax4.set_xticklabels(LabelList)
            
            LabelList = ax4.get_yticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax4.set_yticklabels(LabelList)
            
            LabelList = ax4.get_zticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax4.set_zticklabels(LabelList)
            
            fig2.tight_layout()            
            
            # The puncture plots
            fig3 = plt.figure(3)
            fig3.set_size_inches(5, 10, forward = True)
            fig3.suptitle('$Puncture$' + ' ' + '$Plots$' + '\n' + '$Started$' + ' ' + '(' + str(Xinit*SIMds) + ',' + str(Yinit*SIMds) + ',' + str(Zinit*SIMds) + ')' + '\n' + '$ds$' + ' ' + '$=$' + ' ' + str(ds*SIMds) + '$,$' + ' ' + ' ' + '$Passes$' + ' ' + '$=$' + ' ' + str(float(Steps)*ds/B[0].shape[0]), fontsize = 18, y = .99)
            plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
            fig3.patch.set_facecolor('lightgrey')
            
            ax5 = fig3.add_subplot(311)
            ax5.pcolormesh(xx, yy, Bm5.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
            ax5.scatter(Xz, Yz, c = 'k', s = 1)
            ax5.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax5.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax5.set_aspect('equal')
            ax5.set_title('$Slice$' + ' ' + '$Z$' + '$=$' + str(Zinit*SIMds), fontsize = 14)
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
            ax6.set_title('$Slice$' + ' ' + '$Y$' + '$=$' + str(Yinit*SIMds), fontsize = 14)
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
            ax7.set_title('$Slice$' + ' ' + '$X$' + '$=$' + str(Xinit*SIMds), fontsize = 14)
            ax7.locator_params(nbins = 6, axis = 'y')
            ax7.locator_params(nbins = 8, axis = 'x')
            ax7.set_ylim([0, Zsize])
            ax7.set_xlim([0, Ysize])
            ax7.set_ylabel('Z', rotation = 0)
            ax7.set_xlabel('Y')
            ax7.yaxis.set_label_position('right')
            ax7.yaxis.labelpad = 10
            
            # Forcing matplot to generate the axes tick labels
            fig3.canvas.draw()
            
            # Returns an list of tick label objects
            LabelList = ax5.get_xticklabels()
            # Adjusting the tick labels
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            # Resetting the tick labels
            ax5.set_xticklabels(LabelList)
            
            LabelList = ax6.get_xticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax6.set_xticklabels(LabelList)
            
            LabelList = ax7.get_xticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax7.set_xticklabels(LabelList)
            
            LabelList = ax5.get_yticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax5.set_yticklabels(LabelList)
            
            LabelList = ax6.get_yticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax6.set_yticklabels(LabelList)
            
            LabelList = ax7.get_yticklabels()
            for LabelObj in LabelList:
                try: 
                    LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                except:
                    1
            ax7.set_yticklabels(LabelList)            
            
            break
        elif ((ans == 'N') or (ans == 'n')):
            break
        else:
            print('invalid input try again \n')
        
    plt.show()


# 2D tracing wrapper and plotter
def FieldLine_2D(SIMds, Xinit, Yinit, B1, B2, Xsize, Ysize, ds, Steps):
    # Identifying the C tracing function
    func          = _lib.RK4_2D
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
    # Arrays for the traced line data points              
    Line_X = np.zeros(int(Steps))
    Line_Y = np.zeros(int(Steps))
    
    # Arrays for restructuring B
    Bx = np.zeros(Xsize * Ysize)
    By = np.zeros(Xsize * Ysize)

    # Restructuring Bx, and By into 1D arrays so they can be passed as
    # double pointers to the C function
    Bx = B1.reshape(Xsize * Ysize, order='F')
    By = B2.reshape(Xsize * Ysize, order='F')
    
    print('\n' + 'Tracing...')
    # Callling the C 2D tracing method which returns the number of steps
    # for the line to bite its own tail so the data arrays can be trimmed
    Length = func(Line_X, Line_Y, Xinit, Yinit, Bx, By, Xsize, Ysize, ds, Steps)
    # Trimming the data arryas to remove unused points
    Line_X = Line_X[:Length]
    Line_Y = Line_Y[:Length]
    
    # The figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(10,10, forward = True)
    fig1.patch.set_facecolor('lightgrey')
    
    # The plot
    ax = fig1.add_subplot(111)
    ax.set_ylim([0, Ysize-1])
    ax.set_xlim([0, Xsize-1])
    ax.set_ylabel('Y', rotation = 0)
    ax.set_xlabel('X')
    ax.yaxis.set_label_position('right')
    ax.yaxis.labelpad = 10
    
    print('\n' + 'Calculating |B|...')
    Bm = np.sqrt(B1**2 + B2**2)

    print('\n' + 'Plotting...')
    # Plotting |B| colormesh
    ax.pcolormesh(Bm.T, vmin = -.1)            
    ax.plot(Line_X,Line_Y, linestyle = 'none', marker = '.', markersize = .01, color='black')
    ax.set_title('$2D$' + ' ' + '$Line$' + ' ' + '$Trace$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$' + '\n' + '$Started$' + ' ' + '$at$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(Xinit*SIMds) + ', ' + '$Y$' + ' ' + '$=$' + ' ' + str(Yinit*SIMds), fontsize=20)
    
    # Forcing matplot to generate the axes tick labels
    fig1.canvas.draw()
    
    # Returns an list of tick label objects
    LabelList = ax.get_xticklabels()
    # Adjusting the tick labels
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    # Resetting the tick labels
    ax.set_xticklabels(LabelList)
    
    LabelList = ax.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax.set_yticklabels(LabelList)
    
    plt.show()                
        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#   
# PUNCTURE METHODS

# Puncture wrapper
# To make a puncture plot of Z = 100, PunctAxisData should be the SepZ data,
# OtherAxisData1 and 2 would be SepX and SepY
def Punct(PunctAxisData, Val, ds, Steps, OtherAxisData1, OtherAxisData2):
    # Identifying the C puncturing function
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
    
    # Arrays to hold puncture point coordinates       
    Points1 = np.zeros(Steps)
    Points2 = np.zeros(Steps)
    # Callling the C puncturing method which returns the total number of
    # puncture points to the arrays can be trimmed
    Total = func(PunctAxisData, Val, ds, Steps, OtherAxisData1, OtherAxisData2, Points1, Points2)
    # Trimming the data arrays
    Points1 = Points1[:Total]
    Points2 = Points2[:Total]
    
    # returns the puncture points
    return Points1, Points2

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# SEPARATOR MAPPING METHODS

# Master separator mapping method
def MapSeparator(SIMds, Saves, B, Ystart, UPorLOW, N = 1):
    # ds and passes are set to provide the fastest possible execution
    # By trial and error it was determined that ds must be <= 3 to avoid line
    # drifitng which will cause significant erros near a separator sheet, this
    # may relate to the debye length of the simulation
    # ds IS IN GRID SPACE
    ds = 3
    # Any field line which reconnects will do so within 1.5 passes
    passes = 1.5
    
    # Checking validity of arguments
    try:
        Saves[0] == 'String'
        Saves[1] == 'String'
        Saves[2] == 'String'
        B[0][0,0,0] * B[1][0,0,0] * B[2][0,0,0]
        assert(B[0].shape == B[1].shape == B[2].shape)
        assert(int(N))
        assert(float(Ystart))
        assert(1 <= N < (B[0].shape[0] / 30))
        assert(0 < Ystart/SIMds < B[0].shape[1])
        assert((UPorLOW == 'Upper') or (UPorLOW == 'Lower'))
        
    except:
        print('invalid arguments')
        return 0
        
    try:
        np.load(Saves[0])
        print('OVERWRITE WARNING, Save name already exists')
        return 0
    except:
        1
    try:
        np.load(Saves[1])
        print('OVERWRITE WARNING, Save name already exists')
        return 0
    except:
        1
    try:
        np.load(Saves[2])
        print('OVERWRITE WARNING, Save name already exists')
        return 0
    except:
        1
    
    # Converting to grid space
    Ystart = Ystart/SIMds
    
    # Arrays for separator data points
    SeparatorX = np.zeros((B[0].shape[0] / N) * (B[0].shape[2] / N))
    SeparatorY = np.zeros((B[0].shape[0] / N) * (B[0].shape[2] / N))
    SeparatorZ = np.zeros((B[0].shape[0] / N) * (B[0].shape[2] / N))
    
    # Calling the separator mapping wrapper
    SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ, Ystart, UPorLOW)    
    
    print('Saving...')
    # Saving the data
    np.save(Saves[0], SeparatorX)
    np.save(Saves[1], SeparatorY)
    np.save(Saves[2], SeparatorZ)
    print('Saved.')


# Separator mapping method, for each point in the X Z plane calls a C 
# method which searches for individual points on the separator surface
def SepPoints(N, B, ds, passes, SeparatorX, SeparatorY, SeparatorZ, Ystart, UPorLOWstr):
    # Identifying the C separator point function
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
                     ndpointer(ctypes.c_double),
                     ctypes.c_uint]
    
    # Determining if the program should search for an upper or lower separator
    if (UPorLOWstr == 'Upper'):
        UPorLOWint = 1
    else:
        UPorLOWint = 0
    
    # Setting the steps
    Steps = passes * (B[0].shape[0] / ds)
    
    # Setting the simulation size
    Xsize = B[0].shape[0]
    Ysize = B[0].shape[1]
    Zsize = B[0].shape[2]
    
    # Arrays for restructuring B
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    # Restructuring Bx, By, and Bz into 1D arrays so they can be passed as
    # double pointers to the C function
    Bx = B[0].reshape(Xsize * Ysize * Zsize, order='F')
    By = B[1].reshape(Xsize * Ysize * Zsize, order='F')
    Bz = B[2].reshape(Xsize * Ysize * Zsize, order='F')
    
    # A starting point for the separator mapping
    Start = np.zeros(3)
    Start[0] = 0
    Start[1] = Ystart
    Start[2] = 0
    
    # Finding the first point on the separator by calling the C function which
    # returns the Y value of the separator point at X = 0 and Z = 0
    Yval = func(int(N), Start, Bx, By, Bz, Xsize, Ysize, Zsize, float(ds), int(Steps), SeparatorX, SeparatorY, SeparatorZ, UPorLOWint)
    
    # Adding this point to the data arrays    
    SeparatorX[0] = 0
    SeparatorY[0] = Yval
    SeparatorZ[0] = 0

    # Loop to search for remaining separator points, this could have been
    # implemented in C, but this offfered a minimal speed increase and would
    # make it impossible to output status updates to the user which are
    # very helpful
    for i in range(0, Xsize/N):
        # Status update
        print(str(i+1) + '/' + str(Xsize/N))
        for j in range(0, Zsize/N):
            # Skip the first point which was already found
            if ((i == 0) and (j == 0)):
                continue
            # Identifying the next point to be searched for
            Start[0] = (N * i)
            # Use the previous point as a starting Y value for the next point,
            # greatly reducing the time required to find points
            Start[1] = SeparatorY[(j + (Zsize/N)*(i))-1]
            Start[2] = (N * j)
            # Calling the C separator point searching method
            Yval = func(int(N), Start, Bx, By, Bz, Xsize, Ysize, Zsize, float(ds), int(Steps), SeparatorX, SeparatorY, SeparatorZ, UPorLOWint)
            # Adding data points            
            SeparatorX[j + ((Zsize / N) * i)] = N * i
            SeparatorY[j + ((Zsize / N) * i)] = Yval
            SeparatorZ[j + ((Zsize / N) * i)] = N * j

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# SEPARATOR ANALYSIS METHODS

# Separator slice method
def SeparatorSlice(SIMds, PathSepX, PathSepY, PathSepZ, Xsize, Ysize, Zsize, B = None):
    # Checking validiity of arguments
    try:
        Bx = B[0]
        By = B[1]
        Bz = B[2]
        Bx[0,0,0] * By[0,0,0] * Bz[0,0,0]
        assert(Bx.shape == By.shape == Bz.shape)
        assert(Bx.shape == (int(Xsize/SIMds), int(Ysize/SIMds), int(Zsize/SIMds)))
        assert(int(Xsize))
        assert(int(Ysize))
        assert(int(Zsize))
    except:
        try:
            assert(B == None)
            assert(int(Xsize/SIMds))
            assert(int(Ysize/SIMds))
            assert(int(Zsize/SIMds))
        except:
            print('invalid arguments')
            return 0
    
    # Ensuring these are integers and converting to grid space
    Xsize = int(Xsize/SIMds)
    Ysize = int(Ysize/SIMds)
    Zsize = int(Zsize/SIMds)
    
    # Checking if the user would like to plot slices along the X or the Z axis
    # A slice along the Y axis is not useful
    while True:
        print(' ')
        Axis = str(raw_input('Slice X or Z value: \n'))
        # X slices
        if ((Axis == 'X') or (Axis == 'x')):
            while True:
                while True:
                    try:
                        # The value to slice
                        X = int(raw_input('Enter X value for slice: \n'))
                        assert(0 <= int(X/SIMds) < Xsize)
                        break;
                    except:
                        print('invalid X try again. \n')
                        continue
                
                X = int(X / SIMds)                
                
                print('\n' + 'Loading')
                # loading the separator data
                SepX = np.load(PathSepX)
                SepY = np.load(PathSepY)
                SepZ = np.load(PathSepZ)
                print('Loaded' + '\n')
                
                print('Slicing...')
                
                # Arrays for the separator X slice data
                SepYSlice = np.zeros(Zsize)        
                SepZSlice = np.zeros(Zsize)             
                
                # Identifying the points in the separator data corresponding to
                # the the X value to be sliced
                First = Zsize * (X)
                SepYSlice = SepY[First:(First + Zsize - 1)]
                SepZSlice = SepZ[First:(First + Zsize - 1)]
                
                print('\n' + 'Generating Figures...')
                # The figure
                fig = plt.figure(1)
                fig.set_size_inches(30,6, forward = True)
                fig.patch.set_facecolor('lightgrey')
                plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                
                # The plot
                ax = fig.add_subplot(111)
                ax.plot(SepZSlice, SepYSlice, linewidth = .5, color = 'g')
                ax.set_ylabel('Y', rotation = 0)
                ax.set_xlabel('Z')
                ax.yaxis.set_label_position('right')
                ax.yaxis.labelpad = 10
                ax.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(X*SIMds), fontsize=20)
                
                fig.canvas.draw()
        
                # Returns an list of tick label objects
                LabelList = ax.get_xticklabels()
                # Adjusting the tick labels
                for LabelObj in LabelList:
                    try: 
                        LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                    except:
                        1
                # Resetting the tick labels
                ax.set_xticklabels(LabelList)
                
                LabelList = ax.get_yticklabels()
                for LabelObj in LabelList:
                    try:
                        LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                    except:
                        1
                ax.set_yticklabels(LabelList)
                
                # If B is included, plots a puncture plot of the slice by
                # tracing a line started right above and another right below
                # the separator
                try:
                    Bz[1]
                    
                    # The slice / puncture plot
                    fig2 = plt.figure(2)
                    fig2.set_size_inches(30,6, forward = True)
                    fig2.patch.set_facecolor('lightgrey')
                    plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                    
                    ax2 = fig2.add_subplot(111)
                    ax2.plot(SepZSlice, SepYSlice, linewidth = 2, color = 'g', label = 'Separator Sheet')
                    ax2.set_ylabel('Y', rotation = 0)
                    ax2.set_xlabel('Z')
                    ax2.yaxis.set_label_position('right')
                    ax2.yaxis.labelpad = 10
                    ax2.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Puncture$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(X*SIMds), fontsize=20)
                    
                    # Both lines are traced for 100 passes
                    Steps = 1000*Bx.shape[0]
                    # Tracing the line above
                    Upper = FieldLine_3D_SepSlice(0, SepY[0] + .3, 0, Bx, By, Bz, Xsize, Ysize, Zsize, .1, Steps)
                    # And below
                    Lower = FieldLine_3D_SepSlice(0, SepY[0] - .3, 0, Bx, By, Bz, Xsize, Ysize, Zsize, .1, Steps)
                    # Finding puncture points of both lines
                    UpperPunct = Punct(Upper[0], X, .1, Steps, Upper[2], Upper[1])
                    LowerPunct = Punct(Lower[0], X, .1, Steps, Lower[2], Lower[1])
                    ax2.scatter(UpperPunct[0], UpperPunct[1], s = 3, c = 'r', label = 'Line Started Above')
                    ax2.scatter(LowerPunct[0], LowerPunct[1], s = 3, c = 'b', label = 'Line Started Below')
                    ax2.legend(loc='center right', fontsize = 'large')
                    
                    fig2.canvas.draw()
        
                    # Returns an list of tick label objects
                    LabelList = ax2.get_xticklabels()
                    # Adjusting the tick labels
                    for LabelObj in LabelList:
                        try: 
                            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                        except:
                            1
                    # Resetting the tick labels
                    ax2.set_xticklabels(LabelList)
                    
                    LabelList = ax2.get_yticklabels()
                    for LabelObj in LabelList:
                        try: 
                            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                        except:
                            1
                    ax2.set_yticklabels(LabelList)
                                   
                except:
                    1
                
                print('\n' + 'Finding XLine...')
                # Finding and plotting the X line of the separator
                XLine(SIMds, SepX, SepY, SepZ, Xsize, Ysize, Zsize)                
                
                print('\n' + 'Plotting...')
                plt.show()
                
                # Checks if the user would like to slice another point
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
        # Z slices
        elif ((Axis == 'Z') or (Axis == 'z')):
            while True:
                while True:
                    try:
                        Z = int(raw_input('Enter Z value for slice: \n'))
                        assert(0 <= int(Z/SIMds) < Zsize)
                        break
                    except:
                        print('invalid Z try again. \n')
                        continue
                
                Z = int(Z / SIMds)                
                
                # Loading the separator data
                print('\n' + 'Loading')
                SepX = np.load(PathSepX)
                SepY = np.load(PathSepY)
                SepZ = np.load(PathSepZ)
                print('Loaded' + '\n')
                
                print('Slicing...')
                
                # Arrays for the separator X slice data
                SepYSlice = np.zeros(Xsize)        
                SepXSlice = np.zeros(Xsize)        
                
                # Identifying the points in the separator data corresponding to
                # the the Z value to be sliced
                for i in range(0, Xsize):
                    SepYSlice[i] = SepY[Zsize * i + Z]
                    SepXSlice[i] = SepX[Zsize * i + Z]
                
                print('\n' + 'Generating Figures...')            
                # The figure
                fig = plt.figure(1)
                fig.set_size_inches(30,6, forward = True)
                fig.patch.set_facecolor('lightgrey')
                plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
    
                # The plot
                ax = fig.add_subplot(111)
                ax.plot(SepXSlice, SepYSlice, linewidth = .5, color = 'g')
                ax.set_ylabel('Y', rotation = 0)
                ax.set_xlabel('X')
                ax.yaxis.set_label_position('right')
                ax.yaxis.labelpad = 10
                ax.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(Z*SIMds), fontsize=20)
                
                fig.canvas.draw()
        
                # Returns an list of tick label objects
                LabelList = ax.get_xticklabels()
                # Adjusting the tick labels
                for LabelObj in LabelList:
                    try: 
                        LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                    except:
                        1
                # Resetting the tick labels
                ax.set_xticklabels(LabelList)
                
                LabelList = ax.get_yticklabels()
                for LabelObj in LabelList:
                    try: 
                        LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                    except:
                        1
                ax.set_yticklabels(LabelList)                
                
                # If B is included, plots a puncture plot of the slice by
                # tracing a line started right above and another right below
                # the separator
                try:
                    Bz[1]
                    
                    # The slice / puncture plot
                    fig2 = plt.figure(2)
                    fig2.set_size_inches(30,6, forward = True)
                    fig2.patch.set_facecolor('lightgrey')
                    plt.subplots_adjust(left = .05, bottom = .1, right = .95, top = .9)
                    
                    ax2 = fig2.add_subplot(111)
                    ax2.plot(SepXSlice, SepYSlice, linewidth = 2, color = 'g', label = 'Separator Sheet')
                    ax2.set_ylabel('Y', rotation = 0)
                    ax2.set_xlabel('X')
                    ax2.yaxis.set_label_position('right')
                    ax2.yaxis.labelpad = 10
                    ax2.set_title('$Separator$' + ' ' + '$Sheet$' + ' ' + '$Slice$' + ' ' + '$Puncture$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(Z*SIMds), fontsize=20)
                    
                    # Both lines are traced for 100 passes
                    Steps = 1000*Bx.shape[0]
                    # Tracing the line above
                    Upper = FieldLine_3D_SepSlice(0, SepY[0] + .3, 0, Bx, By, Bz, Xsize, Ysize, Zsize, .1, Steps)
                    # And below
                    Lower = FieldLine_3D_SepSlice(0, SepY[0] - .3, 0, Bx, By, Bz, Xsize, Ysize, Zsize, .1, Steps)
                    # Finding puncture points of both lines
                    UpperPunct = Punct(Upper[2], Z, .1, Steps, Upper[0], Upper[1])
                    LowerPunct = Punct(Lower[2], Z, .1, Steps, Lower[0], Lower[1])
                    ax2.scatter(UpperPunct[0], UpperPunct[1], s = 3, c = 'r', label = 'Line Started Above')
                    ax2.scatter(LowerPunct[0], LowerPunct[1], s = 3, c = 'b', label = 'Line Started Below')
                    ax2.legend(loc='center right', fontsize = 'large')
                    
                    fig2.canvas.draw()
        
                    # Returns an list of tick label objects
                    LabelList = ax2.get_xticklabels()
                    # Adjusting the tick labels
                    for LabelObj in LabelList:
                        try: 
                            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                        except:
                            1
                    # Resetting the tick labels
                    ax2.set_xticklabels(LabelList)
                    
                    LabelList = ax2.get_yticklabels()
                    for LabelObj in LabelList:
                        try: 
                            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
                        except:
                            1
                    ax2.set_yticklabels(LabelList)                    
                    
                except:
                    1
                
                print('\n' + 'Finding XLine...')
                # Finding and plotting the X line of the separator
                XLine(SIMds, SepX, SepY, SepZ, Xsize, Ysize, Zsize)
                
                print('\n' + 'Plotting...')
                plt.show()            
                
                # Checks if the user would like to slice another point
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


# 3D tracing wrapper and plotter for separator slice method
# same as FieldLine_3D accept it does not plot the lines, returning them to
# the separator slice method to be plotted instead
def FieldLine_3D_SepSlice(Xinit, Yinit, Zinit, B1, B2, B3, Xsize, Ysize, Zsize, ds, Steps):
    # Identifying the C tracing function
    func          = _lib.RK4_3D
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
                     ctypes.c_uint,
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double), 
                     ndpointer(ctypes.c_double),
                     ndpointer(ctypes.c_double)]
          
    # Arrays for the traced line data points           
    Line_X = np.zeros(Steps)
    Line_Y = np.zeros(Steps)
    Line_Z = np.zeros(Steps)
    
    # Arrays for restructuring B
    Bx = np.zeros(Xsize * Ysize * Zsize)
    By = np.zeros(Xsize * Ysize * Zsize)
    Bz = np.zeros(Xsize * Ysize * Zsize)

    # Restructuring Bx, By, and Bz into 1D arrays so they can be passed as
    # double pointers to the C function
    Bx = B1.reshape(Xsize * Ysize * Zsize, order='F')
    By = B2.reshape(Xsize * Ysize * Zsize, order='F')
    Bz = B3.reshape(Xsize * Ysize * Zsize, order='F')
    
    # E is a required argument for the C tracing function so here is set to 0
    Ex = np.zeros(3)
    Ey = np.zeros(3)
    Ez = np.zeros(3)
    EI = np.zeros(3)
    
    print('\n' + 'Tracing...')
    # Callling the C 3D tracing method
    func(Line_X, Line_Y, Line_Z, Xinit, Yinit, Zinit, Bx, By, Bz, Xsize, Ysize, Zsize, ds, Steps, Ex, Ey, Ez, EI)
    print('Done Tracing')
    
    return Line_X, Line_Y, Line_Z
    
def XLine(SIMds, SepX, SepY, SepZ, Xsize, Ysize, Zsize):
    # Arrays for the X line data points
    XLineX = np.zeros(Zsize)
    XLineY = np.zeros(Zsize)
    XLineZ = np.zeros(Zsize)
    
    # Checking if the X line is above or below the separator by checking if
    # various points in the middle of the separator are lower or higher than
    # the first point
    Linepos = ' '
    while True:
        # A fraction of the size of the separator
        num = 40
        # Comparing the first point, and the point at the fraction of the Sep.
        if SepY[0] > SepY[Zsize * (int(Xsize / num))]:
            Linepos = 'Lower'
            break
        elif SepY[0] < SepY[Zsize * (int(Xsize / num))]:
            Linepos = 'Upper'
            break
        else:
            # Adjusting the size of the fraction
            if num < 5:
                print('Could not find Xline')
                return
            num = num / 2
    
    # If a lower X line, program looks for all the lowest Y points for each
    # value of Z
    if (Linepos == 'Lower'):
        for k in range(0, Zsize):
            # Min value for search
            MinY = Ysize + 1
            MinX = 0
            MinZ = 0
            # Checks each X value for each Z
            for i in range(0, Xsize):
                if SepY[(Zsize*i) + k] < MinY:
                    MinY = SepY[(Zsize*i) + k]
                    MinX = i
                    MinZ = k
            # Adding data point
            XLineX[k] = MinX
            XLineY[k] = MinY
            XLineZ[k] = MinZ
    # If an upper X line, program looks for all the highest Y points for each
    # value of Z
    else:
        for k in range(0, Zsize):
            # Max value for search
            MaxY = -1
            MaxX = 0
            MaxZ = 0
            # Checks each X value for each Z
            for i in range(0, Xsize):
                if SepY[(Zsize*i) + k] > MaxY:
                    MaxY = SepY[(Zsize*i) + k]
                    MaxX = i
                    MaxZ = k
            # Adding data point
            XLineX[k] = MaxX
            XLineY[k] = MaxY
            XLineZ[k] = MaxZ
    
    # The X line plots
    fig3 = plt.figure(3)
    fig3.set_size_inches(14,12, forward = True)
    fig3.patch.set_facecolor('lightgrey')
    fig3.suptitle('$XLine$' + ' ' + '$Projections$', fontsize = 20)
    plt.subplots_adjust(left = .08, bottom = .1, right = .92, top = .9, hspace = .3)
    
    ax3 = fig3.add_subplot(211)
    ax3.plot(XLineZ, XLineY, linewidth = .5, color = 'b', label = 'Separator Sheet')
    ax3.set_ylabel('Y', rotation = 0)
    ax3.set_xlabel('Z')
    ax3.yaxis.set_label_position('right')
    ax3.yaxis.labelpad = 10
    ax3.set_title('$XLine$' + ' ' + '$Viewed$' + ' ' + '$Down$' + ' ' + '$X$' + ' ' + '$Axis$')
    
    ax4 = fig3.add_subplot(212)
    ax4.plot(XLineZ, XLineX, linewidth = .5, color = 'b', label = 'Separator Sheet')
    ax4.set_ylabel('X', rotation = 0)
    ax4.set_xlabel('Z')
    ax4.yaxis.set_label_position('right')
    ax4.yaxis.labelpad = 10
    ax4.set_title('$XLine$' + ' ' + '$Viewed$' + ' ' + '$Down$' + ' ' + '$Y$' + ' ' + '$Axis$')
    
    fig3.canvas.draw()
        
    # Returns an list of tick label objects
    LabelList = ax3.get_xticklabels()
    # Adjusting the tick labels
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    # Resetting the tick labels
    ax3.set_xticklabels(LabelList)
    
    LabelList = ax3.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax3.set_yticklabels(LabelList)
    
    LabelList = ax4.get_xticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax4.set_xticklabels(LabelList)
    
    LabelList = ax4.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax4.set_yticklabels(LabelList)
    

def SeparatorLoader(SIMds, PathSepX, PathSepY, PathSepZ, Xsize, Ysize, Zsize, B = None):
    # Checking validity of arguments
    try:
        Bx = B[0]
        By = B[1]
        Bz = B[2]
        Bx[0,0,0] * By[0,0,0] * Bz[0,0,0]
        assert(Bx.shape == By.shape == Bz.shape)
        assert(Bx.shape == (int(Xsize/SIMds), int(Ysize/SIMds), int(Zsize/SIMds)))
        assert(int(Xsize/SIMds))
        assert(int(Ysize/SIMds))
        assert(int(Zsize/SIMds))
    except:
        # If B = None the first check will throw an error
        try:
            assert(B == None)
            assert(int(Xsize/SIMds))
            assert(int(Ysize/SIMds))
            assert(int(Zsize/SIMds))
        except:
            print('invalid arguments')
            return 0
            
    # Ensuring these values are integers and converting to grid space
    Xsize = int(Xsize/SIMds)
    Ysize = int(Ysize/SIMds)
    Zsize = int(Zsize/SIMds)
    
    print('Loading Separator...')
    # Loading the separator data
    X = np.load(PathSepX)
    Y = np.load(PathSepY)
    Z = np.load(PathSepZ)
    print('Loaded \n')
    
    # Plotting the basic separator projections which can be done quickly 
    # since no colormeshes are involved
    print('Generating Figures...')
    # The figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(5, 10, forward = True)
    plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
    fig1.suptitle('$Separator$' + ' ' + '$Projections$', fontsize = 20)
    fig1.patch.set_facecolor('lightgrey')
    
    # The plots
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
    
    # Forcing matplot to generate the axes tick labels
    fig1.canvas.draw()
    
    # Returns an list of tick label objects
    LabelList = ax.get_xticklabels()
    # Adjusting the tick labels
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    # Resetting the tick labels
    ax.set_xticklabels(LabelList)
    
    LabelList = ax2.get_xticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax2.set_xticklabels(LabelList)
    
    LabelList = ax3.get_xticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax3.set_xticklabels(LabelList)
    
    LabelList = ax.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax.set_yticklabels(LabelList)
    
    LabelList = ax2.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax2.set_yticklabels(LabelList)
    
    LabelList = ax3.get_yticklabels()
    for LabelObj in LabelList:
        try: 
            LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
        except:
            1
    ax3.set_yticklabels(LabelList)
    
    # If B is included, plots the separator in 3D over colormeshes of |B|as 
    # well as a large projection of the separator down the Z axis over a 
    # colormesh of |B|, this takes extra time
    try:
        Bz[0,0,0]
        
        print('\n'+'Calculating |B|...')
        
        # Only calculates |B| for the 6 planes which are plotted to
        # avoid excessive time required to find |B| for the entire sim.
        Bm0 = np.sqrt(Bx[0,:,:]**2 + By[0,:,:]**2 + Bz[0,:,:]**2)
        Bm1 = np.sqrt(Bx[:,0,:]**2 + By[:,0,:]**2 + Bz[:,0,:]**2)
        Bm2 = np.sqrt(Bx[:,:,0]**2 + By[:,:,0]**2 + Bz[:,:,0]**2)
        Bm3 = np.sqrt(Bx[:,:,int(Zsize / 2)]**2 + By[:,:,int(Zsize / 2)]**2 + Bz[:,:,int(Zsize / 2)]**2)
        
        # Min and Max values for normalizing colormaps
        MIN = Bm2.min()
        MAX = Bm2.max()
        
        print('\n'+'Making Colormeshs...')
        # Arrays for making meshgrids
        Xval = np.linspace(0, Xsize - 1, Xsize)
        Yval = np.linspace(0, Ysize - 1, Ysize)
        Zval = np.linspace(0, Zsize - 1, Zsize)
    
        # Meshgrids to be plotted in 3D after applying the |B| colormap
        # allows plotting 2D colormesh in 3D
        XZ, YZ = np.meshgrid(Xval, Yval)
        # Field lines actually follow the transpose of the data
        XZ = XZ.T
        YZ = YZ.T
        # This meshgrid will be on the plane Z = 0
        ZZ = np.zeros((Xsize, Ysize))
        cmp = plt.cm.bwr
        # Normalizing the colomap
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsZ = cmp(norm(Bm2))
        
        XY, ZY = np.meshgrid(Xval, Zval)
        XY = XY.T
        ZY = ZY.T
        # This meshgrid will be on the plane Y = 0
        YY = np.zeros((Xsize, Zsize))
        cmp = plt.cm.bwr
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsY = cmp(norm(Bm1))
        
        YX, ZX = np.meshgrid(Yval, Zval)
        ZX = ZX.T
        YX = YX.T
        # This meshgrid will be on the plane X = 0
        XX = np.zeros((Zsize, Ysize))
        cmp = plt.cm.bwr
        norm = mpt.colors.Normalize(vmin = MIN, vmax = MAX)
        colorsX = cmp(norm(Bm0))     
        
        print('\n' + 'Plotting...')        
        
        # The 3D plot
        fig2 = plt.figure(2)
        fig2.set_size_inches(10, 10, forward = True)
        fig2.patch.set_facecolor('lightgrey')
        fig2.suptitle('$Separator$' + ' ' + '$3D$' + ' ' + '$Over$' + ' ' + '$|$' + '$B$' + '$|$', fontsize = 20, y = .95)

        ax4 = fig2.add_subplot(111, projection = '3d')        
        ax4.plot(Z, X, Y, linestyle = 'none', marker = '.', markersize = .01, color = 'k')
        ax4.view_init(elev = 10, azim = 10)
        # rstride and cstride may be adjusted to change the resolution of
        # the 3D colormeshs, and the plotting speed, 10 means every 10th point,
        # 5 every 5th etc.
        ax4.plot_surface(ZZ, XZ, YZ, facecolors = colorsZ, shade = False, rstride = 10, cstride = 10)
        ax4.plot_surface(ZY, XY, YY, facecolors = colorsY, shade = False, rstride = 10, cstride = 10)
        ax4.plot_surface(ZX, XX, YX, facecolors = colorsX, shade = False, rstride = 10, cstride = 10) 
        ax4.set_zlim([0, Ysize-1])
        ax4.set_ylim([0, Xsize-1])
        ax4.set_xlim([0, Zsize-1])
        ax4.set_xlabel('Z')
        ax4.set_ylabel('X')
        ax4.set_zlabel('Y') 
        
        # Forcing matplot to generate the axes tick labels
        fig2.canvas.draw()
        
        # Returns an list of tick label objects
        LabelList = ax4.get_xticklabels()
        # Adjusting the tick labels
        for LabelObj in LabelList:
            try: 
                LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
            except:
                1
        # Resetting the tick labels
        ax4.set_xticklabels(LabelList)
        
        LabelList = ax4.get_yticklabels()
        for LabelObj in LabelList:
            try: 
                LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
            except:
                1
        ax4.set_yticklabels(LabelList)
        
        LabelList = ax4.get_zticklabels()
        for LabelObj in LabelList:
            try: 
                LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
            except:
                1
        ax4.set_zticklabels(LabelList)
        
        fig2.tight_layout()
        
        # The colormesh projection plot
        fig3 = plt.figure(3)
        fig3.set_size_inches(18, 9, forward = True)
        plt.subplots_adjust(left = .1, bottom = .06, right = .9, top = .86, wspace = None, hspace = .45)
        fig3.suptitle('$Separator$' + ' ' + '$Projection$', fontsize = 20)
        fig3.patch.set_facecolor('lightgrey')
        
        # Arrays for plotting
        xx = np.linspace(0, Xsize - 1, Xsize)
        yy = np.linspace(0, Ysize - 1, Ysize)
        
        ax5 = fig3.add_subplot(111)
        ax5.set_title('$Over$' + ' ' + '$|$' + '$B$' + '$|$' + ' ' + '$with$' + ' ' + '$Z$' + ' ' + '$=$' + ' ' + str(int(Zsize/2)*SIMds), fontsize = 14)
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
        
        fig3.canvas.draw()
        
        # Returns an list of tick label objects
        LabelList = ax5.get_xticklabels()
        # Adjusting the tick labels
        for LabelObj in LabelList:
            try: 
                LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
            except:
                1
        # Resetting the tick labels
        ax5.set_xticklabels(LabelList)
        
        LabelList = ax5.get_yticklabels()
        for LabelObj in LabelList:
            try: 
                LabelObj.set_text(str(float(LabelObj.get_text())*SIMds))
            except:
                1
        ax5.set_yticklabels(LabelList)
        
    except:
        1

    plt.show()
    
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# TESTS

# 2D

#d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by'], 0)   
#Bx = d['by']
#By = d['bx']

#TraceField(.05, [Bx, By], [50, 50])


# 3D

print('Loading')
BX =  np.load('/scratch-fast/asym030/bx.npy')
BY =  np.load('/scratch-fast/asym030/by.npy')
BZ =  np.load('/scratch-fast/asym030/bz.npy')

#EX =  np.load('/scratch-fast/asym030/ex.npy')
#EY =  np.load('/scratch-fast/asym030/ey.npy')
#EZ =  np.load('/scratch-fast/asym030/ez.npy')
print('Loaded')

#SepY0 = np.load('LowerSepY.npy')[0]
#TraceField(.025, [BX, BY, BZ], [0, SepY0, 0], .075, 1.5)
#TraceField(.025, [BX, BY, BZ], [12.5, 5, 12.5], .0125, 100, ['TraceX6.npy', 'TraceY6.npy', 'TraceZ6.npy'])
#TraceField(.025, [BX, BY, BZ, EX, EY, EZ], [25.6, 12.8, 12.8], .0125, 100)
TraceField(.025, [BX, BY, BZ], [25.6, 12.8, 12.8], .0025, 100)

#SeparatorSlice(.025, 'LowerSepX.npy', 'LowerSepY.npy', 'LowerSepZ.npy', 51.2, 25.6, 25.6, [BX, BY, BZ])
#SeparatorSlice(.025, 'UpperSepX.npy', 'UpperSepY.npy', 'UpperSepZ.npy', 51.2, 25.6, 25.6)

#SeparatorLoader(.025, 'LowerSepX.npy','LowerSepY.npy','LowerSepZ.npy', 51.2, 25.6, 25.6, [BX, BY, BZ])


# 318*.025 for Upper
# 150*.025 for Lower
#MapSeparator(.025, ['SepX10.npy', 'SepY10.npy', 'SepZ10.npy'], [BX, BY, BZ], 8, 'Upper', 16)
#MapSeparator(.025, ['SepX9.npy', 'SepY9.npy', 'SepZ9.npy'], [BX, BY, BZ], 3.75, 'Lower', 8)

#SeparatorLoader(.025, 'SepX10.npy','SepY10.npy','SepZ10.npy', 51.2, 25.6, 25.6)

