
# 2D Field Line Tracing Program RK4

import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie
from scipy.ndimage.filters import gaussian_filter


# line tracing methods
#------------------------------------------------------------------------------

def CalcSlopes(x, y, Bx, By):
    # distance between field line point (X,Y) and lower left data grid 
    #point (i,j), such that Wx = X - i, Wy = Y - j
    if x < 0 and y < 0:
        i = -1
        Wx = 1 + x
        j = -1
        Wy = 1 + y
    elif x < 0:
        i = -1
        Wx = 1 + x
        Wy = y % 1
        j = int(y - Wy)
    elif y < 0:
        j = -1
        Wy = 1 + y
        Wx = x % 1
        i = int(x - Wx)
    else:
        Wx = x % 1
        Wy = y % 1
        # corresponding closest lower left data grid point (i,j) to field line
        # point (X,Y)
        i = int(x - Wx)
        j = int(y - Wy)
    # other three nearest neighbors identified with i+1 and j+1
    i1 = i + 1
    j1 = j + 1
    
    # check i and j for edge conditions
    #if (i != iStore or j != jStore):
    if i < 0:
        i = 8191
        i1 = 0
    if i > 8191:
        i = 0
        i1 = 1
        
    if j < 0:
        j = 8191
        j1 = 0
    if j > 8191:
        j = 0
        j1 = 1
        
    if i == 8191:
        i1 = 0
    if j == 8191:
        j1 = 0
    
    #if (i != iStore or j != jStore):
    # identifying Bx, By at closest 4 data grid points
    Bx_ij = Bx[i,j]
    Bx_i1j = Bx[i1,j]
    Bx_ij1 = Bx[i,j1]
    Bx_i1j1 = Bx[i1,j1]
    
    By_ij = By[i,j]
    By_i1j = By[i1,j]
    By_ij1 = By[i,j1]
    By_i1j1 = By[i1,j1]
            
            # stores coordinates of lower left neighbor to test if neighbors
            # have changed between steps so they don't need to be recalculated
            # saves considerable run time
            #iStore = i
            #jStore = j
    
    # finding average Bx, By and Bm from 4 nearest neighboring data points
    # (i,j), (i+1,j), (i,j+1) and (i+1,j+1) at field line point (X,Y)
    B_Wx = (1-Wx)*(1-Wy)*Bx_ij + (1-Wx)*Wy*Bx_ij1 + Wx*(1-Wy)*Bx_i1j + Wx*Wy*Bx_i1j1
    B_Wy = (1-Wx)*(1-Wy)*By_ij + (1-Wx)*Wy*By_ij1 + Wx*(1-Wy)*By_i1j + Wx*Wy*By_i1j1
    B_Wm = np.sqrt(B_Wx**2 + B_Wy**2)
    
    return B_Wx/B_Wm, B_Wy/B_Wm

def Line(InitX, InitY, Bx, By): 
    
    # maximum number of steps along line
    MaxSteps = 3000000
    
    K1x = 0
    K2x = 0
    K3x = 0
    K4x = 0
    
    K1y = 0
    K2y = 0
    K3y = 0
    K4y = 0
    
    # arrays to hold X and Y coordinates of field line points
    Line_X = np.zeros(MaxSteps)
    Line_Y = np.zeros(MaxSteps)
    
    # initial point for a field line
    X = InitX
    Y = InitY
    
    # position holders to avoid recalculating nearest neighbors if new (X,Y)
    # field line coordinate lies between the same data grid spaces as the
    # previous step
    #iStore = -1
    #jStore = -1
    
    # differential step towards next point on field line
    # dx has little effect on accuracy of tracing, accuracy is most strongly
    # affected by the curl of the field line
    dx = .02
    # additional dx variable incase of variable dx
    dx_Original = dx
    
    # loop to step forward line from initial point
    for step in range(0,MaxSteps):
        
        # checking if X or Y has moved outside range (-.5 to 8191.5) and if so
        # switches to other side of data grid
        if X < -.5:
            X = X + 8192
        if Y < -.5:
            Y = Y + 8192
        if X > 8191.5:
            X = X - 8192
        if Y > 8191.5:
            Y = Y - 8192
        
        # just a status update, sometimes it takes a long time
        #if (step % 1000000) == 0 and step != 0:
        #    print("Step = ",step)
        
        # adding points (X,Y) to field line
        Line_X[step] = X/80
        Line_Y[step] = Y/80
        
        # breaks out of loop when X and Y are within dx of their values at 
        # the second step (first step not chosen since it has higher error)
        if step > 1000 and abs(X - (Line_X[2]*80)) < 25*dx_Original and abs(Y - (Line_Y[2]*80)) < 25*dx_Original:
            print('Line Completed')
            # trimming the arrays to proper length after completion
            Line_X = Line_X[:step + 1]
            Line_Y = Line_Y[:step + 1]
            break
        
        Slopes = CalcSlopes(X, Y, Bx, By)
        SlopeX = Slopes[0]
        SlopeY = Slopes[1]
        
        K1x = SlopeX
        K1y = SlopeY
        
        Slopes2 = CalcSlopes(X + (dx/2)*SlopeX , Y + (dx/2)*SlopeY, Bx, By)
        SlopeX2 = Slopes2[0]
        SlopeY2 = Slopes2[1]
        
        K2x = SlopeX2
        K2y = SlopeY2
        
        Slopes3 = CalcSlopes(X + (dx/2)*SlopeX2 , Y + (dx/2)*SlopeY2, Bx, By)
        SlopeX3 = Slopes3[0]
        SlopeY3 = Slopes3[1]
        
        K3x = SlopeX3
        K3y = SlopeY3
        
        Slopes4 = CalcSlopes(X + dx*SlopeX3 , Y + dx*SlopeY3, Bx, By)
        SlopeX4 = Slopes4[0]
        SlopeY4 = Slopes4[1]
        
        K4x = SlopeX4
        K4y = SlopeY4
        
        
        X = X + (dx/6)*(K1x + 2*K2x + 2*K3x + K4x)
        Y = Y + (dx/6)*(K1y + 2*K2y + 2*K3y + K4y)

    # returns line as two arrays of X and Y points
    return Line_X, Line_Y
    
    
# data loading and plotting method
#------------------------------------------------------------------------------    

def Draw_2D_Field_Lines():
    
    # loading data
    d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by', 'jz'], 0)
    # for loading just a piece of the data     
    # d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by', 'jz'], 0, slc=np.s_[0,0,:100,:100])
    
    # subsampling may be necessary if plotting over Jz
    #def SubSample(Par_Orig, rate):
    #    print('subsampling....')
    #    Par_Sub = np.zeros((8192/rate, 8192/rate))
    #    for i in range(0,8192/rate):
    #        for n in range(0,8192/rate):
    #            Par_Sub[i,n] = Par_Orig[i*rate,n*rate]
    #    return Par_Sub
    
    # Subsampling has been removed, to reimplement must divide indexes
    # by SS_Rate where appropriate
    # SS_Rate may be set to 1 to avoid subsampling without disturbing code
    #SS_Rate = 1
    # subsampling calls if necessary
    #d['jz'] = SubSample(d['jz'], SS_Rate)
    #d['bx'] = SubSample(d['bx'], SS_Rate)
    #d['by'] = SubSample(d['by'], SS_Rate)    
    
    # the magnetic field components, gaussian filter is applied to imporve
    # tracing accuracy
    # Bx and By are switched intentionally
    print('filtering...')
    By = gaussian_filter(d['bx'], sigma = 1, mode = 'wrap')
    Bx = gaussian_filter(d['by'], sigma = 1, mode = 'wrap')
    
    # the figure
    fig1 = plt.figure(1)
    fig1.set_size_inches(10,10, forward = True)
    ax = fig1.add_subplot(111)
    
    # Plotting single field line
    # most troublesome line at 4596,4596, movie time 199
    #print('tracing')
    #Contour = Line(4596, 4596, Bx, By)
    #print('plotting')
    #ax.plot(Contour[0],Contour[1], linestyle = 'none', marker = '.')
    
    # dictionary is used to conveiniently calculate 64 lines starting at points in
    # an evenly spaced grid across the 2-D field
    Lines = dict()
    for t in range(1,9):
        for q in range(1,9):
            print("Tracing Line # " + str(((t-1)*8)+q) + " of 64")
            Lines[((t-1)*8)+q] = Line(((1000)*t)-(404),((1000)*q)-(404), Bx, By)
#    
#    # plotting the field lines
    for i in range(1,65):
        ax.plot(Lines[i][0],Lines[i][1], linestyle = 'none', marker = '.', markersize = .1, color='black')
        ax.set_ylim([0, 102.4])
        ax.set_xlim([0, 102.4])
        
    # plotting Jz if desired
#   #ax.pcolormesh(d['jz'].T)
#    
    plt.show()
    

# user interface
#-------------------------------------------------------------------------------

i = 0
while(i == 0):
    # catch for non-string entry
    while True:
        try:
            X = str(raw_input('Load Movie? Y or N \n '))
            break
        except ValueError:
            print("invalid input try again...\n")
    # handling string entries
    if ((X == 'Y') or (X == 'y')):
        i = 1        
        Draw_2D_Field_Lines()
    elif((X == 'N') or (X == 'n')):
        i = 1
        print("Load Data Cancelled")
    else:
        print("invalid input try again...\n")
