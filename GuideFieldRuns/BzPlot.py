import numpy as np
import matplotlib as mpt
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import sys
sys.path.append('/global/homes/r/ransom/Py3D')
from py3d.sub import *


run = raw_input('Which run? \n')
movies = '/global/cscratch1/sd/ransom/' + str(run) + '/staging'
paramFile = 'param_'+ str(run)
runNum = run.replace('R', '')


print(' ')
print('Loading...')
d = load_movie(0, paramFile, movies, ['bx', 'by', 'bz', 'ex', 'ey', 'ez', 'ni', 'ne', 'jex', 'jey', 'jez', 'jix', 'jiy', 'jiz', 'pexx', 'peyy', 'pezz', 'pexy', 'pexz', 'peyz', 'pixx', 'piyy', 'pizz', 'pixy', 'pixz', 'piyz'])
print('Loaded')
print(' ')
print('Calculating params...')
dx = d['xx'][1] - d['xx'][0]

Bx = d['bx']
By = d['by']
Bz = d['bz']
Bm = np.sqrt(Bx**2 + By**2 + Bz**2)

Ex = d['ex']
Ey = d['ey']
Ez = d['ez']
Em = np.sqrt(Ex**2 + Ey**2 + Ez**2)

EdotB = Bx*Ex + By*Ey + Bz*Ez
EcrBx = (Ey*Bz) - (Ez*By)
EcrBy = (Ez*Bx) - (Ex*Bz)
EcrBz = (Ex*By) - (Ey*Bx)
EcrBm = np.sqrt(EcrBx**2 + EcrBy**2 + EcrBz**2)

Ni = d['ni']
Ne = d['ne']

# Velocity = current/number density
Vix = d['jix']/Ni
Viy = d['jiy']/Ni
Viz = d['jiz']/Ni

Vex = d['jex']/Ne
Vey = d['jey']/Ne
Vez = d['jez']/Ne

Jx = d['jix'] + d['jex']
Jy = d['jiy'] + d['jey']
Jz = d['jiz'] + d['jez']

# Temperature (energy (F*d)) = pressure (F/d^2) / number density (1/d^3)
Texx = d['pexx']/Ne
Teyy = d['peyy']/Ne
Tezz = d['pezz']/Ne

Texy = d['pexy']/Ne
Texz = d['pexz']/Ne
Teyz = d['peyz']/Ne

Tixx = d['pixx']/Ni
Tiyy = d['piyy']/Ni
Tizz = d['pizz']/Ni

Tixy = d['pixy']/Ni
Tixz = d['pixz']/Ni
Tiyz = d['piyz']/Ni

# The temperature tensor is a strange concept, but it related to the pressure as stated above.
# Imagine a box, Pexx is the pressure on the X-normal face of the box acting in the X direction, 
# so Pexx, Peyy, and Pezz can be imagined as the force on each side acting along the three basis vecotors,
# normal to each side, while Pexy is the pressure on an X-normal side in the Y direction, and therefore
# can be pictured as a torque acting to spin the box.  Luckily, Pexy = Peyx, Pexz = Pezx and Peyz = Pezy,
# so the torques balance and no rotation occurs.
# The standard basis is {<1, 0, 0>, <0, 1, 0>, <0, 0, 1>}
# normalized to {<1/sqrt(3), 0, 0>, <0, 1/sqrt(3), 0>, <0, 0, 1/sqrt(3)>}
TeTens = np.array([ [Texx, Texy, Texz], 
                    [Texy, Teyy, Teyz], 
                    [Texz, Teyz, Tezz] ])
                    
TiTens = np.array([ [Tixx, Tixy, Tixz],
                    [Tixy, Tiyy, Tiyz],
                    [Tixz, Tiyz, Tizz] ])

# {u1, u2, u3} form a new basis where u1 is parallel to the magnetic field, and u2 and u3 are perependicular 
# to B.  u3 = Bx(ExB)
u1 = np.array([Bx, By, Bz])/Bm
u2 = np.array([EcrBx, EcrBy, EcrBz])/EcrBm
u3 = np.cross(u1, u2, axis=0)

# the next step is to rotate the Temperature tensor from the standard basis to the new B par/perp basis.
# this is done by defining a rotation matrix Q such that
# Q = [ [u1]   =  [[bx           by           bz]
#       [u2]       [EcrBx        EcrBy        EcrBz]
#       [u3] ]     [Bcr(EcrB)x   Bcr(EcrB)y   Bcr(EcrB)z]]
# this matrix is found by taking dot products of the old basis and the new basis such that
# Q11 = <bx, by, bz> * <1, 0, 0> = bx and
# Q12 = <bx, by, bz> * <0, 1, 0> = by etc.
# these dot products represent the cosine of the angle between the old basis vectors
# and the new.

# let the standard basis be denoted e, and the new basis b,
# then TeTens{b} = Q x TeTens{e} x Q.T
# and TeTens{b} will have Tpar, Tperp1, Tperp2 along its diagonal like this:
# TeTens{b} = [ [Tpar   Tprod   Tprod]
#               [Tprod  Tperp1  Tprod]
#               [Tprod  Tprod   Tperp2] ]
# where Tprod are the various products of the tensor which are unimportant
TeU1 = np.array([np.sum(u1*TeTens[0,:,:,:], axis=0), np.sum(u1*TeTens[1,:,:,:], axis=0), np.sum(u1*TeTens[2,:,:,:], axis=0)])
TeU2 = np.array([np.sum(u2*TeTens[0,:,:,:], axis=0), np.sum(u2*TeTens[1,:,:,:], axis=0), np.sum(u2*TeTens[2,:,:,:], axis=0)])
TeU3 = np.array([np.sum(u3*TeTens[0,:,:,:], axis=0), np.sum(u3*TeTens[1,:,:,:], axis=0), np.sum(u3*TeTens[2,:,:,:], axis=0)])

TiU1 = np.array([np.sum(u1*TiTens[0,:,:,:], axis=0), np.sum(u1*TiTens[1,:,:,:], axis=0), np.sum(u1*TiTens[2,:,:,:], axis=0)])
TiU2 = np.array([np.sum(u2*TiTens[0,:,:,:], axis=0), np.sum(u2*TiTens[1,:,:,:], axis=0), np.sum(u2*TiTens[2,:,:,:], axis=0)])
TiU3 = np.array([np.sum(u3*TiTens[0,:,:,:], axis=0), np.sum(u3*TiTens[1,:,:,:], axis=0), np.sum(u3*TiTens[2,:,:,:], axis=0)])

TePar = np.sum(TeU1*u1, axis=0)
TiPar = np.sum(TiU1*u1, axis=0)

TePerp1 = np.sum(TeU2*u2, axis=0)
TiPerp1 = np.sum(TiU2*u2, axis=0)

TePerp2 = np.sum(TeU3*u3, axis=0)
TiPerp2 = np.sum(TiU3*u3, axis=0)

# The trace, (sum of diagonal) of a tensor is invariant, so Texx + Teyy + Tezz = Tpar + Tperp1 + Tperp2
# Plasmas are generally gyrotropic meaning that Tperp1 = Tperp2 != Tpar
print(' ')
print('Tpar + Tperp1 + Tperp2 / Txx + Tyy + Tzz = ') 
print(np.average((TePar + TePerp1 + TePerp2)/(Texx + Teyy + Tezz)))

dBx = (np.roll(Bx, 1, 0) - np.roll(Bx, -1, 0))/dx
dBy = (np.roll(By, 1, 1) - np.roll(By, -1, 1))/dx
divB = dBx + dBy

print('DivB Average = ')
print(np.average(divB))

print(' ')


print('Done')
print(' ')

#---------------------------------------------------#
def showTestFig(bz, Xpoints):
    MIN = bz.min()
    MAX = bz.max()
    
    XSize = bz.shape[0]
    YSize = bz.shape[1]
    
    fig1 = plt.figure(1)
    fig1.set_size_inches(12, 6, forward = True)
    fig1.suptitle('$Test$' + ' ' + '$Colormesh$', fontsize = 18, y = .99)
    plt.subplots_adjust(left = .08, bottom = .08, right = .92, top = .92) #wspace = None, hspace = .45)
    fig1.patch.set_facecolor('lightgrey')
    
    ax1 = fig1.add_subplot(111)
    ax1.pcolormesh(bz.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
    ax1.plot([Xpoints[0], Xpoints[2]], [Xpoints[1], Xpoints[3]], 'x', color = 'k')
    ax1.set_aspect('equal')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.locator_params(nbins = 8, axis = 'x')
    ax1.locator_params(nbins = 6, axis = 'y')
    ax1.set_xlim([0, XSize])
    ax1.set_ylim([0, YSize])
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y', rotation = 0)
    ax1.yaxis.set_label_position('right')
    ax1.yaxis.labelpad = 10
    
    plt.show()

def findBackground(Slice):
    # starts at max and min values and shrinks the gap between the upper and
    # lower bounds until 50% of the data is between the upper and lower bound.
    # the average of the upper and lower bound then represents the background value of B
    minimum = Slice.min()
    maximum = Slice.max()
    DataInRange = 1
    while DataInRange > .5:
        minimum = minimum + .001
        pointsInRange = 0;
        for i in range(0, len(Slice)):
            if minimum < Slice[i]:
                pointsInRange = pointsInRange + 1
        DataInRange = float(pointsInRange)/float(len(Slice))
    minimum = minimum - .1
    DataInRange = 1
    while DataInRange > .5:
        maximum = maximum - .001
        pointsInRange = 0;
        for i in range(0, len(Slice)):
            if Slice[i] < maximum:
                pointsInRange = pointsInRange + 1
        DataInRange = float(pointsInRange)/float(len(Slice))
    maximum = maximum + .1
    back = (minimum + maximum)/2
    
    return back

# returns enhance width, suppress width and sorted list of min, middle and max values demarcating widths
def findWidths(Slice, back, ds, UpOrLow):
    # starts at max and min values, walking left and right from each until hitting
    # the backgournd level of B, marking the edge of the enhanced and suppressed 
    # regions of B.
    maxSlice = Slice.max()
    minSlice = Slice.min()
    minYVal = 0
    maxYVal = 0
    for i in range(0,len(Slice)):
        if Slice[i] == maxSlice:
             maxYVal = i
        if Slice[i] == minSlice:
             minYVal = i

    if maxYVal > minYVal:
        i = maxYVal
        while True:
            i = i + 1
            if Slice[i] <= back + .01:
                maxUpBound = i
                break
            
        i = maxYVal
        while True:
            i = i - 1
            if Slice[i] <= back:
                maxLowBound = i
                break
            
        i = minYVal
        while True:
            i = i + 1
            if Slice[i] >= back:
                minUpBound = i
                break
            
        i = minYVal
        while True:
            i = i - 1
            if Slice[i] >= back - .01:
                minLowBound = i
                break
            
    if maxYVal < minYVal:
        i = maxYVal
        while True:
            i = i + 1
            if Slice[i] <= back:
                maxUpBound = i
                break
            
        i = maxYVal
        while True:
            i = i - 1
            if Slice[i] <= back + .01:
                maxLowBound = i
                break
            
        i = minYVal
        while True:
            i = i + 1
            if Slice[i] >= back - .01:
                minUpBound = i
                break
            
        i = minYVal
        while True:
            i = i - 1
            if Slice[i] >= back:
                minLowBound = i
                break
        
    EnhancWidth = (maxUpBound - maxLowBound)*ds
    SuppressWidth = (minUpBound - minLowBound)*ds
    
    widths = np.array([minUpBound, maxUpBound, minLowBound, maxLowBound])  
    if UpOrLow == 'Upper':
        widths = widths + len(Slice)
    Sorted = False
    while Sorted == False:
        Sorted = True
        for i in range(0, 3):
            if widths[i] > widths[i + 1]:
                holder = widths[i + 1]
                widths[i + 1] = widths[i]
                widths[i] = holder
                Sorted = False
    
    widths[1] = (widths[1] + widths[2])/2
    widths = np.array([widths[0], widths[1], widths[3]])*ds
    
    return [EnhancWidth, SuppressWidth, widths]
        
# B = Zhat cross - grad(Psi)
# <Bx, By, Bz> = <0, 0, 1> X <-dPsi/dx, -dPsi/dy, -dPsi/dz>
# therefore Bx = dPsi/dy, and By = -dPsi/dx
# so Psi[i, j] = Psi[i, j] - Psi[0, 0] = lineIntegral(B) along path from [0, 0] to [i, j]
# path has 2 parts, path from [0, 0] to [0, j] and then [0, j] to [i, j]
# so Psi = sum(Bx * dy) - sum(By * dx)
# this is the magnetic scalar potential
def calcPsi(bx, by):
    Psi = 0.0*bx
    Psi[0,1:] = np.cumsum(Bx[0,1:])*dx
    Psi[1:,:] = (Psi[0,:] - np.cumsum(by[1:,:], axis=0)*dx)
    
    return Psi
    
# by observing a plot of Psi, it is clear that the upper X Point
# lies along the Line Y = j where j is the index such that the sum of Psi[x, j]
# along all values of x is a maximum.  The upper X point also lies on the line 
# X = i where i is the index such that the sum of Psi[i, y] for values of y between
# Ysize/2 and Ysize is a minimum (for lower X Point, sum y values from 0 to Ysize/2.)
# the upper X point is at the intersection of the maximum y line, and minimum upper half x line.
# the lower X point is at the intersection of the minimum y line, and maximum lower half x line.
def findXpt(Psi):
    rowTotals = np.zeros(YSize)
    for j in range(0, YSize):
        rowTotals[j] = np.sum(Psi[:, j])
        
    rowMin = rowTotals.min()
    rowMax = rowTotals.max()
    for j in range(0, YSize):
        if rowTotals[j] == rowMin:
            jLower = j
        if rowTotals[j] == rowMax:
            jUpper = j
    
    
    colTotals = np.zeros(XSize)
    for i in range(0, XSize):
        colTotals[i] = np.sum(Psi[i,:(YSize/2)])
    
    lowColMax = colTotals.max()
    for i in range(0, XSize):
        if colTotals[i] == lowColMax:
            iLower = i
            
    for i in range(0, XSize):
        colTotals[i] = np.sum(Psi[i,(YSize/2):])
    
    uppColMin = colTotals.min()
    for i in range(0, XSize):
        if colTotals[i] == uppColMin:
            iUpper = i
    
    return [iUpper, jUpper, iLower, jLower]
    
# gets energy array from stdout file
def getEnergy(stdoutPath = None):
    
        if stdoutPath is None:
            stdoutPath = raw_input('Enter path for p3d.stdout file: ')
            
        foundData = 0
        try:
            stdout = open(stdoutPath, 'r')
            #fort = open(fortPath, 'r')
            eTot = []
            eEM = []
            eK = []
            eKe = []
            eKi = []
            for line in stdout:
                if line.find('ENERGY') > -1 and line.find('ENERGY:') < 0:              
                    eTot.append(float(line.split()[1:2][0]))
                    eEM.append(float(line.split()[2:3][0]))
                    eK.append(float(line.split()[3:4][0]))
                    foundData = 1
            stdout.close()
            if foundData == 0:
                raise Exception
        except: 
            print('Could not load energies')
    
        return np.array(eTot), np.array(eEM), np.array(eK), np.array(eKi), np.array(eKe)

# returns the 100% cyclo dE and total cyclo times ran
def EConPerc(Run):
    for i in range(1, 10):
            try:
                test = open('/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + str(i), 'r')
                num = i
            except:
                pass
    if int(num) == 0: 
        pathStdout = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + num
        pathToMovieData = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging'
        #pathExp = '/global/cscratch1/sd/ransom/' + str(Run) + '/exp3d_' + str(Run)
        paramFile = 'param_'+ str(Run)
        paramPath = pathToMovieData + '/' + paramFile
        
        param = open(paramPath, 'r')
        for line in param:
            if line.find('define') > 0 and line.find('dt') > 0:              
                dt = float(line.split()[2:3][0])        
        
        
        Energies = getEnergy(pathStdout)
        
        return [(100*(Energies[0][-1] - Energies[0][0])/Energies[0][0])*(100/(dt*Energies[0].shape[0])), len(Energies[0])*dt]
    else:
        pathStdoutRestart = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + str(num)
        pathStdoutOriginal = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging_00/staging/p3d.stdout.000'
        EnergiesLength = 0;
        for i in range(0,int(num)):
            pathStdout = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging_0' + str(i) + '/staging/p3d.stdout.00' + str(i)
            Energy = getEnergy(pathStdout)
            EnergiesLength = EnergiesLength + len(Energy[0])
        
        pathToMovieData = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging'
        #pathExp = '/global/cscratch1/sd/ransom/' + str(Run) + '/exp3d_' + str(Run)
        paramFile = 'param_'+ str(Run)
        paramPath = pathToMovieData + '/' + paramFile
        
        param = open(paramPath, 'r')
        for line in param:
            if line.find('define') > 0 and line.find('dt') > 0:              
                dt = float(line.split()[2:3][0])        
        
        
        EnergiesRestart = getEnergy(pathStdoutRestart)
        EnergiesOriginal = getEnergy(pathStdoutOriginal)
        
        EnergiesLength = EnergiesLength + len(EnergiesRestart[0])     
        
        return [(100*(EnergiesRestart[0][-1] - EnergiesOriginal[0][0])/EnergiesOriginal[0][0])*(100/(dt*EnergiesLength)), EnergiesLength*dt]


        
#--------------------------------------------------#
        
XSize = Bz.shape[0]
YSize = Bz.shape[1]

psi = calcPsi(Bx, By)
xpts = findXpt(psi)
showTestFig(Bz, xpts)

while True:
    UpLow = raw_input('Upper or Lower X point? U or L: \n')
    if UpLow == 'U' or UpLow == 'u':
        UpLow = 'Upper'
        break
    elif UpLow == 'L' or UpLow == 'l':
        UpLow = 'Lower'
        break
    else:
        continue
    
fifteenDi = int(15/dx)
while True:
    leftRight = raw_input('Left or Right of X point? L or R: \n')
    if leftRight == 'R' or leftRight == 'r':
        if UpLow == "Upper":
            SliceXVal = xpts[0] + fifteenDi
        if UpLow == "Lower":
            SliceXVal = xpts[2] + fifteenDi
        break
    elif leftRight == 'L' or leftRight == 'l':
        if UpLow == "Upper":
            SliceXVal = xpts[0] - fifteenDi
        if UpLow == "Lower":
            SliceXVal = xpts[2] - fifteenDi
        break
    else:
        continue
    try:
        SliceXVal = int(raw_input('Enter X Slice Value: \n'))
        assert(0 <= SliceXVal <= XSize - 1)
        break
    except:
        continue

if UpLow == 'Upper':
    BzSlice = Bz[SliceXVal, (YSize/2):]
    BxSlice = Bx[SliceXVal, (YSize/2):]
    BySlice = By[SliceXVal, (YSize/2):]
    BmSlice = Bm[SliceXVal, (YSize/2):]
    
    ExSlice = Ex[SliceXVal, (YSize/2):]
    EySlice = Ey[SliceXVal, (YSize/2):]
    EzSlice = Ez[SliceXVal, (YSize/2):]
    EmSlice = Em[SliceXVal, (YSize/2):]
    
    EdotBSlice = EdotB[SliceXVal, (YSize/2):]
    EcrBxSlice = EcrBx[SliceXVal, (YSize/2):]
    EcrBySlice = EcrBy[SliceXVal, (YSize/2):]
    EcrBzSlice = EcrBz[SliceXVal, (YSize/2):]
    EcrBmSlice = EcrBm[SliceXVal, (YSize/2):]
    
    NiSlice = Ni[SliceXVal, (YSize/2):]
    NeSlice = Ne[SliceXVal, (YSize/2):]
    
    VixSlice = Vix[SliceXVal, (YSize/2):]
    ViySlice = Viy[SliceXVal, (YSize/2):]
    VizSlice = Viz[SliceXVal, (YSize/2):]
    
    VexSlice = Vex[SliceXVal, (YSize/2):]
    VeySlice = Vey[SliceXVal, (YSize/2):]
    VezSlice = Vez[SliceXVal, (YSize/2):]
    
    JxSlice = Jx[SliceXVal, (YSize/2):]
    JySlice = Jy[SliceXVal, (YSize/2):]
    JzSlice = Jz[SliceXVal, (YSize/2):]
    
    TeParSlice = TePar[SliceXVal, (YSize/2):]
    TiParSlice = TiPar[SliceXVal, (YSize/2):]
    
    TePerpSlice = TePerp1[SliceXVal, (YSize/2):]
    TiPerpSlice = TiPerp1[SliceXVal, (YSize/2):]
    
    XVals = np.linspace(dx + (YSize/2)*dx, YSize*dx, YSize/2)
    LowerBound = YSize/2
elif UpLow == 'Lower':
    BzSlice = Bz[SliceXVal, :(YSize/2)]
    BxSlice = Bx[SliceXVal, :(YSize/2)]
    BySlice = By[SliceXVal, :(YSize/2)]
    BmSlice = Bm[SliceXVal, :(YSize/2)]
    
    ExSlice = Ex[SliceXVal, :(YSize/2)]
    EySlice = Ey[SliceXVal, :(YSize/2)]
    EzSlice = Ez[SliceXVal, :(YSize/2)]
    EmSlice = Em[SliceXVal, :(YSize/2)]
    
    EdotBSlice = EdotB[SliceXVal, :(YSize/2)]
    EcrBxSlice = EcrBx[SliceXVal, :(YSize/2)]
    EcrBySlice = EcrBy[SliceXVal, :(YSize/2)]
    EcrBzSlice = EcrBz[SliceXVal, :(YSize/2)]
    EcrBmSlice = EcrBm[SliceXVal, :(YSize/2)]
    
    NiSlice = Ni[SliceXVal, :(YSize/2)]
    NeSlice = Ne[SliceXVal, :(YSize/2)]
    
    VixSlice = Vix[SliceXVal, :(YSize/2)]
    ViySlice = Viy[SliceXVal, :(YSize/2)]
    VizSlice = Viz[SliceXVal, :(YSize/2)]
    
    VexSlice = Vex[SliceXVal, :(YSize/2)]
    VeySlice = Vey[SliceXVal, :(YSize/2)]
    VezSlice = Vez[SliceXVal, :(YSize/2)]
    
    JxSlice = Jx[SliceXVal, :(YSize/2)]
    JySlice = Jy[SliceXVal, :(YSize/2)]
    JzSlice = Jz[SliceXVal, :(YSize/2)]
    
    TeParSlice = TePar[SliceXVal, :(YSize/2)]
    TiParSlice = TiPar[SliceXVal, :(YSize/2)]
    
    TePerpSlice = TePerp1[SliceXVal, :(YSize/2)]
    TiPerpSlice = TiPerp1[SliceXVal, :(YSize/2)]    
    
    XVals = np.linspace(dx, (YSize/2)*dx, YSize/2)
    LowerBound = 0    


plotSliceY = np.linspace(LowerBound, LowerBound - 1 + YSize/2, YSize/2)
plotSliceX = np.full(YSize/2, SliceXVal)


BackGround = findBackground(BzSlice)
print(' ')
print(' ')
print('Run ' + str(runNum) + ' Data:')
print('--------------------')
print('Background = ' + str(round(BackGround, 4)))
print('Max = ' + str(round(BzSlice.max(), 4)))
print('Min = ' + str(round(BzSlice.min(), 4)))

BzWidths = findWidths(BzSlice, BackGround, dx, UpLow)
print('Enhancement Width = ' + str(BzWidths[0]))
print('Suppression Width = ' + str(BzWidths[1]))
print(' ')
print(' ')

iUpp = xpts[0]
jUpp = xpts[1]
iLow = xpts[2]
jLow = xpts[3] 

MIN = Bz.min()
MAX = Bz.max()
    
paramPath = '/global/cscratch1/sd/ransom/'+ str(run) + '/staging/param_'+ str(run)
param = open(paramPath, 'r')
for line in param:
    if line.find('define') > 0 and line.find('T_e1') > 0:              
        Te = float(line.split()[2:3][0])
    if line.find('define') > 0 and line.find('T_i1') > 0:              
        Ti = float(line.split()[2:3][0])
    if line.find('define') > 0 and line.find('b0') > 0:              
        B = float(line.split()[2:3][0])
    if line.find('define') > 0 and line.find('c_2') > 0:              
        C = int(np.sqrt(float(line.split()[2:3][0])))
    if line.find('define') > 0 and line.find('ppg') > 0:              
        PPG = int(line.split()[2:3][0])
    if line.find('define') > 0 and line.find('dt') > 0:              
        dt = float(line.split()[2:3][0])
    if line.find('define') > 0 and line.find('substeps') > 0:              
        sub = int(line.split()[2:3][0])


# Initial params: Te, Ti, B
#  Energy Cons Params: C, dx, dt, sub, PPG 
#  Results: tot Cyclos, 100 Cyc dE, SliceVal, Background, min, max, Enhance, Suppress
ECon = EConPerc(run)
if UpLow == 'Upper':
    UL = 2
if UpLow == 'Lower':
    UL = 1

saveArray = np.array([int(runNum), Te, Ti, B, C, dx, dt, sub, PPG, ECon[1], round(ECon[0], 4), SliceXVal, UL, round(BackGround, 4), round(BzSlice.min(), 4), round(BzSlice.max(), 4), BzWidths[0], BzWidths[1]])
print('Saving data...')
np.save(run, saveArray)
print('Saved.')


Page1 = plt.figure(1)
Page1.set_size_inches(8.25, 9.25, forward = True)
Page1.subplots_adjust(left = .08, bottom = .02, right = .92, top = .91, wspace = .25, hspace = 2)
Page1.suptitle('$Run$' + ' ' + runNum + ':  ' + '$T_e$' + '$=$' + str(Te) + ', ' + '$T_i$' + '$=$' + str(Ti) + ', ' + '$B_z$' + '$=$' + str(B), fontsize = 18, y = .995)
Page1.patch.set_facecolor('lightgrey')

sp1 = plt.subplot2grid((15, 11), (0, 0), colspan = 5, rowspan = 2)
sp1.plot(XVals, BzSlice, color = 'g')
sp1.xaxis.set_minor_locator(AutoMinorLocator(5))
sp1.yaxis.set_minor_locator(AutoMinorLocator(5))
sp1.locator_params(nbins = 8, axis = 'x')
sp1.locator_params(nbins = 4, axis = 'y')
sp1.set_xlim([XVals[0], XVals[-1]])
sp1.set_xlabel('Y', fontsize = 10)
sp1.set_ylabel('Bz', rotation = 0, fontsize = 10)
sp1.xaxis.set_label_coords(1.03, .05)
sp1.yaxis.set_label_coords(0, 1.05)
sp1.legend(['$B_z$'], ncol = 1, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin1, ymax1 = sp1.get_ylim()
sp1.plot([BzWidths[2][0], BzWidths[2][0]],[ymin1, ymax1], '--', color = 'k', linewidth = 1)
sp1.plot([BzWidths[2][1], BzWidths[2][1]],[ymin1, ymax1], '--', color = 'k', linewidth = 1)
sp1.plot([BzWidths[2][2], BzWidths[2][2]],[ymin1, ymax1], '--', color = 'k', linewidth = 1)

sp2 = plt.subplot2grid((15, 11), (0, 6), colspan = 5, rowspan = 3)
sp2.pcolormesh(Bz.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
sp2.plot(plotSliceX, plotSliceY, color = 'k')
sp2.plot([iUpp, iLow], [jUpp, jLow], 'x', color = 'k')
subTitle = sp2.set_title('$B_z$' + ' ' + '$at$' + ' ' + str(int(ECon[1])) + ' ' + '$Cyclo$' + ' ' + '$Times$', fontsize = 11)
subTitle.set_position([.5, 1.155])
sp2.set_aspect('equal')
sp2.xaxis.set_minor_locator(AutoMinorLocator(5))
sp2.yaxis.set_minor_locator(AutoMinorLocator(5))
sp2.locator_params(nbins = 8, axis = 'x')
sp2.locator_params(nbins = 6, axis = 'y')
sp2.set_xlim([0, XSize])
sp2.set_ylim([0, YSize])
sp2.set_xlabel('X', fontsize = 10)
sp2.set_ylabel('Y', rotation = 0, fontsize = 10)
sp2.xaxis.set_label_coords(1.05, .05)
sp2.yaxis.set_label_coords(-.1, 1.05)
if UpLow == 'Upper':
    sp2.legend(['$Upper$' + ' ' + '$Slice$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(dx*SliceXVal), '$X-point$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (.5, 1.23), fontsize = 8)
if UpLow == 'Lower':
    sp2.legend(['$Lower$' + ' ' + '$Slice$' + ' ' + '$X$' + ' ' + '$=$' + ' ' + str(dx*SliceXVal), '$X-point$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (.5, 1.23), fontsize = 8)

sp3 = plt.subplot2grid((15, 11), (3, 0), colspan = 5, rowspan = 2)
sp3.plot(XVals, BxSlice)
sp3.plot(XVals, BySlice)
sp3.plot(XVals, BzSlice)
sp3.plot(XVals, BmSlice)
sp3.xaxis.set_minor_locator(AutoMinorLocator(5))
sp3.yaxis.set_minor_locator(AutoMinorLocator(5))
sp3.locator_params(nbins = 8, axis = 'x')
sp3.locator_params(nbins = 4, axis = 'y')
sp3.set_xlim([XVals[0], XVals[-1]])
sp3.set_xlabel('Y', fontsize = 10)
sp3.set_ylabel('B', rotation = 0, fontsize = 10)
sp3.xaxis.set_label_coords(1.03, .05)
sp3.yaxis.set_label_coords(0, 1.05)
sp3.legend(['$B_x$', '$B_y$', '$B_z$', '$\mid B\mid$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin3, ymax3 = sp3.get_ylim()
sp3.plot([BzWidths[2][0], BzWidths[2][0]],[ymin3, ymax3], '--', color = 'k', linewidth = 1)
sp3.plot([BzWidths[2][1], BzWidths[2][1]],[ymin3, ymax3], '--', color = 'k', linewidth = 1)
sp3.plot([BzWidths[2][2], BzWidths[2][2]],[ymin3, ymax3], '--', color = 'k', linewidth = 1)

sp4 = plt.subplot2grid((15, 11), (3, 6), colspan = 5, rowspan = 2)
sp4.plot(XVals, ExSlice)
sp4.plot(XVals, EySlice)
sp4.plot(XVals, EzSlice)
sp4.plot(XVals, EmSlice)
sp4.xaxis.set_minor_locator(AutoMinorLocator(5))
sp4.yaxis.set_minor_locator(AutoMinorLocator(5))
sp4.locator_params(nbins = 8, axis = 'x')
sp4.locator_params(nbins = 4, axis = 'y')
sp4.set_xlim([XVals[0], XVals[-1]])
sp4.set_xlabel('Y', fontsize = 10)
sp4.set_ylabel('E', rotation = 0, fontsize = 10)
sp4.xaxis.set_label_coords(1.03, .05)
sp4.yaxis.set_label_coords(0, 1.05)
sp4.legend(['$E_x$', '$E_y$', '$E_z$', '$\mid E\mid$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin4, ymax4 = sp4.get_ylim()
sp4.plot([BzWidths[2][0], BzWidths[2][0]],[ymin4, ymax4], '--', color = 'k', linewidth = 1)
sp4.plot([BzWidths[2][1], BzWidths[2][1]],[ymin4, ymax4], '--', color = 'k', linewidth = 1)
sp4.plot([BzWidths[2][2], BzWidths[2][2]],[ymin4, ymax4], '--', color = 'k', linewidth = 1)

sp5 = plt.subplot2grid((15, 11), (6, 0), colspan = 5, rowspan = 2)
sp5.plot(XVals, EcrBxSlice)
sp5.plot(XVals, EcrBySlice)
sp5.plot(XVals, EcrBzSlice)
sp5.plot(XVals, EcrBmSlice)
sp5.plot(XVals, EdotBSlice)
sp5.xaxis.set_minor_locator(AutoMinorLocator(5))
sp5.yaxis.set_minor_locator(AutoMinorLocator(5))
sp5.locator_params(nbins = 8, axis = 'x')
sp5.locator_params(nbins = 4, axis = 'y')
sp5.set_xlim([XVals[0], XVals[-1]])
sp5.set_xlabel('Y', fontsize = 10)
sp5.set_ylabel('ExB', rotation = 0, fontsize = 10)
sp5.xaxis.set_label_coords(1.03, .05)
sp5.yaxis.set_label_coords(0, 1.05)
sp5.legend([r'$(E\times B)_x$', r'$(E\times B)_y$', r'$(E\times B)_z$', r'$\mid E\times B\mid$', r'$E\cdot B$'], ncol = 3, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin5, ymax5 = sp5.get_ylim()
sp5.plot([BzWidths[2][0], BzWidths[2][0]],[ymin5, ymax5], '--', color = 'k', linewidth = 1)
sp5.plot([BzWidths[2][1], BzWidths[2][1]],[ymin5, ymax5], '--', color = 'k', linewidth = 1)
sp5.plot([BzWidths[2][2], BzWidths[2][2]],[ymin5, ymax5], '--', color = 'k', linewidth = 1)

sp6 = plt.subplot2grid((15, 11), (6, 6), colspan = 5, rowspan = 2)
sp6.plot(XVals, TeParSlice)
sp6.plot(XVals, TePerpSlice)
sp6.plot(XVals, TiParSlice)
sp6.plot(XVals, TiPerpSlice)
sp6.xaxis.set_minor_locator(AutoMinorLocator(5))
sp6.yaxis.set_minor_locator(AutoMinorLocator(5))
sp6.locator_params(nbins = 8, axis = 'x')
sp6.locator_params(nbins = 4, axis = 'y')
sp6.set_xlim([XVals[0], XVals[-1]])
sp6.set_xlabel('Y', fontsize = 10)
sp6.set_ylabel('T', rotation = 0, fontsize = 10)
sp6.xaxis.set_label_coords(1.03, .05)
sp6.yaxis.set_label_coords(0, 1.05)
sp6.legend(['$T_{epar}$', '$T_{eperp}$', '$T_{ipar}$', '$T_{iperp}$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin6, ymax6 = sp6.get_ylim()
sp6.plot([BzWidths[2][0], BzWidths[2][0]],[ymin6, ymax6], '--', color = 'k', linewidth = 1)
sp6.plot([BzWidths[2][1], BzWidths[2][1]],[ymin6, ymax6], '--', color = 'k', linewidth = 1)
sp6.plot([BzWidths[2][2], BzWidths[2][2]],[ymin6, ymax6], '--', color = 'k', linewidth = 1)

sp7 = plt.subplot2grid((15, 11), (9, 0), colspan = 5, rowspan = 2)
sp7.plot(XVals, NiSlice)
sp7.plot(XVals, NeSlice)
sp7.xaxis.set_minor_locator(AutoMinorLocator(5))
sp7.yaxis.set_minor_locator(AutoMinorLocator(5))
sp7.locator_params(nbins = 8, axis = 'x')
sp7.locator_params(nbins = 4, axis = 'y')
sp7.set_xlim([XVals[0], XVals[-1]])
sp7.set_xlabel('Y', fontsize = 10)
sp7.set_ylabel('N', rotation = 0, fontsize = 10)
sp7.xaxis.set_label_coords(1.03, .05)
sp7.yaxis.set_label_coords(0, 1.05)
sp7.legend(['$N_i$', '$N_e$'], ncol = 1, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin7, ymax7 = sp7.get_ylim()
sp7.plot([BzWidths[2][0], BzWidths[2][0]],[ymin7, ymax7], '--', color = 'k', linewidth = 1)
sp7.plot([BzWidths[2][1], BzWidths[2][1]],[ymin7, ymax7], '--', color = 'k', linewidth = 1)
sp7.plot([BzWidths[2][2], BzWidths[2][2]],[ymin7, ymax7], '--', color = 'k', linewidth = 1)

sp8 = plt.subplot2grid((15, 11), (9, 6), colspan = 5, rowspan = 2)
sp8.plot(XVals, JxSlice)
sp8.plot(XVals, JySlice)
sp8.plot(XVals, JzSlice)
sp8.xaxis.set_minor_locator(AutoMinorLocator(5))
sp8.yaxis.set_minor_locator(AutoMinorLocator(5))
sp8.locator_params(nbins = 8, axis = 'x')
sp8.locator_params(nbins = 4, axis = 'y')
sp8.set_xlim([XVals[0], XVals[-1]])
sp8.set_xlabel('Y', fontsize = 10)
sp8.set_ylabel('J', rotation = 0, fontsize = 10)
sp8.xaxis.set_label_coords(1.03, .05)
sp8.yaxis.set_label_coords(0, 1.05)
sp8.legend(['$J_x$', '$J_y$', '$J_z$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin8, ymax8 = sp8.get_ylim()
sp8.plot([BzWidths[2][0], BzWidths[2][0]],[ymin8, ymax8], '--', color = 'k', linewidth = 1)
sp8.plot([BzWidths[2][1], BzWidths[2][1]],[ymin8, ymax8], '--', color = 'k', linewidth = 1)
sp8.plot([BzWidths[2][2], BzWidths[2][2]],[ymin8, ymax8], '--', color = 'k', linewidth = 1)

sp9 = plt.subplot2grid((15, 11), (12, 0), colspan = 5, rowspan = 2)
sp9.plot(XVals, VexSlice)
sp9.plot(XVals, VeySlice)
sp9.plot(XVals, VezSlice)
sp9.xaxis.set_minor_locator(AutoMinorLocator(5))
sp9.yaxis.set_minor_locator(AutoMinorLocator(5))
sp9.locator_params(nbins = 8, axis = 'x')
sp9.locator_params(nbins = 4, axis = 'y')
sp9.set_xlim([XVals[0], XVals[-1]])
sp9.set_xlabel('Y', fontsize = 10)
sp9.set_ylabel('V', rotation = 0, fontsize = 10)
sp9.xaxis.set_label_coords(1.03, .05)
sp9.yaxis.set_label_coords(0, 1.05)
sp9.legend(['$V_{ex}$', '$V_{ey}$', '$V_{ez}$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin9, ymax9 = sp9.get_ylim()
sp9.plot([BzWidths[2][0], BzWidths[2][0]],[ymin9, ymax9], '--', color = 'k', linewidth = 1)
sp9.plot([BzWidths[2][1], BzWidths[2][1]],[ymin9, ymax9], '--', color = 'k', linewidth = 1)
sp9.plot([BzWidths[2][2], BzWidths[2][2]],[ymin9, ymax9], '--', color = 'k', linewidth = 1)

sp10 = plt.subplot2grid((15, 11), (12, 6), colspan = 5, rowspan = 2)
sp10.plot(XVals, VixSlice)
sp10.plot(XVals, ViySlice)
sp10.plot(XVals, VizSlice)
sp10.xaxis.set_minor_locator(AutoMinorLocator(5))
sp10.yaxis.set_minor_locator(AutoMinorLocator(5))
sp10.locator_params(nbins = 8, axis = 'x')
sp10.locator_params(nbins = 4, axis = 'y')
sp10.set_xlim([XVals[0], XVals[-1]])
sp10.set_xlabel('Y', fontsize = 10)
sp10.set_ylabel('V', rotation = 0, fontsize = 10)
sp10.xaxis.set_label_coords(1.03, .05)
sp10.yaxis.set_label_coords(0, 1.05)
sp10.legend(['$V_{ix}$', '$V_{iy}$', '$V_{iz}$'], ncol = 2, loc = 'upper center', bbox_to_anchor = (0.5, -0.25), fontsize = 8)

ymin10, ymax10 = sp10.get_ylim()
sp10.plot([BzWidths[2][0], BzWidths[2][0]],[ymin10, ymax10], '--', color = 'k', linewidth = 1)
sp10.plot([BzWidths[2][1], BzWidths[2][1]],[ymin10, ymax10], '--', color = 'k', linewidth = 1)
sp10.plot([BzWidths[2][2], BzWidths[2][2]],[ymin10, ymax10], '--', color = 'k', linewidth = 1)


Page1.canvas.draw()

LabelListX = sp2.get_xticklabels()
for LabelObj in LabelListX:
    try:
        if float(LabelObj.get_text()) == 0:
            LabelObj.set_text('')
        LabelObj.set_text(str(float(LabelObj.get_text())*dx))
    except:
        pass
sp2.set_xticklabels(LabelListX)

LabelListY = sp2.get_yticklabels()
for LabelObj in LabelListY:
    try: 
        if float(LabelObj.get_text()) == 0:
            LabelObj.set_text('')
        LabelObj.set_text(str(float(LabelObj.get_text())*dx))
    except:
        pass
sp2.set_yticklabels(LabelListY)

plt.show()
   
    
#fig1 = plt.figure(1)
#fig1.set_size_inches(10, 6, forward = True)
#fig1.suptitle('$Bz$', fontsize = 18, y = .99)
#plt.subplots_adjust(left = .08, bottom = .08, right = .92, top = .92) #wspace = None, hspace = .45)
#fig1.patch.set_facecolor('lightgrey')
#
#ax1 = fig1.add_subplot(111)
#ax1.pcolormesh(Bz.T, vmin = MIN, vmax = MAX, cmap=plt.cm.bwr)
#ax1.plot(plotSliceX, plotSliceY, color = 'k')
#ax1.set_title('$Slice$'+ ' ' + '$X$' + ' ' + '$=$' + ' ' + str(SliceXVal), fontsize = 14)
#ax1.set_aspect('equal')
#ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
#ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
#ax1.locator_params(nbins = 8, axis = 'x')
#ax1.locator_params(nbins = 6, axis = 'y')
#ax1.set_xlim([0, XSize])
#ax1.set_ylim([0, YSize])
#ax1.set_xlabel('X')
#ax1.set_ylabel('Y', rotation = 0)
#ax1.yaxis.set_label_position('right')
#ax1.yaxis.labelpad = 10
#
#fig2 = plt.figure(2)
#fig2.set_size_inches(10, 6, forward = True)
#fig2.suptitle('$Run$' + ' ' + runNum + '\n' + '$Bz$' + ' ' + '$Slice$'+ ' ' + '$X$' + ' ' + '$=$' + ' ' + str(SliceXVal), fontsize = 18, y = .99)
#plt.subplots_adjust(left = .1, bottom = .1, right = .9, top = .9) #wspace = None, hspace = .45)
#fig2.patch.set_facecolor('lightgrey')
#
#ax2 = fig2.add_subplot(111)
#ax2.plot(XVals, BzSlice)
#ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
#ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
#ax2.locator_params(nbins = 8, axis = 'x')
#ax2.locator_params(nbins = 6, axis = 'y')
#ax2.set_xlim([XVals[0], XVals[-1]])
#ax2.set_ylim([BzSlice.min() - .5, BzSlice.max() + .5])
#ax2.set_xlabel('Y', fontsize = 12)
#ax2.set_ylabel('Bz', rotation = 0, fontsize = 12)
#ax2.yaxis.set_label_position('right')
#ax2.yaxis.labelpad = 20

        
        
#sortBz = BzSlice[:]
#print(id(sortBz))
#print(id(BzSlice))
#BzSorted = False
#while BzSorted == False:
#    BzSorted = True
#    for i in range(0, len(sortBz) - 1):
#        if sortBz[i] > sortBz[i + 1]:
#            holder = sortBz[i + 1]
#            sortBz[i + 1] = sortBz[i]
#            sortBz[i] = holder
#            BzSorted = False
#
#median = (sortBz[len(sortBz)/2 - 1] + sortBz[len(sortBz)/2])/2
#print('Median = ' + str(median))
