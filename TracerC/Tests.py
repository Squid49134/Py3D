# FIELD TRACER TESTS
# William Ransom, University of Delaware 2017

# To run these tests, place TracerWrapper.py and TracerFunctions.c in the
#  same directory as this file, then use:
#  gcc -o TracerFunctions.so -shared -fPIC TracerFunctions.c
#  to compile TracerFunctions.so, the shared object file.
#  This test file can then be ran using python or ipython
# These test can be run remotely but only cover the TraceField() method in 2D and 3D.
#  Testing MapSeparator(), SeparatorSlice() and SeparatorLoader() requires access
#  to the large data files from a real reconnection simulation, and therefore can
#  only be done locally, but are provided below in the commented out section.

import numpy as np
import matplotlib.pyplot as plt
from TracerWrapper import *

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 2D Test Field

print('Creating 2D Field...')
# dx of the test simulation
dx = .05
# creating array corresponding to number of grid spaces in test simulation
# from dx/2 to (length of simulation)-dx/2 each dx in length
x = np.arange(2048)*dx + dx/2.
# length of the test simulation
lx = x[-1] + x[0]

# creating 2D grid of x and y values
yy,xx = np.meshgrid(x,x)
# creating 2D array containing values equal to the length of the space from
# the center of the grid (forms concentric circles)
rr = np.sqrt((xx - lx/2.)**2 + (yy - lx/2.)**2)

# setting magnitude of magnetic field
psi = np.exp(-rr**1/4.)

# creating circular magnetic field components - picture shifting 2 circles
# veritcally to create X component such that it is greatest at the top and 
# bottom of circle and zero on either side, and shifting 2 circles horizontally
# to create Y component such that it is greatest at left and right of circle 
# and zero at top and bottom
Bx = (np.roll(psi,-1,axis=1) - np.roll(psi,1,axis=1))/(2.*dx)
By = -(np.roll(psi,-1,axis=0) - np.roll(psi,1,axis=0))/(2.*dx)

# 2D fields must be float32 type
Bx = np.float32(Bx)
By = np.float32(By)

# checking DivB = dBx/dx + dBy/dy = 0
DivX = (np.roll(Bx, -1, axis=0) - np.roll(Bx, 1, axis=0))/(2.*dx)
DivY = (np.roll(By, -1, axis=1) - np.roll(By, 1, axis=1))/(2.*dx)

DivB = DivX + DivY

DivB_Tot = np.sum(DivB)

print('DivB = ' + str(DivB_Tot))

# 2D trace call
TraceField(.05, [Bx, By], [32, 32])

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# 3D Test Fields

print('\n' + 'Creating 3D Field...')
dx = .05
x = np.arange(256)*dx + dx/2.
lx = x[-1] + x[0]

yy,xx,zz = np.meshgrid(x,x,x)
rr = np.sqrt((xx - lx/2.)**2 + (yy - lx/2.)**2)

psi = np.exp(-rr**1/4.)

Bx = (np.roll(psi,-1,axis=1) - np.roll(psi,1,axis=1))/(2.*dx)
By = -(np.roll(psi,-1,axis=0) - np.roll(psi,1,axis=0))/(2.*dx)
Bm2 = np.sqrt(Bx**2 + By**2)
# constant magnetic field in the Z direction (spiral magneic field)
Bz = Bm2

DivX = (np.roll(Bx, -1, axis=0) - np.roll(Bx, 1, axis=0))/(2.*dx)
DivY = (np.roll(By, -1, axis=1) - np.roll(By, 1, axis=1))/(2.*dx)
DivZ = (np.roll(Bz, -1, axis=2) - np.roll(Bz, 1, axis=2))/(2.*dx)

DivB = DivX + DivY + DivZ

DivB_Tot = np.sum(DivB)

print('DivB = ' + str(DivB_Tot))

# 3D trace call
TraceField(.05, [Bx, By, Bz], [8, 8, 0], .001, 12.5)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# LOCAL DATA TESTS
    
#import sys
#sys.path.append('/home/ransom/Py3D')
#from py3d.sub import *
#from TracerWrapper import *
#
## 2D
#
#d = load_movie( 6, 'param_turb8192r1', '/scratch-fast/ransom/turb_data', ['bx', 'by'], 0)   
#Bx = d['by']
#By = d['bx']
#
## 3D
#
#print('Loading')
#BX =  np.load('/scratch-fast/asym030/bx.npy')
#BY =  np.load('/scratch-fast/asym030/by.npy')
#BZ =  np.load('/scratch-fast/asym030/bz.npy')
#
##EX =  np.load('/scratch-fast/asym030/ex.npy')
##EY =  np.load('/scratch-fast/asym030/ey.npy')
##EZ =  np.load('/scratch-fast/asym030/ez.npy')
#print('Loaded')
#
##TraceField(.05, [Bx, By], [50, 50])
##
##SepY0 = np.load('LowerSepY.npy')[0]
##raw_input("Continue?")
##TraceField(.025, [BX, BY, BZ], [0, SepY0, 0], .075, 1.5)
##raw_input("Continue?")
##TraceField(.025, [BX, BY, BZ], [12.5, 5, 12.5], .0125, 100, ['TraceX6.npy', 'TraceY6.npy', 'TraceZ6.npy'])
##raw_input("Continue?")
##TraceField(.025, [BX, BY, BZ, EX, EY, EZ], [25.6, 12.8, 12.8], .0125, 100)
##raw_input("Continue?")
#TraceField(.025, [BX, BY, BZ], [25.6, 12.8, 12.8], .0025, 100)
#
##raw_input("Continue?")
##SeparatorSlice(.025, 'LowerSepX.npy', 'LowerSepY.npy', 'LowerSepZ.npy', 51.2, 25.6, 25.6, [BX, BY, BZ])
##raw_input("Continue?")
##SeparatorSlice(.025, 'UpperSepX.npy', 'UpperSepY.npy', 'UpperSepZ.npy', 51.2, 25.6, 25.6)
##
##
### 318*.025 for Upper
### 150*.025 for Lower
##raw_input("Continue?")
##MapSeparator(.025, ['SepX11.npy', 'SepY11.npy', 'SepZ11.npy'], [BX, BY, BZ], 8, 'Upper', 16)
##raw_input("Continue?")
##MapSeparator(.025, ['SepX9.npy', 'SepY9.npy', 'SepZ9.npy'], [BX, BY, BZ], 3.75, 'Lower', 8)
##
##raw_input("Continue?")
##SeparatorLoader(.025, 'LowerSepX.npy','LowerSepY.npy','LowerSepZ.npy', 51.2, 25.6, 25.6)
##raw_input("Continue?")
##SeparatorLoader(.025, 'LowerSepX.npy','LowerSepY.npy','LowerSepZ.npy', 51.2, 25.6, 25.6, [BX, BY, BZ])
