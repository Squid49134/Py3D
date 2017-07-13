#-----------------------------------------------------------------------------#
# 2D and 3D magnetic field line tracing program

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from types import StringType
from Py3D.sub import load_movie




def TraceField(X = None, Y = None, Z = None, Bx = None, By = None, Bz = None, ds = .1, passes = 10):
    if X == None:
        while True:        
            try:
                X = abs(int(raw_input('\n Enter X value of starting position: \n')))
                break
            except:
                print('invalid input try again \n')
                continue
    if Y == None:
        while True:        
            try:
                Y = abs(int(raw_input('Enter Y value of starting position: \n')))
                break
            except:
                print('invalid input try again \n')
                continue
    if Z == None:
        while True:        
            try:
                Z = abs(int(raw_input('Enter Z value of starting position: \n')))
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
            MaxSteps = passes * (Bx.shape[0] / ds)
        elif Bx.shape[1] >= Bx.shape[0] and Bx.shape[1] >= Bx.shape[2]:
            MaxSteps = passes * (Bx.shape[1] / ds)
        else:
            MaxSteps = passes * (Bx.shape[2] / ds)
    elif len(Bx.shape) == 2:
        Dimension = 2
        if Bx.shape[0] >= Bx.shape[1]:
            MaxSteps = passes * (Bx.shape[0] / ds)
        else:
            MaxSteps = passes * (Bx.shape[1] / ds)
            
    if Dimension == 3:
        3DTrace
        Generate Puncture Plots?
    elif Dimension == 2:
        2DTrace
        Generate Multiline plot?
    else:
        print('invalid data file, aborting...')
    
TraceField()
    
