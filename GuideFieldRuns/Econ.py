import numpy as np
import os
import sys
sys.path.append('/global/homes/r/ransom/Py3D')
from py3d.sub import *

def get_energy(stdoutPath = None):
    
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

def ECons(Run = None):
    if Run == None:
        Run = str(raw_input('Which Run? \n'))
    num = 0
    for i in range(1, 10):
            try:
                test = open('/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + str(i), 'r')
                num = i
            except:
                pass
    if num == 0: 
        pathStdout = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + num
        pathToMovieData = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging'
        #pathExp = '/global/cscratch1/sd/ransom/' + str(Run) + '/exp3d_' + str(Run)
        paramFile = 'param_'+ str(Run)
        paramPath = pathToMovieData + '/' + paramFile
        
        param = open(paramPath, 'r')
        for line in param:
            if line.find('define') > 0 and line.find('dt') > 0:              
                dt = float(line.split()[2:3][0])        
        
        
        Energies = get_energy(pathStdout)
        print('Total Cyclo Times = ' + str(len(Energies[0])*dt))
        
        return (100*(Energies[0][-1] - Energies[0][0])/Energies[0][0])*(100/(dt*Energies[0].shape[0]))
    else:
        pathStdoutRestart = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging/p3d.stdout.00' + str(num)
        pathStdoutOriginal = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging_00/staging/p3d.stdout.000'
        EnergiesLength = 0;
        for i in range(0,int(num)):
            pathStdout = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging_0' + str(i) + '/staging/p3d.stdout.00' + str(i)
            Energy = get_energy(pathStdout)
            EnergiesLength = EnergiesLength + len(Energy[0])
        
        pathToMovieData = '/global/cscratch1/sd/ransom/'+ str(Run) + '/staging'
        #pathExp = '/global/cscratch1/sd/ransom/' + str(Run) + '/exp3d_' + str(Run)
        paramFile = 'param_'+ str(Run)
        paramPath = pathToMovieData + '/' + paramFile
        
        param = open(paramPath, 'r')
        for line in param:
            if line.find('define') > 0 and line.find('dt') > 0:              
                dt = float(line.split()[2:3][0])        
        
        
        EnergiesRestart = get_energy(pathStdoutRestart)
        EnergiesOriginal = get_energy(pathStdoutOriginal)
        
        EnergiesLength = EnergiesLength + len(EnergiesRestart[0])
        print('Total Cyclo Times = ' + str(EnergiesLength*dt))        
        
        return (100*(EnergiesRestart[0][-1] - EnergiesOriginal[0][0])/EnergiesOriginal[0][0])*(100/(dt*EnergiesLength))
    
print('100 Cyclo dEs')
print('-------------------')
for i in range(1, 13):
    print('R' + str(i))
    print(str(ECons('R' + str(i))))
    print(' ')
