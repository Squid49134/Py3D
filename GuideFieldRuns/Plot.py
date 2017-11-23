import numpy as np
import matplotlib as mpt
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import sys
sys.path.append('/global/homes/r/ransom/Py3D')
from py3d.sub import *

def EPlot():
    run = raw_input('Which run? \n')
    
    for i in range(1, 10):
            try:
                test = open('/global/cscratch1/sd/ransom/'+ str(run) + '/staging/p3d.stdout.00' + str(i), 'r')
                restartNum = i
            except:
                pass
            
    print('Sim has been resarted ' + str(restartNum) + ' times,')
    
    num = raw_input('Which restart would you like to plot from? \n')
    
    if restarts == 0 or num == restarts:
        movies = '/global/cscratch1/sd/ransom/' + str(run) + '/staging'
        paramFile = 'param_'+ str(run)   
        
        d = load_movie(0, paramFile, movies, ['bz'])   
    else:
        movies = '/global/cscratch1/sd/ransom/' + str(run) + '/staging_0' + str(num) + '/staging'
        paramFile = 'param_'+ str(run)   
        
        d = load_movie(0, paramFile, movies, ['bz'])
    
    ims(d, 'bz')

#    Energies = get_energy()
#    file = open(param, 'r')
#    for line in file:
#        if line.find('define') > 0 and line.find('dt') > 0:              
#            dt = float(line.split()[2:3][0])
#    print('Avg Ki + Ke / K = ' + str(np.average((Energies[3] + Energies[4])/Energies[2] )))
#    times = np.linspace(0, (Energies[0].shape[0] - 1)*dt, Energies[0].shape[0])
#    mins = [(100*(Energies[0] - Energies[0][0])/Energies[0][0]).min(), (100*(Energies[1] - Energies[1][0])/Energies[1][0]).min(), (100*(Energies[2] - Energies[2][0])/Energies[2][0]).min(), (100*(Energies[3] - Energies[3][0])/Energies[3][0]).min(), (100*(Energies[4] - Energies[4][0])/Energies[4][0]).min()]
#    maxs = [(100*(Energies[0] - Energies[0][0])/Energies[0][0]).max(), (100*(Energies[1] - Energies[1][0])/Energies[1][0]).max(), (100*(Energies[2] - Energies[2][0])/Energies[2][0]).max(), (100*(Energies[3] - Energies[3][0])/Energies[3][0]).max(), (100*(Energies[4] - Energies[4][0])/Energies[4][0]).max()]
#    mins = np.array(mins)
#    maxs = np.array(maxs)
#    minTot = mins.min()
#    maxTot = maxs.max()
#    
#    mn = minTot - .1*(maxTot - minTot)
#    mx = maxTot + .1*(maxTot - minTot) 
#    
#    print('E Change at 100 Cyclo Times = ' + str((100*(Energies[0][Energies[0].shape[0] - 1] - Energies[0][0])/Energies[0][0])*(100/(dt*Energies[0].shape[0]))))
#    
#    fig1 = plt.figure(1)
#    fig1.set_size_inches(12, 5, forward = True)
#    ax1 = fig1.add_subplot(111)
#    ax1.plot(times, Energies[0], label = 'E Tot')
#    ax1.plot(times, Energies[1], label = 'E EM')
#    ax1.plot(times, Energies[2], label = 'E Kin')
#    ax1.legend(loc='upper right', fontsize = 'large')
#    ax1.set_xlim([0, Energies[0].shape[0]])
#    ax1.set_ylim([mins.min() - abs(mins[0] / 4), maxs.max() + abs(mins[0] / 4)])
#    ax1.set_title('$Energy$' + ' ' + '$Conservation$', fontsize = 20)
#    ax1.set_ylabel('$Energy$', fontsize = 16)
#    ax1.set_xlabel('$Movie$' + ' ' + '$Time$', fontsize = 16)
#    ax1.yaxis.set_label_position('right')
#    
#    fig2 = plt.figure(2)
#    fig2.set_size_inches(12, 5, forward = True)
#    ax2 = fig2.add_subplot(111)
#    ax2.plot(times, 100*(Energies[0] - Energies[0][0])/Energies[0][0], 'k', label = 'E Tot', )
#    ax2.plot(times, 100*(Energies[1] - Energies[1][0])/Energies[1][0], 'r', label = 'E EM')
#    ax2.plot(times, 100*(Energies[3] - Energies[3][0])/Energies[3][0], 'b', label = 'E Kin i')
#    ax2.plot(times, 100*(Energies[4] - Energies[4][0])/Energies[4][0], 'g', label = 'E Kin e')
#    ax2.legend(loc='upper right', fontsize = 'large')
#    ax2.set_xlim([0, times[times.shape[0] - 1]])
#    ax2.set_ylim([mn, mx])
#    ax2.set_title('$Energy$' + ' ' + '$Conservation$', fontsize = 20)
#    ax2.set_ylabel('$Energy$' + ' ' + '$Percent$' + ' ' + '$Change$', fontsize = 16)
#    ax2.set_xlabel('$Cyclotron$' + ' ' + '$Times$', fontsize = 16)
#    ax2.yaxis.set_label_position('right')
#    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))    
    
    plt.show()
    
EPlot()
