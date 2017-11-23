import numpy as np
import math
import numpy as np
import math

Runs = [[0.25, 0.25, 1], [0.25, 0.25, 2], [0.25, 0.25, 4], [0.25, 1, 1], [0.25, 1, 4], [0.25, 4, 1], [0.25, 4, 2], [0.25, 4, 4], [1, 0.25, 1], [1, 0.25, 4], [1, 4, 1], [1, 4, 4], [4, 0.25, 1], [4, 0.25, 2], [4, 0.25, 4], [4, 1, 1], [4, 1, 4], [4, 4, 1], [4, 4, 2], [4, 4, 4]]

def parameters(Te, Ti, Bz, i):
    N = 1
    Me = .04
    Bx = 1
    By = 0
    
    # |B|
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    
    # Velocities
    velTherm = np.sqrt(Te/Me)
    velAlfven = B/np.sqrt(N*Me)
    if Te == .25:
        C = 10
    if Te == 1:
        C = 20
    if Te == 4:
        C = 40
    #C = int(max(velTherm, velAlfven)) + 1
    
    # Lengths
    lenDebye = np.sqrt(Te/N)/C
    lenGyro = np.sqrt(2*Me*Te)/B
    lenInert = np.sqrt(Me/N)
    
    # Times
    timCyclo = Me*2*np.pi/B
    timPlasma = 2*np.pi*np.sqrt(Me/N)/C
    
    # DeltaX
    l = min(lenDebye, lenGyro, lenInert)
    DeltaX = .5*(l - (l % .025))
    
    
    # DeltaT
    DeltaT = .5*DeltaX/velTherm
    if DeltaT > (timPlasma):
        print('dt adjusted for Plasma Time')
        DeltaT = (timPlasma*.9)
    if DeltaT > timCyclo:
        print('dt adjusted for Cyclo Time')
        DeltaT = (timCyclo*.9)
    DeltaT = round(DeltaT, -int(math.floor(math.log10(abs(DeltaT))) - 2))
    
    # Substeps
    substeps = int(C*DeltaT/(2*DeltaX))
    if substeps < 1:
        substeps = 1
    
    if i == 1:
        PPG = 100
    if i == 2:
        PPG = 100
    if i == 3:
        PPG = 100
    if i == 4:
        PPG = 100
    if i == 5:
        PPG = 100
    if i == 6:
        PPG = 100
    if i == 7:
        PPG = 100
    if i == 8:
        PPG = 100
    if i == 9:
        PPG = 300
    if i == 10:
        PPG = 100
    if i == 11:
        PPG = 150
    if i == 12:
        PPG = 100
    if i == 13:
        PPG = 900
    if i == 14:
        PPG = 650
    if i == 15:
        PPG = 400
    if i == 16:
        PPG = 800
    if i == 17:
        PPG = 400
    if i == 18:
        PPG = 600
    if i == 19:
        PPG = 500
    if i == 20:
        PPG = 350
    
    return DeltaX, DeltaT, substeps, PPG, C
    

params = np.zeros((20, 5))
TotHours = 0
for i in range(0, 20):
    #print(i + 1)
    params[i] = parameters(Runs[i][0], Runs[i][1], Runs[i][2], i + 1)
    hours = (100/params[i][1])*(.00002*params[i][3]/(params[i][0]**2))/4
    TotHours = TotHours + hours
print('Total core hours = ' + str(TotHours))



run = 'R'
lsimX = 102.4
lsimY = 51.2
nex = 64
ney = 64

while True:
    i = int(raw_input("Which run 1 - 20 with no R? \n"))
    i = i - 1

    dx = params[i][0]
    
    GroupNum = 1
    C = params[i][4]
    dt = params[i][1]
    sub = params[i][2]
    ppg = params[i][3]
    Te = Runs[i][0]
    Ti = Runs[i][1]
    B = Runs[i][2]
       
    try:
        existingFile = open('param_' + run + str(i+1 +(20*(GroupNum - 1))), 'r')
        print('File name already exists, aborting...')
    except:    
        file = open('param_' + run + str(i+1 +(20*(GroupNum - 1))), 'w')
        file.write('! Precision of Real Numbers\n')
        file.write('! default:\n')
        file.write('#define drk 8\n')
        file.write('! Particles:\n')
        file.write('#define prk 4\n')
        file.write('\n')
        file.write('! box dimensions\n')
        file.write('! processors in each space directions\n')
        file.write('#define pex ' + str(int(round(lsimX/(dx*nex)))) + '\n')
        file.write('#define pey ' + str(int(round(lsimY/(dx*ney)))) + '\n')
        file.write('#define pez 1\n')
        file.write('! grid points per processor (power of two in each direction highly recommended)\n')
        file.write('! nz must be power of two, nx and should be powers of two for multigrid\n')
        file.write('! if nz=1 set pez=1 (2D case)\n')
        file.write('#define nx ' + str(int(nex)) + '\n')
        file.write('#define ny ' + str(int(ney)) + '\n')
        file.write('#define nz 1\n')
        file.write('! box dimensions\n')
        file.write('!#define lx 25.6\n')
        file.write('!#define ly 12.8\n')
        file.write('!#define lz 12.8\n')
        file.write('\n')
        file.write('#define lx ' + str(float(lsimX)) + '\n')
        file.write('#define ly ' + str(float(lsimY)) + '\n')
        file.write('#define lz ' + str(float(1.0)) + '\n')
        file.write('\n')
        file.write('! physical and numerical parameters\n')
        file.write('#define c_2 ' + str(float(C**2)) + '\n')
        file.write('#define m_e 0.04\n')
        file.write('#define dt ' + str(float(dt)) + '\n')
        file.write('#define substeps ' + str(int(sub)) + '\n')
        file.write('#define boundary_condition periodic\n')
        file.write('#define diagout 0.2\n')
        file.write('#define subtract_average_rho\n')
        file.write('#define relativistic\n')
        file.write('#define smooth_e_j_rho\n')
        file.write('#define steps_per_energy_out 1\n')
        file.write('##define steps_per_partsort 10\n')
        file.write('\n')
        file.write('! movie diagnostics\n')
        file.write('#define movie_header "movie_pic3.0.h"\n')
        file.write('!#define double_byte\n')
        file.write('!#define movieout 0.99999999999\n')
        file.write('#define n_movieout ' + str(int(1/dt)) + '\n')
        file.write('\n')
        file.write('! other definitions\n')
        file.write('#define head \'p3d\'\n')
        file.write('#define bufsize 10000\n')
        file.write('#define maxparticles 18000000\n')
        file.write('#define nchannels pex\n')
        file.write('#define maxruntime 6000\n')
        file.write('#define USE_IO_V2\n')
        file.write('\n')
        file.write('! initialization\n')
        file.write('#define init_scheme initasymreconn\n')
        file.write('#define ppg ' + str(int(ppg)) + '\n')
        file.write('! Dont change these\n')
        file.write('#define b1 1.0\n')
        file.write('#define b2 1.0\n')
        file.write('#define n1 1.0\n')
        file.write('#define n2 1.0\n')
        file.write('! Do change these\n')
        file.write('#define T_e1 ' + str(float(Te)) + '\n')
        file.write('#define T_i1 ' + str(float(Ti)) + '\n')
        file.write('#define T_i2 ' + str(float(Ti)) + '\n')
        file.write('#define b0 ' + str(float(B)) + '\n')
        file.write('\n')
        file.write('#define w0 0.6\n')
        file.write('#define psi0 0.15\n')
        file.write('!#define alt_initial_current\n')
        file.write('!#define smooth_sheet\n')
        file.write('!#define n_frac 0.\n')
        file.write('! Default value for randseedadd is 54321\n')
        file.write('#define RANDSEEDADD 54321\n')
        file.write('\n')
        file.write('! multigrid parameters\n')
        file.write('#define nu1 4\n')
        file.write('#define nu2 4\n')
        file.write('#define nu3 20\n')
        file.write('#define maxiter 2000\n')
        file.write('#define eps 1e-8\n')
        file.write('#define norm_res\n')
        file.write('#define mgverbose\n')
        file.write('!#define mgverbose2\n')
        file.write('#define alt_restrict_prolong\n')
        file.write('#define skip_poisson 1\n')
        file.write('#define coarsegrid_sor\n')
        file.write('#define miniter 4\n')
        file.write('#define JACOBI\n')
        file.write('!#define GS_RB\n')
        file.write('!#define WCYCLE\n')
        file.write('\n')
        file.write('! domain for post processing\n')
        file.write('#define pexv 1\n')
        file.write('#define peyv 1\n')
        file.write('#define pezv 1\n')
        file.write('#define nxv 256\n')
        file.write('#define nyv 256\n')
        file.write('#define nzv 1\n')
        file.write('\n')
        
    cont = raw_input("make another param file? Y/N \n")
    if cont == 'Y' or cont == 'y':
        continue
    else:
        break




