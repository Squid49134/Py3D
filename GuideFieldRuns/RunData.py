import numpy as np

dataFile = open('RunData', 'w')
dataFile.write('Guide Field Run Data \n')
dataFile.write('-------------------- \n')
dataFile.write('-------------------- \n')
dataFile.write('\n')

for i in range(1, 21):
    try:
        # 0     , 1 , 2 , 3, 4, 5 , 6 , 7  , 8  , 9      , 10, 11    , 12, 13  , 14 , 15 , 16     , 17
        # runNum, Te, Ti, B, C, dx, dt, sub, PPG, totCycl, dE, SliceX, UL, Back, min, max, Enhance, suppress
        data = np.load('R' + str(i) + '.npy')
        dataFile.write('Run ' + str(int(data[0])) + ':\n')
        dataFile.write('------- \n')
        dataFile.write('Te = ' + str(data[1]) + ', Ti = ' + str(data[2]) + ', B = ' + str(data[3]) + '\n')
        dataFile.write('  Initial Parameters: \n')
        dataFile.write('     C = ' + str(int(data[4])) + ', dx = ' + str(data[5]) + ', dt = ' + str(data[6]) + ', sub = ' + str(int(data[7])) + ', PPG = ' + str(int(data[8])) + '\n')
        dataFile.write('  Energy Conservation / General Data: \n')
        dataFile.write('     100 Cyclo dE % = ' + str(data[10]) + '\n')
        dataFile.write('     Data Cyclo Time = ' + str(int(data[9])) + '\n')
        if data[12] == 2:
            dataFile.write('     Slice X = ' + str(int(data[11])) + ' Upper \n')
        if data[12] == 1:
            dataFile.write('     Slice X = ' + str(int(data[11])) + ' Lower \n')
        dataFile.write('  Bz: \n')
        dataFile.write('     Background = ' + str(data[13]) + '\n')
        dataFile.write('     Min = ' + str(data[14]) + '\n')
        dataFile.write('     Max = ' + str(data[15]) + '\n')
        dataFile.write('     Enhance Width = ' + str(data[16]) + '\n')
        dataFile.write('     Suppress Width = ' + str(data[17]) + '\n')
        dataFile.write('\n')
        
    except:
        pass
    
dataFile.close()
