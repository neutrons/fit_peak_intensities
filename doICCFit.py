import numpy as np
import ICCFitTools as ICCFT
reload(ICCFT)
import sys
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import ICConvoluted as ICC
reload(ICC)
import os 
import sys
from parsePeaksFile import getPeaks
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)


# Some parameters
dtSpread = 0.05 #how far we look on either side of the nominal peak
gridBox = 201 #Number of points for peak MDBoxes
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
'''
#Scolecite - 2016A
sampleRuns = range(15629,  15644)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
crystalSystem = 'monoclinic'
latticeConstants = [6.5175,18.9722,9.7936,90.0000,108.9985,90.0000]
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
descriptor = 'scolecite_0p04' #Does not end with '/'
'''
#Si - 2016A
sampleRuns = range(15647,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'
crystalSystem ='cubic'
latticeConstants = [5.43071] #Since it's cubic, this we only need a (in angstrom)
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
descriptor = 'si_detCal_0p05' #Does not end with '/'




if os.path.isdir(workDir + descriptor):
    inp = raw_input('!!!!!! WARNING: PATH %s ALREADY EXIST!!!!! CONTINUE? (Y/<n>):'%(workDir + descriptor))
    print inp
    if inp is 'y' or inp is 'Y':
        print 'Will overwrite files in %s...'%workDir + descriptor
    else:
        print 'exiting...'
        sys.exit()
else:
    os.mkdir(workDir + descriptor)
    os.mkdir(workDir + descriptor + '/figs/')

figsFormat = workDir + descriptor+'/figs/mantid_%i_%i.png'
peaks_ws = LoadIsawPeaks(Filename = peaksFile)

#Get detector bank for each peak
peakList = getPeaks(peaksFile)
detBankList = np.zeros(len(peakList))
for i, peak in enumerate(peakList):
    detBankList[i] = peak.detnum


LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
for sampleRun in sampleRuns:
    paramList = list()
    MDdata = ICCFT.getSample(sampleRun, UBFile, DetCalFile, workDir, loadDir)

    peaks_ws,paramList= ICCFT.integrateSample(sampleRun, MDdata, latticeConstants,crystalSystem, gridBox, peaks_ws,paramList,detBankList, UBMatrix, figsFormat=figsFormat,dtSpread=dtSpread)
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename=workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
    np.savetxt(workDir+descriptor+'/params_%i_%s.dat'%(sampleRun, descriptor), np.array(paramList))

    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
