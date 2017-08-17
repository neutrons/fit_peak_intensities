import numpy as np
import ICCFitTools as ICCFT
reload(ICCFT)
import sys
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *


gridBox = 201 #Number of points for peak MDBoxes
'''
#Scolecite - 2016A
sampleRuns = range(15629,  15644)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
'''
#Si - 2016A
sampleRuns = range(15647,15670)
#sampleRuns = range(15666,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
#peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.integrate'
#peaksFile = '/SNS/users/vel/Dropbox (ORNL)/first62.peaks'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.mat'
crystalSystem ='cubic'
latticeConstants = [5.43071] #Since it's cubic, this we only need a (in angstrom)

DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
figsFormat = 'testfigs/mantid_%i_%i.png'
peaks_ws = LoadIsawPeaks(Filename = peaksFile)
for sampleRun in sampleRuns:
    paramList = list()
    MDdata = ICCFT.getSample(sampleRun, UBFile, DetCalFile, workDir, loadDir)
    peaks_ws,paramList= ICCFT.integrateSample(sampleRun, MDdata, latticeConstants,crystalSystem, gridBox, peaks_ws,paramList,figsFormat=figsFormat)
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename='peaks_%i_dynamic_binning.integrate'%(sampleRun))
    np.savetxt('params_%i_dynamic_binning.dat'%sampleRun, np.array(paramList))
    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
