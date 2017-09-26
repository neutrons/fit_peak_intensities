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
import ICFitLog
reload(ICFitLog)
import getEdgePixels as EdgeTools
reload(EdgeTools)
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)


# Some parameters
dtSpread = 0.03 #how far we look on either side of the nominal peak
dtBinWidth = 4
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
doVolumeNormalization = True #True if you want to normalize TOF profiles by volume
refineCenter = False
removeEdges = True 
fracHKL = 0.8 #Fraction of HKL to look on either side
fracStop = 0.01 #Fraction of max counts to include in peak selection

'''
#Scolecite - 2016A
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = range(15629,  15644)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/%i_Niggli.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
descriptor = 'scolecite_removeEdges_0p8hkl' #Does not end with '/'
'''

'''
#Natrolite - 2016 - MANDI
loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
sampleRuns = [8401]
peaksFile=None#'/SNS/MANDI/IPTS-8776/shared/Natrolite/New/8041_Niggli.integrate'
DetCalFile = '/SNS/MANDI/shared/calibration/MANDI_500.DetCal'
descriptor = 'natrolite' #Does not end with '/'
'''

#Si - 2016A
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = range(15647,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'

UBFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/%i_Niggli.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
descriptor = 'si_newUB_0p8hkl' #Does not end with '/'


figsFormat = workDir + descriptor+'/figs/mantid_%i_%i.png'

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

#Load our peaks files
if peaksFile is not None:
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)

#Write the log
logFile = workDir + descriptor + '/log.log'
ICFitLog.writeLog(logFile, workDir, loadDir, nxsTemplate, figsFormat, sampleRuns, dtSpread, dtBinWidth, fracHKL, fracStop, refineCenter, removeEdges, doVolumeNormalization, peaksFile, UBFormat, DetCalFile, descriptor)

for sampleRun in sampleRuns:
    
    #Set up a few things for the run
    paramList = list()
    fileName = nxsTemplate%sampleRun

    #If we want to remove edges, we rebuild the panel dictionary every run
    # TODO this can be reformulated in QLab and apply R each box.
    if removeEdges:
        instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
        panelDict = EdgeTools.getInstrumentDict(instrumentFile, peaks_ws, sampleRun, fitOrder=2)
    else:
        panelDict = EdgeTools.getPanelDictionary(instrumentFile)

    #Conver the sample to reciprocal space
    MDdata = ICCFT.getSample(sampleRun, DetCalFile, workDir, fileName)
    
    #Load the new UB and find peaks in this run if we need to.

    if peaksFile is None:
        peaks_ws = FindPeaksMD(InputWorkspace='MDdata', PeakDistanceThreshold=1.1304, MaxPeaks=1000, DensityThresholdFactor=30, OutputWorkspace='peaks_ws')
        LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFormat%sampleRun)
        UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
    else:
        LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFormat%sampleRun)
        UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()

    #Do the actual integration
    peaks_ws,paramList= ICCFT.integrateSample(sampleRun, MDdata, peaks_ws, paramList, panelDict, UBMatrix, figsFormat=figsFormat,dtBinWidth = dtBinWidth, dtSpread=dtSpread, fracHKL = fracHKL, refineCenter=refineCenter, doVolumeNormalization=doVolumeNormalization, minFracPixels=0.0075, fracStop=fracStop, removeEdges=removeEdges)

    #Save the results and delete the leftovers
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename=workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
    np.savetxt(workDir+descriptor+'/params_%i_%s.dat'%(sampleRun, descriptor), np.array(paramList))
    
    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
