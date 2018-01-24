import numpy as np
import ICCFitTools as ICCFT
reload(ICCFT)
import sys
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import ICConvoluted as ICC
reload(ICC)
import os 
import pickle
import ICFitLog
reload(ICFitLog)
import getEdgePixels as EdgeTools
reload(EdgeTools)
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)


# Some parameters
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
bgPolyOrder = 1 #int- bg order for IC fit
doVolumeNormalization = False #True if you want to normalize TOF profiles by volume
refineCenter = False #True if you want to determine new centers - still not very good
removeEdges = False #True if you want to not consider q-pixels that are off detector faces
calcTOFPerPixel = False #True to calculate TOF for each pixel in a qBox - uses interpolation for volNorm (if doVolumeNormalization is True)
fracHKL = 0.5 #Fraction of HKL to look on either side
fracStop = 0.01 #Fraction of max counts to include in peak selection
moderatorCoefficientsFile = 'franz_coefficients_2017.dat'
#moderatorCoefficientsFile = 'franz_coefficients_2010.dat'
calibrationDictFile = 'det_calibration/calibration_dictionary.pkl'
neigh_length_m = 3 #Will average over a (neigh_length_m)**3 box
zBG = 1.96 #z score to keep this with

'''
#Scolecite - 2016A
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = range(15629,  15644)
peaksFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/%i_Niggli.integrate'
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/%i_Niggli.mat'
UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
qLow = -15.0; qHigh = 15.0
dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 4 #Width (in us) in TOF profile bins
dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
descriptor = 'scol_dynamicQuadBG_constraints' #Does not end with '/'
dQMax = 0.4
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
'''

'''
#IPTS22331
loadDir = '/SNS/TOPAZ/IPTS-16903/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = [22331]
peaksFormat='/SNS/users/ntv/peaks_%i.integrate'
peaksFile = '/SNS/users/ntv/peaks_22331.integrate'
UBFormat = None
UBFile = '/SNS/TOPAZ/IPTS-16903/shared/Chunruo/ep_100K_MoTe2_adpQ_find/22331_Niggli.mat'
DetCalFile = '/SNS/TOPAZ/IPTS-16903/shared/calibration/TOPAZ_2016B.DetCal'
qLow = -25.0; qHigh = 25.0
dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 4 #Width (in us) in TOF profile bins
dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.4
doIterativeBackgroundFitting = False
descriptor = 'run22331' #Does not end with '/'
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
'''
'''
#NaK - MANDI
loadDir = '/SNS/MANDI/IPTS-17495/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
sampleRuns = range(8275,8282+1)
peaksFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_NAK_UB.mat'
peaksFormat = peaksFile
UBFormat = UBFile
DetCalFile = None
qLow = -5.0; qHigh = 5.0
dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 30 #Width (in us) in TOF profile bins
dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.15 #tune this
descriptor = 'nak_predpws5_lab' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
'''


#PsbO - 2016 - MANDI
loadDir = '/SNS/MANDI/IPTS-16286/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
sampleRuns = range(6154,6165+1)
peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.mat'
peaksFormat = peaksFile
UBFormat = UBFile
DetCalFile = None
qLow = -5.0; qHigh = 5.0
dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 30 #Width (in us) in TOF profile bins
dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.15 #tune this
descriptor = 'psbo_lab' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'



'''
#Beta lactamase - 2016 - MANDI
loadDir = '/SNS/MANDI/IPTS-15000/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
sampleRuns = range(4999,5003+1)
peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
#peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_3.integrate'
#UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
peaksFormat = peaksFile
UBFormat = UBFile
DetCalFile = None
qLow = -10.0; qHigh = 10.0
dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 30 #Width (in us) in TOF profile bins
dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.15 #tune this
descriptor = 'beta_lac_lab' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'
'''

'''
#Natrolite - 2016 - MANDI
loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
sampleRuns = range(8733,8750+1)
peaksFile = '/SNS/users/ntv/integrate/mandi_natrolite/mandi_natrolite_peaks.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_natrolite/mandi_natrolite_niggli_ub.mat'
peaksFormat = peaksFile
UBFormat = UBFile
DetCalFile = None
qLow = -5.0; qHigh = 5.0
dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 30 #Width (in us) in TOF profile bins
dQPixel = [0.005,0.005] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.15 #tune this
descriptor = 'beta_lac_predpws6' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
'''

'''
#Natrolite - 2017 - MANDI
loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
sampleRuns = [8679]
peaksFile='/SNS/users/ntv/integrate/mandi_natrolite/peaks_8679.integrate'
UBFile='/SNS/users/ntv/integrate/mandi_natrolite/mandi_8679_niggli.mat'
#DetCalFile = '/SNS/MANDI/shared/calibration/MANDI_500.DetCal'
DetCalFile = None
peaksFormat = peaksFile
UBFormat = UBFile
qLow = -10.0; qHigh = 10.0
dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 4 #Width (in us) in TOF profile bins
dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.2
doIterativeBackgroundFitting = False
descriptor = 'natrolite' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame = 'sample'
'''
'''
#Si - 2016A
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = range(15647,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
peaksFormat = peaksFile
UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'
UBFormat = UBFile
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
qLow = -25.0; qHigh = 25.0
dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 4 #Width (in us) in TOF profile bins
dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
dQMax = 0.4
descriptor = 'si_0p015_removeBG' #Does not end with '/'
doIterativeBackgroundFitting = False
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
'''
#==================WORK STARTS HERE==========================================
figsFormat = None# workDir + descriptor+'/figs/mantid_%i_%i.png'
'''
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
peaks_ws = LoadIsawPeaks(Filename = peaksFile)
LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
'''


#Load our peaks files and detector fitting parameters

if peaksFile is not None:
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
    UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
    dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=fracHKL))
    dQ[dQ>dQMax] = dQMax
    qMask = list()
    for dQP in dQPixel:
        print 'Getting qMask for dQPixel=%f'%dQP
        qMask.append(ICCFT.getHKLMask(UBMatrix, frac=fracHKL, dQPixel=dQP,dQ=dQ))

padeCoefficients = ICCFT.getModeratorCoefficients(moderatorCoefficientsFile)
calibrationDict = pickle.load(open(calibrationDictFile, 'rb'))

#Write the log
logFile = workDir + descriptor + '/log.log'
ICFitLog.writeLog(logFile, workDir, loadDir, nxsTemplate, figsFormat, sampleRuns, dtSpread, dtBinWidth, fracHKL, fracStop, refineCenter, removeEdges, doVolumeNormalization, peaksFormat, UBFormat, DetCalFile, moderatorCoefficientsFile, calibrationDictFile, descriptor,zBG,neigh_length_m)
for sampleRun in sampleRuns[1::2]:
    #Set up a few things for the run
    paramList = list()
    fileName = nxsTemplate%sampleRun

    #If we want to remove edges, we rebuild the panel dictionary every run
    # TODO this can be reformulated in QLab and apply R each box.
    instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
    if removeEdges:
        panelDict = EdgeTools.getInstrumentDict(instrumentFile, peaks_ws, sampleRun, fitOrder=2)
    else:
        panelDict = EdgeTools.getPanelDictionary(instrumentFile)

    #Conver the sample to reciprocal space
    MDdata = ICCFT.getSample(sampleRun, DetCalFile, workDir, fileName, qLow=qLow, qHigh=qHigh, q_frame=q_frame)
    
    #Do the actual integration
    peaks_ws,paramList,fitDict = ICCFT.integrateSample(sampleRun, MDdata, peaks_ws, paramList, panelDict, UBMatrix, dQ, qMask, padeCoefficients,parameterDict, figsFormat=figsFormat,dtBinWidth = dtBinWidth, dtSpread=dtSpread, fracHKL = fracHKL, refineCenter=refineCenter, doVolumeNormalization=doVolumeNormalization, minFracPixels=0.01, fracStop=fracStop, removeEdges=removeEdges, calibrationDict=calibrationDict,dQPixel=dQPixel, calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG, bgPolyOrder=bgPolyOrder, nBG=nBG, doIterativeBackgroundFitting=doIterativeBackgroundFitting,predCoefficients=predpplCoefficients, q_frame=q_frame)

    #Save the results and delete the leftovers
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename=workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
    np.savetxt(workDir+descriptor+'/params_%i_%s.dat'%(sampleRun, descriptor), np.array(paramList))
    pickle.dump(fitDict, open(workDir+descriptor+'/fitDict_%i_%s.pkl'%(sampleRun,descriptor),'wb'))
    
    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
