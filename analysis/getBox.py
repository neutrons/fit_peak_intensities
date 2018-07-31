import sys
#remove the original mantid path

popList = []
for i in range(len(sys.path))[::-1]:
    if 'antid' in sys.path[i]:
        sys.path.pop(i)
import socket
if 'sns' in socket.gethostname():
    sys.path.append('/SNS/users/ntv/mantid/mantid/release/bin')
    #sys.path.append('/SNS/users/ntv/workspace/mantid/release/bin')
else: 
    #sys.path.append('/home/ntv/mantid/mantid/bin/')
    sys.path.append('/home/ntv/workspace/mantid/release/bin/')

import matplotlib.pyplot as plt
plt.ion()
#if '../' not in sys.path: sys.path.append('../')
import numpy as np
from scipy.optimize import curve_fit
from mantid.simpleapi import *
from mantid.kernel import V3D
import ICCFitTools as ICCFT
import BVGFitTools as BVGFT
import pickle
from timeit import default_timer as timer
reload(ICCFT)
reload(BVGFT)

print "Which peak?"
try: print 'Current peak is %i'%peakToGet
except: pass
peakToGet = int(input())
'''
#Scolecite
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
#peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/15629_Niggli.integrate'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
#UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/15629_Niggli.mat'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
dQPixel=0.006#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=2; maxdtBinWidth=15
'''

#Si 2016
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
q_frame = 'lab'
pplmin_frac=0.9; pplmax_frac=1.1;
'''
#pth
peaksFile = '/SNS/users/ntv/integrate/pth/852_pred.integrate'
UBFile = '/SNS/users/ntv/integrate/pth/combined_3.mat'
peaksFile = '/home/ntv/Desktop/pth/peaks_combined.integrate'
UBFile = '/home/ntv/Desktop/pth/UB_combined.mat'
DetCalFile = '/home/ntv/Desktop/runReduction/MaNDi2015.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
nxsTemplate = '/data/pth/MANDI_{1}_event.nxs'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([ 6.12383767,  8.8677518 , -0.02761688]) #Go with ICCFT.oldScatFun
#predpplCoefficients = np.array([[23.2736324 ,  10.10909695,   0.6229528 ]]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
pplmin_frac=0.7; pplmax_frac=1.5; mindtBinWidth=15
'''
'''
#gfp
peaksFile = '/SNS/users/ntv/integrate/gfp/combined.integrate'
UBFile = '/SNS/users/ntv/integrate/gfp/test.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = 'SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
nxsTemplate = '/SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
#predpplCoefficients = np.array([[23.2736324 ,  10.10909695,   0.6229528 ]]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
pplmin_frac=0.7; pplmax_frac=1.5; mindtBinWidth=15
'''
'''
#PsbO 2016
#peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.integrate'
peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal_highres_2p0A.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-16286/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
#predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
#pplmin_frac=0.8; pplmax_frac=1.5; mindtBinWidth=15
pplmin_frac=0.99; pplmax_frac=1.01; mindtBinWidth=15
'''

'''
#NaK 2017
peaksFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_NAK_UB.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-17495/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
pplmin_frac=0.8; pplmax_frac=1.5; mindtBinWidth=4
'''

'''
#MaNDI natrolite
peaksFile = '/SNS/users//ntv/integrate/mandi_natrolite/peaks_8679.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_natrolite/mandi_8679_niggli.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
predpplCoefficients = None##np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
'''
'''
#CORELLI - natrolite 42357 July 2017
loadDir = '/data/corelli_natrolite/'
peaksFile = '/data/corelli_natrolite/peaks_42357.integrate'
UBFile =  '/data/corelli_natrolite/peaks_42357.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
nxsTemplate = loadDir+'CORELLI_%i.nxs.h5'
dQPixel=0.007#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=1; maxdtBinWidth=60
'''

'''
#CORELLI - beryl
loadDir = '/data/corelli_beryl/IPTS-20302/'
peaksFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal_indexedonly.integrate'
UBFile =  '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
nxsTemplate = loadDir+'CORELLI_%i.nxs.h5'
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=4; maxdtBinWidth=60
'''
'''
#lpmo
peaksFile = '/SNS/users/ntv/integrate/mandi_lpmo/lpmo_combined.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_lpmo/lpmo_combined.mat'
DetCalFile = '/SNS/users/ntv/integrate/mandi_lpmo/MANDI_June2018.DetCal'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-21379/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
q_frame='lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=15; maxdtBinWidth=50;
'''
'''
#DNA
peaksFile = '/SNS/users/ntv/integrate/mandi_dna2/combined_1p5A.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_dna2/combined_1p5A.mat'
DetCalFile = '/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/data/dna/IPTS-18552/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.008#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.8; pplmax_frac=5; mindtBinWidth=25
'''

'''
#cryo
peaksFile = '/SNS/users/ntv/integrate/mandi_cryo/cryo_combined_2.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_cryo/cryo_combined_2.mat'
DetCalFile = None#'/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-19172/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients =  np.array([28.73949834,  13.04192586,   0.41210929]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.4; pplmax_frac=1.5; mindtBinWidth=15
'''

'''
#secondDNA
peaksFile = '/SNS/users/ntv/integrate/mandi_secondDNA/combined.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_secondDNA/combined.mat'
DetCalFile = None#'/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-15151/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.6; pplmax_frac=1.5; mindtBinWidth=25
'''
'''
#Beta Lac
#peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
#UBFile =  '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
#peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
#UBFile =  '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.mat'
peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.mat'

DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-15000/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.8; pplmax_frac=2.0; mindtBinWidth=15
   #---mutant
peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.integrate'
#peaksFile = '/SNS/MANDI/shared/ProfileFitting/demo_5921.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.mat'
UBFile =  '/SNS/MANDI/shared/ProfileFitting/demo_5921.mat'
nxsTemplate = '/SNS/MANDI/IPTS-8776/data/MANDI_%i_event.nxs'
predpplCoefficients = np.array([ 3.56405187,  8.34071842,  0.14134522])
#predpplCoefficients = np.array([  4.88049788,  9.29823399,  0.14255074]) #Go with ICCFT.oldScatFun
#pplmin_frac=0.4; pplmax_frac=1.5; mindtBinWidth=15
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=15
'''

'''
#MnSOD
peaksFile = '/home/ntv/mandi_preprocessing/mnsod/8298_Niggli.integrate'
UBFile =  '/home/ntv/mandi_preprocessing/mnsod/8298_Niggli.mat'

DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-17425/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.9; pplmax_frac=1.1; mindtBinWidth=15
'''

'''
#beta_lac_cryo
peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase_cryo/combined.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_beta_lactamase_cryo/8799.mat'
DetCalFile = None#'/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients =  np.array([28.73949834,  13.04192586,   0.41210929]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.0; pplmax_frac=1.0; mindtBinWidth=15
'''

# Some parameters
removeEdges = False 
importPeaks = True
for ws in mtd.getObjectNames():
    if mtd[ws].getComment() == '%s'%peaksFile:
        print '    using already loaded peaks file'
        importPeaks = False
        peaks_ws = mtd[ws]
if importPeaks:
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    peaks_ws.setComment(peaksFile)
peak = peaks_ws.getPeak(peakToGet)
LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()

importFlag = True
for ws in mtd.getObjectNames():
    if mtd[ws].getComment() == 'BSGETBOX%i'%peak.getRunNumber():
        print '   Using already loaded MDdata'
        MDdata = mtd[ws]
        importFlag = False
        break
if importFlag:
    try:
        fileName = nxsTemplate%peak.getRunNumber()
    except:
        fileName = nxsTemplate.format(0, peak.getRunNumber())
    MDdata = ICCFT.getSample(peak.getRunNumber(), DetCalFile, workDir, fileName, q_frame=q_frame)
    MDdata.setComment('BSGETBOX%i'%peak.getRunNumber())

figNumber =1 

fracHKL = 0.5
dtSpread = 0.002
dtSpreadToPlot = [0.01]
wavelength = peak.getWavelength() #in Angstrom
energy = 81.804 / wavelength**2 / 1000.0 #in eV
flightPath = peak.getL1() + peak.getL2() #in m
scatteringHalfAngle = 0.5*peak.getScattering()
dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=fracHKL))
dQ[dQ>0.25]=0.25

dQPixel = peaks_ws.getInstrument().getNumberParameter("DQPixel")[0]
Box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakToGet, dQ, fracHKL=fracHKL,dQPixel=dQPixel,  q_frame=q_frame)
box = Box
Box.setTitle('Box for peak %i'%peakToGet)
#SaveMD(InputWorkspace=Box, Filename = 'Box.nxs')

n_events = Box.getNumEventsArray()

qMask = ICCFT.getHKLMask(UBMatrix, frac=0.4, dQPixel=dQPixel, dQ=dQ)
mask = np.ones_like(n_events)

xaxis = Box.getXDimension()
qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
yaxis = Box.getYDimension()
qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
zaxis = Box.getZDimension()
qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    
QX, QY, QZ = np.meshgrid(qx, qy, qz,indexing='ij')
qSq = QX**2 + QY**2 + QZ**2
tof = 3176.507 * (peak.getL1()+peak.getL2())*np.sin(peak.getScattering()*0.5) * 1/np.sqrt(QX**2 + QY**2 + QZ**2)
tBins = np.arange(np.min(tof), np.max(tof), 10)
tofSpread = [0.95*peak.getTOF(), 1.05*peak.getTOF()]

theta = np.arctan(QY/QX)[qMask]
phi = np.arccos(QZ/np.sqrt(QX**2 + QY**2 + QZ**2))[qMask]
print '====================================****'

padeCoefficients = ICCFT.getModeratorCoefficients('/SNS/users/ntv/integrate/franz_coefficients_2017.dat')
strongPeakParams = None
instrumentName = peaks_ws.getInstrument().getFullName()
mindtBinWidth = peaks_ws.getInstrument().getNumberParameter("minDTBinWidth")[0]
maxdtBinWidth = peaks_ws.getInstrument().getNumberParameter("maxDTBinWidth")[0]
nTheta = peaks_ws.getInstrument().getIntParameter("numBinsTheta")[0]
nPhi   = peaks_ws.getInstrument().getIntParameter("numBinsPhi")[0]
iccFitDict = ICCFT.parseConstraints(peaks_ws)
Y3D, gIDX, pp_lambda, params = BVGFT.get3DPeak(peak, peaks_ws, box, padeCoefficients,qMask,nTheta=nTheta, nPhi=nPhi, plotResults=True, zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1, strongPeakParams=strongPeakParams, q_frame=q_frame, mindtBinWidth=mindtBinWidth, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac,forceCutoff=-10,edgeCutoff=0,maxdtBinWidth=maxdtBinWidth, iccFitDict=iccFitDict)

peakIDX = Y3D/Y3D.max()>0.05
plt.pause(0.01)
print 'ell: %4.4f; new: %4.4f'%(peak.getIntensity(), np.sum(Y3D[peakIDX]))
