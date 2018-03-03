import sys
#remove the original mantid path
popList = []
for i in range(len(sys.path))[::-1]:
    if 'antid' in sys.path[i]:
        sys.path.pop(i)
import socket
if 'sns' in socket.gethostname():
    sys.path.append('/SNS/users/ntv/mantid/mantid/release/bin')
else:
    sys.path.append('/home/ntv/mantid/mantid/bin/')


import matplotlib.pyplot as plt
plt.ion()
if '../' not in sys.path: sys.path.append('../')
import numpy as np
from scipy.optimize import curve_fit
from mantid.simpleapi import *
from mantid.kernel import V3D
import ICCFitTools as ICCFT
import BVGFitTools as BVGFT
import pickle
import getEdgePixels as EdgeTools
from timeit import default_timer as timer
reload(EdgeTools)
reload(ICCFT)


print "Which peak?"
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
dtBinWidth = 4 
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='sample'
pplmin_frac=3.8; pplmax_frac=4.0; mindtBinWidth=4
'''


'''
#Si 2016
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
'''
#gfp
peaksFile = '/SNS/users/ntv/integrate/mandi_gfp/combined_short.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_gfp/combined.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = 'SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
nxsTemplate = '/SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
pplmin_frac=0.0; pplmax_frac=100.5; mindtBinWidth=15

'''
#PsbO 2016
peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-16286/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
#predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
q_frame = 'lab'
pplmin_frac=0.8; pplmax_frac=1.5; mindtBinWidth=15
'''

'''
#NaK 2017
peaksFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
UBFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_NAK_UB.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-17495/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
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
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.005#np.array([0.003, 0.003, 0.003])
predpplCoefficients = None##np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
'''
'''
#CORELLI - beryl
loadDir = '/data/corelli_beryl/IPTS-20302/'
peaksFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal_indexedonly.integrate'
UBFile =  '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
nxsTemplate = loadDir+'CORELLI_%i.nxs.h5'
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.02#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.; pplmax_frac=4.5; mindtBinWidth=10
'''


'''
#DNA
peaksFile = '/SNS/users/ntv/integrate/mandi_dna/combined_orthorhombic.integrate'
UBFile =  '/SNS/users/ntv/integrate/mandi_dna/combined_orthorhombic.mat'
DetCalFile = None
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
loadDir = '/SNS/MANDI/IPTS-18552/nexus/'
nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.007#np.array([0.003, 0.003, 0.003])
#predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.8; pplmax_frac=2.0; mindtBinWidth=15
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
#panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
dtBinWidth = 25 
dQPixel=0.003#np.array([0.003, 0.003, 0.003])
predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
q_frame='lab'
pplmin_frac=0.8; pplmax_frac=2.0; mindtBinWidth=15
   #---mutant
#peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.integrate'
#UBFile =  '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.mat'
#nxsTemplate = '/SNS/MANDI/IPTS-8776/data/MANDI_%i_event.nxs'
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
#dQPixel = ICCFT.getPixelStep(peak)
dtSpread = 0.015
dtSpreadToPlot = [0.01]
wavelength = peak.getWavelength() #in Angstrom
energy = 81.804 / wavelength**2 / 1000.0 #in eV
flightPath = peak.getL1() + peak.getL2() #in m
scatteringHalfAngle = 0.5*peak.getScattering()
dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=fracHKL))
dQ[dQ>0.4]=0.4


print dQPixel
print peak.getQSampleFrame()
Box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakToGet, dQ, fracHKL=fracHKL,dQPixel=dQPixel, refineCenter=False, q_frame=q_frame)
box = Box
Box.setTitle('Box for peak %i'%peakToGet)
#SaveMD(InputWorkspace=Box, Filename = 'Box.nxs')



instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
if removeEdges:
    try:
        panelDict = pickle.load(open('panelDict_%i.pkl'%peak.getRunNumber(),'rb'))
    except:
        panelDict = EdgeTools.getInstrumentDict(instrumentFile, peaks_ws, peaks_ws.getPeak(peakToGet).getRunNumber(), fitOrder=2)
        pickle.dump(panelDict,open('panelDict_%i.pkl'%peak.getRunNumber(),'wb'))

else:
    panelDict = EdgeTools.getPanelDictionary(instrumentFile)

n_events = Box.getNumEventsArray()

qMask = ICCFT.getHKLMask(UBMatrix, frac=0.25, dQPixel=dQPixel, dQ=dQ)
if not removeEdges:
    mask = np.ones_like(n_events)
else:
    mask = EdgeTools.getMask(peak, Box, panelDict, qMask)

xaxis = Box.getXDimension()
qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
yaxis = Box.getYDimension()
qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
zaxis = Box.getZDimension()
qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
if False:
    ctsX = np.sum(n_events, axis=(1,2))
    ctsY = np.sum(n_events, axis=(0,2))
    ctsZ = np.sum(n_events, axis=(0,1))
    plt.figure(figNumber)
    plt.clf()
    plt.plot(qx, np.sum(n_events, axis=(1,2)),label='qx')
    plt.plot(qy, np.sum(n_events, axis=(0,2)),label='qy')
    plt.plot(qz, np.sum(n_events, axis=(0,1)),label='qz')
    plt.legend(loc='best')

if True:
    #qMask = np.ones_like(n_events).astype(np.bool)
    QX, QY, QZ = np.meshgrid(qx, qy, qz,indexing='ij')
    qSq = QX**2 + QY**2 + QZ**2
    #tList = 1/np.sqrt(qSq)
    #tof = tList * 3176.507 * (peak.getL1() + peak.getL2()) + np.sin(0.5*peak.getScattering())
    tof = 3176.507 * (peak.getL1()+peak.getL2())*np.sin(peak.getScattering()*0.5) * 1/np.sqrt(QX**2 + QY**2 + QZ**2)
    tBins = np.arange(np.min(tof), np.max(tof), 10)
    tofSpread = [0.95*peak.getTOF(), 1.05*peak.getTOF()]
    plt.figure(figNumber)
    plt.clf()
    plt.subplot(2,1,1)
    if removeEdges:
        plt.hist((tof*qMask*mask).flatten(), tBins, edgecolor='none')
    else:
        plt.hist((tof*qMask).flatten(), tBins,edgecolor='none')
    for i, dtS in enumerate(dtSpreadToPlot):
        tofSpread = [(1.0-dtS)*peak.getTOF(), (1.0+dtS)*peak.getTOF()]
        plt.plot(tofSpread, (0.9+0.025*i)*plt.ylim()[1]*np.ones(2),lw=3, label='%4.3f'%dtS)
    plt.legend(loc='best')
    plt.title('Peak %i'%peakToGet)
    plt.xlabel('TOF (us)')
    plt.ylabel('Num Pixels in Box')
    topXLim = plt.xlim()
    peakTOF = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    edgesToCheck = []#EdgeTools.needsEdgeRemoval(Box,panelDict,peak)
    for dtS in dtSpreadToPlot:
        tofWS,pp_lambda = ICCFT.getTOFWS(Box,flightPath, scatteringHalfAngle, peakTOF, peak, panelDict, qMask, dtBinWidth=dtBinWidth,dtSpread=dtS,doVolumeNormalization=False, minFracPixels=0.015, removeEdges=removeEdges, edgesToCheck=edgesToCheck, calcTOFPerPixel=False, zBG=-1., calc_pp_lambda=False, mindtBinWidth=mindtBinWidth)
        YDATA1 = tofWS.readY(0).copy()
        plt.subplot(2,1,2)
        plt.plot(tofWS.readX(0), tofWS.readY(0),'o',label='%2.3f'%dtS)
        tofWS,pp_lambda = ICCFT.getTOFWS(Box,flightPath, scatteringHalfAngle, peakTOF, peak, panelDict, qMask, dtBinWidth=dtBinWidth,dtSpread=dtS,doVolumeNormalization=False, minFracPixels=0.005, removeEdges=removeEdges, edgesToCheck=edgesToCheck, calcTOFPerPixel=False, zBG=1.96,neigh_length_m=3, padeCoefficients=ICCFT.getModeratorCoefficients('franz_coefficients_2017.dat'), predCoefficients=predpplCoefficients, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac,mindtBinWidth=mindtBinWidth)
        YDATA2 = tofWS.readY(0).copy()
        print pp_lambda
        plt.subplot(2,1,2)
        plt.plot(tofWS.readX(0), tofWS.readY(0),'-o',label='%2.3f'%dtS)
    plt.xlabel('TOF (us)')
    plt.ylabel('Counts')
    plt.legend(loc='best')
    plt.xlim(topXLim)
    plt.plot([peakTOF, peakTOF],plt.ylim(),'k',lw=2)
    plt.subplot(2,1,1)
    plt.plot([peakTOF, peakTOF],plt.ylim(),'k',lw=2)
    theta = np.arctan(QY/QX)[qMask]
    phi = np.arccos(QZ/np.sqrt(QX**2 + QY**2 + QZ**2))[qMask]
    '''
    def gaussian(x, amp, cen, wid, BG):
        return amp * np.exp(-(x-cen)**2 /wid) + BG

    numPixelsInBin = plt.hist((tof*qMask*mask).flatten(), tBins, edgecolor='none')[0]
    t = 0.5*(tBins[1:] + tBins[:-1])
    plt.figure(2); plt.clf();
    plt.plot(t, numPixelsInBin)
    gp = curve_fit(gaussian, t, numPixelsInBin,p0=[np.max(numPixelsInBin),np.mean(t),50.0,0])[0]
    yG = gaussian(t,gp[0], gp[1],gp[2],gp[3])
    plt.plot(t,yG,'r--')

    tofWS = ICCFT.getTOFWS(Box,flightPath, scatteringHalfAngle, peakTOF, peak, panelDict, qMask, dtBinWidth=dtBinWidth,dtSpread=0.05,doVolumeNormalization=False, minFracPixels=0.015, removeEdges=removeEdges,tListMode=1)
    y = tofWS.readY(0) 
    t = tofWS.readX(0)
    bgIDX = np.logical_or(t<peak.getTOF()*0.995, t>peak.getTOF()*1.005)
    #bgIDX = np.ones_like(t).astype(np.bool)
    plt.figure(2); plt.clf();
    plt.plot(t, y,'bo')
    
    gpD = curve_fit(gaussian, t[bgIDX], y[bgIDX],p0=[np.mean(y),gp[1], gp[2], 0])[0]
    yGD = gaussian(t,gpD[0], gpD[1],gpD[2],gpD[3])
    plt.plot(t,yGD,'r--')
    
    gP = np.polyfit(t[bgIDX],y[bgIDX],2)
    yP = np.polyval(gP,t)
    #plt.plot(t,yP,'m--')
    '''
