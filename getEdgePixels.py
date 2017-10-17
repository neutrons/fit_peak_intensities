import xml.etree.ElementTree
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from mpl_toolkits.mplot3d import Axes3D
import sys
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import ICCFitTools as ICCFT
import itertools
import re
import pickle
from timeit import default_timer as timer
#reload(ICCFT)

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, residual, _, _ = np.linalg.lstsq(G, z)
    return m, residual

#n.b. this is hardcoded at 2nd order right now
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    z = np.zeros_like(x)
    XSQ = x*x
    YSQ = y*y
    z = m[0] + m[1]*y + m[2]*YSQ + m[3]*x + m[4]*x*y + m[5]*x*YSQ + m[6]*XSQ + m[7]*x*x*y + m[8]*XSQ*YSQ
    #ij = itertools.product(range(order+1), range(order+1))
    #z = np.zeros_like(x)
    #for a, (i,j) in zip(m, ij):
    #    z += a * x**i * y**j
    return z


def getDetectorBank(panelDict, detID):
    for key in panelDict.keys():
        panel = panelDict[key]
        if panel['idstart'] <= detID <= panel['idend']:
            return panel


def getBaseEdge(N, offset = 0):
    edgePixels =list()
    #First the top
    for i in range(N):
        edgePixels.append(i)
    #Now the left side
    for i in range(1,N):
        edgePixels.append(i*N)
    #Now the right side
    for i in range(1,N):
        if i != N-1:
            edgePixels.append((i+1)*N-1)
    #Now the sides
    for i in range(N*(N-1),N*N):
        edgePixels.append(i)
    return np.array(edgePixels) + offset

def getPanelDictionary(instrumentFile):
    et = xml.etree.ElementTree.parse(instrumentFile)
    root = et.getroot()
    panelDict = {} 
    for i, it in enumerate(root.getiterator()):
        if 'type' in it.attrib:
            if it.attrib['type'] == 'panel':
                panel = {}
                panel['idstart'] = int(it.attrib['idstart'])
                panel['idend'] = int(it.attrib['idstart']) + int(it.attrib['idstepbyrow'])**2 - 1
                panel['idstepbyrow'] = int(it.attrib['idstepbyrow'])
                panel['idfillbyfirst'] = it.attrib['idfillbyfirst']
                panel['name'] = it[0].attrib['name']
                panel['bankNumber'] = int(it[0].attrib['name'].strip('bank'))
                #panelList.append(panel)
                panelDict[panel['name']] = panel
    return panelDict

def getEdgePixels(panel):
    pixelList = np.empty(0)
    try:
            nSide = int(panel['idstepbyrow'])
            offSet = int(panel['idstart'])
            tmpPxl = getBaseEdge(nSide,offSet)
            pixelList = np.append(pixelList, tmpPxl)
            return pixelList
    except:
            print 'Error with panel starting with ', panel['idstart']

def getQList(peak, pixelList, WLList, edgeNumbers):
    qList = list()
    for wl in WLList:
        peak.setWavelength(wl)
        for edgeNumber in edgeNumbers:
            for pxl in pixelList[edgeNumber*1021//4+1:(edgeNumber+1)*1021//4]:
                try:
                    peak.setDetectorID(int(pxl))
                    qList.append(peak.getQSampleFrame())
                except:
                    print "Cannot set detector ID %i"%pxl
    qList = np.asarray(qList)
    return qList

def fitSurface(qList, order=3):
    surfaceCoefficients, resid = polyfit2d(qList[:,0],qList[:,1], qList[:,2],order=order)
    switchXZ = False
    if resid > 5:
        switchXZ = True
        surfaceCoefficients, resid = polyfit2d(qList[:,2], qList[:,1], qList[:,0],order=order)
        surf = polyval2d(qList[:,2], qList[:,1], surfaceCoefficients) 

    else:
        surf = polyval2d(qList[:,0], qList[:,1], surfaceCoefficients) 
    return surfaceCoefficients, surf, resid, switchXZ


def makeFigure(peak, qList, surf,UBMatrix, fracHKL=0.8, figNumber=2, showQBoxOnly=False):
    plt.clf()
    fig = plt.figure(figNumber)
    ax = Axes3D(fig)
    ax.scatter3D(qList[:,0],qList[:,1],qList[:,2],color='b',linewidth=0)
    ax.scatter3D(qList[:,0],qList[:,1],surf,color='g',linewidth=0)
    ax.set_xlabel('qx')
    ax.set_ylabel('qy')
    ax.set_zlabel('qz')

    if showQBoxOnly: 
        qS = peak.getQSampleFrame()
        dQ = np.abs(ICCFT.getDQFracHKL(peak, UBMatrix, frac = fracHKL))
        badIDX = (np.sum(qList < qS-dQ[:,0],axis=1) + np.sum(qList>qS+dQ[:,1],axis=1) > 0) > 0
        qList = qList[~badIDX,:] 
        ax.set_xbound(qS[0] - dQ[0,0], qS[0] + dQ[0,1])
        ax.set_ybound(qS[1] - dQ[1,0], qS[1] + dQ[1,1])
        ax.set_zbound(qS[2] - dQ[2,0], qS[2] + dQ[2,1])

def findPeakInRun(peaks_ws, runNumber):
    for i in range(peaks_ws.getNumberPeaks()):
        if peaks_ws.getPeak(i).getRunNumber() == runNumber:
            print 'Building panelDict on peak %i, run numbers should match: %i, %i'%(i, runNumber, peaks_ws.getPeak(i).getRunNumber())
            return peaks_ws.getPeak(i)
    print 'getEdgePixels.py - findPeakInRunNumber cannot find a peak matching the run number.'
    return None

def getInstrumentDict(instrumentFile, peaks_ws, runNumber, edgeNumbers=[0,1,2,3], fitOrder=3):

    print 'Setting Q boundaries based on instrument file: %s'%instrumentFile
    panelDict = getPanelDictionary(instrumentFile)
    # Get the edge pixels of our detector
    #  TODO: Right now, everything is in QSample - so we need to recalculate Q coverage for each orientation
    #  to account for the rotation of the detctors w.r.t. as seen by the sample.  In principle (I think)
    #  once for each instrument configuration and convert to the current orientation using the goniometer
    #  matrix.  To do that, we need to expose peak.getGoniometer() to python (should be easy).
    peak = findPeakInRun(peaks_ws, runNumber)
    wavelength_original = peak.getWavelength()
    detID_original = peak.getDetectorID()
    panel = getDetectorBank(panelDict, peak.getDetectorID())
    goodFitNames = []
    failedFitNames = []
    for pID, key in enumerate(panelDict.keys()):
        panel = panelDict[key]
        print 'Fitting detector %s - detector number %i of %i'%(panel['name'], pID+1, len(panelDict))
        pixelList = getEdgePixels(panel)
        #WLList = np.linspace(0.95*wavelength_original, 1.05*wavelength_original, 10)
        energyList = np.arange(0.001,0.5,0.005)
        WLList = np.sqrt(81.804/1000/energyList)
        try:
            for edge in edgeNumbers:
                qList = getQList(peak, pixelList, WLList, [edge])
                # Append our dictionary
                surfaceCoefficients, surf, resid, switchXZ = fitSurface(qList, order=fitOrder)
                panel['surf_coefficients_%i'%edge] = surfaceCoefficients
                panel['switch_xz_%i'%edge] = switchXZ
                panel['resid_%i'%edge] = resid
                goodFitNames.append(panel['name'])
        except KeyboardInterrupt:
            print "KEYBOARDINTERRUPT - Exiting"
            sys.exit()
        except:
            print 'ERROR WITH DETECTOR %s - skipping!'%panel['name']
            failedFitNames.append(panel['name'])
    #Reset peak properties
    peak.setWavelength(wavelength_original)
    peak.setDetectorID(detID_original)

    if len(failedFitNames) > 0:
        print 'Failed to fit detectors:', ' '.join(name for name in failedFitNames)
    else:
        print 'Successfully fit all detectors!'
    return panelDict


def needsEdgeRemoval(Box, panelDict, peak):
    panel = getDetectorBank(panelDict, peak.getDetectorID())
    cx = [Box.getDimension(0).getMinimum(), Box.getDimension(0).getMaximum()]
    cy = [Box.getDimension(1).getMinimum(), Box.getDimension(1).getMaximum()]
    cz = [Box.getDimension(2).getMinimum(), Box.getDimension(2).getMaximum()]
    CX, CY, CZ = np.meshgrid(cx,cy,cz,copy=False) 
    edgesToCheck = list()
    for i in range(4):
        if not panel['switch_xz_%i'%i]:
            gtVect = CZ > polyval2d(CX, CY, panel['surf_coefficients_%i'%i])
        else:#We had a poorly formed polynomial, so we fit in X instead of Z
            gtVect = CX > polyval2d(CZ, CY, panel['surf_coefficients_%i'%i])
        if not ((gtVect == False).all() or (gtVect == True).all()):
            edgesToCheck.append(i)
    return edgesToCheck


def getMask(peak, Box, panelDict, qMask, edgesToCheck=[0,1,2,3]):
    maskList = list()
    qS = peak.getQSampleFrame()
    panel = getDetectorBank(panelDict, peak.getDetectorID())

    xaxis = Box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = Box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = Box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    QX,QY,QZ = np.meshgrid(qx,qy,qz,indexing='ij',copy=False)

    t1 = timer()
    # -- Prototype for mask generation
    for i in edgesToCheck:
        tmpMask = np.zeros_like(QX).astype(np.bool)
        if not panel['switch_xz_%i'%i]:
            gtVect = qS > polyval2d(qS[0], qS[1], panel['surf_coefficients_%i'%i])
            peakGreaterZ = gtVect[2]
            if peakGreaterZ:
                tmpMask[qMask] += QZ[qMask] > polyval2d(QX[qMask], QY[qMask], panel['surf_coefficients_%i'%i])
            else:
                tmpMask[qMask] += QZ[qMask] < polyval2d(QX[qMask], QY[qMask], panel['surf_coefficients_%i'%i]) 
        else:#We had a poorly formed polynomial, so we fit in X instead of Z
            gtVect = qS > polyval2d(qS[2], qS[1], panel['surf_coefficients_%i'%i])
            peakGreaterX = gtVect[0]
            if peakGreaterX:
                tmpMask[qMask] += QX[qMask] > polyval2d(QZ[qMask], QY[qMask], panel['surf_coefficients_%i'%i])
            else:
                tmpMask[qMask] += QX[qMask] < polyval2d(QZ[qMask], QY[qMask], panel['surf_coefficients_%i'%i]) 
        maskList.append(tmpMask)
    maskList = np.asarray(maskList)
    mask = reduce(np.logical_and, maskList)
    print 'Made mask in %f s'%(timer()-t1)
    return mask

def getPeaksWS(peaksFile, UBFile):
    importPeaks = True
    for ws in mtd.getObjectNames():
        if mtd[ws].getComment() == '%s'%peaksFile:
            print '    using already loaded peaks file'
            importPeaks = False
            peaks_ws = mtd[ws]
    if importPeaks:
        peaks_ws = LoadIsawPeaks(Filename = peaksFile)
        peaks_ws.setComment(peaksFile)
        LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
    return peaks_ws, peaks_ws.sample().getOrientedLattice().getUB()

def getMDData(peak,loadDir,DetCalFile ):
    importFlag = True
    for ws in mtd.getObjectNames():
        if mtd[ws].getComment() == 'BSGETBOX%i'%peak.getRunNumber():
            print '   Using already loaded MDdata'
            MDdata = mtd[ws]
            importFlag = False
            break
    if importFlag:
        nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
        fileName = nxsTemplate%peak.getRunNumber()
        MDdata = ICCFT.getSample(peak.getRunNumber(), DetCalFile, '', fileName)
        MDdata.setComment('BSGETBOX%i'%peak.getRunNumber())
    return MDdata

def getInstrumentFile(peaks_ws, peaksFile):
    instrument = peaks_ws.getInstrument().getName()
    with open(peaksFile, 'r') as f:
        peaksFileFirstLine = f.readline()
    runStart = re.findall('\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d', peaksFileFirstLine)[0]
    instrumentFile = '/opt/Mantid/instrument/TOPAZ_Definition_2017-02-24.xml'
    instrumentFile = peaks_ws.getInstrumentFilename(instrument, runStart)
    return instrumentFile



if __name__ == '__main__':

    # --- import the peaks if they're not loaded
    peakToGet = 195
    peaksFile ='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
    UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'

    #peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
    #UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'

    loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
    DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'

    peaks_ws, UBMatrix = getPeaksWS(peaksFile, UBFile)
    #UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()

    instrumentFile = getInstrumentFile(peaks_ws, peaksFile)
    #panelDict = pickle.load(open('panelDict_15647.pkl','rb'))
    panelDict = getInstrumentDict(instrumentFile, peaks_ws,fitOrder=2)

    peak = peaks_ws.getPeak(peakToGet)
    MDdata = getMDData(peak, loadDir, DetCalFile)
    Box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakToGet, fracHKL=0.5, refineCenter=False)
    n_events = Box.getNumEventsArray()
    mask = getMask(peak,Box, panelDict )

    iiid = mask.shape[1]//2
    q = mask[:,mask.shape[1]//2,:]
    r = n_events[:,n_events.shape[1]//2,:]
    plt.clf()
    plt.imshow(r)
    plt.hold('on')
    plt.imshow(q,cmap='gray',alpha=0.2)
