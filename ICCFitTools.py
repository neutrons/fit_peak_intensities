import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import sys
import os
from scipy.interpolate import interp1d
from scipy.misc import factorial
from scipy.optimize import curve_fit
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
from mantid.kernel import V3D
import pickle
from scipy import interpolate
import ICConvoluted as ICC
reload(ICC)
import getEdgePixels as EdgeTools
reload(EdgeTools)
import itertools
from scipy.interpolate import LinearNDInterpolator
from timeit import default_timer as timer
from scipy.ndimage.filters import convolve

def scatFun(x, A, bg):
    return A/x+bg

def oldScatFun(x,A,k,bg):
    return 1.0*A*np.exp(-k*x) + bg

def calcSomeTOF(box, peak, refitIDX = None):
    xaxis = box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    QX, QY, QZ = getQXQYQZ(box)

    if refitIDX  is None:
        refitIDX = np.ones_like(QX).astype(np.bool)

    from mantid.kernel import V3D
    qS0 = peak.getQSampleFrame()
    PIXELFACTOR = np.ones_like(QX)*(peak.getL1() + peak.getL2())*np.sin(0.5*peak.getScattering())
    for i, x in enumerate(qx):
        print i
        for j, y in enumerate(qy):
            for k, z in enumerate(qz):
                if refitIDX[i,j,k]:
                    qNew = V3D(x,y,z)
                    peak.setQSampleFrame(qNew)
                    L = peak.getL1() + peak.getL2()
                    HALFSCAT = 0.5*peak.getScattering()
                    PIXELFACTOR[i,j,k] = L*np.sin(HALFSCAT)
    peak.setQSampleFrame(qS0)


    tofBox = 3176.507 * PIXELFACTOR * 1.0/np.sqrt(QX**2 + QY**2 + QZ**2)
    return tofBox


# (x,y,z) -> (r,phi,theta)
def cart2sph(x,y,z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return r, az, el

def getQXQYQZ(box):
    xaxis = box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    QX, QY, QZ = np.meshgrid(qx, qy, qz,indexing='ij',copy=False)
    return QX, QY, QZ

def getQuickTOFWS(box, peak, padeCoefficients, goodIDX=None, dtSpread=0.03, dtBinWidth=30, qMask=None, pp_lambda=None, nBG=15):
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    detNumber = 0#EdgeTools.getDetectorBank(panelDict, peak.getDetectorID())['bankNumber']
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool)
  
    calc_pp_lambda=False
    if pp_lambda is None:
        calc_pp_lambda=True

    tofWS,ppl = getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, None, 0, qMask, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=False, minFracPixels=0.01, removeEdges=False,calcTOFPerPixel=False,neigh_length_m=3,zBG=1.96,pp_lambda=pp_lambda,calc_pp_lambda=calc_pp_lambda)
    print 'ySum:', np.sum(tofWS.readY(0))
    fitResults,fICC = doICCFit(tofWS, energy, flightPath, padeCoefficients, 0, None,nBG=nBG,fitOrder=1,constraintScheme=2)
    h = [tofWS.readY(0), tofWS.readX(0)]
    chiSq = fitResults.OutputChi2overDoF
    
    r = mtd['fit_Workspace']
    param = mtd['fit_Parameters']
    n_events = box.getNumEventsArray()
     
 
    iii = fICC.numParams() - 1
    fitBG = [param.row(int(iii+i+1))['Value'] for i in range(1+1)]


    #Set the intensity before moving on to the next peak
    icProfile = r.readY(1)
    bgCoefficients = fitBG[::-1]

    #peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))i
    t0 = param.row(3)['Value']
    intensity, sigma, xStart, xStop = integratePeak(r.readX(0), icProfile,r.readY(0), np.polyval(bgCoefficients, r.readX(1)), pp_lambda=pp_lambda, fracStop=0.01,totEvents=np.sum(n_events[goodIDX*qMask]), bgEvents=np.sum(goodIDX*qMask)*pp_lambda)

 
    return chiSq, h, intensity, sigma


def getPoissionGoodIDX(n_events, zBG=1.96, neigh_length_m=3):
    hasEventsIDX = n_events>0
    #Set up some things to only consider good pixels
    N = np.shape(n_events)[0]
    neigh_length_m = neigh_length_m #Set to zero for "this pixel only" mode - performance is optimized for neigh_length_m=0 
    maxBin = np.shape(n_events)

    pp_lambda = get_pp_lambda(n_events,hasEventsIDX) #Get the most probably number of events
    found_pp_lambda = False
    convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
    conv_n_events = convolve(n_events,convBox)
    allEvents = np.sum(n_events[hasEventsIDX])
    if allEvents > 0:
        while not found_pp_lambda and pp_lambda < 3.0:
            goodIDX = np.logical_and(hasEventsIDX, conv_n_events > pp_lambda+zBG*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
            boxMean = n_events[goodIDX]
            boxMeanIDX = np.where(goodIDX)
            if allEvents > np.sum(boxMean):
                found_pp_lambda = True
            else:
                pp_lambda *= 1.05
    return goodIDX, pp_lambda

def getOptimizedGoodIDX(n_events, padeCoefficients, zBG=1.96, neigh_length_m=3,dtBinWidth=4, qMask=None, peak=None, box=None, pp_lambda=None,peakNumber=-1,nBG=15):
    #Set up some things to only consider good pixels
    hasEventsIDX = n_events>0
    N = np.shape(n_events)[0]
    neigh_length_m = neigh_length_m #Set to zero for "this pixel only" mode - performance is optimized for neigh_length_m=0 
    found_pp_lambda = False
    convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
    conv_n_events = convolve(n_events,convBox)
    pp_lambda = get_pp_lambda(n_events,hasEventsIDX) #Get the most probable number of events
    print pp_lambda, conv_n_events.max(), '~~~~~~'
    pp_lambda_toCheck = np.unique(conv_n_events)
    pp_lambda_toCheck = pp_lambda_toCheck[1:][np.diff(pp_lambda_toCheck)>0.001]
    if peak is not None: #TODO: This MUST be parameterized, keep it hard coded ONLY for testing
        pred_ppl = scatFun(np.sin(0.5*peak.getScattering())**2/peak.getWavelength()**4, 0.00122958,  0.29769245)
        pred_ppl = oldScatFun(peak.getScattering()/peak.getWavelength(),5.24730283,  7.23719321,  0.27449887) 
        minppl = 0.8*pred_ppl
        maxppl = 1.5*pred_ppl 
    else:
        minppl=0
        maxppl = pp_lambda_toCheck.max() + 0.5 #add some just to make sure we don't skip any
    pp_lambda_toCheck = pp_lambda_toCheck[pp_lambda_toCheck > minppl]
    pp_lambda_toCheck = pp_lambda_toCheck[pp_lambda_toCheck < maxppl]
#    zBG = 1.000

    chiSqList = 1.0e30*np.ones_like(pp_lambda_toCheck)
    ISIGList = 1.0e-30*np.ones_like(pp_lambda_toCheck)
    IList = 1.0e-30*np.ones_like(pp_lambda_toCheck)
    #hList = []
    oldGoodIDXSum = -1.0
    for i, pp_lambda in enumerate(pp_lambda_toCheck):
        try:
            goodIDX = np.logical_and(hasEventsIDX, conv_n_events > pp_lambda+zBG*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
            if np.sum(goodIDX) == oldGoodIDXSum: #No new points removed, we skip this
                #print '#############skipping pp_lambda=%4.4f because no new entries'%pp_lambda
                continue
            else:
                oldGoodIDXSum = np.sum(goodIDX)
            try: 
                chiSq, h, intens, sigma = getQuickTOFWS(box, peak, padeCoefficients, goodIDX=goodIDX,qMask=qMask,pp_lambda=pp_lambda,dtBinWidth=dtBinWidth,nBG=nBG)
            except:
                break
            chiSqList[i] = chiSq
            ISIGList[i] = intens/sigma
            IList[i] = intens
            #hList.append((pp_lambda, chiSq, h))
            if len(h[0])<10:#or np.sum(h[0])<10: #or (chiSq > 100 and np.min(chiSqList)<5):
                 break
        except RuntimeError:
            #This is caused by there being fewer datapoints remaining than parameters.  For now, we just hope
            # we found a satisfactory answer.  TODO: we can rebin and try that, though it may not help much.
            break
        except KeyboardInterrupt:
            0/0
    #pickle.dump(hList, open('/home/ntv/analysis/data/hList_beta_lac_peak2.pkl','wb'))
    #print chiSqList[:i+1], 'is chiSqList'
    #print ISIGList[:i+1], 'is ISIG'
    #print IList[:i+1], 'is Intens'
    print '\n'.join([str(v) for v in zip(chiSqList[:i+1], ISIGList[:i+1], IList[:i+1])])
    chiSqConsider = np.logical_and(chiSqList < 1.4, chiSqList>0.9)
    if np.sum(chiSqConsider) > 1.0:
        use_ppl = np.argmax(ISIGList[chiSqConsider])
        pp_lambda = pp_lambda_toCheck[chiSqConsider][use_ppl]
        print 'USING PP_LAMBDA', pp_lambda, 'WITH CHISQ:', chiSqList[chiSqConsider][use_ppl]
    else:
        use_ppl = np.argmin(np.abs(chiSqList[:i+1]-1.0))
        pp_lambda = pp_lambda_toCheck[use_ppl]
        print 'USING PP_LAMBDA', pp_lambda, 'WITH CHISQ:', chiSqList[use_ppl]
    #goodIDX = np.logical_and(hasEventsIDX, conv_n_events > pp_lambda+zBG*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
    goodIDX, _ = getBGRemovedIndices(n_events, pp_lambda=pp_lambda)

    chiSq, h, intens, sigma = getQuickTOFWS(box, peak, padeCoefficients, goodIDX=goodIDX,qMask=qMask,pp_lambda=pp_lambda,dtBinWidth=dtBinWidth,nBG=nBG)
    return goodIDX, pp_lambda
 


#Must give this a peak, box, and qMask to do iterative pp_lambda
def getBGRemovedIndices(n_events,zBG=1.96,calc_pp_lambda=False, neigh_length_m=3,dtBinWidth=4, qMask=None, 
                        peak=None, box=None, pp_lambda=None,peakNumber=-1, padeCoefficients=None,nBG=15):

    if calc_pp_lambda is True and pp_lambda is not None:
        import sys
        sys.exit('Error in ICCFT:getBGRemovedIndices: You should not calculate and specify pp_lambda.')

    if calc_pp_lambda is True and padeCoefficients is None:
        import sys
        sys.exit('Error in ICCFT:getBGRemovedIndices: calc_pp_lambda is True, but no moderator coefficients are provided.')


    if pp_lambda is not None:
        #Set up some things to only consider good pixels
        hasEventsIDX = n_events>0
        N = np.shape(n_events)[0]
        neigh_length_m = neigh_length_m #Set to zero for "this pixel only" mode - performance is optimized for neigh_length_m=0 
        convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
        conv_n_events = convolve(n_events,convBox)
        goodIDX = np.logical_and(hasEventsIDX, conv_n_events > pp_lambda+zBG*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
        return goodIDX, pp_lambda


    if calc_pp_lambda is False:
        return getPoissionGoodIDX(n_events, zBG=zBG, neigh_length_m=neigh_length_m) 
    
    if peak is not None and box is not None and padeCoefficients is not None:
        return getOptimizedGoodIDX(n_events, padeCoefficients, zBG=1.96, neigh_length_m=neigh_length_m,
            dtBinWidth=dtBinWidth, qMask=qMask, peak=peak, box=box, pp_lambda=pp_lambda,peakNumber=peakNumber,nBG=nBG)
    print 'ERROR WITH ICCFT:getBGRemovedIndices!' 

def getDQTOF(peak, dtSpread=0.03, maxDQ=0.5):
    dQ=np.zeros(3)
    dtTarget = dtSpread*peak.getTOF()
    gamma = 3176.507*(peak.getL1()+peak.getL2())*np.sin(peak.getScattering()*0.5)
    for i in range(3):
        gradientVector=np.zeros(3)
        gradientVector[i]= np.sign(peak.getQSampleFrame()[i]) 
        for qStep in np.linspace(0,maxDQ, 100):
            dt = np.abs((gamma /np.linalg.norm(peak.getQSampleFrame()+gradientVector*qStep)) - peak.getTOF())
            if dt >dtTarget:
                dQ[i] = qStep
                break;
        if dQ[i] == 0:
            dQ[i] = maxDQ
    dQ2d = np.array([[dQ[0],dQ[0]],[dQ[1],dQ[1]],[dQ[2],dQ[2]]])
    return dQ2d

def getPixelStep(peak, dtBin=4):
    gamma = 3176.507*(peak.getL1()+peak.getL2())*np.sin(peak.getScattering()*0.5)
    gradientVector = -1.0*np.sign(peak.getQSampleFrame())
    gradientVector *= 1.0/np.sqrt(3.0)
    
    for qStep in np.linspace(0,0.005, 100):
            dt = np.abs((gamma /np.linalg.norm(peak.getQSampleFrame()+gradientVector*qStep)) - peak.getTOF())
            if dt >dtBin/np.sqrt(3):
                return qStep

    return 0.005
    

# UB = UBmatrix as loaded by LoadIsawUB().  Only works in the
#    Qsample frame right now
#   TODO - calculate this once per run, not every peak
def getDQFracHKL(UB, frac=0.5):
    dQ = np.zeros((3,2))

    q = [2*np.pi*frac*UB.dot(v) for v in [seq for seq in itertools.product([-1.0,1.0],repeat=3)]]
    dQ[:,0] = np.max(q,axis=0)#TODO THIS CAN BE 1D since it's symmetric
    dQ[:,1] = np.min(q,axis=0)
    return dQ

def getHKLMask(UB, frac=0.5,dQPixel=0.005, dQ=None):
    if dQ is None:
        dQ = np.abs(getDQFracHKL(UB, frac=frac))
        dQ[dQ>0.5] = 0.5
    nPtsQ = np.round(np.sum(dQ/dQPixel,axis=1)).astype(int)
    h0 = 1.0; k0 = 27.0; l0=7.0
    qDummy = 2*np.pi*UB.dot(np.asarray([h0, k0, l0]))
    qx = np.linspace(qDummy[0]-dQ[0,0], qDummy[0]+dQ[0,1], nPtsQ[0])
    qy = np.linspace(qDummy[1]-dQ[1,0], qDummy[1]+dQ[1,1], nPtsQ[1])
    qz = np.linspace(qDummy[2]-dQ[2,0], qDummy[2]+dQ[2,1], nPtsQ[2])
    QX,QY,QZ = np.meshgrid(qx,qy,qz,indexing='ij',copy=False)
    UBinv = np.linalg.inv(UB)
    tHKL = UBinv.dot([QX.ravel(),QY.ravel(),QZ.ravel()])/2/np.pi
    H = np.reshape(tHKL[0,:], QX.shape)
    K = np.reshape(tHKL[1,:], QX.shape)
    L = np.reshape(tHKL[2,:], QX.shape)
    mask = reduce(np.logical_and,  [H>h0-frac, H<h0+frac, K>k0-frac, K<k0+frac, L>l0-frac, L<l0+frac])
    return mask 


#Wrapper for pade with parametres as separate arguments for compatability with scipy.optimize.curve_fit
def padeWrapper(x,a,b,c,d,f,g,h,i,j,k):
    #print a,b,c,d,f,g,h,i,j,k
    pArr = np.zeros(10)
    pArr[0] = a; pArr[1] = b; pArr[2] = c; pArr[3] = d; pArr[4] = f; 
    pArr[5] = g; pArr[6] = h; pArr[7] = i; pArr[8] = j; pArr[9] = k;
    return pade(pArr,x)

#Standard Pade function used for getting initial guesses
def pade(c,x): #c are coefficients, x is the energy in eV
    return c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9])

def integratePeak(x, yFit, yData, bg, pp_lambda=0, fracStop = 0.01,totEvents=1, bgEvents=1):
    #Find out start/stop point
    yScaled = (yFit-bg) / np.max(yFit-bg)
    goodIDX = yScaled > fracStop
    if np.sum(goodIDX) > 0:
        iStart = np.min(np.where(goodIDX))
        iStop = np.max(np.where(goodIDX))
        xStart = x[iStart]
        xStop = x[iStop]
    else:
        print 'ICCFITTOOLS:integratePeak - NO GOOD START/STOP POINT!!'
        return 0.0, 1.0, x[0], x[-1]
 
    #Do the integration
    intensity = np.sum(yFit[iStart:iStop] - bg[iStart:iStop])

    #Calculate the background sigma = sqrt(var(Fit) + sum(BG))
    yFitSum = np.sum(yFit[iStart:iStop])
    bgSum = np.abs(np.sum(bg[iStart:iStop]))
    #varFit = np.average((yData-yFit)**2,weights=(yData-bg))   
    #sigma = np.sqrt(varFit + bgSum)
    sigma = np.sqrt(totEvents + bgEvents) 
    #sigma = np.sqrt(totEvents)
    print 'Intensity: ', intensity, 'Sigma: ', sigma, 'pp_lambda:', pp_lambda
    return intensity, sigma, xStart, xStop

#Poission distribution
def poisson(k,lam):
    return (lam**k)*(np.exp(-lam)/factorial(k))

#Normalized Poission distribution
def normPoiss(k,lam,eventHist):
    pp_dist = poisson(k,lam)
    pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
    return pp_dist*pp_n


#Estimate the most likely number of events based on a Poission
#distribution of the box.  n_events an ND array (N=3 for Qx,Qy,Qz)
#containing the number of events in each pixel.  Returns pp_lambda,
#the most likely number of events.
def get_pp_lambda(n_events, hasEventsIDX ):

    eventValues = n_events[hasEventsIDX]
    numEvents = np.sum(eventValues)
    eventHist = np.bincount(eventValues.astype('int'))[1:]
    eventX = np.array(range(len(eventHist)))+1
    normPoissL = lambda k,lam: normPoiss(k,lam,eventHist)
    pp_lambda,cov = curve_fit(normPoissL, eventX, eventHist,p0=0.1)
    return pp_lambda


#Builds a TOF profile from the data in box which is nominally centered around a peak.
#Input:
#    box: in a binned MD box.
#    flightPath = L1+L2 (units: m)
#    scatteringHalfAngle: the scattering half angle (units: rad)
#    tofPeak: the nominal TOF of the peak (units: us)
#    dtBinWidth: the spacing between points in the resulting TOF histogram (units: us)
#    zBG: the z score at which we will accept pixels (i.e. 1.96 for 95% CI).  If zBG<0
#         then we will not remove background and will instead just consider each pixel
#         individually
#    dtSpread: the fraction of t around tofPeak we will consider
#Output:
#    tofWS: a Workspace2D containing the TOF profile.  X-axis is TOF (units: us) and
#           Y-axis is the number of events.
def getTOFWS(box, flightPath, scatteringHalfAngle, tofPeak, peak, panelDict, peakNumber, qMask, dtBinWidth=2, zBG=-1.0, dtSpread = 0.02, doVolumeNormalization=False, minFracPixels = 0.005, removeEdges=False, edgesToCheck=None, calcTOFPerPixel=False, workspaceNumber=None,neigh_length_m=0, pp_lambda=None, calc_pp_lambda=False, padeCoefficients=None):
    #Find the qVoxels to use
    n_events = box.getNumEventsArray()
    hasEventsIDX = np.logical_and(n_events>0,qMask)
    #print '~~~ ', np.sum(n_events), np.sum(n_events[qMask]), np.sum(n_events[hasEventsIDX])

    #Set up some things to only consider good pixels
    N = np.shape(n_events)[0]
    maxBin = np.shape(n_events)
    if zBG >= 0:
        if pp_lambda is None:
            calc_pp_lambda=True
        goodIDX, pp_lambda= getBGRemovedIndices(n_events,box=box, qMask=qMask, peak=peak, pp_lambda=pp_lambda, peakNumber=peakNumber,dtBinWidth=dtBinWidth, calc_pp_lambda=calc_pp_lambda, padeCoefficients=padeCoefficients)
        hasEventsIDX = np.logical_and(goodIDX, qMask) #TODO bad naming, but a lot of the naming in this function assumes it
        boxMean = n_events[hasEventsIDX]
        boxMeanIDX = np.where(hasEventsIDX)
    else: #don't do background removal - just consider one pixel at a time
        pp_lambda = 0 
        boxMean = n_events[hasEventsIDX]
        boxMeanIDX = np.where(hasEventsIDX)
    boxMeanIDX = np.asarray(boxMeanIDX) 
    useIDX = boxMeanIDX.transpose()
 
    #Setup our axes -- ask if there is a way to just get this
    xaxis = box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    QX, QY, QZ = np.meshgrid(qx, qy, qz,indexing='ij',copy=False)

    #Create our TOF distribution from bg corrected data
    if calcTOFPerPixel == False:
        tList = 1.0/np.sqrt(QX[hasEventsIDX]**2 + QY[hasEventsIDX]**2 + QZ[hasEventsIDX]**2)
        tList = 3176.507 * flightPath * np.sin(scatteringHalfAngle) * tList #convert to microseconds
        #refitIDX = hasEventsIDX
        #tList = calcSomeTOF(box, peak, refitIDX=refitIDX)
        #tList=tList[hasEventsIDX] 
    if calcTOFPerPixel == True:
        origFlightPath = flightPath
        origScatteringHalfAngle = scatteringHalfAngle
        def getTList(peak,qx,qy,qz,boxMeanIDX):
            origQS = peak.getQSampleFrame()
            tList = []
            for idx in boxMeanIDX.transpose():
                newQ = V3D(qx[idx[0]],qy[idx[1]],qz[idx[2]])
                peak.setQSampleFrame(newQ)
                flightPath = peak.getL1() + peak.getL2()
                scatteringHalfAngle=0.5*peak.getScattering()
                tList.append(3176.507 * flightPath * np.sin(scatteringHalfAngle) / np.linalg.norm(newQ)) #convert to microseconds)
            tList = np.asarray(tList)
            peak.setQSampleFrame(origQS)
            return tList
            
        tList = getTList(peak, qx, qy, qz, boxMeanIDX)
    #Set up our bins for histogramming
    tMin = np.min(tList)
    tMax = np.max(tList)
    dt = tofPeak*dtSpread #time in us on either side of the peak position to consider
    dt = max(dt, 400)
    tMin = min(tMin, tofPeak - dt)
    tMax = max(tMax, tofPeak + dt)
    
    tBins = np.arange(tMin, tMax, dtBinWidth)
    weightList = n_events[hasEventsIDX] #- pp_lambda
    if removeEdges:
        mask = EdgeTools.getMask(peak, box, panelDict,qMask, edgesToCheck=edgesToCheck)
        print np.shape(mask), np.shape(useIDX), np.shape(useIDX[0])
        h = np.histogram(tList,tBins,weights=weightList*mask[hasEventsIDX]);
        '''
        plt.figure(2); plt.clf()
        q = mask[:,mask.shape[1]//2,:]
        r = n_events[:,n_events.shape[1]//2,:]
        plt.imshow(r)
        plt.hold('on')
        plt.imshow(q,cmap='gray',alpha=0.2)
        plt.savefig('/SNS/users/ntv/dropbox/si_removeEdges2/maskfigs/maskfig_%i.png'%peakNumber)
        '''
    else:
        h = np.histogram(tList,tBins,weights=weightList);
    #For and plot the TOF distribution
    tPoints = 0.5*(h[1][1:] + h[1][:-1])
    yPoints = h[0]

    if doVolumeNormalization:
        if calcTOFPerPixel == False:
            tofBox = 3176.507 * flightPath *np.sin(scatteringHalfAngle) * 1/np.sqrt(QX[qMask]**2 + QY[qMask]**2 + QZ[qMask]**2)
            #tofMask = np.ones_like(tofBox).astype(np.bool)
        if calcTOFPerPixel == True:
            qS0 = peak.getQSampleFrame()
            qqx = qx[::20] #Points we'll interpolate from for speed
            qqy = qy[::20]
            qqz = qz[::20]
            LLIST = np.zeros([qqx.size, qqy.size, qqz.size])
            for i, x in enumerate(qqx):
                for j, y in enumerate(qqy):
                    for k, z in enumerate(qqz):
                        qNew = V3D(x,y,z)
                        peak.setQSampleFrame(qNew)
                        L = peak.getL1() + peak.getL2()
                        HALFSCAT = 0.5*peak.getScattering()
                        LLIST[i,j,k] = L*np.sin(HALFSCAT)
            peak.setQSampleFrame(qS0)
            LLIST = np.asarray(LLIST)
            LLIST[LLIST>19.0] = np.nan
            QQX,QQY,QQZ = np.meshgrid(qqx,qqy,qqz,indexing='ij',copy=False)
            cartcoord = list(zip(QQX.flatten(), QQY.flatten(), QQZ.flatten()))
            nn = LinearNDInterpolator(cartcoord, LLIST.flatten())

            PIXELFACTOR = np.array(nn(QX,QY,QZ)) #=(L1*L2) + sin(theta)

            tofBox = 3176.507 * PIXELFACTOR * 1.0/np.sqrt(QX[qMask]**2 + QY[qMask]**2 + QZ[qMask]**2)
            tofMask = ~np.isnan(tofBox) 
        if removeEdges:
            print 'REMOVING EDGES'
            if not calcTOFPerPixel:
                numPixels = np.histogram((tofBox*mask[qMask]), tBins)[0]
            else:
                numPixels = np.histogram((tofBox*mask[qMask])[tofMask], tBins)[0]
        else:
            if not calcTOFPerPixel:
                numPixels = np.histogram(tofBox, tBins)[0]
            else:
                numPixels = np.histogram(tofBox[tofMask], tBins)[0]
            
        yPoints = 1.0*yPoints / numPixels
        useIDX = 1.0*numPixels/np.sum(numPixels) > minFracPixels
        if np.sum(useIDX < 1): #Bad threshold, we'll juse use it all
            useIDX = np.ones(np.size(yPoints)).astype(np.bool)
        tPoints = tPoints[useIDX]
        yPoints = yPoints[useIDX] * np.mean(numPixels)
        #if dtSpread > 0.04:
        #    plt.figure(8); plt.clf()
        #    plt.plot(tBins[1:], 1.0*numPixels/np.sum(numPixels))
        #    plt.figure(2)
    if workspaceNumber is None: 
        tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=yPoints, DataE=np.sqrt(yPoints))
    else:
        tofWS = CreateWorkspace(OutputWorkspace='tofWS%i'%workspaceNumber, DataX=tPoints, DataY=yPoints, DataE=np.sqrt(yPoints))
    return tofWS, float(pp_lambda)

#Determines the T0 shift for comparing moderator simulations (done at L=0)
# to our data (measured at L=L1+L2).  E is the neutron energy (units: eV)
# and L has units of m.  The returned value t0Shift has units of us.
def getT0Shift(E,L):
    #E is energy in eV, L is flight path in m
    mn = 1.674929e-27 #neutron mass, kg
    E = E * 1.60218e-19 #convert eV to J
    t0Shift = L*np.sqrt(mn/2/E) #units = s
    t0Shift = t0Shift * 1.0e6 #units =us
    return t0Shift

def getModeratorCoefficients(fileName):
    r = np.loadtxt(fileName)
    d = dict()
    d['A'] = r[0]
    d['B'] = r[1]
    d['R'] = r[2]
    d['T0'] = r[3]
    return d

def oneOverXSquared(x, A, bg):
    return A/np.sqrt(x) + bg

# d = calibration dictionary
def getInitialGuessByDetector(tofWS, paramNames, energy, flightPath, detNumber, calibrationDict):
    x0 = np.zeros(len(paramNames))
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    d = calibrationDict
    x0[0] = np.polyval(d['det_%i'%detNumber]['A'],energy)
    x0[1] = np.polyval(d['det_%i'%detNumber]['B'],energy)
    x0[2] = np.polyval(d['det_%i'%detNumber]['R'],energy)
    x0[3] = oneOverXSquared(energy, d['det_%i'%detNumber]['T0'][0], d['det_%i'%detNumber]['T0'][1])
    x0[4] = 1.0*(np.max(y)) 
    x0[5] = 0.5
    x0[6] = 200
    return x0

#Returns intial parameters for fitting based on a few quickly derived TOF
# profile parameters.  tofWS is a worskapce containng the TOF profile,
# paramNames is the list of parameter names
# energy is the energy of the peak (units: eV)
# flightPath is L = L1 + L2 (units: m)
def getInitialGuess(tofWS, paramNames, energy, flightPath, padeCoefficients,detNumber, calibDict):
    x0 = np.zeros(len(paramNames))
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    x0[0] = pade(padeCoefficients['A'], energy)
    x0[1] = pade(padeCoefficients['B'], energy)
    x0[2] = pade(padeCoefficients['R'], energy)
    #x0[3] = oneOverXSquared(energy, calibDict['det_%i'%detNumber]['T0'][0], calibDict['det_%i'%detNumber]['T0'][1])
    x0[3] = pade(padeCoefficients['T0'], energy) + getT0Shift(energy, flightPath) #extra is for ~18m downstream we are

    #These are still phenomenological
    x0[0] /= 1.2
    x0[2] += 0.05
    x0[3] -= 10 #This is lazy - we can do it detector-by-detector
    x0[3] = x[np.argmax(y)]
    x0[4] = (np.max(y))/x0[0]*2*2.5  #Amplitude
    x0[5] = 0.5 #hat width in IDX units
    x0[6] = 120.0 #Exponential decay rate for convolution
    return x0

#Get sample loads the NeXus evnts file and converts from detector space to
# reciprocal space.
# run is the run number.
# DetCalFile is a string for the file containng the detector calibration
# workDir is not used
# loadDir is the directory to extract the data from
def getSample(run, DetCalFile,  workDir, fileName, qLow=-25, qHigh=25):
    #data
    print 'Loading file', fileName
    data = Load(Filename = fileName)
    if DetCalFile is not None:
        LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic',
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1',
      MinValues = '%f, %f, %f'%(qLow, qLow, qLow), Maxvalues = '%f, %f, %f'%(qHigh, qHigh, qHigh), MaxRecursionDepth=10)
    return MDdata

def plotFitPresentation(filenameFormat, r,tofWS,fICC,runNumber, peakNumber, energy, chiSq,bgFinal, xStart, xStop, bgx0=None):
    plt.figure(1); plt.clf();

    import matplotlib
    import matplotlib.ticker as plticker
    ax = plt.gca()
    loc = plticker.MultipleLocator(base=50.0)
    ax.xaxis.set_major_locator(loc)

    font={'size':18,'family':'sans-serif'}
    matplotlib.rc('font',**font)
    plt.plot(r.readX(0),r.readY(0),'o',color=[30./256,104./256,255./256],ms=8,mew=2,zorder=10000)

    if bgx0 is not None:
        plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0))+np.polyval(bgx0, tofWS.readX(0)),'g',label='Moderator Prediction',lw=3)
    else:
        plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0)),color=[30./256,118./256,64./256],label='Initial Guess',lw=3)
        
    #plt.plot(r.readX(1),r.readY(1),'-',color=[0.3,0.3,0.7],label='Fit',lw=3)
    plt.plot(r.readX(1), np.polyval(bgFinal, r.readX(1)),color=[0.3,0.3,0.4],label='Background',lw=3)
    yLims = plt.ylim()
    #plt.plot([xStart, xStart], yLims, 'k')
    #plt.plot([xStop, xStop], yLims, 'k')
    #plt.title('E0=%4.4f meV, redChiSq=%4.4e'%(energy*1000,chiSq))
    plt.legend(loc='best')
    plt.xlabel('TOF ($\mu$s)')
    plt.ylabel('Counts')
    plt.savefig(filenameFormat%(runNumber, peakNumber))


#Function to make and save plots of the fits. bgx0=polynomial coefficients (polyfit order)
# for the initial guess
def plotFit(filenameFormat, r,tofWS,fICC,runNumber, peakNumber, energy, chiSq,bgFinal, xStart, xStop, bgx0=None):
    plt.figure(1); plt.clf()
    plt.plot(r.readX(0),r.readY(0),'o',label='Data')
    if bgx0 is not None:
        plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0))+np.polyval(bgx0, tofWS.readX(0)),'b',label='Initial Guess')
    else:
        plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0)),'b',label='Initial Guess')
        
    plt.plot(r.readX(1),r.readY(1),'.-',label='Fit')
    plt.plot(r.readX(1), np.polyval(bgFinal, r.readX(1)),'r',label='Background')
    yLims = plt.ylim()
    plt.plot([xStart, xStart], yLims, 'k')
    plt.plot([xStop, xStop], yLims, 'k')
    plt.title('E0=%4.4f meV, redChiSq=%4.4e'%(energy*1000,chiSq))
    plt.legend(loc='best')
    plt.savefig(filenameFormat%(runNumber, peakNumber))


def getRefinedCenter(peak, MDdata, UBMatrix, dQPixel,nPtsQ, neigh_length_m = 5, fracHKLSearch = 0.2):
    from mantid.kernel._kernel import V3D
    QSample = peak.getQSampleFrame()
    print QSample
    Qx = QSample[0]
    Qy = QSample[1]
    Qz = QSample[2]
    #dQ = np.abs(getDQFracHKL(UBMatrix,frac=fracHKLSearch))
    dQ = np.abs(getDQTOF(peak))
    oldWavelength = peak.getWavelength()
    Box = BinMD(InputWorkspace = 'MDdata',
        AlignedDim0='Q_sample_x,'+str(Qx-dQ[0,0])+','+str(Qx+dQ[0,1])+','+str(nPtsQ[0]),
        AlignedDim1='Q_sample_y,'+str(Qy-dQ[1,0])+','+str(Qy+dQ[1,1])+','+str(nPtsQ[1]),
        AlignedDim2='Q_sample_z,'+str(Qz-dQ[2,0])+','+str(Qz+dQ[2,1])+','+str(nPtsQ[2]),
        OutputWorkspace = 'MDbox')
    n_events = Box.getNumEventsArray()
    xaxis = Box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = Box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = Box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    meanList = list()
    neigh_length_m = 3
    maxBin = np.shape(n_events)
    hasEventsIDX = np.array(np.where((n_events>0) & ~np.isnan(n_events))).transpose()
    for i,idx in enumerate(hasEventsIDX):
        dataBox = n_events[max(idx[0] - neigh_length_m,0):min(idx[0] + neigh_length_m+1, maxBin[0]),
               max(idx[1] - neigh_length_m,0):min(idx[1] + neigh_length_m+1, maxBin[1]),
               max(idx[2] - neigh_length_m,0):min(idx[2] + neigh_length_m+1, maxBin[2])]
        meanList.append(np.mean(dataBox))
    peakIDX = hasEventsIDX[np.argmax(meanList)]
    qS = V3D()
    qS[0] = qx[peakIDX[0]]
    qS[1] = qy[peakIDX[1]]
    qS[2] = qz[peakIDX[2]]
    print qS
    peak.setQSampleFrame(qS)
    newWavelength = 4*np.pi*np.sin(peak.getScattering()/2.0)/np.sqrt(np.sum(np.power(peak.getQSampleFrame(),2)))
    peak.setWavelength(newWavelength)
    print 'Wavelength %4.4f --> %4.4f'%(oldWavelength, newWavelength)
    return np.array([qx[peakIDX[0]], qy[peakIDX[1]], qz[peakIDX[2]]])


#getBoxFracHKL returns the binned MDbox ranging from (hkl-0.5)-(hkl+0.5) (i.e. half integers 
# in hkl space) in q space.  
# Inputs:
#    peak: peak object to be analyzed.  HKL and peak centers must be defined
#    MDdata: the MD events workspace to be binned
#    latticeConstants: an array-type object containing lattice constants - note that
#            only the minmium number of constants is necessary for simpler systems (e.g.
#            one value for cubic crystals). 
#    crystalSystem: a string containng the crystal type.  Currently only works for 'cubic'
#    peakNumber: integer peak number within a  dataset - basically an index
#  Returns:
#  Box, an MDWorkspace with histogrammed events around the peak
def getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, dQPixel=0.005,fracHKL = 0.5, refineCenter=False, fracHKLRefine = 0.2):
    runNumber = peak.getRunNumber()
    QSample = peak.getQSampleFrame()
    Qx = QSample[0]
    Qy = QSample[1]
    Qz = QSample[2]
    #dQ = np.abs(getDQFracHKL(UBMatrix, frac = fracHKL))
    #print dQ
    #dQ = np.abs(getDQTOF(peak))
    #dQPixel = getPixelStep(peak)
    dQ = np.abs(dQ)
    dQ[dQ>0.5] = 0.5
    nPtsQ = np.round(np.sum(dQ/dQPixel,axis=1)).astype(int)
    print dQ, nPtsQ
    if refineCenter: #Find better center by flattining the cube in each direction and fitting a Gaussian

        #Get the new centers and new Box
        Qxr,Qyr,Qzr = getRefinedCenter(peak, MDdata, UBMatrix, dQPixel,nPtsQ, neigh_length_m = 5, fracHKLSearch = fracHKLRefine)

        Box = BinMD(InputWorkspace = 'MDdata',
            AlignedDim0='Q_sample_x,'+str(Qxr-dQ[0,0])+','+str(Qxr+dQ[0,1])+','+str(nPtsQ[0]),
            AlignedDim1='Q_sample_y,'+str(Qyr-dQ[1,0])+','+str(Qyr+dQ[1,1])+','+str(nPtsQ[1]),
            AlignedDim2='Q_sample_z,'+str(Qzr-dQ[2,0])+','+str(Qzr+dQ[2,1])+','+str(nPtsQ[2]),
            OutputWorkspace = 'MDbox')
    
    else: #We'll juse use the given center
        Box = BinMD(InputWorkspace = 'MDdata',
            AlignedDim0='Q_sample_x,'+str(Qx-dQ[0,0])+','+str(Qx+dQ[0,1])+','+str(nPtsQ[0]),
            AlignedDim1='Q_sample_y,'+str(Qy-dQ[1,0])+','+str(Qy+dQ[1,1])+','+str(nPtsQ[1]),
            AlignedDim2='Q_sample_z,'+str(Qz-dQ[2,0])+','+str(Qz+dQ[2,1])+','+str(nPtsQ[2]),
            OutputWorkspace = 'MDbox')
            #OutputWorkspace = 'MDbox_'+str(runNumber)+'_'+str(peakNumber))

    return Box


def doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,constraintScheme=None,nBG=15, outputWSName='fit', fitOrder=1):
    #Set up our inital guess
    fICC = ICC.IkedaCarpenterConvoluted()
    fICC.init()
    paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
    x0 = getInitialGuess(tofWS,paramNames,energy,flightPath,padeCoefficients,detNumber,calibrationDict)
    [fICC.setParameter(iii,v) for iii,v in enumerate(x0[:fICC.numParams()])]
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    if len(y)//2 < nBG: nBG = len(y)//2
    bgx0 = np.polyfit(x[np.r_[0:nBG,-nBG:0]], y[np.r_[0:nBG,-nBG:0]], fitOrder)

    nPts = x.size                
    scaleFactor = np.max((y-np.polyval(bgx0,x))[nPts//3:2*nPts//3])/np.max(fICC.function1D(x)[nPts//3:2*nPts//3])
    x0[4] = x0[4]*scaleFactor
    fICC.setParameter(4,x0[4])
    #fICC.setPenalizedConstraints(A0=[0.01, 1.0], B0=[0.005, 1.5], R0=[0.01, 1.0], T00=[0,1.0e10], k_conv0=[10,500],penalty=1.0e20)
    if constraintScheme == 1:
        try:
            fICC.setPenalizedConstraints(A0=[0.5*x0[0], 1.5*x0[0]], B0=[0.5*x0[1], 1.5*x0[1]], R0=[0.5*x0[2], 1.5*x0[2]], T00=[0,1.0e10],k_conv0=[5,500],penalty=1.0e20)
        except:
            fICC.setPenalizedConstraints(A0=[0.5*x0[0], 1.5*x0[0]], B0=[0.5*x0[1], 1.5*x0[1]], R0=[0.5*x0[2], 1.5*x0[2]], T00=[0,1.0e10],k_conv0=[5,500],penalty=None)
    if constraintScheme == 2:
        try:
            fICC.setPenalizedConstraints(A0=[0.01, 1.0], B0=[0.005, 1.5], R0=[0.00, 1.], scale0=[0.0, 1.0e10],T00=[0,1.0e10], k_conv0=[5.,500], penalty=1.0e20)
        except:
            fICC.setPenalizedConstraints(A0=[0.01, 1.0], B0=[0.005, 1.5], R0=[0.00, 1.], scale0=[0.0, 1.0e10], T00=[0,1.0e10], k_conv0=[5.,500], penalty=None)
    f = FunctionWrapper(fICC)
    bg = Polynomial(n=fitOrder)
    
    for i in range(fitOrder+1):
        bg['A'+str(fitOrder-i)] = bgx0[i]
    bg.constrain('-1.0 < A%i < 1.0'%fitOrder)
    fitFun = f + bg
    fitResults = Fit(Function=fitFun, InputWorkspace='tofWS', Output=outputWSName)
    return fitResults, fICC

#Does the actual integration and modifies the peaks_ws to have correct intensities.
def integrateSample(run, MDdata, peaks_ws, paramList, panelDict, UBMatrix, dQ, qMask, padeCoefficients, parameterDict, figsFormat=None, dtBinWidth = 4, nBG=15, dtSpread=0.02, fracHKL = 0.5, refineCenter=False, doVolumeNormalization=False, minFracPixels=0.0000, fracStop = 0.01, removeEdges=False, calibrationDict=None,dQPixel=0.005,calcTOFPerPixel=False, p=None,neigh_length_m=0,zBG=-1.0,bgPolyOrder=1, doIterativeBackgroundFitting=False):
    if removeEdges is True and panelDict is None:
        import sys
        sys.exit('ICCFT:integrateSample - trying to remove edges without a panelDict, this is impossible!')
    
    if p is None:
        p = range(peaks_ws.getNumberPeaks())
    fitDict = {}
    for i in p:
        peak = peaks_ws.getPeak(i)
        if peak.getRunNumber() == run:
            try:#for ppppp in [3]:#try:
                Box = getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel[0])
                tof = peak.getTOF() #in us
                wavelength = peak.getWavelength() #in Angstrom
                energy = 81.804 / wavelength**2 / 1000.0 #in eV
                flightPath = peak.getL1() + peak.getL2() #in m
                scatteringHalfAngle = 0.5*peak.getScattering()
                detNumber = 0#EdgeTools.getDetectorBank(panelDict, peak.getDetectorID())['bankNumber']
                print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
                if Box.getNEvents() < 1 or np.all(np.abs(peak.getHKL())==0):
                    print "Peak %i has 0 events or is HKL=000. Skipping!"%i
                    peak.setIntensity(0)
                    peak.setSigmaIntensity(1)
                    paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())]+[0])

                    mtd.remove('MDbox_'+str(run)+'_'+str(i))
                    continue
                n_events = Box.getNumEventsArray()
                goodIDX, pp_lambda = getBGRemovedIndices(n_events, peak=peak, box=Box,qMask=qMask[0], calc_pp_lambda=True, padeCoefficients=padeCoefficients, dtBinWidth=dtBinWidth,nBG=nBG)
                #Do background removal (optionally) and construct the TOF workspace for fitting
                if removeEdges:
                    edgesToCheck = EdgeTools.needsEdgeRemoval(Box,panelDict,peak) 
                    if edgesToCheck != []: #At least one plane intersects so we have to fit
                        tofWS,ppl = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[0], dtBinWidth=dtBinWidth,dtSpread=dtSpread[0], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=removeEdges, edgesToCheck=edgesToCheck, calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG,pp_lambda=pp_lambda)
                    else:
                        tofWS,ppl = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[0], dtBinWidth=dtBinWidth,dtSpread=dtSpread[0], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG,pp_lambda=pp_lambda)
                else:
                    #tofWS,pp_lambda = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[0], dtBinWidth=dtBinWidth,dtSpread=dtSpread[0], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG,pp_lambda=pp_lambda)
                    tofWS = mtd['tofWS'] # --IN PRINCIPLE!!! WE CALCULATE THIS BEFORE GETTING HERE
                    #TODO: Make sure we calculate it here - it seems to not work well for scolecute, but works for beta lac?

                print pp_lambda, 'is ppl'
                if doIterativeBackgroundFitting:
                    nBGToTry = range(2,tofWS.readX(0).size,4)
                    lowChiSq = 1.0e99
                    lowNBG = 0
                    for nBG in nBGToTry:
                        fitResults,fICC = doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,nBG=nBG,fitOrder=bgPolyOrder,constraintScheme=1) 
                        fitStatus = fitResults.OutputStatus
                        chiSq = fitResults.OutputChi2overDoF
                        if chiSq<lowChiSq:
                            lowChiSq = chiSq
                            lownBG = nBG
                        if chiSq < 2.0:
                            break
                else:
                   lownBG = nBG
                fitResults,fICC = doICCFit(tofWS, energy, flightPath, padeCoefficients, 0, None,nBG=lownBG,fitOrder=bgPolyOrder,constraintScheme=2)
                fitStatus = fitResults.OutputStatus
                chiSq = fitResults.OutputChi2overDoF

                #plt.close('all')
                #plt.figure(1); plt.clf()
                #plt.plot(mtd['fit_Workspace'].readX(0), mtd['fit_Workspace'].readY(0))
                #plt.plot(mtd['fit_Workspace'].readX(0), mtd['fit_Workspace'].readY(1))
                #plt.title('Chi Sq: %f, Peak Number: '%chiSq + str(i))
                #plt.pause(0.01)
            


 
                chiSq2  = 1.0e99
                if (chiSq > 1.0e99) and (tofWS.readY(0).max() > 10): #The initial fit isn't great - let's see if we can do better
                    Box = getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel[1])

                    if removeEdges:
                        if edgesToCheck != []: #At least one plane intersects so we have to fit
                            tofWS2,ppl2 = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[1], dtBinWidth=dtBinWidth,dtSpread=dtSpread[1], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=removeEdges, edgesToCheck=edgesToCheck, calcTOFPerPixel=calcTOFPerPixel, workspaceNumber=2,neigh_length_m=neigh_length_m,zBG=zBG)
                        else:
                            tofWS2,ppl2 = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[1], dtBinWidth=dtBinWidth,dtSpread=dtSpread[1], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel, workspaceNumber=2,neigh_length_m=neigh_length_m,zBG=zBG)
                    else:
                        tofWS2,ppl2 = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, qMask[1], dtBinWidth=dtBinWidth,dtSpread=dtSpread[1], doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,workspaceNumber=2,neigh_length_m=neigh_length_m,zBG=zBG)



                    print '############REFITTING########### on %4.4f'%chiSq
                    try:
                        fitResults2,fICC2 = doICCFit(tofWS2, energy, flightPath, padeCoefficients, detNumber, calibrationDict,nBG=nBG,outputWSName='fit2',fitOrder=bgPolyOrder,constraintScheme=2) 
                        fitStatus2 = fitResults2.OutputStatus
                        chiSq2 = fitResults2.OutputChi2overDoF
                    except: 
                        print 'CANNOT DO SECOND FIT, GOING BACK TO FIRST!!'
                if(chiSq < chiSq2):
                    r = mtd['fit_Workspace']
                    param = mtd['fit_Parameters']
                    tofWS = mtd['tofWS']
                else:
                    print 'USING SECOND FIT'
                    r = mtd['fit2_Workspace']
                    param = mtd['fit2_Parameters']
                    chiSq = chiSq2
                    fICC = fICC2
                    tofWS = mtd['tofWS2']
                    ppl = ppl2

                iii = fICC.numParams() - 1
                fitBG = [param.row(int(iii+bgIDX+1))['Value'] for bgIDX in range(bgPolyOrder+1)]


                #Set the intensity before moving on to the next peak
                icProfile = r.readY(1)
                bgCoefficients = fitBG[::-1]

                #peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))i
                t0 = param.row(3)['Value']
                
                #from scipy.ndimage import label
                #g = label(goodIDX*qMask)
                #bgPixels = np.sort(np.bincount(g[0].ravel()))[-2] #num bg pixels - -1 is 0s
                #bgIDX = np.argsort(np.bincount(g[0].ravel()))[-2] #num bg pixels - -1 is 0s

                 
                intensity, sigma, xStart, xStop = integratePeak(r.readX(0), icProfile,r.readY(0), np.polyval(bgCoefficients, r.readX(1)), pp_lambda=pp_lambda, fracStop=fracStop,totEvents=np.sum(n_events[goodIDX*qMask[0]]), bgEvents=np.sum(goodIDX*qMask[0])*pp_lambda*(neigh_length_m)**3/8)
                #print '~~~ ', intensity, sigma
                icProfile = icProfile - np.polyval(bgCoefficients, r.readX(1)) #subtract background
                peak.setIntensity(intensity)
                peak.setSigmaIntensity(sigma)
                if figsFormat is not None:
                    plotFit(figsFormat, r,tofWS,fICC,peak.getRunNumber(), i, energy, chiSq,fitBG, xStart, xStop, bgx0=None)
                    #plotFitPresentation('/SNS/users/ntv/med_peak.pdf', r, tofWS,fICC,peak.getRunNumber(), i, energy, chiSq,fitBG, xStart, xStop, bgx0)
                fitDict[i] = np.array([r.readX(0),r.readY(0), r.readY(1), r.readY(2)])
                paramList.append([i, energy, np.sum(icProfile), 0.0,chiSq] + [param.row(i)['Value'] for i in range(param.rowCount())]+[pp_lambda])
                if param.row(2)['Value'] < 0:
                    print i, [param.row(i)['Value'] for i in range(param.rowCount())]
                mtd.remove('MDbox_'+str(run)+'_'+str(i))
                
            except KeyboardInterrupt:
                print 'KeyboardInterrupt: Exiting Program!!!!!!!'
                sys.exit()
            except: #Error with fitting
                #raise
                import sys
                peak.setIntensity(0)
                peak.setSigmaIntensity(1)
                print 'Error with peak ' + str(i)
                paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(10)]+[0])
                #paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())]+[0])
                continue
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
        mtd.remove('MDbox_'+str(run)+'_'+str(i))
    return peaks_ws, paramList, fitDict


