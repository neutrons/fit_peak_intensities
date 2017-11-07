import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import ICCFitTools as ICCFT
from mantid.simpleapi import *
from scipy.interpolate import interp1d
from scipy.misc import factorial
from scipy.optimize import curve_fit
from scipy.ndimage.filters import convolve
from scipy.stats import multivariate_normal
import ICConvoluted as ICC
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)

#Example usage for bvg
# run getBox.py (or get the box object from ICCFT)
# params,h,th,ph = BVGFT.doBVGFit(box)
# BVGFT.compareBVGFitData(box,params[0])


def boxToTOFThetaPhi(box,peak):
    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    R, THETA, PHI = ICCFT.cart2sph(QX,QY,QZ)
    flightPath = peak.getL1() + peak.getL2()
    scatteringHalfAngle = 0.5*peak.getScattering()
    TOF = 3176.507 * flightPath *np.sin(scatteringHalfAngle) * 1/np.abs(R)
    X = np.empty(TOF.shape+(3,))
    X[:,:,:,0] = TOF
    X[:,:,:,1] = THETA
    X[:,:,:,2] = PHI
    return X

def fitScaling(n_events, YTOF, YBVG):
    YJOINT = 1.0*YTOF * YBVG
    YJOINT /= 1.0*YJOINT.sum()
#    goodIDX = n_events > 0.5*n_events.max()
    goodIDX = ICCFT.getBGRemovedIndices(n_events)[0]
    p0 = np.array([np.sum(goodIDX), 0.])
    p, cov = curve_fit(fitScalingFunction,YJOINT[goodIDX],n_events[goodIDX],sigma=np.sqrt(n_events[goodIDX]),p0=p0)
    #p = [0,0]
    #p[1] = np.mean(YJOINT[~goodIDX])
    #p[0] = np.mean(n_events[goodIDX]) / np.mean(YJOINT[goodIDX]-p[1] )
    print p
    YRET = p[0]*YJOINT  
    
    weights = 1.0*n_events.copy()
    weights[weights<1] = 1.
    chiSq = np.sum((YRET[goodIDX]-n_events[goodIDX])**2 / weights[goodIDX])
    chiSqRed = chiSq / (np.sum(goodIDX) - 2)
    print chiSqRed, 'is chiSqRed' 
    return YRET, chiSqRed
   
def fitScalingFunction(x,a,bg):
    return a*x #+ bg
 

def fitTOFCoordinate(box,peak, padeCoefficients,dtBinWidth=4,dtSpread=0.03,doVolumeNormalization=False,minFracPixels=0.01,removeEdges=False,calcTOFPerPixel=False,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,panelDict=None,qMask=None,calibrationDict=None,nBG=15):
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    detNumber = 0#EdgeTools.getDetectorBank(panelDict, peak.getDetectorID())['bankNumber']
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool) 
    tofWS,ppl = ICCFT.getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, panelDict, 0, qMask, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG)

    fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,nBG=nBG,fitOrder=bgPolyOrder)
    for i, param in enumerate(['A','B','R','T0','scale', 'hatWidth', 'k_conv']):
        fICC[param] = mtd['fit_Parameters'].row(i)['Value']
    print fICC, mtd['fit_Parameters'].row(2)
    tofxx = np.linspace(tofWS.readX(0).min(), tofWS.readX(0).max(),1000)
    tofyy = fICC.function1D(tofxx)
    plt.figure(1); plt.clf(); plt.plot(tofxx,tofyy)
    ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
    XTOF = boxToTOFThetaPhi(box,peak)[:,:,:,0]
    YTOF = ftof(XTOF)
    return YTOF
    

def fitPeak3D(box, X, n_events, peak,goodIDX):
    paramNames = [0 for i in range(7)]
    padeCoefficients = ICCFT.getModeratorCoefficients('franz_coefficients_2017.dat')
    flightPath = peak.getL1() + peak.getL2()
    halfScatteringAngle = 0.5*peak.getScattering()
    paramNames = [0 for i in range(7)]
    energy = peak.getInitialEnergy()/1000.0
    tofWS = CreateWorkspace(DataX = np.ones(50), DataY=np.ones(50))
    x0 = ICCFT.getInitialGuess(tofWS, paramNames, energy, flightPath, padeCoefficients,None, None)     
    alpha = x0[0]
    beta = x0[1]
    R = x0[2]
    T0 = x0[3]
    k_conv = 120
    params,h,t,p = doBVGFit(box,nTheta=400,nPhi=400)
    params, cov = params
    p0 = [np.max(n_events), np.max(n_events), np.max(n_events), params[1], params[2], params[3], params[4], params[5], alpha, beta, R, T0, k_conv, 0.0 ]
    bounds = ([0,-np.inf, -np.inf, -np.inf,-np.inf,0,0,-np.inf] + [0.5*v for v in x0[:4]] + [10, 0],
                [np.inf for i in range(8)] + [1.5*v for v in x0[:4]] + [500, np.inf]  )
    params= curve_fit(peak3DFitFunction, X[goodIDX], n_events[goodIDX],  p0,maxfev=1000, sigma=np.sqrt(n_events[goodIDX]),bounds=bounds)
    return params

def peak3DFitFunction(X, A, ATOF, ABVG, mu0, mu1, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg):
    return peak3D(X, A, ATOF, ABVG, mu0, mu1, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg)[0].ravel()
    
def peak3DFromParams(X,params):
    return peak3D(X, params[0],params[1],params[2],params[3],params[4],params[5],
                    params[6],params[7],params[8],params[9],params[10],params[11],
                    params[12],params[13])

def peak3D(X, A, ATOF, ABVG, mu0, mu1, sigX, sigY, p, alpha, beta, R, T0, k_conv, bg):
    #Axis 0 = TOF, axis1 = theta, axis2 = phi
    #First we calculate the TOF distribution to estimate sigma_TOF

    if X.ndim == 4:
        XTOF = X[:,:,:,0]
        XTHETA = X[:,:,:,1]
        XPHI = X[:,:,:,2]
        XANGLE = X[:,:,:,1:]
    elif X.ndim == 2:
        XTOF = X[:,0]
        XTHETA = X[:,1]
        XPHI = X[:,2]
        XANGLE = X[:,1:]
        
    else: 0/0 

    #Do the TOF Fit
    fICC = ICC.IkedaCarpenterConvoluted()
    fICC.init()
    fICC['A'] = alpha
    fICC['B'] = beta
    fICC['R'] = R
    fICC['T0'] = T0
    fICC['hatWidth'] = 0.5
    fICC['scale'] = 1/alpha
    fICC['k_conv'] = k_conv
    tofMin = np.min(XTOF)
    tofMax = np.max(XTOF)
    tofxx = np.linspace(tofMin, tofMax, 1000)
    tofyy = fICC.function1D(tofxx.ravel())
    ftof = interp1d(tofxx, tofyy)

    YTOF = ftof(XTOF)
    YTOF /= np.max(YTOF)


    #Do the bivariate normal for the angles
    while XANGLE.ndim < 0:
        XANGLE = np.expand_dims(XANGLE,axis=0)
    
    sigma = np.array([[sigX**2,p*sigX*sigY], [p*sigX*sigY,sigY**2]])
    mu = np.array([mu0,mu1])
    YBVG= bvg(ABVG, mu,sigma,XTHETA,XPHI,0)
    YBVG /= np.max(YBVG)
    #combine the results
    return A*YTOF*YBVG+bg, YTOF, YBVG


def getAngularHistogram(box, useIDX=None, nTheta=200, nPhi=200,zBG=1.96,neigh_length_m=3):
    n_events = box.getNumEventsArray()
    hasEventsIDX = n_events>0
    if zBG >=0:
        goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events)
    else: 
        goodIDX = hasEventsIDX

    useIDX = goodIDX

    #Setup our coordinates
    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    R, THETA, PHI = ICCFT.cart2sph(QX,QY,QZ)
    thetaMin = np.min(THETA); thetaMax = np.max(THETA)
    phiMin = np.min(PHI); phiMax = np.max(PHI) 
    thetaBins = np.linspace(thetaMin, thetaMax, nTheta)
    phiBins = np.linspace(phiMin, phiMax, nPhi)
    thetaVect = THETA[useIDX]
    phiVect = PHI[useIDX]
    nVect = n_events[useIDX]
   
    #Do the histogram 
    h, thBins, phBins = np.histogram2d(thetaVect, phiVect, weights=nVect, bins=[thetaBins,phiBins])
    return h,thBins, phBins 

def integrateBVGFit(Y,params):
    bg = params[0][-1]
    pIDX = Y > bg
    fitSum = np.sum(Y[pIDX])
    bgSum = bg*np.sum(pIDX) 
    intensity = fitSum - bgSum 
    sigma = np.sqrt(fitSum + bgSum)
    return intensity, sigma

def integrateBVGPeak(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, dQPixel=0.005,fracHKL = 0.5, refineCenter=False, fracHKLRefine = 0.2,nTheta=400,nPhi=400):
    box = getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel)
    params,h,t,p = doBVGFit(box,nTheta=nTheta,nPhi=nPhi)
    Y = getBVGResult(box, params[0],nTheta=nTheta,nPhi=nPhi)
    intens, sigma = integrateBVGFit(Y,params)
    print peak.getIntensity(), intens, sigma
    peak.setIntensity(intens)
    peak.setSigmaIntensity(sigma)

def getBVGResult(box, params,nTheta=200,nPhi=200):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi)
    thCenters = 0.5*(thBins[1:] + thBins[:-1])
    phCenters = 0.5*(phBins[1:] + phBins[:-1])
    TH, PH = np.meshgrid(thCenters, phCenters,indexing='ij',copy=False)
    Y = bvgFitFun([TH,PH],params[0],params[1],params[2],params[3],params[4],params[5],params[6])
    Y = Y.reshape([nTheta-1,nPhi-1])
    return Y

def compareBVGFitData(box,params,nTheta=200,nPhi=200,figNumber=2):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi)
    Y = getBVGResult(box,params,nTheta=nTheta,nPhi=nPhi)
    pLow = 0.3; pHigh = 0.7
    nX, nY = Y.shape
    plt.figure(figNumber); plt.clf()
    plt.subplot(2,2,1);
    plt.imshow(h,vmin=0,vmax=np.max(h),interpolation='None')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    plt.title('Measured Peak')
    plt.colorbar() 
    plt.subplot(2,2,2);
    plt.imshow(Y,vmin=0,vmax=np.max(h),interpolation='None' )
    plt.title('Modeled Peak')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    plt.colorbar() 
    plt.subplot(2,2,3);
    plt.imshow(h-Y,interpolation='None')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    plt.xlabel('Difference')
    plt.colorbar() 

def doBVGFit(box,nTheta=200, nPhi=200, zBG=1.96):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,zBG=zBG)
    dtH = np.mean(np.diff(thBins))
    dpH = np.mean(np.diff(phBins))
    thSum = np.sum(h,axis=1)
    phSum = np.sum(h,axis=0)
    thStd = np.std(thSum)*dtH
    phStd = np.std(phSum)*dpH
    thCenters = 0.5*(thBins[1:] + thBins[:-1])
    phCenters = 0.5*(phBins[1:] + phBins[:-1])
    TH, PH = np.meshgrid(thCenters, phCenters,indexing='ij',copy=False)
    
    params= curve_fit(bvgFitFun, [TH, PH], h.ravel(), p0=[np.max(h)/300., thCenters.mean(), phCenters.mean(), 0.005, 0.005, 0.05, 0.0] )
    #params= curve_fit(bvgFitFun, [TH, PH], h.ravel(), p0=[np.max(h)/300., thCenters.mean(), phCenters.mean(), thStd, phStd, 0.05, 0.0] )
    return params, h, thBins, phBins

def bvgFitFun(x, A, mu0, mu1,sigX,sigY,p,bg):
    sigma = np.array([[sigX**2,p*sigX*sigY], [p*sigX*sigY,sigY**2]])
    mu = np.array([mu0,mu1])
    return bvg(A, mu,sigma,x[0],x[1],bg).ravel() 

def bvg(A, mu,sigma,x,y,bg):
    pos = np.empty(x.shape+(2,))
    print pos.shape
    if pos.ndim == 4:
        pos[:,:,:,0] = x; pos[:,:,:,1] = y
    elif pos.ndim == 3:
        pos[:,:,0] = x; pos[:,:,1] = y
    else:
        pos[:,0] = x; pos[:,1] = y

    def is_pos_def(x):
        return np.all(np.linalg.eigvals(x) > 0)
    if is_pos_def(sigma):
        rv = multivariate_normal(mu, sigma)
        return A*rv.pdf(pos)+bg 
    else:
        print 'not PSD Matrix'
        return 0.0*np.ones_like(x)
