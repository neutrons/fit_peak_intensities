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

#Linear interpolation - only works for doubling samples in the range
def interpolateXGrid(X):
    XT = np.zeros(tuple(2*np.array(X.shape[:-1]))+(3,))
    XT[::2,::2,::2,:] = X
    #Now we need to interpolate - linear will be fast

    # -- single +1
    XT[1:-2:2,::2,::2,:] = (XT[:-2:2,::2,::2,:] + XT[2::2,::2,::2,:])*0.5
    XT[::2,1:-2:2,::2,:] = (XT[::2,:-2:2,::2,:] + XT[::2,2::2,::2,:])*0.5
    XT[::2,::2,1:-2:2,:] = (XT[::2,::2,:-2:2,:] + XT[::2,::2,2::2,:])*0.5

    # -- double +1
    XT[1:-2:2,1:-2:2,::2,:] = (XT[:-2:2,:-2:2,::2,:] + XT[2::2,2::2,::2,:])*0.5
    XT[1:-2:2,::2,1:-2:2] = (XT[:-2:2,::2,:-2:2,:] + XT[2::2,::2,2::2,:])*0.5
    XT[::2,1:-2:2,1:-2:2,:] = (XT[::2,:-2:2,:-2:2,:] + XT[::2,2::2,2::2,:])*0.5

    # -- triple +1
    XT[1:-1:2,1:-1:2,1:-1:2,:] = (XT[:-2:2,:-2:2,:-2:2,:] + XT[2::2,2::2,2::2,:])*0.5

    # --the last plane in each dimension - not interpolated, just a copy
    # TODO: real interpolation, but note that if you're integrating your last
    # column, you probably want to draw a bigger box
    XT[-1,:,:,:] = XT[-2,:,:,:]  
    XT[:,-1,:,:] = XT[:,-2,:,:]  
    XT[:,:,-1,:] = XT[:,:,-2,:]  
    return XT

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
    YJOINT /= 1.0*YJOINT.max() #Joint PDF - max is 1, integral is unknown
    #goodIDX = n_events > 0.5*n_events.max()
    #goodIDX = YJOINT > 0.5
    goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events)
    p0 = np.array([4*np.max(n_events), 0.])
    weights = np.sqrt(n_events)
    weights[weights<1] = 1.
    print YJOINT[goodIDX].shape, n_events[goodIDX].shape
    p, cov = curve_fit(fitScalingFunction,YJOINT[goodIDX],n_events[goodIDX],p0=p0)
    
    #highIDX = YJOINT > 0.7
    #p[0] = np.mean(n_events[highIDX] / YJOINT[highIDX])
    print p
    YRET = p[0]*YJOINT #+ p[1] 
    
    weights = 1.0*n_events.copy()
    weights[weights<1] = 1.
    chiSq = np.sum((YRET[goodIDX]-n_events[goodIDX])**2 / weights[goodIDX])
    chiSqRed = chiSq / (np.sum(goodIDX) - 2)
    print chiSqRed, 'is chiSqRed' 
    return YRET, chiSqRed, p[0]
   
def fitScalingFunction(x,a,bg):
    return a*x 

def xyzlinear(X, a0, a1, a2, a3, a4, a5, a6, a7):
    x = X[0]; y = X[1]; z = X[2];
    return a0*x*y*z + a1*x*y + a2*x*z + a3*y*z + a4*x + a5*y + a6*z + a7

def getXTOF(box, peak):
    from mantid.kernel import V3D
    QX,QY,QZ = ICCFT.getQXQYQZ(box) 
    origQS = peak.getQSampleFrame()
    tList = np.zeros_like(QX)
    for i in xrange(QX.shape[0]):
        print i
        for j in xrange(QX.shape[1]):
            for k in xrange(QX.shape[2]):
                newQ = V3D(QX[i,j,k],QY[i,j,k],QZ[i,j,k])
                peak.setQSampleFrame(newQ)
                flightPath = peak.getL1() + peak.getL2()
                scatteringHalfAngle=0.5*peak.getScattering()
                tList[i,j,k] = 3176.507 * flightPath * np.sin(scatteringHalfAngle) / np.linalg.norm(newQ) #convert to microseconds)
    peak.setQSampleFrame(origQS)
    return tList


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
    
    tofxx = np.linspace(tofWS.readX(0).min(), tofWS.readX(0).max(),1000)
    tofyy = fICC.function1D(tofxx)
    plt.figure(1); plt.clf(); plt.plot(tofxx,tofyy)
    plt.plot(tofWS.readX(0), tofWS.readY(0),'o')
    print 'sum:', np.sum(fICC.function1D(tofWS.readX(0)))
    plt.plot(mtd['fit_Workspace'].readX(1), mtd['fit_Workspace'].readY(1))
    ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
    XTOF = boxToTOFThetaPhi(box,peak)[:,:,:,0]
    #XTOF = getXTOF(box,peak)
    YTOF = ftof(XTOF)
    return YTOF, fICC, [tofWS.readX(0).min(), tofWS.readX(0).max()]

def getYTOF(fICC, XTOF, xlims):
    tofxx = np.linspace(xlims[0], xlims[1],10000)
    tofyy = fICC.function1D(tofxx)
    ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
    YTOF = ftof(XTOF)
    return YTOF

def getTOFParameters(box, peak, padeCoefficients,dtBinWidth=4,dtSpread=0.03,doVolumeNormalization=False,minFracPixels=0.01,removeEdges=False,calcTOFPerPixel=False,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,panelDict=None,qMask=None,calibrationDict=None,nBG=15):
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
    return fICC 

def fitPeak3D(box, X, n_events, peak,goodIDX,padeCoefficients):

    fICC = getTOFParameters(box, peak, padeCoefficients)
    alpha = fICC['A']; beta = fICC['B']; R = fICC['R']; T0 = fICC['T0']; k_conv = fICC['k_conv']; ATOF = fICC['scale']
    x0 = [alpha, beta, R, T0, ATOF, 0.5, k_conv]
    params,h,t,p = doBVGFit(box,nTheta=400,nPhi=400)
    params, cov = params
    p0 = [np.max(n_events), ATOF, params[0], params[1], params[2], params[3], params[4], params[5], alpha, beta, R, T0, k_conv, 0.0 ]
    bounds = ([0,0,0, -np.inf,-np.inf,0.,0.,-np.inf] + [0.9*v for v in x0[:4]] + [10, 0],
                [np.inf for i in range(8)] + [1.1*v for v in x0[:4]] + [500, np.inf]  )
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
    fICC['scale'] = ATOF
    fICC['k_conv'] = k_conv
    tofMin = np.min(XTOF)
    tofMax = np.max(XTOF)
    tofxx = np.linspace(tofMin, tofMax, 100000)
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


def getAngularHistogram(box, useIDX=None, nTheta=200, nPhi=200,zBG=1.96,neigh_length_m=3,fracBoxToHistogram=1.0):
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

    thetaMin = np.min(THETA); 
    thetaMax = np.max(THETA)
    dTheta = thetaMax - thetaMin
    thetaMid = 0.5*(thetaMin + thetaMax)
    thetaMin = max(thetaMin, thetaMid-dTheta*fracBoxToHistogram/2.0)
    thetaMax = min(thetaMax, thetaMid+dTheta*fracBoxToHistogram/2.0)


    phiMin = np.min(PHI); 
    phiMax = np.max(PHI) 
    dPhi = phiMax - phiMin
    phiMid = 0.5*(phiMin + phiMax)
    phiMin = max(phiMin, phiMid-dPhi*fracBoxToHistogram/2.0)
    phiMax = min(phiMax, phiMid+dPhi*fracBoxToHistogram/2.0)



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

def getBVGResult(box, params,nTheta=200,nPhi=200,fracBoxToHistogram=1.0):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram)
    thCenters = 0.5*(thBins[1:] + thBins[:-1])
    phCenters = 0.5*(phBins[1:] + phBins[:-1])
    TH, PH = np.meshgrid(thCenters, phCenters,indexing='ij',copy=False)
    Y = bvgFitFun([TH,PH],params[0],params[1],params[2],params[3],params[4],params[5],params[6])
    Y = Y.reshape([nTheta-1,nPhi-1])
    return Y

def compareBVGFitData(box,params,nTheta=200,nPhi=200,figNumber=2,fracBoxToHistogram=1.0):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram)
    Y = getBVGResult(box,params,nTheta=nTheta,nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram)
    pLow = 0.4; pHigh = 0.6
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
    #plt.imshow(Y,interpolation='None' )
    plt.title('Modeled Peak')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    plt.colorbar() 
    plt.subplot(2,2,3);
    plt.imshow(h/h.max()-Y/Y.max(),interpolation='None')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    plt.xlabel('Normed Difference')
    plt.colorbar() 

def doBVGFit(box,nTheta=200, nPhi=200, zBG=1.96, fracBoxToHistogram=1.0):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,zBG=zBG,fracBoxToHistogram=fracBoxToHistogram)
    dtH = np.mean(np.diff(thBins))
    dpH = np.mean(np.diff(phBins))
    thCenters = 0.5*(thBins[1:] + thBins[:-1])
    phCenters = 0.5*(phBins[1:] + phBins[:-1])
    TH, PH = np.meshgrid(thCenters, phCenters,indexing='ij',copy=False)
    fitIDX = np.ones_like(h).astype(np.bool)
    weights = np.sqrt(h)
    weights[weights<1] = 1 
    params= curve_fit(bvgFitFun, [TH[fitIDX], PH[fitIDX]], h[fitIDX].ravel(), p0=[np.max(h)*50., TH.mean(), PH.mean(), 0.005, 0.005, 0.0005,0.0])

    #params= curve_fit(bvgFitFun, [TH, PH], h.ravel(), p0=[np.max(h)/300., thCenters.mean(), phCenters.mean(), thStd, phStd, 0.05, 0.0] )
    return params, h, thBins, phBins

def bvgFitFun(x, A, mu0, mu1,sigX,sigY,p,bg):
    sigma = np.array([[sigX**2,p*sigX*sigY], [p*sigX*sigY,sigY**2]])
    mu = np.array([mu0,mu1])
    return bvg(A, mu,sigma,x[0],x[1],bg).ravel() 

def bvg(A, mu,sigma,x,y,bg):
    pos = np.empty(x.shape+(2,))
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
        return A*rv.pdf(pos) +bg
    else:
        print '   BVGFT:bvg:not PSD Matrix'
        return 0.0*np.ones_like(x)
