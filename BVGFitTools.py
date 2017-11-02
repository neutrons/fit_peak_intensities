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

def fitPeak3D(X, n_events, peak,goodIDX):
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
    p0 = [np.max(n_events), np.mean(X[:,:,:,1]), np.mean(X[:,:,:,2]), 0.005, 0.005, 0.05, alpha, beta, R, T0, k_conv, 0.0 ]
    bounds = ([0,-np.inf,-np.inf,0,0,-np.inf] + [0.5*v for v in x0[:4]] + [10, 0],
                [np.inf for i in range(6)] + [1.5*v for v in x0[:4]] + [500, np.inf]  )
    params= curve_fit(peak3DFitFunction, X[goodIDX], n_events[goodIDX],  p0,maxfev=1000, sigma=np.sqrt(n_events[goodIDX]),bounds=bounds)
    return params

def peak3DFitFunction(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg):
    return peak3D(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg).ravel()
    
def peak3DFromParams(X,params):
    return peak3D(X, params[0],params[1],params[2],params[3],params[4],params[5],
                    params[6],params[7],params[8],params[9],params[10],params[11])

def peak3D(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg):
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
        
    else: return 

    fICC = ICC.IkedaCarpenterConvoluted()
    fICC.init()
    fICC['A'] = alpha
    fICC['B'] = beta



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
    plt.figure(figNumber); plt.clf()
    plt.subplot(2,2,1);
    plt.imshow(h,vmin=0,vmax=np.max(h))
    plt.title('Measured Peak')
    plt.subplot(2,2,2);
    plt.imshow(Y,vmin=0,vmax=np.max(h))
    plt.title('Modeled Peak')
    plt.subplot(2,2,3);
    plt.imshow(h-Y)
    

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
    return bvg(A, mu,sigma,x[0],x[1],bg) 

def bvg(A, mu,sigma,x,y,bg):
    pos = np.empty(x.shape+(2,))
    pos[:,:,0] = x; pos[:,:,1] = y
    rv = multivariate_normal(mu, sigma)
    return A*rv.pdf(pos).ravel() + bg
