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


def get3DPeak(peak, box, padeCoefficients, qMask, nTheta=150, nPhi=150,fracBoxToHistogram=1.0,numTimesToInterpolate=1, plotResults=False, dtBinWidth=4,zBG=1.96,bgPolyOrder=1, fICCParams = None, oldICCFit=None, strongPeakParams=None, forceCutoff=250, edgeCutoff=15, predCoefficients=None, neigh_length_m=3, q_frame = 'sample', dtSpread=0.03, pplmin_frac=0.8, pplmax_frac=1.5, mindtBinWidth=1):
    n_events = box.getNumEventsArray()

    if q_frame == 'lab':
        q0 = peak.getQLabFrame()
    elif q_frame == 'sample':
        q0 = peak.getQSampleFrame()
    else:
        raise ValueError('BVGFT:get3DPeak - q_frame must be either \'lab\' or \'sample\'; %s was provided'%q_frame)


    if fICCParams is None:
        goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events, peak=peak, box=box,qMask=qMask, calc_pp_lambda=True, padeCoefficients=padeCoefficients, dtBinWidth=dtBinWidth, predCoefficients=predCoefficients,neigh_length_m=neigh_length_m, pp_lambda=None, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)
        YTOF, fICC, x_lims = fitTOFCoordinate(box,peak,padeCoefficients,dtSpread=dtSpread,dtBinWidth=dtBinWidth,qMask=qMask,bgPolyOrder=bgPolyOrder,zBG=zBG,plotResults=plotResults, pp_lambda=pp_lambda, neigh_length_m=neigh_length_m, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

    else: #we already did I-C profile, so we'll just read the parameters
        pp_lambda = fICCParams[-1]
        fICC = ICC.IkedaCarpenterConvoluted()
        fICC.init()
        fICC['A'] = fICCParams[5]  
        fICC['B'] = fICCParams[6]  
        fICC['R'] = fICCParams[7]  
        fICC['T0'] = fICCParams[8]  
        fICC['scale'] = fICCParams[9]  
        fICC['hatWidth'] = fICCParams[10]  
        fICC['k_conv'] = fICCParams[11] 
        goodIDX, _ = ICCFT.getBGRemovedIndices(n_events, pp_lambda=pp_lambda, qMask=qMask)

        if oldICCFit is not None:
            x_lims = [np.min(oldICCFit[0]), np.max(oldICCFit[0])]
            tofxx = oldICCFit[0]; tofyy = oldICCFit[2]
        else:
            dtSpread = 0.03
            x_lims = [(1-dtSpread)*peak.getTOF(), (1+dtSpread)*peak.getTOF()]
            tofxx = np.arange(x_lims[0], x_lims[1], 5)
            tofyy = fICC.function1D(tofxx)
        ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
        XTOF = boxToTOFThetaPhi(box,peak)[:,:,:,0]
        YTOF = ftof(XTOF)
    goodIDX *= qMask #TODO: we can do this when we get the goodIDX
    X = boxToTOFThetaPhi(box,peak)
    dEdge = edgeCutoff
    useForceParams = peak.getIntensity() < forceCutoff or peak.getRow() <= dEdge or peak.getRow() >= 255-dEdge or peak.getCol() <= dEdge or peak.getCol() >= 255-dEdge
    if strongPeakParams is not None and useForceParams:
        th = np.arctan2(q0[1],q0[0])
        ph = np.arctan2(q0[2],np.hypot(q0[0],q0[1]))
        thphPeak = np.array([th,ph])
        nnIDX = np.argmin(np.linalg.norm(strongPeakParams[:,:2] - thphPeak,axis=1))
        print 'Using ', strongPeakParams[nnIDX,:2], 'for ', thphPeak
        params,h,t,p = doBVGFit(box,nTheta=nTheta,nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram, goodIDX=goodIDX, forceParams=strongPeakParams[nnIDX])
    else:
        params,h,t,p = doBVGFit(box,nTheta=nTheta,nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram, goodIDX=goodIDX)
    
    if plotResults:
        compareBVGFitData(box,params[0],nTheta, nPhi, fracBoxToHistogram=fracBoxToHistogram, useIDX=goodIDX)

    A = params[0][0]
    mu0 = params[0][1]
    mu1 = params[0][2]
    sigX = params[0][3]
    sigY = params[0][4]
    p = params[0][5]
    bgBVG = params[0][6]
    sigma = np.array([[sigX**2,p*sigX*sigY], [p*sigX*sigY,sigY**2]])
    mu = np.array([mu0,mu1])
    
    retParams = {}
    retParams['Alpha'] = fICC['A']
    retParams['Beta'] = fICC['B']
    retParams['R'] = fICC['R']
    retParams['T0'] = fICC['T0']
    retParams['scale'] = fICC['scale']
    retParams['k_conv'] = fICC['k_conv']
    retParams['muTH'] = mu0
    retParams['muPH'] = mu1
    retParams['sigX'] = sigX
    retParams['sigY'] = sigY
    retParams['sigP'] = p
    retParams['bgBVG'] = bgBVG

    XTOF = X[:,:,:,0]
    XTHETA = X[:,:,:,1]
    XPHI = X[:,:,:,2]
    XANGLE = X[:,:,:,1:]

    YBVG= bvg(1.0, mu,sigma,XTHETA,XPHI,0)
    Y,redChiSq, scaleFactor = fitScaling(n_events,box, YTOF, YBVG)
    retParams['scale3d'] = scaleFactor
    retParams['chiSq3d'] = redChiSq

    for i in range(numTimesToInterpolate):
        X = interpolateXGrid(X)
    
    XTOF = X[:,:,:,0]
    XTHETA = X[:,:,:,1]
    XPHI = X[:,:,:,2]
    XANGLE = X[:,:,:,1:]
    YBVG2 = bvg(1.0, mu,sigma,XTHETA,XPHI,0)
    YTOF2 = getYTOF(fICC, XTOF, x_lims)
    Y2 = YTOF2*YBVG2
    Y2 = scaleFactor*Y2/Y2.max() 
    
    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    fitMaxIDX = tuple(np.array(np.unravel_index(Y2.argmax(), Y2.shape)) // (numTimesToInterpolate+1))
    newCenter = np.array([QX[fitMaxIDX], QY[fitMaxIDX], QZ[fitMaxIDX]])

    retParams['dQ'] =  np.linalg.norm(newCenter - q0)
    retParams['newQ'] =  newCenter 

    return Y2, goodIDX, pp_lambda, retParams

def integrateJointPDF(Y, pp_lambda, scale, YCutoff=0.05):
    YPDF = Y/Y.max()
    goodIDX = YPDF > YCutoff
    intensity = np.sum(Y[goodIDX]-pp_lambda)/(1.0*scale)
    bg = pp_lambda*np.sum(goodIDX)
    sigma = np.sqrt(intensity + 2.0*bg) # = sqrt(OBS + BG)
    return intensity, sigma

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


def calcTOFPerPixel(box, peak, refitIDX = None):
    xaxis = box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    QX, QY, QZ = ICCFT.getQXQYQZ(box)

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

def fitScaling(n_events,box, YTOF, YBVG, goodIDX=None, neigh_length_m=3):

    YJOINT = 1.0*YTOF * YBVG
    YJOINT /= 1.0*YJOINT.max() #Joint PDF - max is 1, integral is unknown


    convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
    convYJOINT = convolve(YJOINT, convBox)
    conv_n_events = convolve(n_events, convBox)

    #goodIDX = n_events > -1.0
    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    dP = 8
    fitMaxIDX = tuple(np.array(np.unravel_index(YJOINT.argmax(), YJOINT.shape)))
    if goodIDX is None:
        goodIDX = np.zeros_like(YJOINT).astype(np.bool)
        goodIDX[max(fitMaxIDX[0]-dP,0):min(fitMaxIDX[0]+dP,goodIDX.shape[0]), 
                max(fitMaxIDX[1]-dP,0):min(fitMaxIDX[1]+dP,goodIDX.shape[1]),
                max(fitMaxIDX[2]-dP,0):min(fitMaxIDX[2]+dP,goodIDX.shape[2])] = True
        goodIDX = np.logical_and(goodIDX, conv_n_events>0)
        #goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events)
        #goodIDX = np.logical_and(goodIDX, n_events>0)

    p0 = np.array([np.max(n_events), np.mean(n_events)])
    weights = np.sqrt(n_events).copy()
    weights[weights<1] = 1.
    bounds = ([0,0],[np.inf,np.inf])
    p, cov = curve_fit(fitScalingFunction,convYJOINT[goodIDX],conv_n_events[goodIDX],p0=p0, bounds=bounds)#, sigma=np.sqrt(weights[goodIDX]))
   
    #highIDX = YJOINT > 0.7
    #p[0] = np.mean(n_events[highIDX] / YJOINT[highIDX])
    YRET = p[0]*YJOINT + p[1] 
    
    chiSq = np.sum((YRET[goodIDX]-n_events[goodIDX])**2 / weights[goodIDX]**2)
    chiSqRed = chiSq / (np.sum(goodIDX) - 2)
    print chiSqRed, 'is chiSqRed'
    return YRET, chiSqRed, p[0]
   
def fitScalingFunction(x,a,bg):
    return a*x + bg

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


def fitTOFCoordinate(box,peak, padeCoefficients,dtBinWidth=4,dtSpread=0.03,doVolumeNormalization=False,minFracPixels=0.01,removeEdges=False,calcTOFPerPixel=False,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,panelDict=None,qMask=None,calibrationDict=None, plotResults=False,fracStop=0.01, pp_lambda=None, pplmin_frac=0.8, pplmax_frac=1.5, mindtBinWidth=1):
    
    #Get info from the peak
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    detNumber = 0#EdgeTools.getDetectorBank(panelDict, peak.getDetectorID())['bankNumber']

    #Set the qMask
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool) 
    
    #Calculate the optimal pp_lambda and 
    tofWS,ppl = ICCFT.getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, panelDict, qMask, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG,pp_lambda=pp_lambda, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

    try:
        fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,fitOrder=bgPolyOrder,constraintScheme=1)
    except:
        fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,fitOrder=bgPolyOrder)
    
    for i, param in enumerate(['A','B','R','T0','scale', 'hatWidth', 'k_conv']):
        fICC[param] = mtd['fit_Parameters'].row(i)['Value']
    bgParamsRows = [7 + i for i in range(bgPolyOrder+1)]
    bgCoeffs = []
    for bgRow in bgParamsRows[::-1]:#reverse for numpy order 
        bgCoeffs.append(mtd['fit_Parameters'].row(bgRow)['Value'])
    #print bgCoeffs
    x = tofWS.readX(0)
    bg = np.polyval(bgCoeffs, x)
    yFit = mtd['fit_Workspace'].readY(1)
    yData = tofWS.readY(0)
    
    yScaled = (yFit-bg) / np.max(yFit-bg)
    goodIDX = yScaled > fracStop
    if np.sum(goodIDX) > 0:
        iStart = np.min(np.where(goodIDX))
        iStop = np.max(np.where(goodIDX))

   
    interpF = interp1d(x, yFit, kind='cubic')
    tofxx = np.linspace(tofWS.readX(0).min(), tofWS.readX(0).max(),1000)
    tofyy = interpF(tofxx)
    if plotResults:
        plt.figure(1); plt.clf(); 
        plt.plot(tofxx,tofyy,label='Interpolated')
        plt.plot(tofWS.readX(0), tofWS.readY(0),'o',label='Data')
        print 'sum:', np.sum(fICC.function1D(tofWS.readX(0)))
        print 'bg: ', np.sum(bg[iStart:iStop])
        plt.plot(mtd['fit_Workspace'].readX(1), mtd['fit_Workspace'].readY(1),label='Fit')
        plt.title(fitResults.OutputChi2overDoF)
        plt.legend(loc='best')
    ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
    XTOF = boxToTOFThetaPhi(box,peak)[:,:,:,0]
    YTOF = ftof(XTOF)
    return YTOF, fICC, [tofWS.readX(0).min(), tofWS.readX(0).max()]

def getYTOF(fICC, XTOF, xlims):
    tofxx = np.linspace(xlims[0], xlims[1],1000)
    tofyy = fICC.function1D(tofxx)
    ftof = interp1d(tofxx, tofyy,bounds_error=False,fill_value=0.0)
    YTOF = ftof(XTOF)
    return YTOF

def getTOFParameters(box, peak, padeCoefficients,dtBinWidth=4,dtSpread=0.03,doVolumeNormalization=False,minFracPixels=0.01,removeEdges=False,calcTOFPerPixel=False,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,panelDict=None,qMask=None,calibrationDict=None, pplmin_frac=0.8, pplmax_frac=1.5):
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    detNumber = 0#EdgeTools.getDetectorBank(panelDict, peak.getDetectorID())['bankNumber']
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool)
    tofWS,ppl = ICCFT.getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, panelDict, qMask, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False,calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac)

    fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, detNumber, calibrationDict,fitOrder=bgPolyOrder)
    for i, param in enumerate(['A','B','R','T0','scale', 'hatWidth', 'k_conv']):
        fICC[param] = mtd['fit_Parameters'].row(i)['Value']
    return fICC 

#This is not a very recent version of this
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
    if useIDX is None:
        if zBG >=0:
            goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events)
            #goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events, peak=peak, box=box,qMask=qMask, calc_pp_lambda=True, padeCoefficients=padeCoefficients, dtBinWidth=dtBinWidth)
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

def integrateBVGPeak(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, dQPixel=0.005,fracHKL = 0.5, refineCenter=False, fracHKLRefine = 0.2,nTheta=400,nPhi=400, q_frame='sample'):
    box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel,q_frame=q_frame)
    params,h,t,p = doBVGFit(box,nTheta=nTheta,nPhi=nPhi)
    Y = getBVGResult(box, params[0],nTheta=nTheta,nPhi=nPhi)
    intens, sigma = integrateBVGFit(Y,params)
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

def compareBVGFitData(box,params,nTheta=200,nPhi=200,figNumber=2,fracBoxToHistogram=1.0, useIDX=None):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram, useIDX=useIDX)
    Y = getBVGResult(box,params,nTheta=nTheta,nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram)
    pLow = 0.0; pHigh = 1.0
    nX, nY = Y.shape
    plt.figure(figNumber); plt.clf()
    plt.subplot(2,2,1);
    plt.imshow(h,vmin=0,vmax=0.7*np.max(h),interpolation='None')
    plt.xlim([pLow*nX, pHigh*nX])
    plt.ylim([pLow*nY, pHigh*nY])
    if useIDX is None: plt.title('Measured Peak')
    else: plt.title('BG Removed Measured Peak')
    plt.colorbar() 
    plt.subplot(2,2,2);
    plt.imshow(Y,vmin=0,vmax=0.7*np.max(h),interpolation='None' )
    #plt.imshow(Y,interpolation='None' )
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

    if useIDX is not None:
        h0, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,fracBoxToHistogram=fracBoxToHistogram, useIDX=None)
        plt.subplot(2,2,4);
        plt.imshow(h0,vmin=0,vmax=1.0*np.max(h0),interpolation='None')
        plt.xlim([pLow*nX, pHigh*nX])
        plt.ylim([pLow*nY, pHigh*nY])
        plt.xlabel('Measured Peak')
        plt.colorbar() 





def doBVGFit(box,nTheta=200, nPhi=200, zBG=1.96, fracBoxToHistogram=1.0, goodIDX=None, forceParams=None, forceTolerance = 0.1, dph=10, dth=10):
    h, thBins, phBins = getAngularHistogram(box, nTheta=nTheta, nPhi=nPhi,zBG=zBG,fracBoxToHistogram=fracBoxToHistogram,useIDX=goodIDX)
    dtH = np.mean(np.diff(thBins))
    dpH = np.mean(np.diff(phBins))
    thCenters = 0.5*(thBins[1:] + thBins[:-1])
    phCenters = 0.5*(phBins[1:] + phBins[:-1])
    TH, PH = np.meshgrid(thCenters, phCenters,indexing='ij',copy=False)

    fitIDX = np.ones_like(h).astype(np.bool)
    weights = np.sqrt(h)
    weights[weights<1] = 1

    def fSigX(x,a,k,x0,b):
        return a*np.exp(-k*(x-x0)) + b
        #return a/(x-x0) + b

    def fSigP(x,a,k,phi,b):
        return a*np.sin(k*x-phi) + b*x 


    if forceParams is None:
            meanTH = TH.mean()
            meanPH = PH.mean()
            #sigX0 = np.polyval([ 0.00173264,  0.0002208 ,  0.00185031, -0.00012078,  0.00189967], meanPH)
            sigX0 = ICCFT.oldScatFun(meanPH, 1.71151521e-02,   6.37218400e+00,   3.39439675e-03)
            sigY0 = 0.002#np.polyval([ 0.00045678, -0.0017618 ,  0.0045013 , -0.00480677,  0.00376619], meanTH)
            sigP0 = fSigP(meanTH,  0.1460775 ,  1.85816592,  0.26850086, -0.00725352) 
            p0=[np.max(h), meanTH, meanPH, sigX0, sigY0, sigP0,0.0]
            print p0 
            bounds = ([0.0, thBins[thBins.size//2 - 2], phBins[phBins.size//2 - 2], 0.7*sigX0, 0.000, -0.4, 0], 
                    [np.inf, thBins[thBins.size//2 + 2], phBins[phBins.size//2 + 2], 1.3*sigX0, 0.007, 0.4, np.inf])
            params= curve_fit(bvgFitFun, [TH[fitIDX], PH[fitIDX]], h[fitIDX].ravel(),p0=p0, bounds=bounds, sigma=np.sqrt(weights[fitIDX]))
            #params= curve_fit(bvgFitFun, [TH[fitIDX], PH[fitIDX]], h[fitIDX].ravel(),p0=p0, sigma=np.sqrt(weights[fitIDX]))
            #params= curve_fit(bvgFitFun, [TH[fitIDX], PH[fitIDX]], h[fitIDX].ravel(),p0=p0)
            print '!!!', params[0]

    elif forceParams is not None:
        p0 = np.zeros(7)
        p0[0] = np.max(h); p0[1] = TH.mean(); p0[2] = PH.mean()
        p0[3] = forceParams[5]; p0[4] = forceParams[6]; p0[5] = forceParams[7];
        isPos = np.sign(p0)
        bounds = ((1.0-isPos*forceTolerance)*p0, (1.0+isPos*forceTolerance)*p0)
        bounds[0][0] = 0.0;  bounds[1][0] = np.inf; #Amplitude
        bounds[0][1] = min(thBins[thBins.size//2 - dth], thBins[thBins.size//2 + dth]);
        bounds[1][1] = max(thBins[thBins.size//2 - dth], thBins[thBins.size//2 + dth]);
        bounds[0][2] = min(phBins[phBins.size//2 - dph], phBins[phBins.size//2 + dph]);
        bounds[1][2] = max(phBins[phBins.size//2 - dph], phBins[phBins.size//2 + dph]);
        bounds[1][-1] = np.inf
        print 'forcing: ', forceParams
        print '~~~ p0:',p0
        params = curve_fit(bvgFitFun, [TH[fitIDX], PH[fitIDX]], h[fitIDX],p0=p0, bounds=bounds)
        print '~~ params:', params[0]
    return params, h, thBins, phBins

def bvgFitFun(x, A, mu0, mu1,sigX,sigY,p,bg):
    sigma = np.array([[sigX**2,p*sigX*sigY], [p*sigX*sigY,sigY**2]])
    mu = np.array([mu0,mu1])
    return bvg(A, mu,sigma,x[0],x[1],bg).ravel() 

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)
 
def bvg(A, mu,sigma,x,y,bg):
    pos = np.empty(x.shape+(2,))
    if pos.ndim == 4:
        pos[:,:,:,0] = x; pos[:,:,:,1] = y
    elif pos.ndim == 3:
        pos[:,:,0] = x; pos[:,:,1] = y
    else:
        pos[:,0] = x; pos[:,1] = y

    if is_pos_def(sigma):
        rv = multivariate_normal(mu, sigma)
        return A*rv.pdf(pos) +bg
    else:
        print '   BVGFT:bvg:not PSD Matrix'
        return 0.0*np.ones_like(x)
