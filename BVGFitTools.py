import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import ICCFitTools as ICCFT
from mantid.simpleapi import *
from scipy.interpolate import interp1d
from scipy.misc import factorial
from scipy.optimize import curve_fit
from scipy.ndimage.filters import convolve
from matplotlib.mlab import bivariate_normal
import ICConvoluted as ICC
import mbvg
reload(mbvg)
FunctionFactory.subscribe(mbvg.MBVG)
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)


class mvn():
    '''
    This class is a wrapper for the matplotlib.mlab.bivariate_normal
    implementation of a bivariate Gaussian.  It is designed to be transparent
    with scipy.stats.multivariate_normal, but can be used when older versions
    of scipy are available (e.g. on the analysis machine).
    '''
    def __init__(self, mu, sigma):
        self.mux = mu[0]
        self.muy = mu[1]
        self.sigx = np.sqrt(sigma[0][0])
        self.sigy = np.sqrt(sigma[1][1])
        self.sigxy = sigma[0][1]
        self.p = 1.0*self.sigxy/self.sigx/self.sigy
    
    def __str__(self):
        return 'bivariate guassian with mu=[%4.4f, %4.4f] and sigma = [%4.4f %4.4f; %4.4f %4.4f]'%(mux, muy, sigx, sigxy, sigxy, sigy)
        
    def __repr__(self):
        return 'bivariate guassian with mu=[%4.4f, %4.4f] and sigma = [%4.4f %4.4f; %4.4f %4.4f]'%(mux, muy, sigx, sigxy, sigxy, sigy)
    def pdf(self, pos):
        X = pos[...,0]
        Y = pos[...,1]
        return bivariate_normal(X,Y, sigmax=self.sigx, sigmay=self.sigy,
                            mux=self.mux,muy=self.muy,sigmaxy=self.sigxy)

multivariate_normal = mvn

def get3DPeak(peak, box, padeCoefficients, qMask, nTheta=150, nPhi=150,fracBoxToHistogram=1.0, plotResults=False, zBG=1.96,bgPolyOrder=1, fICCParams = None, oldICCFit=None, strongPeakParams=None, forceCutoff=250, edgeCutoff=15, predCoefficients=None, neigh_length_m=3, q_frame = 'sample', dtSpread=0.03, pplmin_frac=0.8, pplmax_frac=1.5, mindtBinWidth=1):
    n_events = box.getNumEventsArray()

    if q_frame == 'lab':
        q0 = peak.getQLabFrame()
    elif q_frame == 'sample':
        q0 = peak.getQSampleFrame()
    else:
        raise ValueError('BVGFT:get3DPeak - q_frame must be either \'lab\' or \'sample\'; %s was provided'%q_frame)


    if fICCParams is None:
        goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events, peak=peak, box=box,qMask=qMask, calc_pp_lambda=True, padeCoefficients=padeCoefficients, predCoefficients=predCoefficients,neigh_length_m=neigh_length_m, pp_lambda=None, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)
        YTOF, fICC, x_lims = fitTOFCoordinate(box,peak,padeCoefficients,dtSpread=dtSpread,qMask=qMask,bgPolyOrder=bgPolyOrder,zBG=zBG,plotResults=plotResults, pp_lambda=pp_lambda, neigh_length_m=neigh_length_m, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

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
        ph = np.arctan2(q0[1],q0[0])
        th = np.arctan2(q0[2],np.hypot(q0[0],q0[1])) 
        thphPeak = np.array([ph,th])
        tmp = strongPeakParams[:,:2] - thphPeak
        dist = np.sqrt(tmp[:,0]**2 + tmp[:,1]**2)
        nnIDX = np.argmin(dist)
        print 'Using [ph, th] =', strongPeakParams[nnIDX,:2], 'for ', thphPeak
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

    
    XTOF = X[:,:,:,0]
    XTHETA = X[:,:,:,1]
    XPHI = X[:,:,:,2]
    XANGLE = X[:,:,:,1:]
    YBVG2 = bvg(1.0, mu,sigma,XTHETA,XPHI,0)
    YTOF2 = getYTOF(fICC, XTOF, x_lims)
    Y2 = YTOF2*YBVG2
    Y2 = scaleFactor*Y2/Y2.max() 
    
    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    fitMaxIDX = tuple(np.array(np.unravel_index(Y2.argmax(), Y2.shape)))
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


    scaleLinear = Polynomial(n=1)
    scaleX = YJOINT[goodIDX]
    scaleY = n_events[goodIDX]
    scaleWS = CreateWorkspace(OutputWorkspace='scaleWS', dataX=scaleX, dataY=scaleY)   
    fitResultsScaling = Fit(Function=scaleLinear, InputWorkspace='scaleWS', Output='scalefit')
    scaleFitWorkspace = mtd['scaleFit_Workspace']
    A0 = fitResultsScaling[3].row(0)['Value']
    A1 = fitResultsScaling[3].row(1)['Value']
    YRET = A1*YJOINT + A0
    chiSqRed = fitResultsScaling[1]

    print chiSqRed, 'is chiSqRed'
    return YRET, chiSqRed, A1
   
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
        #print i
        for j in xrange(QX.shape[1]):
            for k in xrange(QX.shape[2]):
                newQ = V3D(QX[i,j,k],QY[i,j,k],QZ[i,j,k])
                peak.setQSampleFrame(newQ)
                flightPath = peak.getL1() + peak.getL2()
                scatteringHalfAngle=0.5*peak.getScattering()
                tList[i,j,k] = 3176.507 * flightPath * np.sin(scatteringHalfAngle) / np.linalg.norm(newQ) #convert to microseconds)
    peak.setQSampleFrame(origQS)
    return tList


def fitTOFCoordinate(box,peak, padeCoefficients,dtSpread=0.03,minFracPixels=0.01,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,qMask=None, plotResults=False,fracStop=0.01, pp_lambda=None, pplmin_frac=0.8, pplmax_frac=1.5, mindtBinWidth=1):
    
    #Get info from the peak
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV

    #Set the qMask
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool) 
    
    #Calculate the optimal pp_lambda and 
    tofWS,ppl = ICCFT.getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, qMask, dtSpread=dtSpread, minFracPixels=minFracPixels, neigh_length_m=neigh_length_m,zBG=zBG,pp_lambda=pp_lambda, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

    try:
        fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, fitOrder=bgPolyOrder,constraintScheme=1)
    except:
        fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, fitOrder=bgPolyOrder)
    
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

def getTOFParameters(box, peak, padeCoefficients,dtSpread=0.03,minFracPixels=0.01,neigh_length_m=3,zBG=1.96,bgPolyOrder=1,qMask=None, pplmin_frac=0.8, pplmax_frac=1.5):
    tof = peak.getTOF() #in us
    wavelength = peak.getWavelength() #in Angstrom
    flightPath = peak.getL1() + peak.getL2() #in m
    scatteringHalfAngle = 0.5*peak.getScattering()
    energy = 81.804 / wavelength**2 / 1000.0 #in eV
    if qMask is None:
        qMask = np.ones_like(box.getNumEventsArray()).astype(np.bool)
    tofWS,ppl = ICCFT.getTOFWS(box,flightPath, scatteringHalfAngle, tof, peak, qMask, dtSpread=dtSpread, minFracPixels=minFracPixels, neigh_length_m=neigh_length_m,zBG=zBG, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, padeCoefficients=ICCFT.getModeratorCoefficients('franz_coefficients_2017.dat'))

    fitResults,fICC = ICCFT.doICCFit(tofWS, energy, flightPath, padeCoefficients, fitOrder=bgPolyOrder)
    for i, param in enumerate(['A','B','R','T0','scale', 'hatWidth', 'k_conv']):
        fICC[param] = mtd['fit_Parameters'].row(i)['Value']
    return fICC 

#Does full 3D fit
def fitPeak3D(box, n_events, peak,goodIDX,padeCoefficients,qMask=None):

    QX, QY, QZ = ICCFT.getQXQYQZ(box)
    XSPH = np.array(ICCFT.cart2sph(QX,QY,QZ))
    X = XSPH.swapaxes(0,1).swapaxes(1,2).swapaxes(2,3)
    X[:,:,:,0] = 3176.507*(peak.getL1()+peak.getL2())*np.sin(0.5*peak.getScattering())/X[:,:,:,0]

    predCoefficients=np.array([28.73949834,  13.04192586,   0.41210929])
    neigh_length_m=3
    pplmin_frac=0.4
    pplmax_frac=1.5
    mindtBinWidth=15
    dtSpread=0.015
    bgPolyOrder=1
    zBG=1.96
    plotResults=True

    goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events, peak=peak, box=box,qMask=qMask, calc_pp_lambda=True, padeCoefficients=padeCoefficients, predCoefficients=predCoefficients,neigh_length_m=neigh_length_m, pp_lambda=None, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

    YTOF, fICC, x_lims = fitTOFCoordinate(box,peak,padeCoefficients,dtSpread=dtSpread,qMask=qMask,zBG=1.96,plotResults=False, pp_lambda=pp_lambda, neigh_length_m=3, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac, mindtBinWidth=mindtBinWidth)

    alpha = fICC['A']; beta = fICC['B']; R = fICC['R']; T0 = fICC['T0']; k_conv = fICC['k_conv']; ATOF = fICC['scale']
    x0 = [alpha, beta, R, T0, ATOF, 0.5, k_conv]
    params,h,t,p = doBVGFit(box,nTheta=70,nPhi=70,goodIDX=goodIDX)
    params, cov = params
  
    p0 = [np.max(n_events), ATOF, params[0], params[1], params[2], params[3], params[4], params[5], alpha, beta, R, T0, k_conv, 0.0 ]
    bounds = ([0,0,0, -np.inf,-np.inf,0.,0.,-np.inf] + [0.5*v for v in x0[:4]] + [100, 0],
                [np.inf for i in range(8)] + [2.*v for v in x0[:4]] + [140, np.inf]  )
    bounds = np.array(bounds)
    swapMe = bounds[0] > bounds[1]
    swapIDX = np.where(swapMe)[0]
    for idx in swapIDX:
        tmp = bounds[1,idx]
        bounds[1,idx] = bounds[0,idx]
        bounds[0,idx] = tmp

    #goodIDX = np.zeros_like(n_events,dtype=np.bool)
    #cIDX = np.array(np.shape(goodIDX))//2
    #dX = 8
    #goodIDX[cIDX[0]-dX:cIDX[0]+dX,cIDX[1]-dX:cIDX[1]+dX,cIDX[2]-dX:cIDX[2]+dX] = True
    try:
        params= curve_fit(peak3DFitFunction, X[goodIDX], n_events[goodIDX],  p0,maxfev=1000, sigma=np.sqrt(n_events[goodIDX]),bounds=bounds)
    except ValueError as e:
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print e
        print 'Returning initial guess (p0, p0)'
        params = (p0,p0)
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
    tofxx = np.linspace(tofMin, tofMax, 500)
    tofyy = fICC.function1D(tofxx.ravel())
    ftof = interp1d(tofxx, tofyy)

    YTOF = ftof(XTOF)
    YTOF -= YTOF.min()
    YTOF /= YTOF.max()


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
            #goodIDX,pp_lambda = ICCFT.getBGRemovedIndices(n_events, peak=peak, box=box,qMask=qMask, calc_pp_lambda=True, padeCoefficients=padeCoefficients )
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

def integrateBVGPeak(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, dQPixel=0.005,fracHKL = 0.5,  fracHKLRefine = 0.2,nTheta=400,nPhi=400, q_frame='sample'):
    box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, dQ, fracHKL = fracHKL, dQPixel=dQPixel,q_frame=q_frame)
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

    pos = np.empty(TH.shape + (2,))
    pos[:,:,0] = TH
    pos[:,:,1] = PH

    H = np.empty(h.shape + (2,))
    H[:,:,0] = h
    H[:,:,1] = h


    if forceParams is None:
        meanTH = TH.mean()
        meanPH = PH.mean()
        sigX0 = 0.0018#ICCFT.oldScatFun(meanPH, 1.71151521e-02,   6.37218400e+00,   3.39439675e-03)
        sigY0 = 0.0018#np.polyval([ 0.00045678, -0.0017618 ,  0.0045013 , -0.00480677,  0.00376619], meanTH)
        sigP0 = 0.0#fSigP(meanTH,  0.1460775 ,  1.85816592,  0.26850086, -0.00725352) 

        #Set some constraints
        boundsDict = {}
        boundsDict['A'] = [0.0, np.inf]
        boundsDict['muX'] = [thBins[thBins.size//2 - 4], thBins[thBins.size//2 + 4]]
        boundsDict['muY'] = [phBins[phBins.size//2 - 4], phBins[phBins.size//2 + 4]]
        boundsDict['sigX'] = [0., 0.02]
        boundsDict['sigY'] = [0., 0.02]
        boundsDict['sigP'] = [-1.0, 1.0]

        # Set our initial guess
        m = mbvg.MBVG() 
        m.init()
        m['A'] = np.max(h)
        m['muX'] = meanTH
        m['muY'] = meanPH
        m['sigX'] = sigX0
        m['sigY'] = sigY0
        m['sigP'] = sigP0
        m.setAttributeValue('nX',h.shape[0])
        m.setAttributeValue('nY',h.shape[1])
        m.setConstraints(boundsDict)
        print(m)

        # Do the fit
        bvgWS = CreateWorkspace(OutputWorkspace='bvgWS',DataX=pos.ravel(),DataY=H.ravel(),dataE=np.sqrt(H.ravel()))
        fitFun = FunctionWrapper(m) + Polynomial(n=0)
        fitResults = Fit(Function=fitFun, InputWorkspace='bvgWS', Output='bvgfit')

    elif forceParams is not None:
        p0 = np.zeros(7)
        p0[0] = np.max(h); p0[1] = TH.mean(); p0[2] = PH.mean()
        p0[3] = forceParams[5]; p0[4] = forceParams[6]; p0[5] = forceParams[7];

        #Set some constraints
        isPos = np.sign(p0)
        bounds = ((1.0-isPos*forceTolerance)*p0, (1.0+isPos*forceTolerance)*p0)
        bounds[0][0] = 0.0;  bounds[1][0] = np.inf; #Amplitude
        bounds[0][1] = min(thBins[thBins.size//2 - dth], thBins[thBins.size//2 + dth]);
        bounds[1][1] = max(thBins[thBins.size//2 - dth], thBins[thBins.size//2 + dth]);
        bounds[0][2] = min(phBins[phBins.size//2 - dph], phBins[phBins.size//2 + dph]);
        bounds[1][2] = max(phBins[phBins.size//2 - dph], phBins[phBins.size//2 + dph]);
        bounds[1][-1] = np.inf

        boundsDict = {}
        boundsDict['A'] = [0.0, np.inf]
        boundsDict['muX'] = [bounds[0][1], bounds[1][1]]
        boundsDict['muY'] = [bounds[0][2], bounds[1][2]]
        boundsDict['sigX'] = [bounds[0][3], bounds[1][3]]
        boundsDict['sigY'] = [bounds[0][4], bounds[1][4]]
        boundsDict['sigP'] = [bounds[0][5], bounds[1][5]]

        # Set our initial guess
        m = mbvg.MBVG()
        m.init()
        m['A'] = np.max(h)
        m['muX'] = TH.mean()
        m['muY'] = PH.mean()
        m['sigX'] = forceParams[5]
        m['sigY'] = forceParams[6]
        m['sigP'] = forceParams[7]
        m.setAttributeValue('nX',h.shape[0])
        m.setAttributeValue('nY',h.shape[1])
        m.setConstraints(boundsDict)

        # Do the fit
        bvgWS = CreateWorkspace(OutputWorkspace='bvgWS',DataX=pos.ravel(),DataY=H.ravel(),dataE=np.sqrt(H.ravel()))
        fitFun = FunctionWrapper(m) + Polynomial(n=0)
        fitResults = Fit(Function=fitFun, InputWorkspace='bvgWS', Output='bvgfit')

    # Recover the result
    m = mbvg.MBVG()
    m.init()
    m['A'] = mtd['bvgfit_Parameters'].row(0)['Value']
    m['muX'] = mtd['bvgfit_Parameters'].row(1)['Value']
    m['muY'] = mtd['bvgfit_Parameters'].row(2)['Value']
    m['sigX'] = mtd['bvgfit_Parameters'].row(3)['Value']
    m['sigY'] = mtd['bvgfit_Parameters'].row(4)['Value']
    m['sigP'] = mtd['bvgfit_Parameters'].row(5)['Value']
    m.setAttributeValue('nX',h.shape[0])
    m.setAttributeValue('nY',h.shape[1])
    chiSq = fitResults[1]
    params = [[m['A'], m['muX'], m['muY'], m['sigX'], m['sigY'], m['sigP'], mtd['bvgfit_Parameters'].row(6)['Value']], chiSq]

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
