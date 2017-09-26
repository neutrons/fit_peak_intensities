import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.misc import factorial
from scipy.optimize import curve_fit
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import pickle
from scipy import interpolate
import ICConvoluted as ICC
reload(ICC)
import getEdgePixels as EdgeTools
reload(EdgeTools)

#getDQ determines the q spacing required to make an MDBox 
# extending from [-hkl+0.5,hkl-0.5].  Inputs are a peak 
# object and a vector containing lattice spacings. Returns
# a 3x2 numpy array containing dQ for hkl on the low 
# and high side.
def getDQ(peak, latticeSpacing, crystalSystem):
    dQ = np.zeros((3,2))
    if crystalSystem == 'cubic':
            q0 = peak.getQSampleFrame()
            hkl = peak.getHKL()
            cubicConstant = 2*np.pi/latticeSpacing[0]
            for i in range(3):
                    dhkl = np.zeros(np.shape(hkl))
                    dhkl[i] = 0.5
                    dQ[i,0] = cubicConstant*(np.linalg.norm(hkl - dhkl) - np.linalg.norm(hkl))
                    dQ[i,1] = cubicConstant*(np.linalg.norm(hkl + dhkl) - np.linalg.norm(hkl))

    if crystalSystem == 'monoclinic':
        a = latticeSpacing[0]; b = latticeSpacing[1]; c = latticeSpacing[2];
        alpha = latticeSpacing[3]; beta = latticeSpacing[4]; 
        hkl = peak.getHKL()
        def monoclinic(hkl):
            tmp = (((hkl[0])/(a*np.sin(beta/180.*np.pi)))**2 + (hkl[1]/b)**2 + (hkl[2]/c*np.sin(beta/180*np.pi))**2 + 
                            (2*hkl[0]*hkl[2]*np.cos(beta/180.*np.pi)/(a*c*np.sin(beta/180.*np.pi)**2))  )
            return 2*np.pi/np.sqrt(tmp)

        for i in range(3):
            dhkl = np.zeros(np.shape(hkl))
            dhkl[i] = 0.5
            dQ[i,0] = monoclinic(hkl-dhkl) - monoclinic(hkl)
            dQ[i,1] = monoclinic(hkl+dhkl) - monoclinic(hkl)
    return dQ

# UB = UBmatrix as loaded by LoadIsawUB().  Only works in the
#    Qsample frame right now
def getDQFracHKL(peak, UB, frac=0.5):
    dQ = np.zeros((3,2))
    hkl = peak.getHKL()
    if np.all(hkl == np.array([0,0,0])):
        return 0.5*np.ones((3,2))
    q0 = peak.getQSampleFrame()
    dhkl = np.array([frac, frac, frac])
    qPlus = UB.dot(hkl+dhkl)*2*np.pi
    qMinus = UB.dot(hkl-dhkl)*2*np.pi
    dQ[:,0] = qMinus - q0
    dQ[:,1] = qPlus - q0    
    '''
    for hklIDX in range(3):
        dhkl = np.zeros(3)
        dhkl[hklIDX]=0.5
        #dhkl = np.array([0.5, 0.5, 0.5])
        qPlus = UB.dot(hkl+dhkl)*2*np.pi
        qMinus = UB.dot(hkl-dhkl)*2*np.pi
        print qMinus[hklIDX] - q0[hklIDX]
        dQ[int(hklIDX),0] = qMinus[hklIDX] - q0[hklIDX]
        dQ[int(hklIDX),1] = qPlus[hklIDX] - q0[hklIDX]
    '''
    return dQ

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

def integratePeak(x, y, bg,t0, fracStop = 0.01):
    yScaled = (y-bg) / np.max(y-bg)
    goodIDX = yScaled > fracStop
    if np.sum(goodIDX) > 0:
        iStart = np.min(np.where(goodIDX))
        iStop = np.max(np.where(goodIDX))
        xStart = x[iStart]
        xStop = x[iStop]
    else:
        print 'THIS IS BAD - NO GOOD START/STOP POINT!!'
        return 0.0, 1.0, x[0], x[-1]
    '''
    try:
        iStart = np.min(np.where((y-bg)>ctsStop))
    except:
        raise
        iStart = 0
    xStart = x[iStart]
    try:
        print np.min(y[x>t0]-bg[x>t0])
        iStop = np.min(np.where(np.logical_and(x>t0,(y-bg)<ctsStop)))
    except:
        raise
        iStop = len(x)-1 #TODO: Probably need to go further out
    xStop = x[iStop]
    '''
    intensity = np.sum(y[iStart:iStop] - bg[iStart:iStop])
    sigma = np.sqrt(np.var(y[iStart:iStop])+np.var(bg[iStart:iStop]))
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
def getTOFWS(box, flightPath, scatteringHalfAngle, tofPeak, peak, panelDict, peakNumber, dtBinWidth=2, zBG=-1.0, dtSpread = 0.02, doVolumeNormalization=False, minFracPixels = 0.005, removeEdges=False, edgesToCheck=None):
    #Pull the number of events
    n_events = box.getNumEventsArray()
    hasEventsIDX = np.where((n_events>0) & ~np.isnan(n_events))
    


    #Set up some things to remove bad pixels
    N = np.shape(n_events)[0]
    neigh_length_m = 0 #Set to zero for "this pixel only" mode - can be made faster if using this mode
    maxBin = np.shape(n_events)
    boxMean = np.zeros(len(hasEventsIDX[0]))
    boxMeanIDX = list()

    if zBG >= 0:
        pp_lambda = get_pp_lambda(n_events,hasEventsIDX) #Get the most probably number of events
        #Determine which pixels we want by considering the surrounding box
        for i,idx in enumerate(np.array(hasEventsIDX).transpose()):
            dataBox = n_events[max(idx[0] - neigh_length_m,0):min(idx[0] + neigh_length_m+1, maxBin[0]),
                                           max(idx[1] - neigh_length_m,0):min(idx[1] + neigh_length_m+1, maxBin[1]),
                                           max(idx[2] - neigh_length_m,0):min(idx[2] + neigh_length_m+1, maxBin[2])]
            boxMean[i] = np.mean(dataBox)
            boxMean[i] = n_events[idx[0],idx[1],idx[2]]
            boxMeanIDX.append(idx)
        signalIDX = np.where(boxMean > pp_lambda+zBG*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
    else: #don't do background removal 
        #Determine which pixels we want by considering the surrounding box
        for i,idx in enumerate(np.array(hasEventsIDX).transpose()):
            boxMean[i] = n_events[idx[0],idx[1],idx[2]]
            boxMeanIDX.append(idx)
        signalIDX = np.where(boxMean > 0)
        
    boxMeanIDX = np.asarray(boxMeanIDX)
    realNeutronIDX = boxMeanIDX[signalIDX].astype(int)
    
    #Setup our axes -- ask if there is a way to just get this
    xaxis = box.getXDimension()
    qx = np.linspace(xaxis.getMinimum(), xaxis.getMaximum(), xaxis.getNBins())
    yaxis = box.getYDimension()
    qy = np.linspace(yaxis.getMinimum(), yaxis.getMaximum(), yaxis.getNBins())
    zaxis = box.getZDimension()
    qz = np.linspace(zaxis.getMinimum(), zaxis.getMaximum(), zaxis.getNBins())
    #qSq = np.vstack([np.power(qx,2), np.power(qy,2), np.power(qz,2)]).transpose()

    #Create our TOF distribution from bg corrected data
    useIDX = realNeutronIDX.transpose()
    tList = 1/np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2)) #1/|q|
    tList = 3176.507 * flightPath * np.sin(scatteringHalfAngle) * tList #convert to microseconds
    tMin = np.min(tList)
    tMax = np.max(tList)

    #Set up our bins for histogramming
    dt = tofPeak*dtSpread #time in us on either side of the peak position to consider
    dt = max(dt, 100)
    tMin = max(tMin, tofPeak - dt)
    tMax = min(tMax, tofPeak + dt)
    #tBins = np.linspace(tMin, tMax, nBins+1)
    tBins = np.arange(tMin, tMax, dtBinWidth)
    weightList = n_events[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]
    if removeEdges:
        mask = EdgeTools.getMask(peak, box, panelDict,edgesToCheck=edgesToCheck)
        h = np.histogram(tList,tBins,weights=weightList*mask[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]);
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
        QX, QY, QZ = np.meshgrid(qx, qy, qz,indexing='ij')
        tofBox = 3176.507 * flightPath *np.sin(scatteringHalfAngle) * 1/np.sqrt(QX**2 + QY**2 + QZ**2)
        if removeEdges:
            print 'REMOVING EDGES'
            numPixels = np.histogram((tofBox*mask).flatten(), tBins)[0]
        else:
            numPixels = np.histogram(tofBox.flatten(), tBins)[0]
        yPoints = yPoints / numPixels*1.0
        useIDX = 1.0*numPixels/np.sum(numPixels) > minFracPixels
        if np.sum(useIDX < 1): #Bad threshold, we'll juse use it all
            useIDX = np.ones(np.size(yPoints)).astype(bool)
        tPoints = tPoints[useIDX]
        yPoints = yPoints[useIDX] * np.mean(numPixels)
        #if dtSpread > 0.04:
        #    plt.figure(8); plt.clf()
        #    plt.plot(tBins[1:], 1.0*numPixels/np.sum(numPixels))
        #    plt.figure(2)
        
    tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=yPoints, DataE=np.sqrt(yPoints))
    return tofWS

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
def getInitialGuessByDetector(tofWS, paramNames, energy, flightPath, detNumber):
    x0 = np.zeros(len(paramNames))
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    d = pickle.load(open('det_calibration/calibration_dictionary.pkl','rb'))
    x0[0] = np.polyval(d['det_%i'%detNumber]['A'],energy)
    x0[1] = np.polyval(d['det_%i'%detNumber]['B'],energy)
    x0[2] = np.polyval(d['det_%i'%detNumber]['R'],energy)
    x0[3] = oneOverXSquared(energy, d['det_%i'%detNumber]['T0'][0], d['det_%i'%detNumber]['T0'][1])
    x0[4] = (np.max(y))/x0[0]*2*2.5 
    x0[5] = 0.5
    x0[6] = 200
    return x0

def getInitialGuessSpline(tofWS, paramNames, energy, flightPath):
    x0 = np.zeros(len(paramNames))
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    splineDict = pickle.load(open('get_franz_coefficients/splineDict.pkl','rb'))
    for i, param in enumerate(['A','B','R','T0']):
        x0[i] = interpolate.splev(energy, splineDict[param])
    x0[3] += getT0Shift(energy, flightPath) #Simulates at L~=0, we are ~18m downstream, this adjusts for that time.
    x0[4] = (np.max(y))/x0[0]*2*2.5  #Amplitude
    x0[5] = 0.5 #hat width in IDX units
    x0[6] = interpolate.splev(energy, splineDict['k_conv']) #Exponential decay rate for convolution

    #Phenomenology
    x0[0] /= 1.2
    x0[2] += 0.05
    x0[3] -= 10 #This is lazy - we can do it detector-by-detector

    return x0 

#Returns intial parameters for fitting based on a few quickly derived TOF
# profile parameters.  tofWS is a worskapce containng the TOF profile,
# paramNames is the list of parameter names
# energy is the energy of the peak (units: eV)
# flightPath is L = L1 + L2 (units: m)
def getInitialGuess(tofWS, paramNames, energy, flightPath):
    x0 = np.zeros(len(paramNames))
    x = tofWS.readX(0)
    y = tofWS.readY(0)
    franzCoeff = getModeratorCoefficients('franz_coefficients_2017.dat')
    x0[0] = pade(franzCoeff['A'], energy)
    x0[1] = pade(franzCoeff['B'], energy)
    x0[2] = pade(franzCoeff['R'], energy)
    x0[3] = pade(franzCoeff['T0'], energy)

    #These are still phenomenological
    x0[0] /= 1.2
    x0[2] += 0.05
    #TODO: This can be calculated once and absorbed into the coefficients.
    x0[3] += getT0Shift(energy, flightPath) #Franz simulates at moderator exit, we are ~18m downstream, this adjusts for that time.
    x0[3] -= 10 #This is lazy - we can do it detector-by-detector
    x0[4] = (np.max(y))/x0[0]*2*2.5  #Amplitude
    x0[5] = 0.5 #hat width in IDX units
    #x0[5] = 3.0 #hat width in us
    x0[6] = 120.0 #Exponential decay rate for convolution
    return x0

#Get sample loads the NeXus evnts file and converts from detector space to
# reciprocal space.
# run is the run number.
# UBFile is a string for the file containng the UB matrix for the crystal
# DetCalFile is a string for the file containng the detector calibration
# workDir is not used
# loadDir is the directory to extract the data from
def getSample(run,  UBFile,  DetCalFile,  workDir, fileName):
    #data
    print 'Loading file', fileName
    data = Load(Filename = fileName)
    LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic',
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1',
      MinValues = '-25, -25, -25', Maxvalues = '25, 25, 25')
    return MDdata

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
    dQ = np.abs(getDQFracHKL(peak, UBMatrix,frac=fracHKLSearch))
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
def getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQPixel=0.005,fracHKL = 0.5, refineCenter=False, fracHKLRefine = 0.2):
    runNumber = peak.getRunNumber()
    QSample = peak.getQSampleFrame()
    Qx = QSample[0]
    Qy = QSample[1]
    Qz = QSample[2]
    dQ = np.abs(getDQFracHKL(peak, UBMatrix, frac = fracHKL))
    dQ[dQ > 0.5] = 0.5
    nPtsQ = np.round(np.sum(dQ/dQPixel,axis=1)).astype(int)
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

#Does the actual integration and modifies the peaks_ws to have correct intensities.
def integrateSample(run, MDdata, peaks_ws, paramList, detBankList, UBMatrix, figsFormat=None, dtBinWidth = 4, nBG=15, dtSpread=0.02, fracHKL = 0.5, refineCenter=False, doVolumeNormalization=False, minFracPixels=0.0, fracStop = 0.01, removeEdges=False, panelDict=None):
    if removeEdges is True and panelDict is None:
        print 'REMOVE EDGES WITHOUT panelDict - IMPOSSIBLE!!'
        0/0
    p = range(peaks_ws.getNumberPeaks())
    for i in p:
        peak = peaks_ws.getPeak(i)
        #TODO: does not work if hkl = (0,0,0) - getDQ returns Inf
        if peak.getRunNumber() == run:
            try:#for ppppp in [3]:#try:
                Box = getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, i, fracHKL = fracHKL, refineCenter = refineCenter)
                tof = peak.getTOF() #in us
                wavelength = peak.getWavelength() #in Angstrom
                energy = 81.804 / wavelength**2 / 1000.0 #in eV
                flightPath = peak.getL1() + peak.getL2() #in m
                scatteringHalfAngle = 0.5*peak.getScattering()
                detNumber = detBankList[i]
                print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
                if Box.getNEvents() < 1 or np.all(np.abs(peak.getHKL())==0):
                    print "Peak %i has 0 events or is HKL=000. Skipping!"%i
                    peak.setIntensity(0)
                    peak.setSigmaIntensity(1)
                    paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())])
                    mtd.remove('MDbox_'+str(run)+'_'+str(i))
                    continue
                #Do background removal (optionally) and construct the TOF workspace for fitting
                if removeEdges: 
                    edgesToCheck = EdgeTools.needsEdgeRemoval(Box,panelDict,peak) 
                    if edgesToCheck != []: #At least one plane intersects so we have to fit
                        tofWS = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=removeEdges, edgesToCheck=edgesToCheck)
                    else:
                        tofWS = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False)
                else:
                    tofWS = getTOFWS(Box,flightPath, scatteringHalfAngle, tof, peak, panelDict, i, dtBinWidth=dtBinWidth,dtSpread=dtSpread, doVolumeNormalization=doVolumeNormalization, minFracPixels=minFracPixels, removeEdges=False)

                #Set up our inital guess
                fICC = ICC.IkedaCarpenterConvoluted()
                fICC.init()
                paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
                #x0 = getInitialGuessByDetector(tofWS,paramNames,energy,flightPath, detNumber)
                #x0 = getInitialGuessSpline(tofWS,paramNames,energy,flightPath)
                x0 = getInitialGuess(tofWS,paramNames,energy,flightPath)
                [fICC.setParameter(iii,v) for iii,v in enumerate(x0[:fICC.numParams()])]
                x = tofWS.readX(0)
                y = tofWS.readY(0)
                if len(y)//2 < nBG: nBG = len(y)//2
                bgx0 = np.polyfit(x[np.r_[0:nBG,-nBG:0]], y[np.r_[0:nBG,-nBG:0]], 1)
                
                scaleFactor = np.max(y-np.polyval(bgx0,x))/np.max(fICC.function1D(x))
                x0[4] = x0[4]*scaleFactor
                fICC.setParameter(4,x0[4])
                
                #Form the strings for fitting and do the fit
                paramString = ''.join(['%s=%4.8f, '%(fICC.getParamName(iii),x0[iii]) for iii in range(fICC.numParams())])
                funcString1 = 'name=IkedaCarpenterConvoluted, ' + paramString
                constraintString1 = ''.join(['%s > 0, '%(fICC.getParamName(iii)) for iii in range(fICC.numParams())])
                constraintString1 += 'R < 1'
                #constraintString1 += '100< k_conv < 500'
                
                bgString= '; name=LinearBackground,A0=%4.8f,A1=%4.8f'%(bgx0[1],bgx0[0]) #A0=const, A1=slope
                constraintString2 = ', constraints=(-1.0 < A1 < 1.0)'
                
                functionString = funcString1 + 'constraints=('+constraintString1+')' + bgString + constraintString2  
                #fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit(Function=functionString, InputWorkspace='tofWS', Output='fit') #This is antiquated as of sept 25 2017
                fitResults = Fit(Function=functionString, InputWorkspace='tofWS', Output='fit')
                fitStatus = fitResults.OutputStatus
                chiSq = fitResults.OutputChi2overDoF
                
    
                chiSq2  = 1.0e99
                if chiSq > 2.0: #The initial fit isn't great - let's see if we can do better
                    print '############REFITTING########### on %4.4f'%chiSq
                    #x0 = getInitialGuess(tofWS,paramNames,energy,flightPath)
                    x0 = getInitialGuessByDetector(tofWS,paramNames,energy,flightPath, detNumber)
                    paramWS = mtd['fit_parameters']
                    paramString = ''.join(['%s=%4.8f, '%(fICC.getParamName(iii),x0[iii]) for iii in range(fICC.numParams())])
                    funcString1 = 'name=IkedaCarpenterConvoluted, ' + paramString[:-2]
                    functionString = funcString1 + ', constraints=('+constraintString1+')' + bgString + constraintString2 
                    try:
                            #fitStatus, chiSq2, covarianceTable, paramTable, fitWorkspace = Fit(Function=functionString, InputWorkspace='tofWS', Output='fit2') #Antiquated, Sept 25 2017
                        fitResults2 = Fit(Function=functionString, InputWorkspace='tofWS', Output='fit2')
                        fitStatus2 = fitResults2.OutputStatus
                        chiSq2 = fitResults2.OutputChi2overDoF
                    except:
                            print 'CANNOT DO SECOND FIT, GOING BACK TO FIRST!'
                
                if(chiSq < chiSq2):
                    r = mtd['fit_Workspace']
                    param = mtd['fit_Parameters']
                else:
                    r = mtd['fit2_Workspace']
                    param = mtd['fit2_Parameters']
                    chiSq = chiSq2

                fitBG = [param.cell(iii+2,1),param.cell(iii+1,1)]
                #Set the intensity before moving on to the next peak
                icProfile = r.readY(1)
                bgCoefficients = fitBG
                #peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))i
                t0 = param.row(3)['Value']
                intensity, sigma, xStart, xStop = integratePeak(r.readX(0), icProfile, np.polyval(bgCoefficients, r.readX(1)),t0, fracStop=fracStop)
                icProfile = icProfile - np.polyval(bgCoefficients, r.readX(1)) #subtract background
                peak.setIntensity(intensity)
                peak.setSigmaIntensity(sigma)
                if figsFormat is not None:
                    plotFit(figsFormat, r,tofWS,fICC,peak.getRunNumber(), i, energy, chiSq,fitBG, xStart, xStop, bgx0)
                paramList.append([i, energy, np.sum(icProfile), 0.0,chiSq] + [mtd['fit_Parameters'].row(i)['Value'] for i in range(mtd['fit_parameters'].rowCount())])
                mtd.remove('MDbox_'+str(run)+'_'+str(i))
            except KeyboardInterrupt:
                print 'KeyboardInterrupt: Exiting Program!!!!!!!'
                sys.exit()
            except: #Error with fitting
                #raise
                peak.setIntensity(0)
                peak.setSigmaIntensity(1)
                print 'Error with peak ' + str(i)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())])
           
        mtd.remove('MDbox_'+str(run)+'_'+str(i))
    return peaks_ws, paramList



