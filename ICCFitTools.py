import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.misc import factorial
from scipy.optimize import curve_fit
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import ICConvoluted as ICC
reload(ICC)


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

    return dQ

#Standard Pade function used for getting initial guesses
def pade(c,x): #c are coefficients, x is the energy in eV
    return c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9])

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
def getTOFWS(box, flightPath, scatteringHalfAngle, tofPeak, dtBinWidth=2, zBG=-1.0, dtSpread = 0.03):
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
    qSq = np.vstack([np.power(qx,2), np.power(qy,2), np.power(qz,2)]).transpose()

    #Create our TOF distribution from bg corrected data
    useIDX = realNeutronIDX.transpose()
    tList = 1/np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2)) #1/|q|
    tList = 3176.507 * flightPath * np.sin(scatteringHalfAngle) * tList #convert to microseconds
    tMin = np.min(tList)
    tMax = np.max(tList)

    #Set up our bins for histogramming
    dt = tofPeak*dtSpread #time in us on either side of the peak position to consider
    tMin = tofPeak - dt
    tMax = tofPeak + dt
    #tBins = np.linspace(tMin, tMax, nBins+1)
    tBins = np.arange(tMin, tMax, dtBinWidth)
    weightList = n_events[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]

    #For and plot the TOF distribution
    h = np.histogram(tList,tBins,weights=weightList);
    tPoints = 0.5*(h[1][1:] + h[1][:-1])
    tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=h[0], DataE=np.sqrt(h[0]))
    #tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=h[0])
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
def getSample(run,  UBFile,  DetCalFile,  workDir,  loadDir):
    #data
    print 'Loading file', loadDir+'TOPAZ_'+str(run)+'_event.nxs'
    data = Load(Filename = loadDir+'TOPAZ_'+str(run)+'_event.nxs')
    LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic',
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1',
      MinValues = '-25, -25, -25', Maxvalues = '25, 25, 25')
    return MDdata

#Function to make and save plots of the fits. bgx0=polynomial coefficients (polyfit order)
# for the initial guess
def plotFit(filenameFormat, r,tofWS,fICC,runNumber, peakNumber, energy, chiSq,bgx0=None):
    plt.figure(1); plt.clf()
    plt.plot(r.readX(0),r.readY(0),'o',label='Data')
    if bgx0 is not None:
    	plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0))+np.polyval(bgx0, tofWS.readX(0)),'b',label='Initial Guess')
    else:
    	plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0)),'b',label='Initial Guess')
        
    plt.plot(r.readX(1),r.readY(1),'.-',label='Fit')

    plt.title('E0=%4.4f meV, redChiSq=%4.4e'%(energy*1000,chiSq))
    plt.legend(loc='best')
    plt.savefig(filenameFormat%(runNumber, peakNumber))

#getBoxHalfHKL returns the binned MDbox ranging from (hkl-0.5)-(hkl+0.5) (i.e. half integers 
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
def getBoxHalfHKL(peak, MDdata, latticeConstants,crystalSystem,gridBox,peakNumber):
    run = peak.getRunNumber()
    QSample = peak.getQSampleFrame()
    Qx = QSample[0]
    Qy = QSample[1]
    Qz = QSample[2]
    dQ = np.abs(getDQ(peak, latticeConstants, crystalSystem))
    Box = BinMD(InputWorkspace = 'MDdata',
        AlignedDim0='Q_sample_x,'+str(Qx-dQ[0,0])+','+str(Qx+dQ[0,1])+','+str(gridBox),
        AlignedDim1='Q_sample_y,'+str(Qy-dQ[1,0])+','+str(Qy+dQ[1,1])+','+str(gridBox),
        AlignedDim2='Q_sample_z,'+str(Qz-dQ[2,0])+','+str(Qz+dQ[2,1])+','+str(gridBox),
        OutputWorkspace = 'MDbox_'+str(run)+'_'+str(peakNumber))
    return Box

#Does the actual integration and modifies the peaks_ws to have correct intensities.
def integrateSample(run, MDdata, latticeConstants,crystalSystem, gridBox, peaks_ws, paramList, figsFormat=None, nBG=15):

    p = range(peaks_ws.getNumberPeaks())
    for i in p:
        peak = peaks_ws.getPeak(i)
        if peak.getRunNumber() == run:
            Box = getBoxHalfHKL(peak, MDdata, latticeConstants, crystalSystem, gridBox, i)
            tof = peak.getTOF() #in us
            wavelength = peak.getWavelength() #in Angstrom
            energy = 81.804 / wavelength**2 / 1000.0 #in eV
            flightPath = peak.getL1() + peak.getL2() #in m
            scatteringHalfAngle = 0.5*peak.getScattering()
            print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
            print energy*1000.0,'meV'
            if Box.getNEvents() < 1:
                print "Peak %i has 0 events. Skipping!"%i
                peak.setIntensity(0)
                peak.setSigmaIntensity(1)
                paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())])
                continue
            for ppppp in [3]:#try:
                #Do background removal (optionally) and construct the TOF workspace for fitting
                tofWS = getTOFWS(Box,flightPath, scatteringHalfAngle, tof,dtBinWidth=2)

                #Set up our inital guess
                fICC = ICC.IkedaCarpenterConvoluted()
                fICC.init()
                paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
                x0 = getInitialGuess(tofWS,paramNames,energy,flightPath)
                [fICC.setParameter(iii,v) for iii,v in enumerate(x0[:fICC.numParams()])]
                x = tofWS.readX(0)
                y = tofWS.readY(0)
		bgx0 = np.polyfit(x[np.r_[0:nBG,-nBG:0]], y[np.r_[0:nBG,-nBG:0]], 1)
                bgx0[0] = 0.0
                bgx0[1] = np.mean(y[np.r_[0:nBG,-nBG:0]])
                #TODO: make this permanent, but this will tweak our background
                scaleFactor = np.max(y)/np.max(fICC.function1D(x)+bgx0[1])
                x0[4] = x0[4]*scaleFactor
                fICC.setParameter(4,x0[4])
		#Form the strings for fitting and do the fit
                #TODO: test if '4.4f' is accurate enough
                paramString = ''.join(['%s=%4.8f, '%(fICC.getParamName(iii),x0[iii]) for iii in range(fICC.numParams())])
                funcString = 'name=IkedaCarpenterConvoluted, ' + paramString
                funcString = funcString[:-2] #Remove last comma so we can append with BG
		bgString = '; name=LinearBackground,A0=%4.8f,A1=%4.8f'%(bgx0[1],bgx0[0]) #A0=const, A1=slope
                constraintString = ''.join(['f0.%s > 0, '%(fICC.getParamName(iii)) for iii in range(fICC.numParams())])
                constraintString += 'f0.R < 1'
                constraintString += ', f1.A1 < 0.01, f1.A0<%4.8f'%np.max(y)
                fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit(Function=funcString+bgString, InputWorkspace='tofWS', Output='fit',Constraints=constraintString)

                if chiSq > 10.0: #The initial fit isn't great - let's see if we can do better
                        print '############REFITTING###########'
                        paramWS = mtd['fit_parameters']
                        paramString = ''.join(['%s=%4.8f, '%(fICC.getParamName(iii),paramWS.cell(iii,1)) for iii in range(fICC.numParams()-1)])
                        funcString = 'name=IkedaCarpenterConvoluted, ' + paramString[:-2]
		        bgString = '; name=LinearBackground,A0=%4.8f,A1=%4.8f'%(paramWS.cell(iii+2,1),paramWS.cell(iii+1,1))
                        try:
                                fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit(Function=funcString+bgString, InputWorkspace='tofWS', Output='fit',Constraints=constraintString,Minimizer='Trust Region')
                        except:
                                print 'CANNOT DO SECOND FIT, GOING BACK TO FIRST!'


                r = mtd['fit_Workspace']
                param = mtd['fit_Parameters']
                if figsFormat is not None:
                    plotFit(figsFormat, r,tofWS,fICC,peak.getRunNumber(), i, energy, chiSq,bgx0)
                #Set the intensity before moving on to the next peak
                icProfile = r.readY(1)
                bgCoefficients = [param.row(7)['Value'], param.row(8)['Value']]
                icProfile = icProfile - np.polyval(bgCoefficients, r.readX(1)) #subtract background
                peak.setIntensity(np.sum(icProfile))
                peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))
                paramList.append([i, energy, np.sum(icProfile), 0.0,chiSq] + [mtd['fit_Parameters'].row(i)['Value'] for i in range(mtd['fit_parameters'].rowCount())])

            '''except: #Error with fitting
                peak.setIntensity(0)
                peak.setSigmaIntensity(1)
                print 'Error with peak ' + str(i)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                paramList.append([i, energy, 0.0, 1.0e10,1.0e10] + [0 for i in range(mtd['fit_parameters'].rowCount())])
           '''
            mtd.remove('MDbox_'+str(run)+'_'+str(i))
    return peaks_ws, paramList



