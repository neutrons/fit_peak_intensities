import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.misc import factorial
from scipy.optimize import curve_fit
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import ICConvoluted as ICC
reload(ICC)
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)

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


def pade(c,x): #c are coefficients, x is the energy in eV
     return c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9])

def poisson(k,lam):
        return (lam**k)*(np.exp(-lam)/factorial(k))

def normPoiss(k,lam,eventHist):
        pp_dist = poisson(k,lam)
        pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
        return pp_dist*pp_n

def getTOFWS(box, flightPath, scatteringHalfAngle, tofPeak, dtBinWidth=2):
	#Pull the number of events
	n_events = box.getNumEventsArray()

	#Estimate a poisson background
	hasEventsIDX = np.where((n_events>0) & ~np.isnan(n_events))
	eventValues = n_events[hasEventsIDX]
	numEvents = np.sum(eventValues)
	eventHist = np.bincount(eventValues.astype('int'))[1:]
	eventX = np.array(range(len(eventHist)))+1
	#rmIDX = eventHist == 0
	#eventX = eventX[~rmIDX]
	#eventHist = eventHist[~rmIDX]
	normPoissL = lambda k,lam: normPoiss(k,lam,eventHist)
	pp_lambda,cov = curve_fit(normPoissL, eventX, eventHist,p0=0.1)
	#pp_lambda = 0.083725 #Testing - this is rick's value

	#Set up some things to remove bad pixels
	N = np.shape(n_events)[0]
	neigh_length_m = 0 #Set to zero for "this pixel only" mode - can be made faster if using this mode
	maxBin = np.shape(n_events)
	boxMean = np.zeros(len(hasEventsIDX[0]))
	boxMeanIDX = list()
	#Determine which pixels we want by considering the surrounding box
	for i,idx in enumerate(np.array(hasEventsIDX).transpose()):
			dataBox = n_events[max(idx[0] - neigh_length_m,0):min(idx[0] + neigh_length_m+1, maxBin[0]),
							   max(idx[1] - neigh_length_m,0):min(idx[1] + neigh_length_m+1, maxBin[1]),
							   max(idx[2] - neigh_length_m,0):min(idx[2] + neigh_length_m+1, maxBin[2])]
			boxMean[i] = np.mean(dataBox)
			boxMean[i] = n_events[idx[0],idx[1],idx[2]]
			boxMeanIDX.append(idx)
	
	boxMeanIDX = np.asarray(boxMeanIDX)
	#signalIDX = np.where(boxMean > pp_lambda+1.65*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
	signalIDX = np.where(boxMean > 0)
	#print some info - note that shape signalIDX != num_events if pixels have more than one event
	print pp_lambda, pp_lambda+1.65*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3), np.shape(signalIDX) 
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
	tList = 3176.507 * flightPath * np.sin(scatteringHalfAngle) * tList
	tMin = np.min(tList)
	tMax = np.max(tList)

	#Set up our bins for histogramming
	dt = tofPeak*0.03 #time in us on either side of the peak position to consider
	tMin = tofPeak - dt 
	tMax = tofPeak + dt
	#tBins = np.linspace(tMin, tMax, nBins+1)
	tBins = np.arange(tMin, tMax, dtBinWidth)
	weightList = n_events[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]

	#For and plot the TOF distribution
	h = plt.hist(tList,tBins,weights=weightList);
	tPoints = 0.5*(h[1][1:] + h[1][:-1])
	dt = np.abs(tPoints[1]-tPoints[0]) #assumes uniform time binning
	tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=h[0])
	return tofWS

def getT0Shift(E,L):
	#E is energy in eV, L is flight path in m
	mn = 1.674929e-27 #neutron mass, kg
	E = E * 1.60218e-19 #convert eV to J
	t0Shift = L*np.sqrt(mn/2/E) #units = s
	t0Shift = t0Shift * 1.0e6 #units =us
	return t0Shift

#Returns intial parameters for fitting based on a few quickly derived TOF profile parameters
def getInitialGuess(tofWS, paramNames, energy, flightPath):
	x0 = np.zeros(len(paramNames))
	x = tofWS.readX(0)
	y = tofWS.readY(0)
	franzCoeff = list()
	franzCoeff.append([121.589384271076, 0.810023477482879, 2.6472470527418, 20278.9841066641, 0.000122668075600656, 1.68157845080583, -84.8598732886308, 15903.4011384698, 5.05509989155325e-05, 1.99732564937083])
	franzCoeff.append([0.0852522484075364, 0.463613094284145, -39.890703844218, 438.080722418418, 0.0200866047920943, -88.5200439215893, -28.2058014872526, 249.628780977559, 0.0218095681648157, -88.0548527180181])
	franzCoeff.append([0.000173294401015949,  0.0, 5502.17688451783,  44311.7900549361,  0.0211087412796991,  -92.12364298736,  3.07133940226443,  34.8657659038993,  0.0217699376583704, -92.12364298736])
	franzCoeff.append([30.838010830493,  0.015616308796326,  -4854.14796435711,  39473.6647918935,  0.00014255433182369,  0.955627865605097, 22311.000613891,  138.895395806703,  0.00155257873919149,  2.48432323742575])
	franzCoeff = np.asarray(franzCoeff)
	for i in range(len(franzCoeff)):
		x0[i] = pade(franzCoeff[i], energy)
	#x0[3] = x[np.argmax(y)] #t0 - this just uses the max value as initial guess
	#TODO: This can be calculated once and absorbed into the coefficients.
	x0[0] /= 1.2
	x0[3] += getT0Shift(energy, flightPath) #Franz simulates at moderator exit, we are ~18m downstream, this adjusts for that time.
	x0[7] = np.mean(y[-2:] + y[:2])/2.0 #Background constant
	x0[4] = (np.max(y)-x0[7])/x0[0]*2*2.5*2  #Amplitude
	x0[5] = 0.5 #hat width in IDX units
	#x0[5] = 3.0 #hat width in us
	x0[6] = 30.0 #Exponential decay rate for convolution
	return x0	

def getSample(run,  UBFile,  DetCalFile,  workDir,  loadDir):
    #data
    print 'Loading file', loadDir+'TOPAZ_'+str(run)+'_event.nxs'
    data = Load(Filename = loadDir+'TOPAZ_'+str(run)+'_event.nxs')
    LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic', 
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1', 
      MinValues = '-25, -25, -25', Maxvalues = '25, 25, 25')
    return MDdata

def integrateSample(run, MDdata, latticeConstants,crystalSystem, gridBox, peaks_ws, paramList):

    #getting box for each peak
    p = range(peaks_ws.getNumberPeaks())
    for i in p:
	peak = peaks_ws.getPeak(i)
        if peak.getRunNumber() == run:
            QSample = peak.getQSampleFrame()
            Qx = QSample[0]
            Qy = QSample[1]
            Qz = QSample[2]
            dQ = np.abs(getDQ(peak, latticeConstants, crystalSystem))
            Box = BinMD(InputWorkspace = 'MDdata', 
                AlignedDim0='Q_sample_x,'+str(Qx-dQ[0,0])+','+str(Qx+dQ[0,1])+','+str(gridBox),
                AlignedDim1='Q_sample_y,'+str(Qy-dQ[1,0])+','+str(Qy+dQ[1,1])+','+str(gridBox),
                AlignedDim2='Q_sample_z,'+str(Qz-dQ[2,0])+','+str(Qz+dQ[2,1])+','+str(gridBox),
                OutputWorkspace = 'MDbox_'+str(run)+'_'+str(i))
                
            tof = peak.getTOF() #in us
            wavelength = peak.getWavelength() #in Angstrom
            energy = 81.804 / wavelength**2 / 1000.0 #in eV
            flightPath = peak.getL1() + peak.getL2() #in m
            scatteringHalfAngle = 0.5*peak.getScattering()
            print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
            print energy*1000.0,'meV'
            if Box.getNEvents() < 1:
                print "Peak %i has 0 events. Skipping!%i"
                continue
            for abcd in ['x']:#try:
                #Do background removal and construct the TOF workspace for fitting
                tofWS = getTOFWS(Box,flightPath, scatteringHalfAngle, tof,dtBinWidth=2)

		# Fitting starts here
                #Integrate the peak
                totalCounts = np.sum(tofWS.readY(0)) 
                tPointsPadded = tofWS.readX(0)
                tPoints = tofWS.readX(0)
                dt = np.abs(tPoints[1]-tPoints[0]) #assumes uniform time binning


                fICC = ICC.IkedaCarpenterConvoluted()
                fICC.init()
                paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
                xRange = [np.min(tPointsPadded)-dt/2.0,np.max(tPointsPadded)+dt/2.0]
                peakRange = [np.min(tPoints)-dt/2.0, np.max(tPoints)+dt/2.0]
                x0 = getInitialGuess(tofWS,paramNames,energy,flightPath)

                [fICC.setParameter(iii,v) for iii,v in enumerate(x0[:fICC.numParams()])]		
                paramString = ''.join(['%s=%4.4f, '%(fICC.getParamName(iii),x0[iii]) for iii in range(fICC.numParams())])

                funcString = 'name=IkedaCarpenterConvoluted, ' + paramString
                constraintString = ''.join(['%s > 0, '%(fICC.getParamName(iii)) for iii in range(fICC.numParams())])
                constraintString += 'R < 1'
                fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit(Function=funcString, InputWorkspace='tofWS', Output='fit',Constraints=constraintString)
                r = mtd['fit_Workspace']
		redChiSq = np.sum(np.power(r.readY(0) - r.readY(1),2) / r.readY(1))/(len(x0)-1) #reduced chisq

		
                if chiSq > 10.0: #I'm not sure why it's stopping yet - so let's just try again using our old solution as a new one
                        print '############REFITTING###########'
                        paramWS = mtd['fit_parameters']
                        paramString = ''.join(['%s=%4.4f, '%(fICC.getParamName(iii),paramWS.cell(iii,1)) for iii in range(paramWS.rowCount()-1)])
                        funcString = 'name=IkedaCarpenterConvoluted, ' + paramString[:-2]
                        print funcString
                        try:
                                fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit(Function=funcString, InputWorkspace='tofWS', Output='fit',Constraints=constraintString,Minimizer='Trust Region')
                        except:
                                print 'CANNOT DO SECOND FIT, GOING BACK TO FIRST!'


                r = mtd['fit_Workspace'] 
		
                redChiSq = np.sum(np.power(r.readY(0) - r.readY(1),2) / r.readY(1))/(len(x0)-1) #reduced chisq
                plt.figure(1); plt.clf()
                plt.plot(r.readX(0),r.readY(0),'o',label='Data')
		plt.plot(tofWS.readX(0), fICC.function1D(tofWS.readX(0)),'b',label='Initial Guess')
                plt.plot(r.readX(1),r.readY(1),'.-',label='Fit')

		#set up a nice, smooth IC plot
		p = mtd['fit_Parameters']
		x = r.readX(0)
		xSmooth = np.linspace(np.min(x), np.max(x),400)
		for parami in range(fICC.numParams()):
			fICC.setParameter(parami, p.row(parami)['Value'])
		#plt.plot(xSmooth,fICC.function1D(xSmooth+1.40*np.mean(np.diff(x))),label='Fit')
                #plt.plot(r.readX(0)[2:-2],fICC.function1D(r.readX(0))[2:-2])
                plt.title('E0=%4.4f meV, redChiSq=%4.4e'%(energy*1000,chiSq))
		plt.legend(loc='best')
                plt.savefig('tof_integration_dynamic_binning/mantid_'+str(peak.getRunNumber())+'_'+str(i)+'.png')
                #plt.savefig('tof_integration_figs_withBG_scolecite/mantid_'+str(peak.getRunNumber())+'_'+str(i)+'.png')

                #Set the intensity before moving on to the next peak
		icProfile = r.readY(1)
		icProfile = icProfile - mtd['fit_Parameters'].row(7)['Value'] #subtract background
                peak.setIntensity(np.sum(icProfile))
                peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))
		paramList.append([i, energy, np.sum(icProfile), redChiSq,chiSq] + [mtd['fit_Parameters'].row(i)['Value'] for i in range(mtd['fit_parameters'].rowCount())])
            '''except:
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


sizeBox = 0.5
gridBox = 201
'''
#Scolecite - 2016A
sampleRuns = range(15629,  15644)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
'''
#Si - 2016A
sampleRuns = range(15647,15670)
#sampleRuns = range(15666,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
#peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.integrate'
#peaksFile = '/SNS/users/vel/Dropbox (ORNL)/first62.peaks'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.mat'
crystalSystem ='cubic'
latticeConstants = [5.43071] #Since it's cubic, this we only need a (in angstrom)

DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
peaks_ws = LoadIsawPeaks(Filename = peaksFile)
for sampleRun in sampleRuns: 
    paramList = list()
    MDdata = getSample(sampleRun, UBFile, DetCalFile, workDir, loadDir)
    peaks_ws,paramList= integrateSample(sampleRun, MDdata, latticeConstants,crystalSystem, gridBox, peaks_ws,paramList)
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename='peaks_%i_dynamic_binning.integrate'%(sampleRun))
    np.savetxt('params_%i_dynamic_binning.dat'%sampleRun, np.array(paramList))
    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
