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

def pade(c,x): #c are coefficients, x is the energy in eV
     return c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9])

def poisson(k,lam):
        return (lam**k)*(np.exp(-lam)/factorial(k))

def normPoiss(k,lam,eventHist):
        pp_dist = poisson(k,lam)
        pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
        return pp_dist*pp_n

def getTOFWS(box, flightPath, scatteringHalfAngle, tofPeak, nBins=20):
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
	neigh_length_m = 1
	maxBin = np.shape(n_events)
	boxMean = np.zeros(len(hasEventsIDX[0]))
	boxMeanIDX = list()
	#A more pythonic version
	for i,idx in enumerate(np.array(hasEventsIDX).transpose()):
			dataBox = n_events[max(idx[0] - neigh_length_m,0):min(idx[0] + neigh_length_m+1, maxBin[0]),
							   max(idx[1] - neigh_length_m,0):min(idx[1] + neigh_length_m+1, maxBin[1]),
							   max(idx[2] - neigh_length_m,0):min(idx[2] + neigh_length_m+1, maxBin[2])]
			boxMean[i] = np.mean(dataBox)
			boxMeanIDX.append(idx)

	boxMeanIDX = np.asarray(boxMeanIDX)
	signalIDX = np.where(boxMean > pp_lambda+1.65*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3))
	backgrIDX = np.where(boxMean <= pp_lambda)
	print pp_lambda, pp_lambda+1.65*np.sqrt(pp_lambda/(2*neigh_length_m+1)**3), np.shape(signalIDX) 
	realNeutronIDX = boxMeanIDX[signalIDX].astype(int)
	backgroundIDX = boxMeanIDX[backgrIDX].astype(int)
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
	bgIDX = backgroundIDX.transpose()
	tList = 1/np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2)) #1/|q|
	tList = 3176.507 * flightPath * np.sin(scatteringHalfAngle) * tList
	dt = 35 #time in us on either side of the peak position to consider
	tList = tList
	bgList = 1/np.sqrt(np.power(qx[bgIDX[0]],2)+np.power(qy[bgIDX[1]],2)+np.power(qz[bgIDX[2]],2)) 
	tMin = np.min(tList)
	tMax = np.max(tList)
	tmpv = (tMax-tMin)/(nBins)

	#This section selects only a small region (~100us) around the center TOF
	tMin = tofPeak - dt 
	tMax = tofPeak + dt
	

	#tMin = tMin - tmpv
	#tMax = tMax + tmpv
	tBins = np.linspace(tMin, tMax, nBins+1)
	weightList = n_events[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]
	weightListBG = n_events[bgIDX.transpose()[:,0],bgIDX.transpose()[:,1],bgIDX.transpose()[:,2]]
	#For and plot the TOF distribution
	h = plt.hist(tList,tBins,weights=weightList);
	#hBG = np.sum(normPoiss(range(max(weightList)),pp_lambda,np.histogram(weightList,np.array(range(max(weightList)))+0.5)))
	hBG = plt.hist(bgList,tBins,weights=weightListBG);
	tPoints = 0.5*(h[1][1:] + h[1][:-1])
	dt = np.abs(tPoints[1]-tPoints[0]) #assumes uniform time binning
	tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPoints, DataY=h[0])
	return tofWS, hBG

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
	x0[3] += getT0Shift(energy, flightPath) #Franz simulates at moderator exit, we are ~18m downstream, this adjusts for that time.
	x0[7] = np.mean(y[-2:] + y[:2])/2.0 #Background constant
	print 'x07',x0[7]
	x0[4] = (np.max(y)-x0[7])/x0[0]*2*2.5  #Amplitude
	x0[5] = 0.5 #hat width in IDX units
	#x0[5] = 3.0 #hat width in us
	x0[6] = 30.0 #Exponential decay rate for convolution
	return x0	

	'''
	maxCounts = np.max(y)
	maxIDX = np.argmax(y)
	#Potential issue - if max value is first or last value
	if ((y[maxIDX-1]/y[maxIDX] > 0.50) and (y[maxIDX+1]/y[maxIDX] > 0.50) and (maxCounts > 8)): #This is a sharp peak
		if maxCounts > 100:
			if maxIDX < 7:x0 = [4.0599e+04,6.2156e+03,2.1305e-01,2.8399e+00,1.7260e+01,1.1451e-03,0.0000e+00,1.2244e+02,5.0000e-01]
			else:         x0 = [2.4472e+04,0.0000e+00,0.0000e+00,7.9248e+00,1.8618e+01,1.3615e-03,0.0000e+00,2.0311e+04,5.0000e-01]
		else: x0 = [5.2327e+02,8.4743e-01,9.0309e-01,1.2439e+00,6.7153e+01,9.2673e-01,0.0000e+00,5.2968e+02,5.0000e-01]
			
	else: #Not a very sharp peak
		if maxCounts > 1000:
			if maxIDX < 7:
				x0 = [4.0599e+04,6.2156e+03,2.1305e-01,2.8399e+00,1.7260e+01,1.1451e-03,0.0000e+00,1.2244e+02,5.0000e-01]
			else:
				x0 = [2.4472e+04,0.0000e+00,0.0000e+00,7.9248e+00,1.8618e+01,1.3615e-03,0.0000e+00,2.0311e+04,5.0000e-01]
		elif maxCounts > 600:
			x0 =[8.3851e+03,1.1386e+02,0.0000e+00,1.0199e+01,4.4895e+01,9.6696e-03,0.0000e+00,6.2012e+00,5.0000e-01]
		elif maxCounts > 200:
			x0 = [1.8172e+04,3.4136e+02,0.0000e+00,9.0390e+00,4.5982e+01,2.2817e-03,0.0000e+00,3.0412e+00,5.0000e-01]
		else:
			x0 =[2.6754e+02,3.6010e+00,1.4239e-15,9.5552e+00,3.6145e+01,2.4258e-01,0.0000e+00,4.0584e+00,5.0000e-01]	
	return x0
	'''
def getSample(run,  UBFile,  DetCalFile,  workDir,  loadDir):
    #data
    print 'Loading file', loadDir+'TOPAZ_'+str(run)+'_event.nxs'
    data = Load(Filename = loadDir+'TOPAZ_'+str(run)+'_event.nxs')
    LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic', 
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1', 
      MinValues = '-25, -25, -25', Maxvalues = '25, 25, 25')
    return MDdata

def integrateSample(run, MDdata, sizeBox, gridBox, peaks_ws, paramList):

    #getting box for each peak
    p = range(peaks_ws.getNumberPeaks())
    for i in p:
	peak = peaks_ws.getPeak(i)
        if peak.getRunNumber() == run:
            QSample = peak.getQSampleFrame()
            Qx = QSample[0]
            Qy = QSample[1]
            Qz = QSample[2]
            Box = BinMD(InputWorkspace = 'MDdata', 
                AlignedDim0='Q_sample_x,'+str(Qx-sizeBox)+','+str(Qx+sizeBox)+','+str(gridBox),
                AlignedDim1='Q_sample_y,'+str(Qy-sizeBox)+','+str(Qy+sizeBox)+','+str(gridBox),
                AlignedDim2='Q_sample_z,'+str(Qz-sizeBox)+','+str(Qz+sizeBox)+','+str(gridBox),
                OutputWorkspace = 'MDbox_'+str(run)+'_'+str(i))
                
            tof = peak.getTOF() #in us
            wavelength = peak.getWavelength() #in Angstrom
            energy = 81.804 / wavelength**2 / 1000.0 #in eV
            flightPath = peak.getL1() + peak.getL2() #in m
            scatteringHalfAngle = 0.5*peak.getScattering()
            print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
            print energy*1000.0,'meV'
            try:#for abcd in ['x']:#try:
                #Do background removal and construct the TOF workspace for fitting
                tofWS, hBG = getTOFWS(Box,flightPath, scatteringHalfAngle, tof,nBins=20)

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

                if redChiSq > 10.0: #I'm not sure why it's stopping yet - so let's just try again using our old solution as a new one
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
		plt.plot(xSmooth,fICC.function1D(xSmooth+1.40*np.mean(np.diff(x))),label='Fit')
                #plt.plot(r.readX(0)[2:-2],fICC.function1D(r.readX(0))[2:-2])
                #plt.plot(r.readX(0),hBG[0])
                plt.title('E0=%4.4f meV, redChiSq=%4.4e'%(energy*1000,redChiSq))
		plt.legend(loc='best')
                plt.savefig('tof_integration_figs/mantid_'+str(peak.getRunNumber())+'_'+str(i)+'.png')

                #Set the intensity before moving on to the next peak
		icProfile = r.readY(1)
		icProfile = icProfile - mtd['fit_Parameters'].row(7)['Value'] #subtract background
                peak.setIntensity(np.sum(icProfile))
                peak.setSigmaIntensity(np.sqrt(np.sum(icProfile)))
		paramList.append([i, energy, np.sum(icProfile), redChiSq,chiSq] + [mtd['fit_Parameters'].row(i)['Value'] for i in range(mtd['fit_parameters'].rowCount())])
            except:
                peak.setIntensity(0)
                peak.setSigmaIntensity(1)
		paramList.append([i, energy, np.sum(icProfile), redChiSq,chiSq] + [0 for i in range(mtd['fit_parameters'].rowCount())])
                print 'Error with peak ' + str(i)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
           
    	    mtd.remove('MDbox_'+str(run)+'_'+str(i))
    return peaks_ws, paramList


sizeBox = 0.5
gridBox = 201


'''#Scolecite - 2016A
sampleRuns = range(15629,  15630)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
'''
#Si - 2016A
sampleRuns = range(15647,15670)
#sampleRuns = range(15657,15670)
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
#peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.integrate'
#peaksFile = '/SNS/users/vel/Dropbox (ORNL)/first62.peaks'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.mat'

DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
peaks_ws = LoadIsawPeaks(Filename = peaksFile)
for sampleRun in sampleRuns: 
    paramList = list()
    MDdata = getSample(sampleRun, UBFile, DetCalFile, workDir, loadDir)
    peaks_ws,paramList= integrateSample(sampleRun, MDdata, sizeBox, gridBox, peaks_ws,paramList)
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename='peaks_%i_removeBG_includeBadPeaks.integrate'%(sampleRun))
    np.savetxt('params_%i_removeBG.dat'%sampleRun, np.array(paramList))
    wsList = mtd.getObjectNames()
    for ws in wsList:
        if 'MDbox_' in ws:
            mtd.remove(ws)
