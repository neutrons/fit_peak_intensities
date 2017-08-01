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

def pade(c,x):
     return np.log10(c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9]))

def pade2011(c,x):
    	return c[0]*x**c[1]*(1+c[2]*x+c[3]*x**2+(x/c[4])**c[5])/(1+c[6]*x+c[7]*x**2+(x/c[8])**c[9])
	#return c[0]*(x**c[1])*(1+c[2]*x+c[3]*x**2+(c[4]/c[5])**c[6])/(1+c[7]*x+c[8]*x**2+(c[9]/c[10])**c[11])

def poisson(k,lam):
        return (lam**k)*(np.exp(-lam)/factorial(k))

def normPoiss(k,lam,eventHist):
        pp_dist = poisson(k,lam)
        pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
        return pp_dist*pp_n

def getTOFWS(box):
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
	tList = 1/np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2))
	tMin = np.min(tList)
	tMax = np.max(tList)
	tmpv = (tMax-tMin)/(20)
	tMin = tMin - tmpv
	tMax = tMax + tmpv
	tBins = np.linspace(tMin, tMax, 21)
	weightList = n_events[useIDX.transpose()[:,0],useIDX.transpose()[:,1],useIDX.transpose()[:,2]]

	#For and plot the TOF distribution
	h = plt.hist(tList,tBins,weights=weightList);
	tPoints = 0.5*(h[1][1:] + h[1][:-1])
	dt = np.abs(tPoints[1]-tPoints[0]) #assumes uniform time binning
	padWidth = 3
	tPointsPadded = np.lib.pad(tPoints,(padWidth,), 'linear_ramp',
        end_values=(min(tPoints)-padWidth*dt,max(tPoints)+padWidth*dt))
	countsPadded = np.lib.pad(h[0],(padWidth,),'constant')
	tofWS = CreateWorkspace(OutputWorkspace='tofWS', DataX=tPointsPadded, DataY=countsPadded)
	return tofWS


def getSample(run,  UBFile,  DetCalFile,  workDir,  loadDir):
    #data
    print 'Loading file', loadDir+'TOPAZ_'+str(run)+'_event.nxs'
    data = Load(Filename = loadDir+'TOPAZ_'+str(run)+'_event.nxs')
    LoadIsawDetCal(InputWorkspace = data, Filename = DetCalFile)
    MDdata = ConvertToMD(InputWorkspace = data, QDimensions = 'Q3D', dEAnalysisMode = 'Elastic', 
      Q3DFrames = 'Q_sample', QConversionScales = 'Q in A^-1', 
      MinValues = '-25, -25, -25', Maxvalues = '25, 25, 25')
    return MDdata

def integrateSample(run, MDdata, sizeBox, gridBox, peaksFile):
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)

    #getting box for each peak
    p = range(100)#range(peaks_ws.getNumberPeaks())
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
	    
	tof = peak.getTOF() #in units?
	wavelength = peak.getWavelength() #in Angstrom
	energy = 81.804 / wavelength**2 #in meV
    	print '---fitting peak ' + str(i) + '  Num events: ' + str(Box.getNEvents()), ' ', peak.getHKL()
	print energy,'meV'
	for abcd in [1]:#try:
		#Do background removal and construct the TOF workspace for fitting
		tofWS = getTOFWS(Box)
		tofWS = mtd['tofWS']		

		#Integrate the peak
		totalCounts = np.sum(tofWS.readY(0)) 
		tPointsPadded = tofWS.readX(0)
		tPoints = tofWS.readX(0)[3:-3]
		dt = np.abs(tPoints[1]-tPoints[0]) #assumes uniform time binning
		
		fICC = ICC.IkedaCarpenterConvoluted()
		fICC.init()
		paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
		xRange = [np.min(tPointsPadded)-dt/2.0,np.max(tPointsPadded)+dt/2.0]
		peakRange = [np.min(tPoints)-dt/2.0, np.max(tPoints)+dt/2.0]
		x0 = np.zeros(len(paramNames))	
		x0[0] = np.polyval([9.2653e-05,8.9562e-01,0.0000e+00],totalCounts)
		x0[1] = np.polyval([9.9709e-06,1.6083e-01,0.0000e+00],totalCounts)
		x0[2] = 0.1
		x0[4] = 43.988*np.exp(-0.0432*totalCounts) + 2.3954*np.exp(-2.519e-5*totalCounts) + 41.937 #c(5)
		x0[5] = 2.8514*np.exp(-0.0689*totalCounts) + 0.4744*np.exp(-0.0041*totalCounts) + 0.0113 #c(6) 
		x0[6] = 0.0
		x0[7] = 50.0
		x0[8] = 0.50
		
		'''#x0[0] = pade2011(alphax0, energy/1000.0)*4.0e3 #alpha (1/s)
		#x0[1] = pade2011(betax0, energy/1000.0)*1.0e3 #beta (1/s)
		#x0[2] = pade2011(rx0, energy/1000.0) #R
		x0[3] = 1.0 #tfacta = First index nonzero - we have removed zeros already
		x0[4] = len(peakRange) #tfactb - last index
		x0[5] = 4.0#pade2011(t0x0, energy/1000.0) #"t0" -where we see our increase
		x0[6] = 0.0 #Background built into our model
		x0[7] = 50.0 #gc_alpha 
		x0[8] = 0.50 #pp_width
		'''
		x0[9] = tofWS.readX(0)[np.argmax(tofWS.readY(0))] #Peak center

		fitWS = FitPeak(InputWorkspace='tofWS', OutputWorkspace='fitWS_'+str(i), ParameterTableWorkspace='parmWS_'+str(i), PeakFunctionType='IkedaCarpenterConvoluted', FitWindow=xRange, PeakRange=peakRange,BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterNames=paramNames, PeakParameterValues=x0,FitBackgroundFirst=False)
		
		plt.figure(1); plt.clf()
		r = mtd['fitWS_'+str(i)] 
		plt.plot(r.readX(0)[3:-3],r.readY(0)[3:-3],'o')
		plt.plot(r.readX(1)[3:-3],r.readY(1)[3:-3],'.-')
		plt.savefig('integration_method_comparison_figs/mantid_'+str(i)+'.png')
		#Set the intensity before moving on to the next peak
'''	except:
		print 'Error with peak ', str(i)
		exc_type, exc_obj, exc_tb = sys.exc_info()
    		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    		print(exc_type, fname, exc_tb.tb_lineno)
'''


sizeBox = 0.5
gridBox = 201

#These are initial Pade parameters from Gallmeier's 2010 MCSTAS simulation, table 3 (viewed by BL11)
alphax0 = [0.07611, 0.1374, -3.188, 20279, 3.844e-5,6.030,-79.76,15903,3.827e-5,5.741]
betax0  = [0.1770, 0.4062, -26.09, 206.4, 0.1212, -20.77, -18.26, 196.1, 0.1103, -20.82]
rx0 =     [0.09278, 0.0, -55.18, 253.3, 0.2709, -1.025, -21.85, 105.5, 0.03883, -1.025]
t0x0 =    [107.2, 0.1411, -4859, 39474, 1.750e-4, 0.9896, 22314, 139, 1.244e-3, 2.594]

'''#Scolecite - 2016A
sampleRuns = range(15629,  15630)
peaksFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFile='/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
'''
#Si - 2016A
sampleRuns = range(15647,15648)
#peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.integrate'
peaksFile = '/SNS/users/vel/Dropbox (ORNL)/first62.peaks'
UBFile =  '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/15647_Niggli.mat'

DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
workDir = '/SNS/users/ntv/dropbox/'
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
for sampleRun in sampleRuns: 
    MDdata = getSample(sampleRun, UBFile, DetCalFile, workDir, loadDir)
    peaks = integrateSample(sampleRun, MDdata, sizeBox, gridBox, peaksFile)

