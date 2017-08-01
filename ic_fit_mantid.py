import sys
import os
sys.path.append(os.environ['MANTIDPATH'])
from mantid.simpleapi import *
from mantid.api._api import IFunction, IPeakFunction
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.misc import factorial
from scipy.optimize import curve_fit
from scipy import integrate
import sys
import math
import ICConvoluted as ICC
reload(ICC)
plt.ion()

doPlotBGModel = False 

'''#Use this for testing evaluating newIC
#Define ICC and guess some initial parameters
r = ICC.IkedaCarpenterConvoluted()
r.init()
params = np.array([3.1415e+05, 6.7163e+03, 1.6690e-01, 4.5778e+00, 4.2383e+01,
   1.4795e-03, 3.5690e-02, 7.7825e+01, 5.0000e-01, 5.7125e-02, 5.0000e-01])
x = np.linspace(0.19646,0.20890,20)

parms = np.array([98611.0351,9469.4850,0.2111,2.7892,38.0725,0.0015,0.0000,145.3966,42.0000,0.0566,-0.0112])
parms = np.array([88262.3366,9353.4956,0.1980,2.6622,38.7615,0.0015,0.0000,143.6135,-35.0000,0.0571,-0.0113])
parms
x = np.linspace(0.1964, 0.21036,20)

for paramIDX, param in enumerate(params):
	r.setParameter(paramIDX, param)

r.setCentre(np.mean(x))
r.setFwhm(np.mean(np.diff(x))*3) #Lazy estimate - can be calculated 
f = r.functionLocal(x)
r.setHeight(np.max(f))
plt.figure(1); plt.clf(); plt.plot(x,f,'.-')
plt.xlabel('1/|q|')
plt.ylabel('Events')
sys.exit()
#End testing
'''
def poisson(k,lam):
	return (lam**k)*(np.exp(-lam)/factorial(k))
def normPoiss(k,lam):
	pp_dist = poisson(eventX,lam)
	pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
	return pp_dist*pp_n

def plotBGModel(eventX, eventHist, lam,figNumber=10):
	plt.figure(figNumber)
	plt.clf()
	plt.semilogy(eventX,eventHist,'o',label='Data')
	vec = normPoiss(eventX,lam)
	vec[vec<1] = 1
	plt.semilogy(eventX,vec, label='Poisson model')
	plt.xlim([0,len(eventHist)])
	tmpy = plt.ylim()
	plt.ylim([1.0,tmpy[1]])
	plt.legend(loc='best')
	plt.xlabel('Counts')
	plt.ylabel('Number of q pixels')


box = Load('/SNS/users/ntv/dropbox/run15647/MDpeak_15647_4.nxs')
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
pp_lambda,cov = curve_fit(normPoiss, eventX, eventHist,p0=0.1)
#pp_lambda = 0.083725 #Testing - this is rick's value
if doPlotBGModel:
	plotBGModel(eventX, eventHist, pp_lambda)

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
realNeutronIDX = boxMeanIDX[signalIDX].astype(int)

#This is to change values - I'm just making a note of it here
#box.setSignalAt(index=0,value=1.0)

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
#tList = np.squeeze(sio.loadmat('tau.mat')['tau']) #This is for testing only - values are slightly different (<1e-5)
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
padWidth = 4
tPointsPadded = np.lib.pad(tPoints,(padWidth,), 'linear_ramp',
	end_values=(min(tPoints)-padWidth*dt,max(tPoints)+padWidth*dt))
countsPadded = np.lib.pad(h[0],(padWidth,),'constant')

plt.figure(1); plt.clf();
plt.plot(tPoints,h[0],'o',label='Data')
plt.plot(tPointsPadded,countsPadded,'ko',label='Padded Data')
plt.xlabel('TOF')
plt.ylabel('Counts')

# Let's try to fit the profile as a mantid workspace
tofWS = CreateWorkspace(DataX=tPointsPadded, DataY=countsPadded)

#det = LoadInstrument(tofWS, '/opt/Mantid/instrument/TOPAZ_Definition_2016-04-06.xml',RewriteSpectraMap=True)
#LoadInstrumentFromNexus(tofWS, Filename='/SNS/TOPAZ/shared/PeakIntegration/data/TOPAZ_15647_event.nxs')

#This section will fit to an IkedaCarpenterPV - which does not fit the observed tails as well as our ICConvoluted does. 
xRange = [np.min(tPointsPadded)-dt/2.0,np.max(tPointsPadded)+dt/2.0]
peakRange = [np.min(tPoints)-dt/2.0, np.max(tPoints)+dt/2.0]
x0 = [1.03,7.4e-4,7.4e-4,2.0e-2,1.0e-2,2.0e-9,2.3e-4,0.153]
x0 = [8.66e+00,5.67e-03,7.68e-03,2.00e-02,1.00e-02,2.92e-05,3.28e-02,2.69e-01]
x0[-1] = tPointsPadded[np.argmax(countsPadded)] #Peak center

fitWS = FitPeak(InputWorkspace='tofWS', ParameterTableWorkspace='peakresult', PeakFunctionType='IkedaCarpenterPV (I, Alpha0, Alpha1, Beta0, Kappa, SigmaSquared, Gamma, X0)', FitWindow=xRange, PeakRange=peakRange,BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterValues=x0,FitBackgroundFirst=False)

#Get the function
fIC = FunctionFactory.createFunction('IkedaCarpenterPV')
fitParms = fitWS.FittedPeakParameterValues
for i, param in enumerate(fitParms):
	fIC.setParameter(i,param)
plt.plot(fitWS.OutputWorkspace.readX(1), fitWS.OutputWorkspace.readY(1),'b',label='IkedaCarpenterPV')


#Fit the profile with an exponential/sq wave modified IC
tofWS = CreateWorkspace(DataX=tPoints, DataY=h[0])
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)
x0 = [3.2415e+05, 6.7163e+03, 1.6690e-01, 4.5778e+00, 4.2383e+01,
   1.4795e-03, 3.5690e-02, 7.7825e+01, 0.5, 5.7125e-02, 5.0000e-01,0.0]
x0 = [8804.7119,1.2173,0.0000,7.1541,45.8515,0.0026,0.0000,9.2165,0.5000,0.0093,0.00240,0.0]
#Here we'll try estimating some values based on the total counts

#Estiamte from total counts
totalCounts = np.sum(h[0])
x0[0] = np.polyval([9.2653e-05,8.9562e-01,0.0000e+00],totalCounts)
x0[1] = np.polyval([9.9709e-06,1.6083e-01,0.0000e+00],totalCounts)
x0[2] = 0.1
x0[3] = 6.2117
x0[4] = 43.988*np.exp(-0.0432*totalCounts) + 2.3954*np.exp(-2.519e-5*totalCounts) + 41.937 #c(5)
x0[5] = 2.8514*np.exp(-0.0689*totalCounts) + 0.4744*np.exp(-0.0041*totalCounts) + 0.0113 #c(6) 
x0[8] = 0.50
#x0[9] = 1316.8*np.exp(-0.0119*totalCounts) + 100.0*np.exp(-0.09986*totalCounts) + 21.93 
x0[-1] = tPoints[np.argmax(h[0])] #Peak center
'''
# Estiamte how rick does
x0[0] = np.max(h[0])
x0[1] = x0[0]/10.0
x0[2] = 0.1
x0[3] = np.min(np.where(h[0]>0)).astype(float)
x0[4] = np.max(np.where(h[0]>0)).astype(float)
x0[5] = np.max(np.where(h[0]/np.max(h[0]) > 0.10)).astype(float)
x0[6] = 0.0
x0[7] = 50.0
x0[8] = 4.0
x0[-1] = tPoints[np.argmax(h[0])] #Peak center
'''
#x0 = [9611.0368,0.0000,0.0296,6.8605,44.2194,0.0024,0.0000,7.8570,0.5000,0.0093,0.5000, 0.230, 0.05, 1400]
#x0 = [98611.0351,9469.4850,0.2111,2.7892,38.0725,0.0015,0.0000,145.3966,42.0000,0.0566,-0.0112,0.1975,0.001,1400]
#x0 = [88262.3366,9353.4956,0.1980,2.6622,38.7615,0.0015,0.0000,143.6135,-35.0000,0.0571,-0.0113,0.1985,0.005,8000]
fICC = ICC.IkedaCarpenterConvoluted()
fICC.init()
paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
xRange = [np.min(tPointsPadded)-dt/2.0,np.max(tPointsPadded)+dt/2.0]
peakRange = [np.min(tPoints)-dt/2.0, np.max(tPoints)+dt/2.0]

fitWS = FitPeak(InputWorkspace='tofWS', ParameterTableWorkspace='peakresultConv', PeakFunctionType='IkedaCarpenterConvoluted', FitWindow=xRange, PeakRange=peakRange,BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterNames=paramNames, PeakParameterValues=x0,FitBackgroundFirst=False)

#Get the function
plt.plot(fitWS.OutputWorkspace.readX(1), fitWS.OutputWorkspace.readY(1),'g.-',label='From fitWS.OutputWorkspace')

#Do it just from parameters
r = ICC.IkedaCarpenterConvoluted()
r.init()
for paramIDX, param in enumerate(x0):
	r.setParameter(paramIDX, param)

f = r.functionLocal(tPoints)
#plt.plot(tPoints,f,'r.-',label='Calculated from X0')
plt.legend(loc='best')

#Finally, print the parameters
q = mtd['peakResultConv']
for i in range(q.rowCount()):
    print q.cell(i,0), str(q.cell(i,1)), str(q.cell(i,2))
