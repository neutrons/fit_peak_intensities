import sys
import os
sys.path.append(os.environ['MANTIDPATH'])
from mantid.simpleapi import *
from mantid.api._api import IFunction, IPeakFunction
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import factorial
from scipy.optimize import curve_fit
from scipy import integrate
import sys
import math
import ICConvoluted as ICC
reload(ICC)
plt.ion()

doPlotBGModel = True

'''#Use this for testing evaluating newIC
#Define ICC and guess some initial parameters
r = ICC.IkedaCarpenterConvoluted()
r.init()
params = np.array([3.1415e+05, 6.7163e+03, 1.6690e-01, 4.5778e+00, 4.2383e+01,
   1.4795e-03, 3.5690e-02, 7.7825e+01, 5.0000e-01, 5.7125e-02, 5.0000e-01])
x = np.linspace(0.19646,0.20890,20)

params = np.array([9611.0368,0.0000,0.0296,6.8605,44.2194,0.0024,0.0000,7.8570,0.5000,0.0093,0.5000])
x = np.linspace(0.2304, 0.36684,20)

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
'''#End testing

def poisson(k,lam):
	return (lam**k/factorial(k))*np.exp(-lam)
def normPoiss(k,lam):
	pp_dist = poisson(eventX,lam)
	pp_n = np.dot(pp_dist,eventHist)/np.dot(pp_dist,pp_dist)
	return pp_dist*pp_n

def plotBGModel(eventX, eventHist, lam,figNumber=10):
	plt.figure(figNumber)
	plt.clf()
	xplot = np.linspace(0,np.max(eventX),100)
	plt.semilogy(eventX,eventHist,'o',label='Data')
	plt.semilogy(eventX, normPoiss(eventX,lam),label='Poisson model')
	plt.xlim([1,len(eventHist)])
	tmpy = plt.ylim()
	plt.ylim([1.0,tmpy[1]])
	plt.legend(loc='best')
	plt.xlabel('Counts')
	plt.ylabel('Number of q pixels')


#box = Load('MDpeak_15629_61.nxs')
box = Load('MDpeak_15647_3.nxs')
#Pull the number of events
n_events = box.getNumEventsArray()

#Estimate a poisson background
hasEventsIDX = np.where((n_events>0) & ~np.isnan(n_events))
eventValues = n_events[hasEventsIDX]
numEvents = np.sum(eventValues)
eventHist = np.bincount(eventValues.astype('int'))
eventX = np.array(range(len(eventHist)))
#rmIDX = eventHist == 0
#eventX = eventX[~rmIDX]
#eventHist = eventHist[~rmIDX]
pp_lambda,cov = curve_fit(normPoiss, eventX, eventHist,p0=0.228)
pp_lambda = 0.10907 #Testing - this is rick's value
if doPlotBGModel:
	plotBGModel(eventX, eventHist, pp_lambda)

N = np.shape(n_events)[0]
x_vals = np.zeros([N,N,N])
y_vals = np.zeros([N,N,N])
z_vals = np.zeros([N,N,N])
xyz = range(N)
for j1 in range(N):
    for j2 in range(N):
        x_vals[:,j1,j2] = xyz
        y_vals[j1,:,j2] = xyz
        z_vals[j1,j2,:] = xyz
xyz_vals = np.vstack([x_vals[hasEventsIDX], y_vals[hasEventsIDX], z_vals[hasEventsIDX]])

neigh_length_m = 3
maxBin = np.shape(n_events)
boxMean = np.zeros(len(hasEventsIDX[0]))
boxMeanIDX = list()
for i,idx in enumerate(np.array(hasEventsIDX).transpose()[:10]):
	dataBox = n_events[max(idx[0] - neigh_length_m,0):min(idx[0] + neigh_length_m, maxBin[0]),
			   max(idx[1] - neigh_length_m,0):min(idx[1] + neigh_length_m, maxBin[1]),
			   max(idx[2] - neigh_length_m,0):min(idx[2] + neigh_length_m, maxBin[2])]
	print i+1, max(idx[0] - neigh_length_m,0), min(idx[0] + neigh_length_m, maxBin[0])
	boxMean[i] = np.mean(dataBox)
	boxMeanIDX.append(idx)
boxMeanIDX = np.asarray(boxMeanIDX)
signalIDX = boxMean > 1.65*np.sqrt(pp_lambda/(2*neigh_length_m)**3)
realNeutronIDX = boxMeanIDX[signalIDX]

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
tMin = 1/np.max(np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2)))
tMax = 1/np.min(np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2)))
tBins = np.linspace(tMin, tMax, 20)
tList = 1/np.sqrt(np.power(qx[useIDX[0]],2)+np.power(qy[useIDX[1]],2)+np.power(qz[useIDX[2]],2))

#Plot the TOF distribution
h = plt.hist(tList,tBins);
tPoints = 0.5*(h[1][1:] + h[1][:-1])
plt.figure(1); plt.clf();
plt.plot(tPoints,h[0],'o',label='Data')
plt.xlabel('TOF')
plt.ylabel('Counts')



# Let's try to do this as a mantid workspace
tofWS = CreateWorkspace(DataX=tPoints, DataY=h[0])

#det = LoadInstrument(tofWS, '/opt/Mantid/instrument/TOPAZ_Definition_2016-04-06.xml',RewriteSpectraMap=True)
#LoadInstrumentFromNexus(tofWS, Filename='/SNS/TOPAZ/shared/PeakIntegration/data/TOPAZ_15647_event.nxs')
'''
fitWS = FitPeak(InputWorkspace='tofWS', ParameterTableWorkspace='peakresult', PeakFunctionType='IkedaCarpenterPV (I, Alpha0, Alpha1, Beta0, Kappa, SigmaSquared, Gamma, X0)', FitWindow=[0.22,0.25], PeakRange=[0.222,0.248],BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterValues=[1.03,7.4e-4,7.4e-4,2.0e-2,1.0e-2,2.0e-9,2.3e-4,0.23],FitBackgroundFirst=True)

#Get the function
fIC = FunctionFactory.createFunction('IkedaCarpenterPV')
fitParms = fitWS.FittedPeakParameterValues
for i, param in enumerate(fitParms):
	fIC.setParameter(i,param)
plt.plot(fitWS.OutputWorkspace.readX(1), fitWS.OutputWorkspace.readY(1))
'''

#Now do it with the exponential/sq wave modified IC
tofWS = CreateWorkspace(DataX=tPoints, DataY=h[0])
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)
x0 = [3.1415e+05, 6.7163e+03, 1.6690e-01, 4.5778e+00, 4.2383e+01,
   1.4795e-03, 3.5690e-02, 7.7825e+01, 0.5, 5.7125e-02, 5.0000e-01,
   0.230, 0.05, 1400 ]
#x0 = [9611.0368,0.0000,0.0296,6.8605,44.2194,0.0024,0.0000,7.8570,0.5000,0.0093,0.5000, 0.230, 0.05, 1400]
fICC = ICC.IkedaCarpenterConvoluted()
fICC.init()
paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]

fitWS = FitPeak(InputWorkspace='tofWS', ParameterTableWorkspace='peakresultConv', PeakFunctionType='IkedaCarpenterConvoluted', FitWindow=[0.218,0.25], PeakRange=[0.222,0.248],BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterNames=paramNames, PeakParameterValues=x0,FitBackgroundFirst=False)

#Get the function
plt.plot(fitWS.OutputWorkspace.readX(1), fitWS.OutputWorkspace.readY(1),'-o',label='From fitWS.OutputWorkspace')

#Do it just from parameters
r = ICC.IkedaCarpenterConvoluted()
r.init()
for paramIDX, param in enumerate(x0):
	r.setParameter(paramIDX, param)

f = r.functionLocal(tofWS.readX(0))
plt.plot(tofWS.readX(0),f,'r.-',label='Calculated from X0')
plt.legend(loc='best')
'''
#Now do it with the sample gaussian
tofWS = CreateWorkspace(DataX=tPoints, DataY=h[0])
FunctionFactory.subscribe(ICC.PyGaussian)
x0 = [1400.0, 0.23, 0.005]

fICC = ICC.PyGaussian()
fICC.init()
paramNames = [fICC.getParamName(x) for x in range(fICC.numParams())]
fitWS = FitPeak(InputWorkspace='tofWS', ParameterTableWorkspace='peakresultConv', PeakFunctionType='PyGaussian', FitWindow=[0.22,0.25], PeakRange=[0.222,0.248],BackgroundType='Flat (A0)',BackgroundParameterValues=[0],PeakParameterNames=paramNames, PeakParameterValues=x0,FitBackgroundFirst=True)

#Get the function
plt.plot(fitWS.OutputWorkspace.readX(1), fitWS.OutputWorkspace.readY(1))
'''
