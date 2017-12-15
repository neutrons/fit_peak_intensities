import numpy as np
import ICCFitTools as ICCFT
reload(ICCFT)
import BVGFitTools as BVGFT
reload(BVGFT)
import sys
sys.path.append("/opt/mantidnightly/bin")
from mantid.simpleapi import *
import os 
import pickle
import ICCAnalysisTools as ICAT
import getEdgePixels as EdgeTools

# Some parameters
workDir = '/SNS/users/ntv/dropbox/' #End with '/'
dQPixel = 0.005 #dQ for each voxel in qBox - recommended to decrease for successive fits
doVolumeNormalization = False #True if you want to normalize TOF profiles by volume
refineCenter = False #True if you want to determine new centers - still not very good
removeEdges = False #True if you want to not consider q-pixels that are off detector faces
calcTOFPerPixel = False #True to calculate TOF for each pixel in a qBox - uses interpolation for volNorm (if doVolumeNormalization is True)
fracHKL = 0.5 #Fraction of HKL to look on either side
fracStop = 0.01 #Fraction of max counts to include in peak selection
neigh_length_m = 3 #Will average over a (neigh_length_m)**3 box
zBG = 1.96 #z score to keep this with

'''
#Scolecite - 2016A
loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
sampleRuns = range(15629,  15644)
peaksFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/%i_Niggli.integrate'
peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
UBFormat = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/%i_Niggli.mat'
UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
descriptor = 'scol_3d' #Does not end with '/'
numTimesToInterpolate=1
dtBinWidth = 4 #Width (in us) in TOF profile bins
qLow = -25; qHigh=25;
'''
#Beta Lac
loadDir = '/SNS/MANDI/IPTS-15000/data/'
nxsTemplate = loadDir+'MANDI_%i_event.nxs'
sampleRuns = range(4999,5003+1)
peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
peaksFormat = peaksFile
UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
UBFormat = UBFile
DetCalFile = None
qLow = -5.0; qHigh = 5.0
dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
dtBinWidth = 40 #Width (in us) in TOF profile bins
dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
descriptor = 'beta_lac_3D_full' #Does not end with '/'
doIterativeBackgroundFitting = False
nBG=5
parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
numTimesToInterpolate=0
workDir = '/SNS/users/ntv/dropbox/'
descriptorRead = 'beta_lac_predpws5'



peaks_ws = LoadIsawPeaks(Filename = peaksFile)
LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()

dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=0.5))
dQ[dQ>0.3] = 0.3
qMask = ICCFT.getHKLMask(UBMatrix, frac=fracHKL, dQPixel=dQPixel,dQ=dQ)

padeCoefficients = ICCFT.getModeratorCoefficients('franz_coefficients_2017.dat')
ICCFitParams = ICAT.getFitParameters(workDir, descriptorRead, sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
ICCFitDict = ICAT.getFitDicts(workDir, descriptorRead,sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
strongPeakParams = pickle.load(open('strongPeakParams.pkl', 'rb'))

from timeit import default_timer as timer

badFits = []
oldNewList = []
for sampleRun in sampleRuns:
    fileName = nxsTemplate%sampleRun
    MDdata = ICCFT.getSample(sampleRun, DetCalFile, workDir, fileName,qLow=qLow, qHigh=qHigh)
    t1 = timer()
    numerrors=0
    numgood=0
    paramList = []
    for peakNumber in range(peaks_ws.getNumberPeaks()):
        TPS = timer()
        peak = peaks_ws.getPeak(peakNumber)
        print peakNumber, peak.getIntensity()
        try:
            if peak.getRunNumber() == sampleRun:
                print 'Integrating peak %i'%peakNumber
                box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel)
                Y3D, goodIDX, pp_lambda, params = BVGFT.get3DPeak(peak, box, padeCoefficients,qMask,nTheta=70, nPhi=70, plotResults=False,nBG=5, dtBinWidth=dtBinWidth,zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1,numTimesToInterpolate=numTimesToInterpolate, fICCParams=ICCFitParams[peakNumber], oldICCFit=ICCFitDict[peakNumber], strongPeakParams=strongPeakParams)

                intensity = np.sum(Y3D[Y3D/Y3D.max() >0.05])/2**(3*numTimesToInterpolate)
                skipIDX = 2**numTimesToInterpolate
                bg = np.sum(goodIDX[Y3D[::skipIDX,::skipIDX,::skipIDX]/Y3D[::skipIDX,::skipIDX,::skipIDX].max()>0.05]*pp_lambda)
                sigma = np.sqrt(intensity + bg)
                oldNewVals = [peaks_ws.getPeak(peakNumber).getIntensity(), peaks_ws.getPeak(peakNumber).getSigmaIntensity(), intensity, sigma]
                print 'original: %4.2f +- %4.2f;  new: %4.2f +- %4.2f'%(oldNewVals[0], oldNewVals[1], oldNewVals[2], oldNewVals[3])
                oldNewList.append(oldNewVals)
                params['peakNumber'] = peakNumber
                params['Intens3d'] = intensity
                params['SigInt3d'] = sigma
                paramList.append(params)
                peak.setIntensity(intensity)
                peak.setSigmaIntensity(sigma)
                numgood += 1
                print 'Finished peak in %f s'%(timer()-TPS)

        except KeyboardInterrupt:
            raise
        except:
            #raise
            numerrors += 1
            print 'error with peak %i'%peakNumber 
            #paramList.append({'Alpha':0.0, 'Beta':0.0, 'R':0.0, 'T0':0.0, 'chiSq3d':1.0e10, 'dQ':0.0, 'k_conv':0.0,
            #                    'muPH':0.0, 'muTH':0.0, 'newQSample':np.array([0.,0.,0.]), 'peakNumber':peakNumber,
            #                    'scale':0.0, 'sigP':0.0, 'sigX':0.0, 'sigY':0.0})
            peak.setIntensity(0.0)
            peak.setSigmaIntensity(1.0)
            badFits.append(peakNumber)
            oldNewList.append([peaks_ws.getPeak(peakNumber).getIntensity(), peaks_ws.getPeak(peakNumber).getSigmaIntensity(), 0.0, 1.0])
    t2 = timer()
    print 'Fit %i peaks in %4.4f s'%(numgood, t2-t1)
    os.system('rm %s'%(workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor)))
    SaveIsawPeaks(InputWorkspace='peaks_ws', Filename=workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
    pickle.dump(paramList, open(workDir+descriptor+'/bvgParams_%i_%s.pkl'%(sampleRun, descriptor),'wb'))