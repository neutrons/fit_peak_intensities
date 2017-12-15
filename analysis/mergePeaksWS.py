import numpy as np
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import ICCAnalysisTools as ICAT
sampleRuns = range(4999,5003+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'beta_lac_predpws5_lorentz'


peaksFile = workDir+descriptor+'/peaks_%i_%s.integrate'

loadDir = '/SNS/MANDI/IPTS-15000/data/'
eventWS = Load(loadDir+'MANDI_%i_event.nxs'%sampleRuns[0])
ws = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws",InstrumentWorkspace=eventWS)
#First one - remove extra peaks
peaks_ws = LoadIsawPeaks(FileName=peaksFile%(sampleRuns[0],descriptor),OutputWorkspace='pws%i'%sampleRuns[0])
for i in range(peaks_ws.getNumberPeaks()):
    if peaks_ws.getPeak(i).getRunNumber() == sampleRuns[0]:
        ws.addPeak(peaks_ws.getPeak(i))


#Now we do the rest
for runNumber in sampleRuns[1:]:
    print 'Loading file ', peaksFile%(runNumber,descriptor)
    pws2 = LoadIsawPeaks(FileName=peaksFile%(runNumber,descriptor),OutputWorkspace='pws%i'%runNumber)
    print ' '
    print 'Merging this run with others...'
    print ' '
    for i in range(pws2.getNumberPeaks()):
        ICAT.print_progress(i, pws2.getNumberPeaks(), prefix = 'Progress:', suffix = 'Complete')
        if pws2.getPeak(i).getRunNumber() == runNumber:
            ws.addPeak(pws2.getPeak(i))
            #pws2.removePeak(i)
    #peaks_ws = CombinePeaksWorkspaces(LHSWorkspace=peaks_ws, RHSWorkspace=pws2, OutputWorkspace=peaks_ws, CombineMatchingPeaks=False)

SaveIsawPeaks(InputWorkspace=ws, Filename=workDir+descriptor+'/peaks_combined_good.integrate')
#SaveIsawPeaks(InputWorkspace=peaks_ws2, Filename='/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.integrate')
#SaveLauenorm(InputWorkspace=peaks_ws,Filename='laue_out/laue_out')
