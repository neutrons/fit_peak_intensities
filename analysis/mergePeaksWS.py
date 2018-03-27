import numpy as np
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import sys
import ICCAnalysisTools as ICAT

#pth
sampleRuns = [870,872,873,874,875,876]
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'pth_tof_secondRun'
loadDir = '/SNS/MANDI/IPTS-10943/0/NeXus/'
#eventWS = Load('/SNS/MANDI/IPTS-10943/0/870/NeXus/MANDI_870_event.nxs')
eventWS = Load('/data/pth/MANDI_870_event.nxs')



'''
#gfp
sampleRuns = range(599,607+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'gfp_tof_goodhkl'
loadDir = '/SNS/MANDI/IPTS-8776/data/'
eventWS = Load('/SNS/MANDI/2013_2_11B_SCI/0/599/NeXus/MANDI_599_event.nxs')
'''

'''
#betalac mut
sampleRuns = range(5921,5931+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'beta_lac_lab_highres_mut2'
loadDir = '/SNS/MANDI/IPTS-8776/data/'
eventWS = Load(loadDir+'MANDI_%i_event.nxs'%sampleRuns[0])
'''

'''
#betalac
sampleRuns = range(4999,5003+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'beta_lac_lab_highres2'
loadDir = '/SNS/MANDI/IPTS-15000/data/'
eventWS = Load(loadDir+'MANDI_%i_event.nxs'%sampleRuns[0])
'''

'''
#PsbO
sampleRuns = range(6154,6165+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'psbo_lab_newpredppl_highres'
loadDir = '/SNS/MANDI/IPTS-16286/data/'
eventWS = Load(loadDir+'MANDI_%i_event.nxs'%sampleRuns[0])
'''

'''
# DNA
sampleRuns = range(8758,8769+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'dna_3D_full_lab_newpredppl_newsigi'
loadDir = '/SNS/MANDI/IPTS-18552/nexus/'
eventWS = Load(loadDir+'MANDI_%i.nxs.h5'%sampleRuns[0])
'''

'''
#CORELLI - beryl
sampleRuns = range(58411,58592+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptor = 'beryl_lab_cs1'
loadDir = '/data/corelli_beryl/IPTS-20302/'
eventWS = Load(loadDir+'CORELLI_%i.nxs.h5'%sampleRuns[0])
'''


peaksFile = workDir+descriptor+'/peaks_%i_%s.integrate'

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
