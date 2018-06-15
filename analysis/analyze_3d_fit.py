import numpy as np
import pandas as pd
import sys
for i in range(len(sys.path))[::-1]:
    if 'antid' in sys.path[i]:
        sys.path.pop(i)
sys.path.append("/home/ntv/mantid/mantid/bin")
sys.path.append("./analysis")
from mantid.simpleapi import *
from mantid.geometry import SpaceGroupFactory, PointGroupFactory
import matplotlib.pyplot as plt
plt.ion()
import convertToPandas as pdTOF
import ICCAnalysisTools as ICAT
from scipy.interpolate import interp1d
import getEdgePixels as EdgeTools
import seaborn as sns
import pickle
from scipy.optimize import curve_fit

#------------------------------Load the bvgFit files
#PsbO
sampleRuns = range(6154, 6165+1)
workDir = '/SNS/users/ntv/dropbox/'
#descriptorBVG = 'psbo_3D_full_lab_newpredppl_highres_newsigi'
descriptorBVG = 'psbo_3D_mbvg_predppl'
descriptorTOF = None#'psbo_lab_newpredppl_highres'
#peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorBVG) #TOF file, BVGxTOF is from fitDict
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorBVG, sampleRuns[-1], descriptorBVG)
#ellipseFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal_highres.integrate'
ellipseFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal_highres_2p0A.integrate'
peaksFile = ellipseFile #This is done to get around corrupted ISAW files; we'll only use peak properties anyway June 15 2018 -BS
sg = SpaceGroupFactory.createSpaceGroup("P 61 2 2")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)

'''
# CORELLI - beryl
sampleRuns = range(58411,58592+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'beryl_3D_full_newsigi'
descriptorTOF = 'beryl_lab_cs1'
#peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal_indexedonly.integrate'
pg = PointGroupFactory.createPointGroup("6/mmm")
'''
'''
#Beta lactamase mutant
sampleRuns = range(5921,5931+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'beta_lac_3D_mbvg_mutant_newppl'
descriptorTOF = 'beta_lac_lab_highres_mut2'
peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
ellipseFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 32 2 1")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''

'''
#pth
sampleRuns = [870,872,873,874,875,876]
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'pth_3d'
descriptorTOF = 'pth_tof_secondRun'
peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/mandi_pth/peaks_combined.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 61 2 2")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''

'''
#gfp
sampleRuns = range(599,607+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'gfp_3d_goodhkl'
descriptorTOF = 'gfp_tof_goodhkl'
peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
ellipseFile = '/SNS/users/ntv/integrate/gfp/combined.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 21 21 2")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''

'''
#Beta lactamase
sampleRuns = range(4999,5004)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'beta_lac_3D_lab_highres2'
descriptorTOF = 'beta_lac_lab_highres2'
peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
#ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
ellipseFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 32 2 1")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''

'''
#DNA
sampleRuns = range(8758,8769+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'dna_3D_highres_mbvg'
descriptorTOF = 'dna_tof_highres'
#peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/mandi_dna2/combined_1p5A.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 21 21 21")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''

'''
#secondDNA
sampleRuns = range(5280,5291+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'secondDNA_3D'
descriptorTOF = 'secondDNA_tof'
#peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/mandi_secondDNA/combined.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 1 21 1")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''



'''
#NaK
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'nak_3D_full_lab'
descriptorTOF = 'nak_predpws5_lab'
sampleRuns = range(8275,8282+1)
ellipseFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
sg = SpaceGroupFactory.createSpaceGroup("I 4") 
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''
'''
#cryo
sampleRuns = range(8785,8791+1)
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'cryo_3d_mbvg'
descriptorTOF = 'cryo_tof_2'
peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
#peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptorTOF, sampleRuns[-1], descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/mandi_cryo/cryo_combined_2.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 21 21 21")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
'''


#--------------------------------Load everything
peaks_ws = ICAT.getPeaksWS(peaksFile, wsName='peaks_ws', forceLoad=True)
peaks_ws2 = ICAT.getPeaksWS(ellipseFile,wsName='peaks_ws2', forceLoad=True)
instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
panelDict = EdgeTools.getPanelDictionary(instrumentFile)
if descriptorTOF is not None:
    fitParams = ICAT.getFitParameters(workDir, descriptorTOF, sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
    fitDict = ICAT.getFitDicts(workDir, descriptorTOF,sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)

    #--------------------------------Create the TOF dataframe 
    dfTOF = pdTOF.getDFForPeaksWS(peaks_ws, fitParams, fitDict, panelDict, pg)
    mandiSpec = np.loadtxt('/SNS/users/ntv/integrate/MANDI_current.dat',skiprows=1,delimiter=',')
    sortedIDX = np.argsort(mandiSpec[:,0])
    intensSpec = interp1d(mandiSpec[sortedIDX[::3],0], mandiSpec[sortedIDX[::3],1],kind='cubic')
    dfTOF['scaledIntens'] = dfTOF['Intens'] / intensSpec(dfTOF['Wavelength']) * np.mean(mandiSpec[:,1])
    dfTOF['scaledIntensEll'] = dfTOF['IntensEll'] / intensSpec(dfTOF['Wavelength']) * np.mean(mandiSpec[:,1])
    dfTOF['lorentzFactor'] = np.sin(0.5*dfTOF['Scattering'])**2 / dfTOF['Wavelength']**4
    dfTOF['lorentzInt'] = dfTOF['Intens']*dfTOF['lorentzFactor']
    dfTOF['lorentzSig'] = dfTOF['SigInt']*dfTOF['lorentzFactor']
    dfTOF['lorentzIntEll'] = dfTOF['IntensEll']*dfTOF['lorentzFactor']
    dfTOF['lorentzSigEll'] = dfTOF['SigEll']*dfTOF['lorentzFactor']
else:
    dfTOF = pdTOF.getDFForPeaksWS(peaks_ws, None, None, panelDict, pg)

pdTOF.addPeaksWS(dfTOF,peaks_ws2)#add ellipsoid I, sigma(I)    
#---------------------------------Create the BVG dataframe
bvgParamFiles = [workDir + descriptorBVG + '/bvgParams_%i_%s.pkl'%(sampleRun, descriptorBVG) for sampleRun in sampleRuns]
bvgParams = []
for bvgFile in bvgParamFiles:
    bvgParams += pickle.load(open(bvgFile,'rb'))
dfBVG = pd.DataFrame(bvgParams)

#--------------------------Join the two dataframes and create some geometric columns
dfBVG['PeakNumber'] = dfBVG['peakNumber']

df = dfTOF.merge(dfBVG, on='PeakNumber')
df['theta'] = df['QLab'].apply(lambda x: np.arctan2(x[2],np.hypot(x[0],x[1])))
df['phi'] = df['QLab'].apply(lambda x: np.arctan2(x[1],x[0]))

#---------------------------Let's get some trends
plt.close('all')
goodIDX = (df['Intens']*5 > df['Intens3d']) & (df['Intens']*1.0/5.0 < df['Intens3d'])
goodIDX = goodIDX & (df['Intens']<1.0e7) & (df['Intens']>250) 
goodIDX = goodIDX & (df['SigX'] < 0.0199) & (df['SigY'] < 0.0199) #& ~(df['PeakNumber'].apply(lambda x: x in badPeaks))
/
def fSigX(x,a,k,x0,b):
    return a*np.exp(-k*(x-x0)) + b
    #return a/(x-x0) + b

def fSigP(x,a,k,phi,b):
    return a*np.sin(k*x-phi) + b*x

graph1 = sns.jointplot(df[goodIDX]['theta'], np.abs(df[goodIDX]['SigX']),s=1)
#pX = np.polyfit(df[goodIDX]['phi'], np.abs(df[goodIDX]['SigX']),4)
x = np.linspace(df[goodIDX]['phi'].min(), df[goodIDX]['phi'].max(), 100)
#y = np.polyval(pX, x)
pX,cov = curve_fit(fSigX, df[goodIDX]['phi'], df[goodIDX]['SigX'],maxfev=10000,p0=[0.02,0.4,0.05,0.005])
graph1.x = x; graph1.y = fSigX(x,pX[0],pX[1],pX[2],pX[3]); graph1.plot_joint(plt.plot)

graph2 = sns.jointplot(df[goodIDX]['phi'], df[goodIDX]['SigY'],s=1)
pY = np.polyfit(df[goodIDX]['theta'], np.abs(df[goodIDX]['SigY']),2)
x = np.linspace(-0.5, 2.5, 100)
y = np.polyval(pY, x)


graph3 = sns.jointplot(df[goodIDX]['theta'], df[goodIDX]['SigP'],s=1)
pX,cov = curve_fit(fSigP, df[goodIDX]['theta'], df[goodIDX]['SigP'],maxfev=10000,p0=[0.25,2.,0.,0.])
x = np.linspace(df[goodIDX]['theta'].min(), df[goodIDX]['theta'].max(), 100)
graph3.x = x; graph3.y = fSigP(x,pX[0],pX[1],pX[2],pX[3]); graph3.plot_joint(plt.plot)


#graph2.x = x; graph2.y = y; graph2.plot_joint(plt.plot)


#====to save
#r = np.array(df[goodIDX][['phi', 'theta', 'scale3d', 'MuTH', 'MuPH', 'SigX', 'SigY', 'SigP', 'PeakNumber']])
#pickle.dump(r, open('strongPeakParams.pkl', 'wb'))

#sns.jointplot(np.abs(df[goodIDX]['muTH']), np.abs(df[goodIDX]['SigX']),s=1)
'''
for bvgParam in bvgParams:
    for key in bvgParam.keys():
        df.iloc[bvgParam['peakNumber']][key] = bvgParam[key]

'''

#------------------------Scaled intensities
mandiSpec = np.loadtxt('MANDI_current.dat',skiprows=1,delimiter=',')
sortedIDX = np.argsort(mandiSpec[:,0])
intensSpec = interp1d(mandiSpec[sortedIDX[::3],0], mandiSpec[sortedIDX[::3],1],kind='cubic')
df['scaledIntens'] = df['Intens'] / intensSpec(df['Wavelength']) * np.mean(mandiSpec[:,1])
df['scaledIntensEll'] = df['IntensEll'] / intensSpec(df['Wavelength']) * np.mean(mandiSpec[:,1])
df['scaledIntens3d'] = df['Intens3d'] / intensSpec(df['Wavelength']) * np.mean(mandiSpec[:,1])
df['lorentzFactor'] = np.sin(0.5*df['Scattering'])**2 / df['Wavelength']**4
df['lorentzInt'] = df['Intens']*df['lorentzFactor']
df['lorentzSig'] = df['SigInt']*df['lorentzFactor']
df['lorentzIntEll'] = df['IntensEll']*df['lorentzFactor']
df['lorentzSigEll'] = df['SigEll']*df['lorentzFactor']
df['lorentzInt3d'] = df['Intens3d']*df['lorentzFactor']
df['lorentzSig3d'] = df['SigInt3d']*df['lorentzFactor']
#------------------------Outlier Detection
print ' '
print 'Detecting outliers in the same reflection families.  This can take some time.'


def isOutlier(intensities):
    if len(intensities) > 1:
        meanI = np.mean(intensities)
        stdI = np.std(intensities)
        isOutlier = (intensities > meanI + 4.0*stdI) + (intensities < meanI - 4.0*stdI) + (intensities > 2.0*meanI) + (intensities < 1.0/2.*meanI)
        return isOutlier
    else: return 0.0

checkIDX = (df['Intens3d'] > 3) & (df['chiSq']<50.0 ) & (df['chiSq3d']<10.) & (df['Intens3d']/df['SigInt3d'] > 1.)
df['isOutlier'] = 0.0

df.loc[checkIDX,'isOutlier'] = df[checkIDX].groupby('hklFam')['scaledIntens3d'].transform(isOutlier)
df['isOutlier'] = df['isOutlier'].astype(bool)
df['notOutlier'] = ~df['isOutlier']



#---------------------Select outputs and save a LaueNorm File
goodIDX = (df['chiSq'] < 50.0) & (df['chiSq3d']<10) & (df['notOutlier']) &  (df['Intens3d'] > -1.*np.inf) 
tooFarIDX = (np.abs(df['Intens3d'] > 100)) & ((np.abs(df['Intens3d']-df['IntensEll']) > 2.*df['Intens3d']) |  (np.abs(df['Intens3d']-df['IntensEll']) > 2.*df['Intens3d']) | (df['Intens3d'] > 5.*df['IntensEll']))
tooFarIDX2 = (np.abs(df['Intens3d'] < 100)) & (np.abs(df['Intens3d']-df['IntensEll']) >150)


#goodIDX = goodIDX & ~tooFarIDX #& ~tooFarIDX2
#goodIDX = goodIDX & (df['IntensEll'] != 0.0)
#goodIDX = goodIDX & ~((df['IntensEll']>100)&(df['Intens3d']<50))

dEdge = 1
edgeIDX = (df['Row'] <= dEdge) | (df['Row'] >= 255-dEdge) | (df['Col'] <= dEdge) | (df['Col'] >= 255-dEdge)
goodIDX = goodIDX & ~edgeIDX 

laueOutput = (df['DSpacing'] > 1.2) & (df['Wavelength'] > 2.0) & (df['Wavelength']<4.0) & (df['Intens3d']/df['SigInt3d'] > -3.0)
plt.figure(3); plt.clf();
plt.plot(df[goodIDX & laueOutput]['IntensEll'], df[goodIDX & laueOutput]['Intens3d'],'.',ms=2)

print ' '
print 'Removing bad peaks from peaks_ws.  This can take some time...'
#events = Load('/data/corelli_beryl/IPTS-20302/CORELLI_58417.nxs.h5')
events = Load('/SNS/MANDI/IPTS-16286/data/MANDI_6154_event.nxs')
#events = Load('/SNS/MANDI/IPTS-10943/0/870/NeXus/MANDI_870_event.nxs')
#events = Load('/SNS/MANDI/IPTS-15151/data/MANDI_5280_event.nxs')
#ws = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws", InstrumentWorkspace='events')
#ws2 = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws2", InstrumentWorkspace='events')
ws = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws")
ws2 = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws2")
peaksAdded = 0
peaks_ws_clone = CloneWorkspace(InputWorkspace=peaks_ws, OutputWorkspace='peaks_ws_clone')
peaks_ws_clone2 = CloneWorkspace(InputWorkspace=peaks_ws, OutputWorkspace='peaks_ws_clone2')

for i in range(len(df)):
    ICAT.print_progress(i,len(df),prefix='Cleaning df: ',suffix='Complete')
    try:
        if goodIDX[i] & laueOutput[i]:
            ws.addPeak(peaks_ws_clone.getPeak(df.iloc[i]['peakNumber']))
            ws.getPeak(peaksAdded).setIntensity(float(df.iloc[i]['Intens3d']))
            ws.getPeak(peaksAdded).setSigmaIntensity(float(df.iloc[i]['SigInt3d']))
            ws2.addPeak(peaks_ws_clone2.getPeak(df.iloc[i]['peakNumber']))
            ws2.getPeak(peaksAdded).setIntensity(float(df.iloc[i]['IntensEll']))
            ws2.getPeak(peaksAdded).setSigmaIntensity(float(df.iloc[i]['SigEll']))
            #ws.getPeak(peaksAdded).setIntensity(float(df.iloc[i]['lorentzInt3d']))
            #ws.getPeak(peaksAdded).setSigmaIntensity(float(df.iloc[i]['lorentzSig3d']))
            peaksAdded += 1
    except KeyError:
        raise
        pass

#print 'Saving LaueNorm Input'
#SaveLauenorm(InputWorkspace=ws, Filename=workDir+descriptorBVG+'/laue/laueNorm', ScalePeaks=3.0, minDSpacing=1.2, minWavelength=2.0, MaxWavelength=4.0, SortFilesBy='RunNumber', MinIsigI=1., MinIntensity=0)
print 'Saving LaueScale'
#TransformHKL(PeaksWorkspace='ws', HKLTransform='0,0,1,0,-1,0,1,0,0')

#SaveLauenorm(InputWorkspace=ws, Filename=workDir+descriptorBVG+'/laue_compare/lscale', ScalePeaks=3.0, minDSpacing=1.2, minWavelength=2.0, MaxWavelength=4.0, SortFilesBy='RunNumber', MinIsigI=1., MinIntensity=0.,LaueScaleFormat=True)
