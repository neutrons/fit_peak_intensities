import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import convertToPandas as pdTOF
import ICCAnalysisTools as ICAT
from scipy.interpolate import interp1d
import getEdgePixels as EdgeTools
import seaborn as sns

#Load the bvgFit files
workDir = '/SNS/users/ntv/dropbox/'
descriptorBVG = 'beta_lac_3D_full'
descriptorTOF = 'beta_lac_predpws5'

peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptorTOF)
ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
sg = SpaceGroupFactory.createSpaceGroup("P 32 2 1")
pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)
sampleRuns = range(4999,5004)

#--------------------------------Load everything
peaks_ws = ICAT.getPeaksWS(peaksFile, wsName='peaks_ws')
peaks_ws2 = ICAT.getPeaksWS(ellipseFile,wsName='peaks_ws2')
fitParams = ICAT.getFitParameters(workDir, descriptorTOF, sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
fitDict = ICAT.getFitDicts(workDir, descriptorTOF,sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
panelDict = EdgeTools.getPanelDictionary(instrumentFile)

#--------------------------------Create the TOF dataframe 
dfTOF = pdTOF.getDFForPeaksWS(peaks_ws, fitParams, fitDict, panelDict, pg)
pdTOF.addPeaksWS(dfTOF,peaks_ws2)#add ellipsoid I, sigma(I)    
mandiSpec = np.loadtxt('MANDI_current.dat',skiprows=1,delimiter=',')
sortedIDX = np.argsort(mandiSpec[:,0])
intensSpec = interp1d(mandiSpec[sortedIDX[::3],0], mandiSpec[sortedIDX[::3],1],kind='cubic')
dfTOF['scaledIntens'] = dfTOF['Intens'] / intensSpec(dfTOF['Wavelength']) * np.mean(mandiSpec[:,1])
dfTOF['scaledIntensEll'] = dfTOF['IntensEll'] / intensSpec(dfTOF['Wavelength']) * np.mean(mandiSpec[:,1])
dfTOF['lorentzFactor'] = np.sin(dfTOF['Scattering'])**2 / dfTOF['Wavelength']**4
dfTOF['lorentzInt'] = dfTOF['Intens']*dfTOF['lorentzFactor']
dfTOF['lorentzSig'] = dfTOF['SigInt']*dfTOF['lorentzFactor']
dfTOF['lorentzIntEll'] = dfTOF['IntensEll']*dfTOF['lorentzFactor']
dfTOF['lorentzSigEll'] = dfTOF['SigEll']*dfTOF['lorentzFactor']

#---------------------------------Create the BVG dataframe
sampleRuns = [4999,5000, 5001, 5002, 5003]
bvgParamFiles = [workDir + descriptorBVG + '/bvgParams_%i_%s.pkl'%(sampleRun, descriptorBVG) for sampleRun in sampleRuns]
bvgParams = []
for bvgFile in bvgParamFiles:
    bvgParams += pickle.load(open(bvgFile,'rb'))
dfBVG = pd.DataFrame(bvgParams)

#--------------------------Join the two dataframes and create some geometric columns
dfBVG['PeakNumber'] = dfBVG['peakNumber']
df = dfTOF.merge(dfBVG, on='PeakNumber')
df['phi'] = df['QSample'].apply(lambda x: np.arctan2(x[2],np.hypot(x[0],x[1])))
df['theta'] = df['QSample'].apply(lambda x: np.arctan2(x[1],x[0]))

#---------------------------Let's get some trends
plt.close('all')
goodIDX = (df['Intens']*5 > df['Intens3d']) & (df['Intens']*1.0/5.0 < df['Intens3d'])
goodIDX = goodIDX & (df['Intens']<1.0e7) & (df['Intens']>1)
goodIDX = goodIDX & (df['sigX'] < 0.0199) & (df['sigY'] < 0.0199) #& ~(df['PeakNumber'].apply(lambda x: x in badPeaks))

graph1 = sns.jointplot(df[goodIDX]['phi'], np.abs(df[goodIDX]['sigX']),s=1)
pX = np.polyfit(df[goodIDX]['phi'], np.abs(df[goodIDX]['sigX']),4)
x = np.linspace(-1.5, 1.5, 100)
y = np.polyval(pX, x)
graph1.x = x; graph1.y = y; graph1.plot_joint(plt.plot)

graph2 = sns.jointplot(df[goodIDX]['theta'], np.abs(df[goodIDX]['sigY']),s=1)
pY = np.polyfit(df[goodIDX]['theta'], np.abs(df[goodIDX]['sigY']),4)
x = np.linspace(-0.5, 2.5, 100)
y = np.polyval(pY, x)
graph2.x = x; graph2.y = y; graph2.plot_joint(plt.plot)

#====to save
#r = np.array(df[goodIDX][['theta', 'phi', 'scale3d', 'muTH', 'muPH', 'sigX', 'sigY', 'sigP']])
#pickle.dump(r, open('strongPeakParams.pkl', 'wb'))

#sns.jointplot(np.abs(df[goodIDX]['muTH']), np.abs(df[goodIDX]['sigX']),s=1)
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
df['lorentzFactor'] = np.sin(df['Scattering'])**2 / df['Wavelength']**4
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
        isOutlier = (intensities > meanI + 4.0*stdI) + (intensities < meanI - 4.0*stdI) + (intensities > 3.0*meanI) + (intensities < 1.0/3.0*meanI)
        return isOutlier
    else: return 0.0

checkIDX = (df['Intens3d'] > 0) & (df['chiSq']<50.0)
df['isOutlier'] = 0.0

df.loc[checkIDX,'isOutlier'] = df[checkIDX].groupby('hklFam')['scaledIntens3d'].transform(isOutlier)
df['isOutlier'] = df['isOutlier'].astype(bool)
df['notOutlier'] = ~df['isOutlier']



#---------------------Select outputs and save a LaueNorm File
goodIDX = (df['chiSq'] < 100.0) & (df['Intens3d'] > 1)  & (df['Intens3d']<3.0e7) & (df['chiSq3d'] < 5.0) & (df['notOutlier'])
tooFarIDX = (np.abs(df['Intens3d'] > 100)) & ((np.abs(df['Intens3d']-df['IntensEll']) > 4.0*df['Intens3d']) |  (np.abs(df['Intens3d']-df['IntensEll']) > 4.0*df['Intens3d']))
#goodIDX = goodIDX #& ~tooFarIDX

dEdge = 3
edgeIDX = (df['Row'] <= dEdge) | (df['Row'] >= 255-dEdge) | (df['Col'] <= dEdge) | (df['Col'] >= 255-dEdge)
goodIDX = goodIDX & ~edgeIDX 

goodIDX = goodIDX & (df['sigX'] < 0.0199)& (df['sigY'] < 0.0199) & (df['theta'] > -1.0)

p = np.array([1.14757325, -746.16871354])

plt.figure(3); plt.clf();
plt.plot(df[goodIDX]['Intens3d'], df[goodIDX]['IntensEll'],'.',ms=3)

laueOutput = (df['DSpacing'] > 2.0) & (df['Wavelength'] > 2.0) & (df['Wavelength']<4.0) & (df['Intens']/df['SigInt'] > 1.0)
print ' '
print 'Removing bad peaks from peaks_ws.  This can take some time...'
ws = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws")
peaksAdded = 0
peaks_ws_clone = CloneWorkspace(InputWorkspace=peaks_ws, OutputWorkspace='peaks_ws_clone')
for i in range(peaks_ws.getNumberPeaks()):
    ICAT.print_progress(i,peaks_ws.getNumberPeaks(),prefix='Cleaning df: ',suffix='Complete')
    try:
        if goodIDX[i]:
            ws.addPeak(peaks_ws_clone.getPeak(i))
            ws.getPeak(peaksAdded).setIntensity(float(df.loc[i]['Intens3d']))
            ws.getPeak(peaksAdded).setSigmaIntensity(float(df.loc[i]['SigInt3d']))
            peaksAdded += 1
    except KeyError:
        #raise
        pass

print 'Saving LaueNorm Input'
SaveLauenorm(InputWorkspace=ws, Filename=workDir+descriptorBVG+'/laue/laueNorm', ScalePeaks=3000.0, minDSpacing=2.0, minWavelength=2.0, MaxWavelength=4.0, SortFilesBy='RunNumber', MinIsigI=1, MinIntensity=0)

