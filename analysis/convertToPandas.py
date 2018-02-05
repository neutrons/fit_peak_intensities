import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import pandas as pd
import sys
sys.path.append('../')
import getEdgePixels as EdgeTools
import ICCAnalysisTools as ICAT
#import xgboost as xgb
#from sklearn.model_selection import train_test_split
import ICConvoluted as ICC
import ICCFitTools as ICCFT
from scipy.stats import pearsonr as corr
from scipy.interpolate import interp1d
from mantid.simpleapi import *
from mantid.geometry import PointGroupFactory, SpaceGroupFactory
FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)
reload(ICCFT)


'''
def getPGFam(df):
pgFam = []
for idx,row in df.iterrows():
    hkl = tuple(row[['h', 'k', 'l']])
    pgFam.append(pg.getReflectionFamily(hkl))
return pgFam
df['hklFam'] = getPGFam(df)
'''
def getDictForPandasPeak(peak, fitParams, fitDict, panelDict, i, pg):
    colFunDict = {
        'PeakNumber':'i',
        'RunNumber':'peak.getRunNumber()',  
        'DetID':'peak.getDetectorID()',  
        'h':'peak.getH()',  
        'k':'peak.getK()',  
        'l':'peak.getL()',  
        'Wavelength':'peak.getWavelength()',  
        'Energy': 'peak.getFinalEnergy()',    
        'TOF':'peak.getTOF()',  
        'DSpacing':'peak.getDSpacing()',  
        'Intens':'peak.getIntensity()',  
        'SigInt':'peak.getSigmaIntensity()',  
        'BinCount':'peak.getBinCount()', 
        'Scattering':'peak.getScattering()', 
        'BankName': 'EdgeTools.getDetectorBank(panelDict,peak.getDetectorID())[\'bankNumber\']', 
        'Row':'peak.getRow()',  
        'Col':'peak.getCol()',  
        'QLab':'np.array(peak.getQLabFrame())',  
        'QSample':'np.array(peak.getQSampleFrame())', 
        'absHKL':'tuple(np.abs(peak.getHKL()).astype(int))',
        'Alpha':'fitParams[i,5]',  
        'Beta':'fitParams[i,6]',  
        'R':'fitParams[i,7]',  
        'T0':'fitParams[i,8]',  
        'scale':'fitParams[i,9]',  
        'hatwidth':'fitParams[i,10]',  
        'k_conv':'fitParams[i,11]',
        'bg_linear':'fitParams[i,13]',  
        'bg_const':'fitParams[i,12]',  
        'chiSq':'fitParams[i,4]',  
        'pp_lambda':'fitParams[i,-1]',
        't':'fitDict[i][0]',  
        'yData':'fitDict[i][1]',  
        'yFit':'fitDict[i][2]',
        'hklFam': 'tuple(np.array(pg.getReflectionFamily(peak.getHKL())))'}

    try:
        if len(fitParams[0]) == 17:  
            colFunDict['bg_quad'] = 'fitParams[i,14]'
    except: pass
    d = {}
    for key in colFunDict.keys():
        command = 'd[\'%s\'] = %s'%(key, colFunDict[key])
        try: exec(command)
        except KeyError:
            pass
        except:
            pass#print command 
            #raise#print command
    d['PeakNumber'] = i
    return d

def getDictForPandasPeaksWS(peaks_ws,fitParams,fitDict,panelDict,pg):
    dictList = []
    print 'Assembling df, this can take a little while:'
    ICAT.print_progress(0,peaks_ws.getNumberPeaks(),prefix='Assumbling df: ',suffix='Complete')
    for i in range(peaks_ws.getNumberPeaks()):
        ICAT.print_progress(i,peaks_ws.getNumberPeaks(),prefix='Assumbling df: ',suffix='Complete')
        peak = peaks_ws.getPeak(i)
        dictList.append(getDictForPandasPeak(peak, fitParams, fitDict,panelDict,i,pg))
    print ' '
    return dictList

def getDFForPeaksWS(peaks_ws,fitParams,fitDict,panelDict,pg):
    dictList = getDictForPandasPeaksWS(peaks_ws,fitParams,fitDict, panelDict,pg)
    return pd.DataFrame(dictList)

 
def addCountsNoFit(df):
    #add raw counts
    ctsList = []
    sigList = []
    for idx,row in df.iterrows():
        x = row['t']
        yF = row['yFit']
        yD = row['yData']
        if np.sum(np.isnan(yF)) == 0:
            bgCoeff = np.array([row['bg_linear'], row['bg_const']])
            bgY = np.polyval(bgCoeff,x)
            t = ICCFT.integratePeak(x,yF,yD,bgY)
            peaks_ws.getPeak(int(idx)).setSigmaIntensity(t[1])
            row['SigInt'] = t[1]
            sumIDX = np.logical_and(x>t[2],x<t[3])
            sumCts = yD-bgY
            sumCts = np.sum(sumCts[sumIDX])
            ctsList.append(sumCts)
            sig = np.sqrt(np.sum(yD[sumIDX]) + np.sum(bgY[sumIDX]))
            sigList.append(sig)
        else:   
            ctsList.append(0)
            sigList.append(1)
    df['IntensCts'] = ctsList
    df['SigCts'] = sigList

#Adds a second workspace's intensities and sigma to the df
#must have the same number of peaks
def addPeaksWS(df, peaks_ws):
    IELL = []
    SELL = []
    for i in range(peaks_ws.getNumberPeaks()):
        IELL.append(peaks_ws.getPeak(i).getIntensity())
        SELL.append(peaks_ws.getPeak(i).getSigmaIntensity())
    df['IntensEll'] = IELL
    df['SigEll'] = SELL


if __name__ == '__main__':
    print '----ARE YOU SURE YOU WANT TO RUN THE WHOLE THING? (<Y>/n)'
    inp = raw_input()
    if inp.upper() != 'Y':
        sys.exit(0)

    '''
    # Parameters for our run
    workDir = '/SNS/users/ntv/dropbox/'
    #descriptor = 'scol_removeBG'
    descriptor = 'scol_dynamicQuadBG'
    sampleRuns = range(15629,  15644)
    peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptor,sampleRuns[-1],descriptor)
    ellipseFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
    #descriptor = 'si_0p015_removeBG'
    #sampleRuns = range(15647,15670)
    '''
    #beta lac
    workDir = '/SNS/users/ntv/dropbox/'
    descriptor = 'beta_lac_predpws5'
    sampleRuns = range(4999,5004)
    ellipseFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
    #peaksFile = '%s%s/peaks_%i_%s.integrate'%(workDir,descriptor,sampleRuns[-1],descriptor)
    peaksFile = '%s%s/peaks_combined_good.integrate'%(workDir,descriptor)
    sg = SpaceGroupFactory.createSpaceGroup("P 32 2 1") 
    pg = PointGroupFactory.createPointGroupFromSpaceGroup(sg)


    #--------------------------------Load everything
    peaks_ws = ICAT.getPeaksWS(peaksFile, wsName='peaks_ws')
    peaks_ws2 = ICAT.getPeaksWS(ellipseFile,wsName='peaks_ws2')
    fitParams = ICAT.getFitParameters(workDir, descriptor, sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
    fitDict = ICAT.getFitDicts(workDir, descriptor,sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
    instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
    panelDict = EdgeTools.getPanelDictionary(instrumentFile)

    #--------------------------------Create the dataframe including 
    #   poor man's wavelength normalization and lorentz corrected intensities
    df = getDFForPeaksWS(peaks_ws, fitParams, fitDict, panelDict, pg)
    addPeaksWS(df,peaks_ws2)#add ellipsoid I, sigma(I)    
    
    mandiSpec = np.loadtxt('MANDI_current.dat',skiprows=1,delimiter=',')
    sortedIDX = np.argsort(mandiSpec[:,0])
    intensSpec = interp1d(mandiSpec[sortedIDX[::3],0], mandiSpec[sortedIDX[::3],1],kind='cubic') 
    df['scaledIntens'] = df['Intens'] / intensSpec(df['Wavelength']) * np.mean(mandiSpec[:,1])
    df['scaledIntensEll'] = df['IntensEll'] / intensSpec(df['Wavelength']) * np.mean(mandiSpec[:,1])
    df['lorentzFactor'] = np.sin(df['Scattering'])**2 / df['Wavelength']**4 
    df['lorentzInt'] = df['Intens']*df['lorentzFactor']
    df['lorentzSig'] = df['SigInt']*df['lorentzFactor']
    df['lorentzIntEll'] = df['IntensEll']*df['lorentzFactor']
    df['lorentzSigEll'] = df['SigEll']*df['lorentzFactor']

    #-----------------Remove outliers in same hkl
    print ' '
    print 'Detecting outliers in the same reflection families.  This can take some time.'
   

    def isOutlier(intensities):
        if len(intensities) > 1:
            meanI = np.mean(intensities)
            stdI = np.std(intensities)
            isOutlier = (intensities > meanI + 4.0*stdI) + (intensities < meanI - 4.0*stdI) + (intensities > 2.5*meanI) + (intensities < 1.0/2.5*meanI)
            return isOutlier
        else: return 0.0    

    checkIDX = (df['Intens'] > 0) & (df['chiSq']<50.0)
    df['isOutlier'] = 0.0
    #df.loc[checkIDX,'isOutlier'] = df[checkIDX].groupby('hklFam')['Intens'].transform(isOutlier)
    #df.loc[checkIDX,'isOutlier'] = df[checkIDX].groupby('hklFam')['lorentzInt'].transform(isOutlier)
    df.loc[checkIDX,'isOutlier'] = df[checkIDX].groupby('hklFam')['scaledIntens'].transform(isOutlier)
    df['isOutlier'] = df['isOutlier'].astype(bool)
    df['notOutlier'] = ~df['isOutlier']

    df.loc[checkIDX,'isOutlierEll'] = df[checkIDX].groupby('hklFam')['IntensEll'].transform(isOutlier)
    df['isOutlierEll'] = df['isOutlierEll'].astype(bool)
    df['notOutlierEll'] = ~df['isOutlierEll']
 
    df['predppl'] = ICCFT.oldScatFun(df['Scattering']/df['Wavelength'],5.24730283,  7.23719321,  0.27449887)
 
    #----------------Set the output conditions and make files 
    goodIDX = (df['chiSq'] < 50.0) & (df['Intens'] > 0) & (df['notOutlier'])   & (df['Intens']<3.0e7) 
    tooFarIDX = (df['Intens'] > 100) & ((np.abs(df['Intens']-df['IntensEll']) > 2.0*df['IntensEll']) |  (np.abs(df['Intens']-df['IntensEll']) > 2*df['Intens']))

    goodIDX = goodIDX & ~tooFarIDX
    dEdge = 3 
    edgeIDX = (df['Row'] <= dEdge) | (df['Row'] >= 255-dEdge) | (df['Col'] <= dEdge) | (df['Col'] >= 255-dEdge) 
    goodIDX = goodIDX & ~edgeIDX
    
    plt.figure(2); plt.clf();
    plt.plot(df[goodIDX]['Intens'], df[goodIDX]['IntensEll'],'.',ms=3)
    #plt.plot(df[~goodIDX]['IntensEll'], 1/0.7*df[~goodIDX]['IntensEll'],'.',ms=3)

    laueOutput = (df['DSpacing'] > 2.0) & (df['Wavelength'] > 2.0) & (df['Wavelength']<4.0) & (df['Intens']/df['SigInt'] > 1.0) 
    print ' '
    print 'Removing bad peaks from peaks_ws.  This can take some time...'
    ws = CreatePeaksWorkspace(NumberOfPeaks=0, OutputWorkspace="ws")
    peaksAdded = 0
    peaks_ws_clone = CloneWorkspace(InputWorkspace=peaks_ws, OutputWorkspace='peaks_ws_clone')
    for i in range(peaks_ws.getNumberPeaks()):
        ICAT.print_progress(peaks_ws.getNumberPeaks()-i,peaks_ws.getNumberPeaks(),prefix='Cleaning df: ',suffix='Complete')
        if goodIDX[i]:
            ws.addPeak(peaks_ws_clone.getPeak(i))
            #ws.getPeak(peaksAdded).setIntensity(float(df.iloc[i]['lorentzInt']))
            #ws.getPeak(peaksAdded).setSigmaIntensity(float(df.iloc[i]['lorentzSig']))
            peaksAdded += 1
   
    print 'Saving LaueNorm Input'
    SaveLauenorm(InputWorkspace=ws, Filename=workDir+descriptor+'/laue/laueNorm', ScalePeaks=3000.0, minDSpacing=2.0, minWavelength=2.0, MaxWavelength=4.0, SortFilesBy='RunNumber', MinIsigI=1, MinIntensity=0) 

