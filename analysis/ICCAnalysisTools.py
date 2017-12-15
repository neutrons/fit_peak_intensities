import numpy as np
from mantid.simpleapi import *
import matplotlib.pyplot as plt
plt.ion()
import glob
import sys
import warnings
import collections
import pickle
import re

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '?' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r %s%s %s' % ( percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def getPeaksWS(peaksFile,forceLoad=False, wsName='peaks_ws'):
    importPeaks = True
    if not forceLoad:
        for ws in mtd.getObjectNames():
            if mtd[ws].getComment() == '%s'%peaksFile and ws == wsName:
                print '    using already loaded peaks file'
                importPeaks = False
                peaks_ws = mtd[ws]
        if importPeaks:
            if wsName in mtd.getObjectNames():
                origWSName = wsName
                counter=1
                while wsName not in mtd.getObjectNames():
                    wsName = origWSName+'_%i'
                    i += 1
            peaks_ws = LoadIsawPeaks(Filename = peaksFile, OutputWorkspace=wsName)
            peaks_ws.setComment(peaksFile)
    else:
        peaks_ws = LoadIsawPeaks(Filename = peaksFile, OutputWorkspace=wsName)
    return peaks_ws 

def getFitDicts(baseDir,descriptor,firstRun,lastRun, sampleRuns=None):
    if sampleRuns is None:
        sampleRuns = range(firstRun, lastRun+1)
    dictFiles = [baseDir + descriptor + '/fitDict_%i_%s.pkl'%(i, descriptor) for i in sampleRuns]
    newDict = {}
    for f in dictFiles:
        td = pickle.load(open(f,'rb'))
        newDict.update(td)
    return newDict

def getFitParameters(baseDir, descriptor, firstRun, lastRun, sampleRuns=None):
    if sampleRuns is None:
        sampleRuns = range(firstRun, lastRun+1)
    paramFiles = [baseDir + descriptor + '/params_%i_%s.dat'%(i, descriptor) for i in sampleRuns]
    for i, fileName in enumerate(paramFiles):
        peakParams = np.loadtxt(fileName)
        if i == 0:
            peakParamsSet = peakParams
        else:
            peakParamsSet = np.vstack([peakParamsSet, peakParams])
    return peakParamsSet

def getImageFileNames(baseDir, descriptor, peakNumbers, peaksWS=None):
    if peaksWS is None:
        fileNames = [baseDir + descriptor + '/figs/*_%i.png'%i for i in peakNumbers]
    else:
        fileNames = [baseDir + descriptor + '/figs/mantid_%i_%i.png'%(peaksWS.getPeak(int(i)).getRunNumber(), i) for i in peakNumbers]
    return fileNames

def showImages(baseDir, descriptor, peakNumbers, peaksWS=None):
    imageFileNames = getImageFileNames(baseDir, descriptor, peakNumbers, peaksWS=peaksWS)
    os.system('display ' + ' '.join([im for im in imageFileNames]))


def gethklDict(peaksWS, peaksToGet=None):
    d = dict()
    dPeak = dict()
    if peaksToGet is None:
        peaksToGet = range(peaksWS.getNumberPeaks())
    for i in peaksToGet:
        p = peaksWS.getPeak(i)
        key = tuple(np.abs(p.getHKL()).astype(int))
        if d.has_key(key): #we have this, well just add 
                d[key] = np.append(d[key], p.getIntensity())
                dPeak[key] = np.append(dPeak[key], i)
        else: #we make the key
                d[key] = np.array(p.getIntensity())
                dPeak[key] = np.array(i)
    return dPeak, d

#def findSameRejectedHKL(peaksWS, rejectedPeaks):
#    for peak in rejectedPeaks

def getREFEDTRejects(peaksWS, lstFile):
    #---------
    minRun = 1.0e10

    for i in range(peaksWS.getNumberPeaks()):
        if  peaksWS.getPeak(i).getRunNumber() < minRun:
            minRun = peaksWS.getPeak(i).getRunNumber()

    #Now read the lst file
    lstData = list()
    with open(lstFile,'r') as f:
        for row in f:
            if len(row.split()) == 12:
                lstData.append(row.split())
                #print row.split()
    lstHeader = lstData[0]
    lstData = np.array(lstData[1:-1]) #First is header, last is text

    rejectedPeakNumbers = list()
    print 'Finding rejectd peaks.  This can take a few minutes...'
    print_progress(0, len(lstData), prefix = 'Progress:', suffix = 'Complete')
    lowI = 0
    for iLst, lst in enumerate(lstData):
        print_progress(iLst, len(lstData), prefix = 'Progress:', suffix = 'Complete')
        try:
            thisRun = int(lst[0].strip('.'))-1+minRun
            h = int(lst[1].strip('.'))
            k = int(lst[2].strip('.'))
            l = int(lst[3].strip('.'))
            hkl = np.array([h,k,l])
            for i in range(lowI, peaksWS.getNumberPeaks()):
                if np.all(np.abs(hkl) == np.abs(peaksWS.getPeak(i).getHKL())) and (thisRun == peaksWS.getPeak(i).getRunNumber()):
                    rejectedPeakNumbers.append(i)
                    lowI = max(lowI,i)
                    break
        except KeyboardInterrupt:
            break
        except:
            print 'Error with ', lst
    print_progress(iLst, len(lstData), prefix = 'Progress:', suffix = 'Complete')
    return rejectedPeakNumbers

def getRejectionReason(lstFile):
    #Now read the lst file
    lstData = list()
    with open(lstFile,'r') as f:
        for row in f:
            if len(row.split()) == 12:
                lstData.append(row.split())
                #print row.split()
    lstHeader = lstData[0]
    lstData = np.array(lstData[1:-1]) #First is header, last is text

    rejectedPeakNumbers = list()
    numF0FC = 0
    numFCF0 = 0
    numF0SIGF0 = 0
    rejectedReason = list()
    for iLst, lst in enumerate(lstData):
        try:
            F0SQ = float(lst[6])
            SIGF0SQ = float(lst[7])
            FCSQ = float(lst[8])
            if 1.0*F0SQ/FCSQ > 4.0:
                numF0FC+=1
                rejectedReason.append(0)
                continue
            elif 1.0*FCSQ/F0SQ > 4.0:
                numFCF0+=1
                rejectedReason.append(1)
                continue
            elif 1.0*np.abs(F0SQ-FCSQ) / SIGF0SQ > 3.0:
                numF0SIGF0+=1
                rejectedReason.append(2)
                continue
            Warning( 'CANNOT FIND WHY THIS PEAK WAS REJECTED')
            rejectedReason.append(-1)
     

        except KeyboardInterrupt:
            break
        except:
            raise
            print 'Error with ', lst
    print 'Peaks rejected for F0SQ/FCSQ>4.0:', numF0FC
    print 'Peaks rejected for FCSQ/F0SQ>4.0:', numFCF0
    print 'Peaks rejected for |FCSQ-F0SQ|/SIGF0SQ < 3:', numF0SIGF0
    return rejectedReason


def getSameHKLRejects(peaks_ws, peakNumber, rejectedPeaks, dP=None):
    peak = peaks_ws.getPeak(peakNumber)
    hkl = peak.getHKL()
    if dP is None:
        dP, dI = gethklDict(peaks_ws)
    sameHKLPeaks = dP[tuple(np.abs(hkl))]
    rejects = []; approved = []
    for peak in sameHKLPeaks:
        if peak in rejectedPeaks:
            rejects.append(peak)
        else:
            approved.append(peak)
    return np.array(rejects), np.array(approved)

def getIntensities(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getIntensity())
    return np.array(I1)

def getSigmas(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getSigmaIntensity())
    return np.array(I1)

def getScattering(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getScattering())
    return np.array(I1)

def getFlightPath(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getL1() + peaks_ws.getPeak(int(peak)).getL2())
    return np.array(I1)

def getRunNumber(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getRunNumber())
    return np.array(I1)

def getWavelength(peaks_ws, peaksToGet):
    I1 = []
    for peak in peaksToGet:
        I1.append(peaks_ws.getPeak(int(peak)).getWavelength())
    return np.array(I1)


class lstPeak():
    def __init__(self):
        self.H = 0
        self.K = 0
        self.L = 0
        self.F0 = 0
        self.FC = 0
        self.sigF0 = 0
        self.histogram = 0

    def setIndex(self, index, val):
        if index == 'H':
            self.H = int(val)
            return
        elif index == 'K':
            self.K = int(val)
            return
        elif index == 'L':
            self.L = int(val)
            return
        print index, index == 'H'

def readHKLSortList(fileName, fastIndex='H', slowIndex='K', thirdIndex='L'):
    currSlowIndex = 0
    currThirdIndex = 0
    peakList = []
    with open(fileName, 'r') as f:
        for line in f:
            #if 'Column 1' in line:
                #fastIndex = line.split(' ')[-1]
            if slowIndex+'=' in line:
                #print line
                tmp = re.findall('[+-]?\d+', line)
                currSlowIndex=int(tmp[0])
                currThirdIndex=int(tmp[1])
            else:
                tmp = re.findall('[+-]?\d+', line)
                if len(tmp) == 6: #its a peak
                    p = lstPeak()
                    currFastIndex = int(tmp[0])
                    p.setIndex(fastIndex, currFastIndex)
                    p.setIndex(slowIndex, currSlowIndex)
                    p.setIndex(thirdIndex, currThirdIndex) 
                    p.F0 = int(tmp[1])
                    p.FC = int(tmp[2])
                    p.sigF0 = float(tmp[3])*(10**(float(tmp[4])))
                    p.histogram = int(tmp[5])
                    peakList.append(p)
    return peakList

def peakListToNumpy(peakList,minRun=None):
    for i, peak in enumerate(peakList):
        tmp = np.array([peak.H, peak.K, peak.L, peak.F0, peak.FC, peak.sigF0, peak.histogram],dtype=np.int)
        if i==0:
            npp = tmp
        else:
            npp = np.vstack([npp, tmp])
    if minRun != None:
        npp[:,-1] = npp[:,-1]+minRun-1
    return npp

def showFit(peakDict, peakNumbers, figNumber=1,paramVects=None):

    if not isinstance(peakNumbers,collections.Iterable):
        peakNumbers = [peakNumbers]
    if paramVects is None:    
        paramVects = [None for i in range(len(peakNumbers))]

    for peakNumber,paramVect in zip(peakNumbers, paramVects):
        try: 
            plt.figure(figNumber)
            plt.clf()
            x = peakDict[peakNumber][0]
            yData = peakDict[peakNumber][1]
            yFit = peakDict[peakNumber][2]
            plt.plot(x,yData,'o')
            plt.plot(x,yFit,'-o')
            plt.xlabel('TOF (us)')
            plt.ylabel('Counts')
            if paramVect is None:
                plt.title('Peak %i'%peakNumber)
            else:
                plt.title('Peak %i - %f eV - chiSQ=%f'%(peakNumber, paramVect[1]*1000, paramVect[-1]))
            if len(peakNumbers) > 1:
                plt.pause(1)
        except KeyboardInterrupt:
            print 'KEYBOARD INTERRUPT!! Breaking loop!'
            break
        except:
            print 'ICAT::showFit - Cannot plot peak %i'%peakNumber
