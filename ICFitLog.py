import numpy as np
from datetime import datetime
import sys

def writeLog(logFile, workDir, loadDir, nxsTemplate, figsFormat, sampleRuns, dtSpread, dtBinWidth, fracHKL, fracStop, removeEdges, doVolumeNormalization, peaksFile, UBFormat, DetCalFile, moderatorCoefficientsFile, descriptor,zBG,neigh_length_m, predpplCoefficients, minppl_frac, maxppl_frac):
    mantidPath = [s for s in sys.path if 'antid' in s]
    with open(logFile,'w') as f:
        f.write('---Log file for %s\n'%descriptor)
        f.write('Log created: ' + datetime.now().strftime('%m/%d/%Y %H:%M:%S')+'\n')
        f.write('workDir: %s\n'%workDir)
        f.write('loadDir: %s\n'%loadDir)
        f.write('nxsTemplate: %s\n'%nxsTemplate)
        f.write('figsFormat: %s\n'%figsFormat)
        f.write('sampleRuns: ' + ' '.join(['%i'%i for i in sampleRuns])+'\n')
        f.write('dtSpread: ' + str(dtSpread) + '\n')
        f.write('dtBinWidth: %i\n'%dtBinWidth)
        f.write('fracHKL: %4.6f\n'%fracHKL)
        f.write('fracStop: %4.6f\n'%fracStop)
        f.write('removeEdges: %i\n'%removeEdges)
        f.write('doVolumeNormalization: %i\n'%doVolumeNormalization)
        f.write('peaksFile: %s\n'%peaksFile)
        f.write('UBFormat: %s\n'%UBFormat)
        f.write('moderatorCoefficientsFile: %s\n'%moderatorCoefficientsFile)
        f.write('DetCalFile: %s\n'%DetCalFile)
        f.write('zBG: %4.6f\n'%zBG)
        f.write('neigh_length_m: %4.6f\n'%neigh_length_m)
        f.write('predpplCoefficients: %s\n'%str(predpplCoefficients))
        f.write('mantidPath: %s\n'%str(mantidPath))
        f.write('minppl_frac: %4.4f\n'%minppl_frac) 
        f.write('maxppl_frac: %4.4f\n'%maxppl_frac) 
