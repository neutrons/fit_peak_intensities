import numpy as np
from datetime import datetime

def writeLog(logFile, workDir, loadDir, nxsTemplate, figsFormat, sampleRuns, dtSpread, dtBinWidth, fracHKL, refineCenter, doVolumeNormalization, peaksFile, UBFile, DetCalFile, descriptor):
    with open(logFile,'w') as f:
        f.write('---Log file for %s\n'%descriptor)
        f.write('Log created: ' + datetime.now().strftime('%m/%d/%Y %H:%M:%S')+'\n')
        f.write('workDir: %s\n'%workDir)
        f.write('loadDir: %s\n'%loadDir)
        f.write('nxsTemplate: %s\n'%nxsTemplate)
        f.write('figsFormat: %s\n'%figsFormat)
        f.write('sampleRuns: ' + ' '.join(['%i'%i for i in sampleRuns])+'\n')
        f.write('dtSpread: %4.6f\n'%dtSpread)
        f.write('dtBinWidth: %i\n'%dtBinWidth)
        f.write('fracHKL: %4.6f\n'%fracHKL)
        f.write('refineCenter: %i\n'%refineCenter)
        f.write('doVolumeNormalization: %i\n'%doVolumeNormalization)
        f.write('peaksFile: %s\n'%peaksFile)
        f.write('UBFile: %s\n'%UBFile)
        f.write('DetCalFile: %s\n'%DetCalFile)

