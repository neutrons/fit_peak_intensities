def doBVGFits(sampleRunsList=None):
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
    peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.integrate'
    UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/295K_predict_2016A/SC295K_Monoclinic_C.mat'
    DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
    descriptor = 'scol_3d' #Does not end with '/'
    numTimesToInterpolate=1
    dtBinWidth = 4 #Width (in us) in TOF profile bins
    qLow = -25; qHigh=25;
    '''
    '''
    # CORELLI - beryl
    loadDir = '/data/corelli_beryl/IPTS-20302/'
    nxsTemplate = loadDir+'CORELLI_%i.nxs.h5'
    sampleRuns = range(58411,58592+1)
    peaksFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal_indexedonly.integrate'
    UBFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal.mat'
    DetCalFile = None
    qLow = -25.0; qHigh = 25.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.02 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'beryl_3D_full_newsigi' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'beryl_lab'
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth=10
    '''

    '''
    #PsbO
    loadDir = '/SNS/MANDI/IPTS-16286/data/'
    nxsTemplate = loadDir+'MANDI_%i_event.nxs'
    sampleRuns = range(6154,6165+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.integrate'#_highres.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.mat'#_highres.mat'
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'psbo_3D_full_lab_newpredppl'#_highres_newsigi' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'psbo_lab_newpredppl_highres'
    #predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 15
    '''

    '''
    #DNA
    loadDir = '/data/dna/IPTS-18552/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = range(8758,8769+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_dna2/combined_1p5A.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_dna2/combined_1p5A.mat'
    DetCalFile = '/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
    qLow = -5.0; qHigh = 5.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.007 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'dna_3D_highres' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = None#'dna_tof_2' 
    #predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 25
    fracHKLQMask = 0.5
    pplmin_frac = 0.7; pplmax_frac = 1.5
    '''
    #secondDNA
    loadDir = '/SNS/MANDI/IPTS-15151/data/'
    nxsTemplate = loadDir+'MANDI_%i_event.nxs'
    sampleRuns = range(5280,5291+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_secondDNA/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_secondDNA/combined.mat'
    DetCalFile = None#'/SNS/users/ntv/integrate/mandi_dna2/mandi_dna.DetCal'
    qLow = -5.0; qHigh = 5.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.005 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'secondDNA_3D' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = None#'dna_tof_2' 
    #predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 25
    fracHKLQMask = 0.5
    pplmin_frac = 0.6; pplmax_frac = 1.5


    '''
    #pth
    loadDir = 'SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFile = '/SNS/users/ntv/integrate/mandi_pth/peaks_combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_pth/UB_combined.mat'
    nxsTemplate = '/SNS/MANDI/IPTS-10943/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    sampleRuns = [870,872,873,874,875,876]
    DetCalFile = '/home/ntv/Desktop/runReduction/MaNDi2015.DetCal'
    qLow = -4.0; qHigh = 4.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'pth_3d_detcal' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'pth_tof_secondTry'
    predpplCoefficients = np.array([ 6.12383767,  8.8677518 , -0.02761688]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 15
    fracHKLQMask = 0.25
    pplmin_frac = 0.7; pplmax_frac = 1.5
    '''

    '''
    #gfp
    loadDir = 'SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFile = '/SNS/users/ntv/integrate/gfp/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/gfp/combined.mat'
    nxsTemplate = '/SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'

    sampleRuns = range(599,607+1)
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'gfp_3d_goodhkl' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = None#'beta_lac_lab_highres2'
    #predpplCoefficients = np.array([49.70856213,18.293623,2.58462655]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([23.2736324,10.10909695,0.6229528]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 15
    fracHKLQMask = 0.25
    pplmin_frac = 0.7; pplmax_frac = 1.5
    '''

    '''
    #Beta Lac Mutant
    loadDir = '/SNS/MANDI/IPTS-8776/data/'
    peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.mat'
    nxsTemplate = '/SNS/MANDI/IPTS-8776/data/MANDI_%i_event.nxs'
    sampleRuns = range(5921,5931+1)
    DetCalFile = None
    qLow = -10.0; qHigh = 10.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'beta_lac_3D_lab_mutant' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'beta_lac_lab_highres2'
    predpplCoefficients = np.array([ 3.56405187,  8.34071842,  0.14134522]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 15
    pplmin_frac = 0.4; pplmax_frac = 1.5
    '''

    '''
    #Beta Lac
    loadDir = '/SNS/MANDI/IPTS-15000/data/'
    nxsTemplate = loadDir+'MANDI_%i_event.nxs'
    sampleRuns = range(4999,5003+1)
    #peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
    #UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
    #peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
    #UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.mat'
    peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.mat'
    DetCalFile = None
    qLow = -10.0; qHigh = 10.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'changeme' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'beta_lac_lab_highres2'
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    mindtBinWidth = 15
    '''

    '''
    #NaK
    loadDir = '/SNS/MANDI/IPTS-17495/nexus/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = range(8275,8282+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_NAK_UB.mat'
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = 0.03 #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 40 #Width (in us) in TOF profile bins
    dQPixel = 0.003 #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'nak_3D_full_lab_2' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    numTimesToInterpolate=0
    workDir = '/SNS/users/ntv/dropbox/'
    descriptorRead = 'nak_predpws5_lab'
    predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    '''

    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
    UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()

    dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=0.5))
    dQ[dQ>0.3] = 0.3

    try: qMask = ICCFT.getHKLMask(UBMatrix, frac=fracHKLQMask, dQPixel=dQPixel,dQ=dQ)
    except: qMask = ICCFT.getHKLMask(UBMatrix, frac=fracHKL, dQPixel=dQPixel,dQ=dQ)

    padeCoefficients = ICCFT.getModeratorCoefficients('franz_coefficients_2017.dat')
    try:
        ICCFitParams = ICAT.getFitParameters(workDir, descriptorRead, sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
    except:
        print 'Cannot read ICCFitParams - will set to None and fit as we go!'
        ICCFitParams = None
    try:
        ICCFitDict = ICAT.getFitDicts(workDir, descriptorRead,sampleRuns[0], sampleRuns[-1], sampleRuns=sampleRuns)
    except:
        print 'Cannot read ICCFitDict.  Will set to None...'
        ICCFitDict = None
    strongPeakParams = pickle.load(open('strongPeakParams_betalac_lab.pkl', 'rb'))
    #strongPeakParams = pickle.load(open('strongPeakParams_dna.pkl', 'rb'))

    from timeit import default_timer as timer

    badFits = []
    oldNewList = []
   
    if sampleRunsList != -1:
        sampleRunsToAnalyze = np.array(sampleRuns).astype(np.int)[sampleRunsList]
    else: sampleRunsToAnalyze = sampleRuns
 
    for sampleRun in sampleRunsToAnalyze:

        if '{' not in nxsTemplate:
            fileName = nxsTemplate%sampleRun
        else:
            fileName = nxsTemplate.format(0, sampleRun)
        
        MDdata = ICCFT.getSample(sampleRun, DetCalFile, workDir, fileName,qLow=qLow, qHigh=qHigh, q_frame=q_frame)
        t1 = timer()
        numerrors=0
        numgood=0
        paramList = []
        for peakNumber in range(peaks_ws.getNumberPeaks()):
            progressFile = workDir+descriptor+'/progress_%i_%s.txt'%(sampleRun, descriptor)
            if progressFile is not None and peakNumber%100==0:
                with open(progressFile, 'w') as f:
                    f.write('%i\n'%(peakNumber))

            TPS = timer()
            peak = peaks_ws.getPeak(peakNumber)
            print peakNumber, peak.getIntensity()
            try:
                if peak.getRunNumber() == sampleRun:
                    print 'Integrating peak %i'%peakNumber
                    box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, fracHKL = fracHKL, refineCenter = refineCenter, dQPixel=dQPixel, q_frame=q_frame)
                    #Will force weak peaks to be fit using a neighboring peak profile
                    if ICCFitParams is not None:
                        iccfp = ICCFitParams[peakNumber]
                    else: iccfp = None
                    Y3D, goodIDX, pp_lambda, params = BVGFT.get3DPeak(peak, box, padeCoefficients,qMask,nTheta=50, nPhi=50, plotResults=False, dtBinWidth=dtBinWidth,zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1,numTimesToInterpolate=numTimesToInterpolate, fICCParams=iccfp, strongPeakParams=strongPeakParams, predCoefficients=predpplCoefficients, q_frame=q_frame, mindtBinWidth=mindtBinWidth, pplmin_frac=pplmin_frac, pplmax_frac=pplmax_frac)
                    #Does not force weak peaks
                    #Y3D, goodIDX, pp_lambda, params = BVGFT.get3DPeak(peak, box, padeCoefficients,qMask,nTheta=50, nPhi=50, plotResults=False, dtBinWidth=dtBinWidth,zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1,numTimesToInterpolate=numTimesToInterpolate, fICCParams=ICCFitParams[peakNumber], oldICCFit=ICCFitDict[peakNumber],  predCoefficients=predpplCoefficients, q_frame=q_frame, mindtBinWidth=mindtBinWidth)


                    # First we get the peak intensity
                    peakIDX = Y3D/Y3D.max() > 0.05
                    intensity = np.sum(Y3D[peakIDX])/2**(3*numTimesToInterpolate)
                    skipIDX = 2**numTimesToInterpolate
                   

                    # Now the number of background counts under the peak assuming a constant bg across the box
                    n_events = box.getNumEventsArray()
                    convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
                    conv_n_events = convolve(n_events,convBox)
                    bgIDX = reduce(np.logical_and,[~goodIDX, qMask, conv_n_events>0])
                    bgEvents = np.mean(n_events[bgIDX])*np.sum(peakIDX)

                    # Now we consider the variation of the fit.  These are done as three independent fits.  So we need to consider
                    # the variance within our fit sig^2 = sum(N*(yFit-yData)) / sum(N) and scale by the number of parameters that go into
                    # the fit.  In total: 10 (removing scale variables)
                    # TODO: It's not clear to me if we should be normalizing by #params - so we'll leave it for now.
                    w_events = n_events.copy()
                    w_events[w_events==0] = 1
                    varFit = np.average((n_events[peakIDX]-Y3D[peakIDX])*(n_events[peakIDX]-Y3D[peakIDX]), weights=(w_events[peakIDX]))

                    # Comapre with the old way
                    bgOld = np.sum(goodIDX[Y3D[::skipIDX,::skipIDX,::skipIDX]/Y3D[::skipIDX,::skipIDX,::skipIDX].max()>0.05]*pp_lambda)
                    sigmaOld = np.sqrt(intensity + bgOld)

                    # Now we add them all together.  Variances add linearly, so we just take the square root at the end.
                    sigma = np.sqrt(intensity + bgEvents + varFit)
                    print sigma, sigmaOld


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

        wsList = mtd.getObjectNames()
        for ws in wsList:
            if 'peaks' not in ws:
                print 'Removing %s'%ws
                mtd.remove(ws)

if __name__ == '__main__':
    import argparse
    import os 
    import sys
    import numpy as np
    import pickle

    parser = argparse.ArgumentParser(description='Does TOF fitting using ICConvoluted')
    parser.add_argument('-r', '--runs', action='store', default=-1, type=int, nargs='+', help='List of integers to run ranging [0,N] where N is number of runs-1. Does all runs if no argument given.')
    parser.add_argument('-p', '--pythonPathForMantid', action='store', default=None, type=str, help='path to mantid python folder which will override the default \'/opt/Mantid/bin\'')
    args = parser.parse_args()


    if args.pythonPathForMantid is None:
        sys.path.append("/opt/mantidnightly/bin")
        from mantid.simpleapi import *
    else:
        #remove the original mantid path
        popList = []
        for i in range(len(sys.path))[::-1]:
            if 'antid' in sys.path[i]:
                sys.path.pop(i)
        sys.path.append(args.pythonPathForMantid)

    from scipy.ndimage import convolve
    import ICCFitTools as ICCFT
    reload(ICCFT)
    import BVGFitTools as BVGFT
    reload(BVGFT)
    from mantid.simpleapi import *
    import ICCAnalysisTools as ICAT
    doBVGFits(sampleRunsList=args.runs)
    

