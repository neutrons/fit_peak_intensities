
def doIntegration(sampleRunsList=None):
    # Some parameters
    workDir = '/SNS/users/ntv/dropbox/' #End with '/'
    bgPolyOrder = 1 #int- bg order for IC fit
    doVolumeNormalization = False #True if you want to normalize TOF profiles by volume
    refineCenter = False #True if you want to determine new centers - still not very good
    removeEdges = False #True if you want to not consider q-pixels that are off detector faces
    calcTOFPerPixel = False #True to calculate TOF for each pixel in a qBox - uses interpolation for volNorm (if doVolumeNormalization is True)
    fracHKL = 0.5 #Fraction of HKL to look on either side
    fracStop = 0.01 #Fraction of max counts to include in peak selection
    moderatorCoefficientsFile = 'franz_coefficients_2017.dat'
    #moderatorCoefficientsFile = 'franz_coefficients_2010.dat'
    calibrationDictFile = 'det_calibration/calibration_dictionary.pkl'
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
    qLow = -15.0; qHigh = 15.0
    dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 4 #Width (in us) in TOF profile bins
    dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    descriptor = 'scol_dynamicQuadBG_constraints' #Does not end with '/'
    dQMax = 0.4
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    '''

    '''
    #IPTS22331
    loadDir = '/SNS/TOPAZ/IPTS-16903/data/'
    nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
    sampleRuns = [22331]
    peaksFormat='/SNS/users/ntv/peaks_%i.integrate'
    peaksFile = '/SNS/users/ntv/peaks_22331.integrate'
    UBFormat = None
    UBFile = '/SNS/TOPAZ/IPTS-16903/shared/Chunruo/ep_100K_MoTe2_adpQ_find/22331_Niggli.mat'
    DetCalFile = '/SNS/TOPAZ/IPTS-16903/shared/calibration/TOPAZ_2016B.DetCal'
    qLow = -25.0; qHigh = 25.0
    dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 4 #Width (in us) in TOF profile bins
    dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.4
    doIterativeBackgroundFitting = False
    descriptor = 'run22331' #Does not end with '/'
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    '''
    '''
    #NaK - MANDI
    loadDir = '/SNS/MANDI/IPTS-17495/nexus/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = range(8275,8282+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_nak_8275_8282.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_nak/MANDI_NAK_UB.mat'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'nak_predpws5_lab' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([12.51275, 13.078622, 0.18924]) #Go with ICCFT.oldScatFun
    q_frame = 'lab'
    '''

    '''
    #PsbO - 2016 - MANDI
    loadDir = '/SNS/MANDI/IPTS-16286/data/'
    nxsTemplate = loadDir+'MANDI_%i_event.nxs'
    sampleRuns = range(6154,6165+1)
    #peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.integrate'
    #UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal.mat'
    peaksFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal_highres.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_psbo/combined_hexagonal_highres.mat'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'psbo_lab_newpredppl_highres' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    #predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([14.36827809, 10.889742, 0.28754095]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    minppl_frac=0.8; maxppl_frac=1.5;
    '''
    '''
    # DNA - 2017 - MANDI
    loadDir = '/SNS/MANDI/IPTS-18552/nexus/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = range(8758,8769+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_dna/combined_orthorhombic.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_dna/combined_orthorhombic.mat'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.007,0.007] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'dna_lab_newpredppl' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    #predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([ 10.46241806,  10.53543448,   0.23630636]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    minppl_frac=0.8; maxppl_frac=1.5
    '''
    '''
    #gfp
    sampleRuns = range(599,607+1)
    loadDir = 'SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFile = '/SNS/users/ntv/integrate/gfp/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/gfp/combined.mat'
    nxsTemplate = '/SNS/MANDI/2013_2_11B_SCI/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'gfp_tof_goodhkl' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    #predpplCoefficients = np.array([ 3.56405187,  8.34071842,  0.14134522]) #Go with ICCFT.oldScatFun
    #predpplCoefficients = np.array([49.70856213,18.293623,2.58462655]) #Go with ICCFT.oldScatFun
    predpplCoefficients = np.array([23.2736324 ,  10.10909695,   0.6229528 ]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    fracHKLQMask = 0.25
    minppl_frac=0.7; maxppl_frac=1.5; mindtBinWidth=15
    '''
    #pth
    sampleRuns = [870,872,873,874,875,876]
    loadDir = '/SNS/MANDI/IPTS-10943/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFile = '/SNS/users/ntv/integrate/mandi_pth/peaks_combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_pth/UB_combined.mat'
    nxsTemplate = '/SNS/MANDI/IPTS-10943/{0}/{1}/NeXus/MANDI_{1}_event.nxs'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -4.0; qHigh = 4.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'pth_tof_secondRun' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([ 6.12383767,  8.8677518 , -0.02761688]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    fracHKLQMask = 0.25
    minppl_frac=0.7; maxppl_frac=1.5; mindtBinWidth=15




    '''
    #Beta lactamase mutant
    sampleRuns = range(5921,5931+1)
    loadDir = '/SNS/MANDI/IPTS-8776/data/'
    peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase3/combined.mat'
    nxsTemplate = '/SNS/MANDI/IPTS-8776/data/MANDI_%i_event.nxs'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -10.0; qHigh = 10.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'beta_lac_lab_highres_mut2' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([ 3.56405187,  8.34071842,  0.14134522]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    minppl_frac=0.4; maxppl_frac=1.5; mindtBinWidth=15
    '''

    ''' 
    #Beta lactamase - 2016 - MANDI
    loadDir = '/SNS/MANDI/IPTS-15000/data/'
    nxsTemplate = loadDir+'MANDI_%i_event.nxs'
    sampleRuns = range(4999,5003+1)
    #peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.integrate'
    #UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/combined_triclinic.mat'
    #peaksFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase_2.integrate'
    #UBFile = '/SNS/users/ntv/integrate/mandi_betalactamase/MANDI_betalactamase.mat'
    peaksFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.integrate'
    UBFile = '/SNS/users/ntv/integrate/mandi_beta_lactamase2/combined.mat'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -10.0; qHigh = 10.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'changeme' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    minppl_frac=0.8; maxppl_frac=1.5; mindtBinWidth=15
    '''
    '''
    # CORELLI - beryl
    loadDir = '/data/corelli_beryl/IPTS-20302/'
    peaksFile = '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal_indexedonly.integrate'
    UBFile =  '/SNS/users/ntv/integrate/corelli_beryl/combined_hexagonal.mat'
    DetCalFile = None
    workDir = '/SNS/users/ntv/dropbox/' #End with '/'
    nxsTemplate = loadDir+'CORELLI_%i.nxs.h5'
    sampleRuns = range(58411,58592+1)
    peaksFormat = peaksFile
    UBFormat = UBFile
    qLow = -25.0; qHigh = 25.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.02,0.02] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'beryl_lab_cs1' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    q_frame='lab'
    minppl_frac=0.0; maxppl_frac=4.5; mindtBinWidth=10
    '''
    '''
    #Natrolite - 2016 - MANDI
    loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = range(8733,8750+1)
    peaksFile = '/SNS/users/ntv/integrate/mandi_natrolite/mandi_natrolite_peaks.integrate'
    UBFile =  '/SNS/users/ntv/integrate/mandi_natrolite/mandi_natrolite_niggli_ub.mat'
    peaksFormat = peaksFile
    UBFormat = UBFile
    DetCalFile = None
    qLow = -5.0; qHigh = 5.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 30 #Width (in us) in TOF profile bins
    dQPixel = [0.005,0.005] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.15 #tune this
    descriptor = 'beta_lac_predpws6' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    '''

    '''
    #Natrolite - 2017 - MANDI
    loadDir = '/SNS/MANDI/IPTS-8776/nexus/'
    nxsTemplate = loadDir+'MANDI_%i.nxs.h5'
    sampleRuns = [8679]
    peaksFile='/SNS/users/ntv/integrate/mandi_natrolite/peaks_8679.integrate'
    UBFile='/SNS/users/ntv/integrate/mandi_natrolite/mandi_8679_niggli.mat'
    #DetCalFile = '/SNS/MANDI/shared/calibration/MANDI_500.DetCal'
    DetCalFile = None
    peaksFormat = peaksFile
    UBFormat = UBFile
    qLow = -10.0; qHigh = 10.0
    dtSpread = [0.03,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 4 #Width (in us) in TOF profile bins
    dQPixel = [0.003,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.2
    doIterativeBackgroundFitting = False
    descriptor = 'natrolite' #Does not end with '/'
    doIterativeBackgroundFitting = False
    nBG=5
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    predpplCoefficients = np.array([5.24730283,  7.23719321,  0.27449887]) #Go with ICCFT.oldScatFun
    q_frame = 'sample'
    '''
    '''
    #Si - 2016A
    loadDir = '/SNS/TOPAZ/shared/PeakIntegration/data/'
    nxsTemplate = loadDir+'TOPAZ_%i_event.nxs'
    sampleRuns = range(15647,15670)
    peaksFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.integrate'
    peaksFormat = peaksFile
    UBFile = '/SNS/TOPAZ/shared/PeakIntegration/DataSet/Si2mm_2016A_15647_15669/Si2mm_Cubic_F.mat'
    UBFormat = UBFile
    DetCalFile = '/SNS/TOPAZ/shared/PeakIntegration/calibration/TOPAZ_2016A.DetCal'
    qLow = -25.0; qHigh = 25.0
    dtSpread = [0.025,0.03] #how far we look on either side of the nominal peak for each fit criteria - recommended to increase
    dtBinWidth = 4 #Width (in us) in TOF profile bins
    dQPixel = [0.005,0.003] #dQ for each voxel in qBox - recommended to decrease for successive fits
    dQMax = 0.4
    descriptor = 'si_0p015_removeBG' #Does not end with '/'
    doIterativeBackgroundFitting = False
    parameterDict = pickle.load(open('det_calibration/calibration_dictionary_scolecite.pkl','rb'))
    '''
    #==================WORK STARTS HERE==========================================
    figsFormat = None# workDir + descriptor+'/figs/mantid_%i_%i.png'
    '''
    if os.path.isdir(workDir + descriptor):
        inp = raw_input('!!!!!! WARNING: PATH %s ALREADY EXIST!!!!! CONTINUE? (Y/<n>):'%(workDir + descriptor))
        print inp
        if inp is 'y' or inp is 'Y':
            print 'Will overwrite files in %s...'%workDir + descriptor
        else:
            print 'exiting...'
            sys.exit()
    else:
        os.mkdir(workDir + descriptor)
        os.mkdir(workDir + descriptor + '/figs/')
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
    UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
    '''


    #Load our peaks files and detector fitting parameters

    if peaksFile is not None:
        peaks_ws = LoadIsawPeaks(Filename = peaksFile)
        LoadIsawUB(InputWorkspace=peaks_ws, FileName=UBFile)
        UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
        dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=fracHKL))
        dQ[dQ>dQMax] = dQMax
        qMask = list()
        for dQP in dQPixel:
            print 'Getting qMask for dQPixel=%f'%dQP
            try:
                qMask.append(ICCFT.getHKLMask(UBMatrix, frac=fracHKLQMask, dQPixel=dQP,dQ=dQ))
            except:
                qMask.append(ICCFT.getHKLMask(UBMatrix, frac=fracHKL, dQPixel=dQP,dQ=dQ))

    padeCoefficients = ICCFT.getModeratorCoefficients(moderatorCoefficientsFile)
    calibrationDict = pickle.load(open(calibrationDictFile, 'rb'))

    #Write the log
    logFile = workDir + descriptor + '/log.log'
    ICFitLog.writeLog(logFile, workDir, loadDir, nxsTemplate, figsFormat, sampleRuns, dtSpread, dtBinWidth, fracHKL, fracStop, refineCenter, removeEdges, doVolumeNormalization, peaksFormat, UBFormat, DetCalFile, moderatorCoefficientsFile, calibrationDictFile, descriptor,zBG,neigh_length_m, predpplCoefficients, minppl_frac, maxppl_frac)
    if sampleRunsList != -1:
        sampleRunsToAnalyze = np.array(sampleRuns).astype(np.int)[sampleRunsList]
    else: sampleRunsToAnalyze = sampleRuns
    for sampleRun in sampleRunsToAnalyze:
        #Set up a few things for the run
        paramList = list()
        if '{' not in nxsTemplate:
            fileName = nxsTemplate%sampleRun
        else:
            fileName = nxsTemplate.format(0, sampleRun)
        #If we want to remove edges, we rebuild the panel dictionary every run
        # TODO this can be reformulated in QLab and apply R each box.
        instrumentFile = EdgeTools.getInstrumentFile(peaks_ws, peaksFile)
        if removeEdges:
            panelDict = EdgeTools.getInstrumentDict(instrumentFile, peaks_ws, sampleRun, fitOrder=2)
        else:
            panelDict = EdgeTools.getPanelDictionary(instrumentFile)

        #Conver the sample to reciprocal space
        MDdata = ICCFT.getSample(sampleRun, DetCalFile, workDir, fileName, qLow=qLow, qHigh=qHigh, q_frame=q_frame)
        
        #Do the actual integration
        peaks_ws,paramList,fitDict = ICCFT.integrateSample(sampleRun, MDdata, peaks_ws, paramList, panelDict, UBMatrix, dQ, qMask, padeCoefficients,parameterDict, figsFormat=figsFormat,dtBinWidth = dtBinWidth, dtSpread=dtSpread, fracHKL = fracHKL, refineCenter=refineCenter, doVolumeNormalization=doVolumeNormalization, minFracPixels=0.01, fracStop=fracStop, removeEdges=removeEdges, calibrationDict=calibrationDict,dQPixel=dQPixel, calcTOFPerPixel=calcTOFPerPixel,neigh_length_m=neigh_length_m,zBG=zBG, bgPolyOrder=bgPolyOrder, nBG=nBG, doIterativeBackgroundFitting=doIterativeBackgroundFitting,predCoefficients=predpplCoefficients, q_frame=q_frame, progressFile=workDir+descriptor+'/progress_%i_%s.txt'%(sampleRun, descriptor), mindtBinWidth=mindtBinWidth,minpplfrac=minppl_frac, maxpplfrac=maxppl_frac)
        #Save the results and delete the leftovers
        os.system('rm ' + workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
        SaveIsawPeaks(InputWorkspace='peaks_ws', Filename=workDir+descriptor+'/peaks_%i_%s.integrate'%(sampleRun,descriptor))
        np.savetxt(workDir+descriptor+'/params_%i_%s.dat'%(sampleRun, descriptor), np.array(paramList))
        pickle.dump(fitDict, open(workDir+descriptor+'/fitDict_%i_%s.pkl'%(sampleRun,descriptor),'wb'))
        
        wsList = mtd.getObjectNames()
        for ws in wsList:
            if 'peaks' not in ws:
                mtd.remove(ws)
                print 'Removing workspace %s'%ws


if __name__=='__main__':
    import argparse
    import sys
    import os 
    import pickle
    import numpy as np
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

    from mantid.simpleapi import *
    import ICCFitTools as ICCFT
    reload(ICCFT)
    import ICConvoluted as ICC
    reload(ICC)
    import ICFitLog
    reload(ICFitLog)
    import getEdgePixels as EdgeTools
    reload(EdgeTools)
    FunctionFactory.subscribe(ICC.IkedaCarpenterConvoluted)


    doIntegration(sampleRunsList=args.runs)
