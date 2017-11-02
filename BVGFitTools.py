
def fitPeak3D(X, n_events, peak,goodIDX):
    paramNames = [0 for i in range(7)]
    padeCoefficients = getModeratorCoefficients('franz_coefficients_2017.dat')
    flightPath = peak.getL1() + peak.getL2()
    halfScatteringAngle = 0.5*peak.getScattering()
    paramNames = [0 for i in range(7)]
    energy = peak.getInitialEnergy()/1000.0
    tofWS = CreateWorkspace(DataX = np.ones(50), DataY=np.ones(50))
    x0 = getInitialGuess(tofWS, paramNames, energy, flightPath, padeCoefficients,None, None)     
    alpha = x0[0]
    beta = x0[1]
    R = x0[2]
    T0 = x0[3]
    k_conv = 120
    p0 = [np.max(n_events), np.mean(X[:,:,:,1]), np.mean(X[:,:,:,2]), 0.005, 0.005, 0.05, alpha, beta, R, T0, k_conv, 0.0 ]
    bounds = ([0,-np.inf,-np.inf,0,0,-np.inf] + [0.5*v for v in x0[:4]] + [10, 0],
                [np.inf for i in range(6)] + [1.5*v for v in x0[:4]] + [500, np.inf]  )
    params= curve_fit(peak3DFitFunction, X[goodIDX], n_events[goodIDX],  p0,maxfev=1000, sigma=np.sqrt(n_events[goodIDX]),bounds=bounds)
    return params

def peak3DFitFunction(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg):
    return peak3D(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg).ravel()
    
def peak3DFromParams(X,params):
    return peak3D(X, params[0],params[1],params[2],params[3],params[4],params[5],
                    params[6],params[7],params[8],params[9],params[10],params[11])

def peak3D(X, A, mu1, mu2, sigmaX, sigmaY, p12, alpha, beta, R, T0, k_conv, bg):
    #Axis 0 = TOF, axis1 = theta, axis2 = phi
    #First we calculate the TOF distribution to estimate sigma_TOF

    if X.ndim == 4:
        XTOF = X[:,:,:,0]
        XTHETA = X[:,:,:,1]
        XPHI = X[:,:,:,2]
        XANGLE = X[:,:,:,1:]
    elif X.ndim == 2:
        XTOF = X[:,0]
        XTHETA = X[:,1]
        XPHI = X[:,2]
        XANGLE = X[:,1:]
        
    else: return 

    fICC = ICC.IkedaCarpenterConvoluted()
    fICC.init()
    fICC['A'] = alpha
    fICC['B'] = beta


