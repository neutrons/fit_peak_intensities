import math
import numpy as np
from mantid.api._api import IFunction1D
from matplotlib.mlab import bivariate_normal


class MBVG(IFunction1D): #MantidBVG
    '''
    MBVG implements a bivariate gaussian (BVG) in Mantid (M) as a 1D function.  This is done so that it can be
    fit in a straightforward fashion using Mantid's Fit() function.  To achieve this, we use the flattened
    version of the 2D profile and fit it as a 1D function.  It is built on matplotlib.mlab.bivariate_normal, which
    is available on SNS analysis machines.

    To make it compatible with fitting, X, Y, and E must all be the same shape.  This is possible if we input
    twice and reconstruct the BVG in the function.

    If h is an n*n 2d profile we're trying to fitt and th, ph are the n*1 arrays containing the x and y coordinates,
    we would fit it as follows:


        TH, PH = np.meshgrid(th, ph,indexing='ij') #Get 2D version
        m = mbvg.MBVG()
        m.init()
        m['A'] = 1.0
        # ... #Set initial parameters
        m['nX'] = len(th) #Length needed for reconstruction
        m['nY'] = len(ph) #Length needed for reconstruction

        pos = np.empty(TH.shape + (2,))
        pos[:,:,0] = TH
        pos[:,:,1] = PH
        h2 = m.function2D(pos)

        H = np.empty(h.shape + (2,))
        H[:,:,0] = h
        H[:,:,1] = h

        bvgWS = CreateWorkspace(OutputWorkspace='bvgWS',DataX=pos.ravel(),DataY=H.ravel(),dataE=np.sqrt(H.ravel()))
        fitResults = Fit(Function=m, InputWorkspace='bvgWS', Output='bvgfit')

        BS - May 18 2018
    '''
    def init(self):
        self.declareParameter("A") #Amplitude
        self.declareParameter("muX") #Mean along the x direction
        self.declareParameter("muY") #Mean along the y direction
        self.declareParameter("sigX") #sigma along the x direction
        self.declareParameter("sigY") #sigma along the y direction
        self.declareParameter("sigP") #interaction term rho
        self.declareParameter("nX") #used for reconstructing 2d profile
        self.declareParameter("nY") #used for reconstruction 2d profile
        self.addConstraints("0 < A") #Require amplitude to be positive
        self.addConstraints("0 < sigX") #standard deviations must be positive
        self.addConstraints("0 < sigY") #standard deviations must be positive
 
    def function1D(self, t):
        '''
        function1D returns the flattened version of the function.
        Input, t, may be in one of two forms:
            1) a 1D array (e.g. pos.ravel() from the example).  
            2) a 3D array (e.g. pos from the example)
        Output
            If input is of type 1, a 1D array matching the size of the 
            input is returned.  This allows fitting to take place.
            If input is of type 2, a 1D array matching the size of 
            pos.ravel() is returned.
        '''
        if t.ndim == 1:
            nX = int(self.getParamValue(6))
            nY = int(self.getParamValue(7))
            pos = t.reshape(nX, nY, 2)
        elif t.ndim == 3:
            pos = t
        X = pos[...,0]
        Y = pos[...,1]

        A = self.getParamValue(0)
        muX = self.getParamValue(1)
        muY = self.getParamValue(2)
        sigX = self.getParamValue(3)
        sigY = self.getParamValue(4)
        sigP = self.getParamValue(5)


        sigXY = sigX*sigY*sigP
        Z = A*bivariate_normal(X,Y, sigmax=sigX, sigmay=sigY,
                            mux=muX,muy=muY,sigmaxy=sigXY)
        if t.ndim == 1:
            zRet = np.empty(Z.shape+(2,))
            zRet[:,:,0] = Z
            zRet[:,:,1] = Z
        elif t.ndim == 3:
            zRet = Z
        return zRet.ravel()

    def function2D(self, t):
        '''
        function2D returns the 2D version of the BVG.
        Input may be in two forms:
            1) 1D array (e.g. pos.ravel()).  This will be reshaped into an nX*nY*2 array, so
                it must contain nX*nY*2 elements.
            2) 3D array of size A*B*2. A and B are arbitrary integers.
        Output:
            a 2D array either size nX*nY (intput type 1) or A*B (input type 2) with intensities
            of the BVG.
        '''
        print t.ndim
        if t.ndim == 1:
            nX = int(self.getParamValue(6))
            nY = int(self.getParamValue(7))
            pos = t.reshape(nX, nY, 2)
        elif t.ndim == 3:
            pos = t
        X = pos[...,0]
        Y = pos[...,1]
        A = self.getParamValue(0)
        muX = self.getParamValue(1)
        muY = self.getParamValue(2)
        sigX = self.getParamValue(3)
        sigY = self.getParamValue(4)
        sigP = self.getParamValue(5)

        sigXY = sigX*sigY*sigP
        Z = A*bivariate_normal(X,Y, sigmax=sigX, sigmay=sigY,
                            mux=muX,muy=muY,sigmaxy=sigXY)

        return Z        
 
