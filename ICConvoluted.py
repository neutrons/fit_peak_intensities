# ICConvoluted.py
#
# Defines two IPeakFunctions.  The first is IkedaCarpenterConvoluted
# which is the standard Ikeda-Carpenter (IC) function convoluted with 
# a square wave and a single exponential.  This implementation is 
# the same as is implemented in Rick's 'est_ic_nn_T.m'.  The second
# is a gaussian which is taken from the Mantid homepage as an example
# custom IPeakFunction.
#
import math
import numpy as np
from mantid.api._api import IFunction1D

class IkedaCarpenterConvoluted(IFunction1D):
    def init(self):
        self.declareParameter("A") #c[0]
        self.declareParameter("B") #c[1]
        self.declareParameter("R") #c[2]
        self.declareParameter("T0") #c[3]
        self.declareParameter("scale") #c[4]
        self.declareParameter("hatWidth") #c[5]
        self.declareParameter("k_conv") #c[6]


    def setPenalizedConstraints(self, A0=None, B0=None, R0=None, T00=None, scale0=None, hatWidth0=None, k_conv0=None, penalty=1.0e10):
        if A0 is not None:
            self.addConstraints("{:4.4e} < A < {:4.4e}".format(A0[0], A0[1]))
            self.setConstraintPenaltyFactor("A", penalty)
        if B0 is not None:
            self.addConstraints("{:4.4e} < B < {:4.4e}".format(B0[0], B0[1]))
            self.setConstraintPenaltyFactor("B", penalty)
        if R0 is not None:
            self.addConstraints("{:4.4e} < R < {:4.4e}".format(R0[0], R0[1]))
            self.setConstraintPenaltyFactor("R", penalty)
        if T00 is not None:
            self.addConstraints("{:4.4e} < T0 < {:4.4e}".format(T00[0], T00[1]))
            self.setConstraintPenaltyFactor("T0", penalty)
        if scale0 is not None:
            self.addConstraints("{:4.4e} < scale < {:4.4e}".format(scale0[0], scale0[1]))
            self.setConstraintPenaltyFactor("scale", penalty)
        if hatWidth0 is not None:
            self.addConstraints("{:4.4e} < hatWidth < {:4.4e}".format(hatWidth0[0], hatWidth0[1]))
            self.setConstraintPenaltyFactor("hatWidth", penalty)
        if k_conv0 is not None:
            self.addConstraints("{:4.4e} < k_conv < {:4.4e}".format(k_conv0[0], k_conv0[1]))
            self.setConstraintPenaltyFactor("k_conv", penalty)
    
    #Evaluate the function
    def function1D(self, t):
        #self.setConstraintPenaltyFactor("A", 1.0e10)
        #self.setConstraintPenaltyFactor("B", 1.0e10)
        #self.setConstraintPenaltyFactor("R", 1.0e10)
        #self.setConstraintPenaltyFactor("T0", 1.0e10)
        A  = self.getParamValue(0) 
        B  = self.getParamValue(1) 
        R  = self.getParamValue(2) 
        T0 = self.getParamValue(3)
        scale = self.getParamValue(4) 
        hatWidth  = self.getParamValue(5) 
        k_conv  = self.getParamValue(6)
        #n.b. A/2 scale factor has been removed to make A more independent 
        f_int = scale*( (1-R)*np.power((A*(t-T0)),2)*
            np.exp(-A*(t-T0))+2*R*A**2*B/np.power((A-B),3) *
            (np.exp(-B*(t-T0))-np.exp(-A*(t-T0))*(1+(A-B)*(t-T0)+0.5*np.power((A-B),2)*np.power((t-T0),2)) ) )
        f_int[t<T0] = 0
                
        mid_point_hat = len(f_int)//2
        gc_x = np.array(range(len(f_int))).astype(float)
        ppd = 0.0*gc_x
        lowIDX  = int(np.floor(np.max([mid_point_hat-np.abs(hatWidth),0])))
        highIDX = int(np.ceil(np.min([mid_point_hat+np.abs(hatWidth),len(gc_x)])))
        #print lowIDX, highIDX
        #print c, lowIDX, highIDX
        
        ppd[lowIDX:highIDX] = 1.0;
        ppd = ppd/sum(ppd);

        gc_x = np.array(range(len(f_int))).astype(float)
        gc_x = 2*(gc_x-np.min(gc_x))/(np.max(gc_x)-np.min(gc_x))-1;
        gc_f = np.exp(-k_conv*np.power(gc_x,2));
        gc_f = gc_f/np.sum(gc_f);
        #print gc_f
        
        npad = len(f_int) - 1
        first = npad - npad//2
        f_int = np.convolve(f_int,ppd,'full')[first:first+len(f_int)];
        f_int = np.convolve(f_int,gc_f,'full')[first:first+len(f_int)];
    

        return f_int
    
    #Evaluate the function for a differnt set of paremeters (trialc)
    def function1DDiffParams(self, xvals, trialc):
        #First, grab the original parameters and set to trialc
        c = np.zeros(self.numParams())
        for i in range(self.numParams()):
            c[i] = self.getParamValue(i)
            self.setParameter(i, trialc[i])
        
        #Get the trial values
        f_trial = self.function1D(xvals)
        
        #Now return to the orignial
        for i in range(self.numParams()):
            self.setParameter(i, c[i])
        return f_trial

    #Construction the Jacobian (df) for the function    
    def functionDeriv1D(self, xvals, jacobian, eps=1.e-3):
        f_int = self.function1D(xvals)
        #Fetch parameters into array c
        c = np.zeros(self.numParams())
        for i in range(self.numParams()):
            c[i] = self.getParamValue(i)
        c_shape = c.shape
        nc = np.prod(np.shape(c))
        f_shape = f_int.shape
        nf = np.prod(f_shape)
        for k in range(nc):
            dc = np.zeros(nc)
            dc[k] = max(eps,eps*c[k])
            f_new = self.function1DDiffParams(xvals,c+dc)
            for i,dF in enumerate(f_new-f_int):
                jacobian.set(i,k,dF/dc[k])





