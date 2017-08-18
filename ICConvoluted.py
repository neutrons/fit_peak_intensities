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
	
	#Evaluate the function
	def function1D(self, t):
		A  = self.getParamValue(0) 
		B  = self.getParamValue(1) 
		R  = self.getParamValue(2) 
		T0 = self.getParamValue(3)
		scale = self.getParamValue(4) 
		hatWidth  = self.getParamValue(5) 
		k_conv  = self.getParamValue(6) 
		f_int = scale*A/2*( (1-R)*np.power((A*(t-T0)),2)*
			np.exp(-A*(t-T0))+2*R*A**2*B/np.power((A-B),3) *
			(np.exp(-B*(t-T0))-np.exp(-A*(t-T0))*(1+(A-B)*(t-T0)+0.5*np.power((A-B),2)*np.power((t-T0),2)) ) )
		f_int[t<T0] = 0
                
		mid_point_hat = len(f_int)//2
                gc_x = np.array(range(len(f_int))).astype(float)
                ppd = 0.0*gc_x
                lowIDX  = np.max([mid_point_hat-np.abs(hatWidth),0])
                highIDX = np.min([mid_point_hat+np.abs(hatWidth),len(gc_x)])
                #print lowIDX, highIDX
		#print c, lowIDX, highIDX
		ppd[lowIDX:highIDX] = 1.0;
		ppd = ppd/sum(ppd);

                gc_x = np.array(range(len(f_int))).astype(float)
                gc_x = 2*(gc_x-np.min(gc_x))/(np.max(gc_x)-np.min(gc_x))-1;
                gc_f = np.exp(-k_conv*np.power(gc_x,2));
                gc_f = gc_f/np.sum(gc_f);
		#print gc_f
                #For compatability with Rick's code, we'll have to do a full
                # convolution and select the elements we want.
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





