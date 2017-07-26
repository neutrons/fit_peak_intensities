import math
import numpy as np
from mantid.api._api import IPeakFunction

class IkedaCarpenterConvoluted(IPeakFunction):
        def init(self):
                self.declareParameter("Alpha") #c[0]//C1
                self.declareParameter("Beta") #c[1]//C2
                self.declareParameter("R") #c[2]//C3
                self.declareParameter("Ta") #c[3]//C4
                self.declareParameter("Tb") #c[4]//C5
                self.declareParameter("Tc") #c[5]//C6
                self.declareParameter("Td") #c[6]//C7
                self.declareParameter("gc_alpha") #c[7]//C8
                self.declareParameter("pp_width") #c[8]//C9
                self.declareParameter("C10") #c[9]//C10
                self.declareParameter("C11") #c[10]//C11
                self.declareParameter("PeakCentre") #IPeak peak property
                self.declareParameter("Sigma") #IPeak FWHM property
                self.declareParameter("Height") #IPeak Height property

        def functionLocal(self, t_prim_n):
		print np.shape(t_prim_n)
		singleFlag = False
		if len(t_prim_n) == 1:
			t_prim_n = np.array(range(20))+1
			singleFlag = True
                c = np.zeros(self.numParams())
                for i in range(self.numParams()):
                        c[i] = self.getParamValue(i)
                n_time = np.prod(np.shape(t_prim_n))
                t_prim = np.array(range(n_time))+1
                #print np.min(t_prim), np.max(t_prim)
                #print ['%4.4f'%i for i in c] 
                t_prim = c[5]*(t_prim-c[3])/(c[4]-c[3])
                #print np.min(t_prim), np.max(t_prim)
                f_int = (1-c[2])*(c[0]/2)*\
                        (np.power(c[0]*t_prim,2))*np.exp(-c[0]*t_prim)+\
                        (c[2]*c[1]*np.power(c[0],3)/np.power(c[0]-c[1],3))*\
                        (np.exp(-c[1]*t_prim)-np.exp(-c[0]*t_prim)*(1+(c[0]-c[1])*\
                        t_prim+0.5*np.power((c[0]-c[1])*t_prim,2)))
                f_int[t_prim<0] = 0
                mid_point_hat = len(f_int)//2
                gc_x = np.array(range(len(f_int))).astype(float)
                ppd = 0.0*gc_x
                lowIDX  = np.max([mid_point_hat-np.abs(c[8]),0])
                highIDX = np.min([mid_point_hat+np.abs(c[8]),len(gc_x)])
                ppd[lowIDX:highIDX] = 1.0;
                #print lowIDX, highIDX, ppd
		ppd = ppd/sum(ppd);

                gc_x = np.array(range(len(f_int))).astype(float)
                gc_x = 2*(gc_x-np.min(gc_x))/(np.max(gc_x)-np.min(gc_x))-1;
                gc_f = np.exp(-c[7]*np.power(gc_x,2));
                gc_f = gc_f/np.sum(gc_f);

                #For compatability with Rick's code, we'll have to do a full
                # convolution and select the elements we want.
                npad = len(f_int) - 1
                first = npad - npad//2
                f_int = np.convolve(f_int,ppd,'full')[first:first+len(f_int)];
                f_int = np.convolve(f_int,gc_f,'full')[first:first+len(f_int)];
                f_int = f_int+c[6];
                if singleFlag:
			#print 'returning %4.4f'%np.max(f_int)
			return np.array(np.max(f_int))
		print 'returning length %i'%len(f_int)
		return f_int

	def functionDerivLocal(self, xvals, jacobian):
		print 'functionDerivLocal not implemented for ICC'


        def centre(self):
          return self.getParameterValue("PeakCentre")

        def height(self):
          return self.getParameterValue("Height")

        def fwhm(self): #This I have to calculte for IC
          return 2.0*math.sqrt(2.0*math.log(2.0))*self.getParameterValue("Sigma")

        def setCentre(self, new_centre):
	  print 'setCentre'
          # User picked point new_centre
          self.setParameter("PeakCentre",new_centre)

        def setHeight(self, new_height):
	  print 'setHeight'
          # User set new height for peak
          self.setParameter("Height", new_height)

        def setFwhm(self, new_fwhm):
	  print 'setFWHM'
          # User had a guess at the width using the range picker bar
          sigma = new_fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
          self.setParameter("Sigma",sigma)




class PyGaussian(IPeakFunction):
    
  def category(self):
    return "Examples"

  def init(self):
    self.declareParameter("Height")
    self.declareParameter("PeakCentre")
    self.declareParameter("Sigma")
       
  def functionLocal(self, xvals):
    height = self.getParameterValue("Height")
    peak_centre = self.getParameterValue("PeakCentre")
    sigma = self.getParameterValue("Sigma")
    weight = math.pow(1./sigma,2);

    offset_sq=np.square(xvals-peak_centre)
    out=height*np.exp(-0.5*offset_sq*weight)
    print '------'
    print out
    print type(out)
    return out
    
  def functionDerivLocal(self, xvals, jacobian):
    height = self.getParameterValue("Height");
    peak_centre = self.getParameterValue("PeakCentre");
    sigma = self.getParameterValue("Sigma")
    weight = math.pow(1./sigma,2);
        
    # X index
    i = 0
    for x in xvals:
      diff = x-peak_centre
      exp_term = math.exp(-0.5*diff*diff*weight)
      jacobian.set(i,0, exp_term)
      jacobian.set(i,1, diff*height*exp_term*weight)
      # derivative with respect to weight not sigma
      jacobian.set(i,2, -0.5*diff*diff*height*exp_term)
      i += 1

  def activeParameter(self, index):
    param_value = self.getParameterValue(index)
    if index == 2: #Sigma. Actually fit to 1/(sigma^2) for stability
      return 1./math.pow(param_value,2)
    else:
      return param_value

  def setActiveParameter(self, index, value):
    param_value = value
    explicit = False
    if index == 2:
      param_value = math.sqrt(math.fabs(1.0/value))
    else:
      param_value = value
      # Final explicit argument is required to be false here
      self.setParameter(index, param_value, False) 

  def centre(self):
    return self.getParameterValue("PeakCentre")

  def height(self):
    return self.getParameterValue("Height")

  def fwhm(self):
    return 2.0*math.sqrt(2.0*math.log(2.0))*self.getParameterValue("Sigma")

  def setCentre(self, new_centre):
    self.setParameter("PeakCentre",new_centre)

  def setHeight(self, new_height):
    self.setParameter("Height", new_height)

  def setFwhm(self, new_fwhm):
    sigma = new_fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
    self.setParameter("Sigma",sigma)

