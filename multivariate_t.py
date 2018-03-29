import numpy as np
from scipy.special import gamma

class multivariate_t():
    '''
    Multivariate student's t distribution.  Modified from http://stackoverflow.com/a/29804411/3521179
    This is designed to be interchangeable with scipy.stats.multivariate_normal.  Currently only pdf is 
    implemented.

    BS 28 March 2018 
    '''

    # constructor inputs:
    #  mu = 1xd numpy array with the means (d is the dimensionality of the distribution)
    #  sigma = dxd numpy array for the sigma matrix (related to covariance)
    #  df = number, degrees of freedom
    def __init__(self, mu, sigma, df):
        if len(mu) != sigma.ndim:
            raise ValueError("Cannot create multivariate t distribution with mu=%s and sigma = %s"%(mu,sigma))
        self.mu = mu  #Center
        self.sigma = sigma #~covariance
        self.d = len(mu) #number of dimensions
        self.df = df #DOF


    '''
    Returns the pdf of the distribution evaluated at points x.  If we are evaluating on a d dimensional grid, then
    x should be d+1 dimensions with d entries in the last dimensions.  For instance on a 100x100 grid (d=2), d would be 
    100x100x2.
    '''
    def pdf(self, x):
        numerator = gamma(1. * (self.d+self.df)/2)
        if self.sigma.ndim > 1:
            invSig = np.linalg.inv(1.*self.sigma)
            detSig = np.linalg.det(1.*self.sigma)
            WeirdBit = np.einsum('...k,kl,...l->...', x-self.mu, invSig , x-self.mu)
        else:
            invSig = 1./self.sigma
            detSig = self.sigma
            WeirdBit = ((x-self.mu)/self.sigma)**2
        denom =  (gamma(1.*self.df/2) * pow(self.df*np.pi,1.*self.d/2) * 
                    pow(detSig,1./2) * np.power(1+1.0/self.df*WeirdBit ,1.* (self.d+self.df)/2))
        d = 1. * numerator / denom 
        return d

# multivariate_lorentzian inherents multivariate_t but sets the degrees of 
# freedom (df) to 1 to get a Lorenztian (Cauchy) distribution. This is done
# for transparency with scipy.stats.multivariate_normal 
class multivariate_lorentzian(multivariate_t):
    def __init__(self, mu, sigma):
        multivariate_t.__init__(self, mu, sigma, 1)



