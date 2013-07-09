# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 15:48:17 2013

"""
from ad import ADF,ADV
import numpy as np
import math
from lhd import lhd
import scipy.stats as ss
import matplotlib.pyplot as plt

__version_info__ = (0, 8, 2)
__version__ = '.'.join(map(str, __version_info__))

__author__ = 'Abraham Lee'

__all__ = ['uv', 'covariance_matrix', 'correlation_matrix']

npts = 10000

CONSTANT_TYPES = (float, int, complex, long)

class NotUpcast(Exception):
    'Raised when an object cannot be converted to a number with uncertainty'

def to_uncertain_func(x):
    """
    Transforms x into a constant automatically differentiated function (ADF),
    unless it is already an ADF (in which case x is returned unchanged).

    Raises an exception unless 'x' belongs to some specific classes of
    objects that are known not to depend on ADF objects
    (which then cannot be considered as constants).
    """
    if isinstance(x, UncertainFunction):
        return x

    #! In Python 2.6+, numbers.Number could be used instead, here:
    elif isinstance(x, CONSTANT_TYPES):
        # No variable => no derivative to define:
        return UncertainFunction(x, {}, {}, {}, x)
    
    raise NotUpcast("%s cannot be converted to a number with"
                    " uncertainty" % type(x))
    
    

class UncertainFunction(ADF):
    """
    UncertainFunction objects represent the uncertainty of a result of 
    calculations with uncertain variables. Nearly all basic mathematical
    operations are supported.
    
    This class is mostly intended for internal use.
    
    
    """
    def __init__(self, x, d, d2, d2c, mcpts):
        self._mcpts = np.atleast_1d(mcpts).flatten()
        ADF.__init__(self, x, d, d2, d2c)

    @property
    def mean(self):
        """
        Mean value as a result of an uncertainty calculation
        """
        mn = np.mean(self._mcpts)
        return mn
    
    @property
    def var(self):
        """
        Variance value as a result of an uncertainty calculation
        """
        mn = self.mean
        vr = np.mean((self._mcpts - mn)**2)
        return vr
        
    @property
    def std(self):
        """
        Standard deviation value as a result of an uncertainty calculation, 
        defined as::
            
                    ________
            std = \/variance
            
        """
        return self.var**0.5
    
    @property
    def skew(self):
        """
        Skewness coefficient value as a result of an uncertainty calculation,
        defined as::
            
              _____     m3
            \/beta1 = ------
                      std**3
        
        where m3 is the third central moment and std is the standard deviation
        """
        mn = self.mean
        sd = self.std
        sk = 0.0 if abs(sd)<=1e-8 else np.mean((self._mcpts - mn)**3)/sd**3
        return sk
    
    @property
    def kurt(self):
        """
        Kurtosis coefficient value as a result of an uncertainty calculation,
        defined as::
            
                      m4
            beta2 = ------
                    std**4

        where m4 is the fourth central moment and std is the standard deviation
        """
        mn = self.mean
        sd = self.std
        kt = 0.0 if abs(sd)<=1e-8 else np.mean((self._mcpts - mn)**4)/sd**4
        return kt
    
    @property
    def stats(self):
        """
        The first four standard moments of a distribution: mean, variance, and
        standardized skewness and kurtosis coefficients.
        """
        mn = self.mean
        vr = self.var
        sk = self.skew
        kt = self.kurt
        return [mn, vr, sk, kt]
        
    def _to_general_representation(self,str_func):
        mn, vr, sk, kt = self.stats
        return ('uv({:}, {:}, {:}, {:})'.format(
            str_func(mn), str_func(vr), str_func(sk), str_func(kt)) 
            if any([vr, sk, kt]) else str_func(mn))

    def __str__(self):
        return self._to_general_representation(str)

    def __repr__(self):
        return self._to_general_representation(repr)

    def describe(self):
        """
        Cleanly show what the distribution moments are:
            - Mean, Variance, Skewness and Kurtosis Coefficients
        """
        mn, vr, sk, kt = self.stats
        s = 'MCERP Uncertain Value:\n'
        s += ' > Mean................... {: }\n'.format(mn)
        s += ' > Variance............... {: }\n'.format(vr)
        s += ' > Skewness Coefficient... {: }\n'.format(sk)
        s += ' > Kurtosis Coefficient... {: }\n'.format(kt)
        print s
        
    def plot(self, **kwargs):
        vals = self._mcpts
        low = min(vals)
        high = max(vals)

#        p = ss.kde.gaussian_kde(vals)
#        xp = np.linspace(low,high,100)

        plt.figure()
#        plt.plot(xp,p.evaluate(xp))
        plt.hist(vals, bins=np.round(np.sqrt(npts)), histtype='stepfilled')
        plt.title(str(self))
        plt.xlim(low - (high - low)*0.1, high + (high - low)*0.1)
        plt.ylim(ymin=0)
        plt.show()

    def __add__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts + uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__add__(self, val), mcpts)

    def __radd__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts + uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__radd__(self, val), mcpts)
        
    def __mul__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts * uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__mul__(self, val), mcpts)

    def __rmul__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts * uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__rmul__(self, val), mcpts)
        
    def __sub__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts - uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__sub__(self, val), mcpts)

    def __rsub__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts - uf[0]._mcpts
        return _make_UF_compatible_object(ADF.__rsub__(self, val), mcpts)
        
    def __div__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts/uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__div__(self, val), mcpts)

    def __rdiv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts/uf[0]._mcpts
        return _make_UF_compatible_object(ADF.__rdiv__(self, val), mcpts)
        
    def __truediv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts/uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__truediv__(self, val), mcpts)

    def __rtruediv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts/uf[0]._mcpts
        return _make_UF_compatible_object(ADF.__rtruediv__(self, val), mcpts)
        
    def __pow__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts**uf[1]._mcpts
        return _make_UF_compatible_object(ADF.__pow__(self, val), mcpts)

    def __rpow__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts**uf[0]._mcpts
        return _make_UF_compatible_object(ADF.__rpow__(self, val), mcpts)
    
    def __neg__(self):
        mcpts = -self._mcpts
        return _make_UF_compatible_object(ADF.__neg__(self), mcpts)
        
    def __pos__(self):
        mcpts = self._mcpts
        return _make_UF_compatible_object(ADF.__pos__(self), mcpts)
    
    def __abs__(self):
        mcpts = np.abs(self._mcpts)
        return _make_UF_compatible_object(ADF.__abs__(self), mcpts)
    
    def __eq__(self,val):
        diff = self - val
        return not (diff.mean or diff.var or diff.skew or diff.kurt)
    
    def __ne__(self,val):
        return not self==val
    
    def __lt__(self,val):
        self,val = map(to_uncertain_func, [self,val])
        return True if float(self.mean - val.mean) < 0 else False
    
    def __le__(self,val):
        return (self==val) or self < val
    
    def __gt__(self,val):
        return not self < val
    
    def __ge__(self,val):
        return (self==val) or self > val

    def __nonzero__(self):
        return self!=0

def _make_UF_compatible_object(tmp, mcpts):
    if isinstance(tmp, ADF):
        return UncertainFunction(tmp.x, tmp.d(), tmp.d2(), tmp.d2c(), mcpts)
    else: # for scalars, etc.
        return UncertainFunction(tmp, {}, {}, {}, [tmp]*npts)

################################################################################

class UncertainVariable(UncertainFunction,ADV):
    """
    UncertainVariable objects track the effects of uncertainty, characterized 
    in terms of the first four standard moments of statistical distributions 
    (mean, variance, skewness and kurtosis coefficients). Monte Carlo simulation,
    in conjunction with Latin-hypercube based sampling performs the calculations.

    Parameters
    ----------
    rv : scipy.stats.distribution
        A distribution to characterize the uncertainty
    
    tag : str, optional
        A string identifier when information about this variable is printed to
        the screen
        
    Notes
    -----
    
    The ``scipy.stats`` module contains many distributions which we can use to
    perform any necessary uncertainty calculation. It is important to follow
    the initialization syntax for creating any kind of distribution object:
        
        - *Location* and *Scale* values must use the kwargs ``loc`` and 
          ``scale``
        - *Shape* values are passed in as non-keyword arguments before the 
          location and scale, (see below for syntax examples)..
        
    The mathematical operations that can be performed on uncertain objects will 
    work for any distribution supplied, but may be misleading if the supplied 
    moments or distribution is not accurately defined. Here are some guidelines 
    for creating UncertainVariable objects using some of the most common 
    statistical distributions:
    
    =======================  ============  ===================  ====  =========
    Distribution             scipy.stats   args                 loc   scale
                             class name    (shape parameters)
    =======================  ============  ===================  ====  =========
    Normal(mu,sigma)         norm                               mu    sigma
    Uniform(a,b)             uniform                            a     b-a
    Exponential(lambda)      expon                                    1/lambda
    Gamma(k,theta)           gamma         k                          theta
    Beta(alpha,beta,[a,b])   beta          alpha,beta           a     b-a
    Log-Normal(mu,sigma)     lognorm       sigma                mu
    Chi-Square(dv)           chi2          dv
    F(df_numer,df_denom)     f             df_numer,df_denom
    Triangular(a,b,peak)     triang        peak                 a     b-a
    Student-T(df)            t             df
    =======================  ============  ===================  ====  =========
    
    Thus, each distribution above would have the same call signature::
        
        >>> import scipy.stats as ss
        >>> ss.your_dist_here(args, loc=loc, scale=scale)
        
    Examples
    --------
    A three-part assembly
        
        >>> import scipy.stats as ss
        >>> x1 = uv(ss.norm(loc=24, scale=1))  # normally distributed
        >>> x2 = uv(ss.norm(loc=37, scale=4))  # normally distributed
        >>> x3 = uv(ss.expon(scale=0.5))  # exponentially distributed
        
        >>> Z = (x1*x2**2)/(15*(1.5 + x3))
        >>> Z
        uv(1161.46231679, 116646.762981, 0.345533974771, 3.00791101068)

    The result shows the mean, variance, and standardized skewness and kurtosis
    of the output variable Z, which will vary from use to use due to the random
    nature of Monte Carlo simulation and latin-hypercube sampling techniques.
    
    Basic math operations may be applied to distributions, where all 
    statistical calculations are performed using latin-hypercube enhanced Monte
    Carlo simulation. Nearly all of the built-in trigonometric-, logarithm-, 
    etc. functions of the ``math`` module have uncertainty-compatible 
    counterparts that should be used when possible since they support both 
    scalar values and uncertain objects. These can be used after importing the 
    ``umath`` module::
        
        >>> from mcerp.umath import * # sin(), sqrt(), etc.
        >>> sqrt(x1)
        uv(4.89791765647, 0.0104291897681, -0.0614940614672, 3.00264937735)
    
    At any time, the standardized statistics can be retrieved using::
        
        >>> x1.mean
        >>> x1.var  # x1.std (standard deviation) is also available
        >>> x1.skew
        >>> x1.kurt
    
    or all four together with::
    
        >>> x1.stats
    
    By default, the Monte Carlo simulation uses 10000 samples, but this can be
    changed at any time with::
        
        >>> mcerp.npts = number_of_samples
    
    Any value from 1,000 to 1,000,000 is recommended (more samples means more
    accurate, but also means more time required to perform the calculations). 
    Although it can be changed, since variables retain their samples from one
    calculation to the next, this parameter should be changed before any 
    calculations are performed to ensure parameter compatibility (this may 
    change to be more dynamic in the future, but for now this is how it is).
    
    Also, to see the underlying distribution of the variable, and if matplotlib
    is installed, simply call its plot method::
        
        >>> x1.plot()
    
    Optional kwargs can be any valid kwarg used by matplotlib.pyplot.plot
        
    """
    
    def __init__(self, rv, tag=None):
        assert hasattr(rv, 'dist'), 'Keyword "rv" must be a distribution ' + \
            'from the scipy.stats module.'

        self.rv = rv
        
        # generate the latin-hypercube points
        self._mcpts = lhd(dist=self.rv, size=npts).flatten()
        
        ADV.__init__(self, self.rv.mean(), tag=tag)
    
    @property
    def mean(self):
        return float(self.rv.stats('m'))
    
    @property
    def var(self):
        return float(self.rv.stats('v'))
		
    @property
    def std(self):
        return self.var**0.5
        
    @property
    def skew(self):
        return float(self.rv.stats('s'))
    
    @property
    def kurt(self):
        return float(self.rv.stats('k')) + 3  # remove the 3 for standardization
    
    @property
    def stats(self):
        return [self.mean, self.var, self.skew, self.kurt]
    
    def plot(self, vals=None, **kwargs):
        """
        Plot the distribution of the UncertainVariable.
        
        Optional
        --------
        vals : array-like
            If no values are put in, default x-values between the 0.01% and 
            99.99% points will be used.
        
        kwargs : any valid matplotlib.pyplot.plot kwarg
        
        """
        if vals is None:
            low = self.rv.ppf(0.0001)
            high = self.rv.ppf(0.9999)
            vals = np.linspace(low, high, 500)
        else:
            low = min(vals)
            high = max(vals)

        plt.plot(vals, self.rv.pdf(vals), **kwargs)
        plt.title(str(self))
        plt.xlim(low - (high - low)*0.1, high + (high - low)*0.1)
        plt.show()
                
        
uv = UncertainVariable # a nicer form for the user

def covariance_matrix(nums_with_uncert):
    """
    Calculate the covariance matrix of uncertain variables, oriented by the
    order of the inputs
    
    Parameters
    ----------
    nums_with_uncert : array-like
        A list of variables that have an associated uncertainty
    
    Returns
    -------
    cov_matrix : 2d-array-like
        A nested list containing covariance values
    
    Example
    -------
    
        >>> x = uv(ss.norm(loc=1, scale=0.1))
        >>> y = uv(ss.norm(loc=10, scale=0.1))
        >>> z = x + 2*y
        >>> covariance_matrix([x,y,z])
        [[  9.99694861e-03   2.54000840e-05   1.00477488e-02]
         [  2.54000840e-05   9.99823207e-03   2.00218642e-02]
         [  1.00477488e-02   2.00218642e-02   5.00914772e-02]]

    """
    ufuncs = map(to_uncertain_func,nums_with_uncert)
    cov_matrix = []
    for (i1, expr1) in enumerate(ufuncs):
        coefs_expr1 = []
        mean1 = expr1.mean
        for (i2, expr2) in enumerate(ufuncs[:i1+1]):
            mean2 = expr2.mean
            coef = np.mean((expr1._mcpts - mean1)*(expr2._mcpts - mean2))
            coefs_expr1.append(coef)
        cov_matrix.append(coefs_expr1)
        
    # We symmetrize the matrix:
    for (i, covariance_coefs) in enumerate(cov_matrix):
        covariance_coefs.extend(cov_matrix[j][i]
                                for j in range(i+1, len(cov_matrix)))

    return cov_matrix

def correlation_matrix(nums_with_uncert):
    """
    Calculate the correlation matrix of uncertain variables, oriented by the
    order of the inputs
    
    Parameters
    ----------
    nums_with_uncert : array-like
        A list of variables that have an associated uncertainty
    
    Returns
    -------
    corr_matrix : 2d-array-like
        A nested list containing covariance values
    
    Example
    -------
    
        >>> x = uv(ss.norm(loc=1, scale=0.1))
        >>> y = uv(ss.norm(loc=10, scale=0.1))
        >>> z = x + 2*y
        >>> correlation_matrix([x,y,z])
        [[ 0.99969486  0.00254001  0.4489385 ]
         [ 0.00254001  0.99982321  0.89458702]
         [ 0.4489385   0.89458702  1.        ]]

    """
    ufuncs = map(to_uncertain_func,nums_with_uncert)
    cov_matrix = covariance_matrix(ufuncs)
    corr_matrix = []
    for (i1, expr1) in enumerate(ufuncs):
        row_data = []
        for (i2, expr2) in enumerate(ufuncs):
            row_data.append(cov_matrix[i1][i2]/expr1.std/expr2.std)
        corr_matrix.append(row_data)
    return corr_matrix
    
if __name__=='__main__':
    
    import mcerp.umath as umath
    
    print '*'*80
    print 'TEST FUNCTIONS USING DERIVED MOMENTS FROM SCIPY DISTRIBUTIONS'
    print '*'*80
    print 'Example of a three part assembly'
    x1 = uv(ss.norm(loc=24, scale=1))   # normally distributed
    x2 = uv(ss.norm(loc=37, scale=4))   # normally distributed
    x3 = uv(ss.expon(scale=0.5))       # exponentially distributed
    
    Z = (x1*x2**2)/(15*(1.5 + x3))
    print Z
    
    print '*'*80
    print 'Example of volumetric gas flow through orifice meter'
    H = uv(ss.norm(loc=64, scale=0.5))
    M = uv(ss.norm(loc=16, scale=0.1))
    P = uv(ss.norm(loc=361, scale=2))
    t = uv(ss.norm(loc=165, scale=0.5))
    C = 38.4
    Q = C*umath.sqrt((520*H*P)/(M*(t + 460)))
    print Q

    print '*'*80
    print 'Example of manufacturing tolerance stackup'
    # for a gamma distribution we need the following conversions:
    # scale = var/mean
    # shape = mean**2/var
    mn = 1.5
    vr = 0.25
    scale = vr/mn
    shape = mn**2/vr
    x = uv(ss.gamma(shape, scale=scale)) 
    y = uv(ss.gamma(shape, scale=scale)) 
    z = uv(ss.gamma(shape, scale=scale))
    w = x + y + z
    print w

    print '*'*80
    print 'Example of scheduling facilities (six stations)'
    s1 = uv(ss.norm(loc=10, scale=1))
    s2 = uv(ss.norm(loc=20, scale=2**0.5))
    mn1 = 1.5
    vr1 = 0.25
    scale1 = vr1/mn1
    shape1 = mn1**2/vr1
    s3 = uv(ss.gamma(shape1, scale=scale1))
    mn2 = 10
    vr2 = 10
    scale2 = vr2/mn2
    shape2 = mn2**2/vr2
    s4 = uv(ss.gamma(shape2, scale=scale2))
    s5 = uv(ss.expon(scale=0.2))
    s6 = uv(ss.chi2(10))
    T = s1 + s2 + s3 + s4 + s5 + s6
    print T

    print '*'*80
    print 'Example of two-bar truss'
    H = uv(ss.norm(loc=30, scale=5/3.), tag='H')
    B = uv(ss.norm(loc=60, scale=0.5/3.), tag='B')
    d = uv(ss.norm(loc=3, scale=0.1/3), tag='d')
    t = uv(ss.norm(loc=0.15, scale=0.01/3), tag='t')
    E = uv(ss.norm(loc=30000, scale=1500/3.), tag='E')
    rho = uv(ss.norm(loc=0.3, scale=0.01/3.), tag='rho')
    P = uv(ss.norm(loc=66, scale=3/3.), tag='P')
    pi = math.pi
    wght = 2*pi*rho*d*t*umath.sqrt((B/2)**2 + H**2)
    strs = (P*umath.sqrt((B/2)**2 + H**2))/(2*pi*d*t*H)
    buck = (pi**2*E*(d**2 + t**2))/(8*((B/2)**2 + H**2))
    defl = (P*((B/2)**2 + H**2)**(1.5))/(2*pi*d*t*H**2*E)
    print 'wght:',wght
    print 'strs:',strs
    print 'buck:',buck
    print 'defl:',defl

    print '** TESTS COMPLETE **'
