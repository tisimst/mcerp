"""
================================================================================
mcerp: Real-time latin-hypercube-sampling-based Monte Carlo Error Propagation
================================================================================

Author: Abraham Lee
Copyright: 2013
"""

import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
from lhd import lhd

__version_info__ = (0, 9, 3)
__version__ = '.'.join(map(str, __version_info__))

__author__ = 'Abraham Lee'

__all__ = [
    # core functions
    'uv', 'covariance_matrix', 'correlation_matrix',
    # continuous distribution constructors
    'N',
    'U',
    'Exp',
    'Gamma',
    'Beta',
    'LogN',
    'X2',
    'F',
    'Tri',
    'T',
    'Weib',
    # discrete distribution constructors
    'Bern',
    'B',
    'G',
    'H',
    'Pois'
    ]

npts = 10000

CONSTANT_TYPES = (float, int, long)

class NotUpcast(Exception):
    'Raised when an object cannot be converted to a number with uncertainty'

def to_uncertain_func(x):
    """
    Transforms x into an UncertainFunction-compatible object,
    unless it is already an UncertainFunction (in which case x is returned 
    unchanged).

    Raises an exception unless 'x' belongs to some specific classes of
    objects that are known not to depend on UncertainFunction objects
    (which then cannot be considered as constants).
    """
    if isinstance(x, UncertainFunction):
        return x

    #! In Python 2.6+, numbers.Number could be used instead, here:
    elif isinstance(x, CONSTANT_TYPES):
        # No variable => no derivative to define:
        return UncertainFunction(x, [x]*npts)
    
    raise NotUpcast("%s cannot be converted to a number with"
                    " uncertainty" % type(x))
    
    

class UncertainFunction(object):
    """
    UncertainFunction objects represent the uncertainty of a result of 
    calculations with uncertain variables. Nearly all basic mathematical
    operations are supported.
    
    This class is mostly intended for internal use.
    
    
    """
    def __init__(self, x, mcpts):
        self._mcpts = np.atleast_1d(mcpts).flatten()
        self.tag = None

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
#        return self._to_general_representation(repr)
        return str(self)

    def describe(self, name=None):
        """
        Cleanly show what the four displayed distribution moments are:
            - Mean
            - Variance
            - Standardized Skewness Coefficient
            - Standardized Kurtosis Coefficient
        
        For a standard Normal distribution, these are [0, 1, 0, 3].
        
        If the object has an associated tag, this is presented. If the optional
        ``name`` kwarg is utilized, this is presented as with the moments.
        Otherwise, no unique name is presented.
        
        Example
        =======
        ::
        
            >>> x = N(0, 1, 'x')
            >>> x.describe()  # print tag since assigned
            MCERP Uncertain Value (x):
            ...

            >>> x.describe('foobar')  # 'name' kwarg takes precedence
            MCERP Uncertain Value (foobar):
            ...
            
            >>> y = x**2
            >>> y.describe('y')  # print name since assigned
            MCERP Uncertain Value (y):
            ...

            >>> y.describe()  # print nothing since no tag
            MCERP Uncertain Value:
            ...

         """
        mn, vr, sk, kt = self.stats
        if name is not None:
            s = 'MCERP Uncertain Value ('+name+'):\n'
        elif self.tag is not None:
            s = 'MCERP Uncertain Value ('+self.tag+'):\n'
        else:
            s = 'MCERP Uncertain Value:\n'
        s += ' > Mean................... {: }\n'.format(mn)
        s += ' > Variance............... {: }\n'.format(vr)
        s += ' > Skewness Coefficient... {: }\n'.format(sk)
        s += ' > Kurtosis Coefficient... {: }\n'.format(kt)
        print s
        
    def plot(self, hist=False, **kwargs):
        """
        Plot the distribution of the UncertainFunction. By default, the
        distribution is shown with a kernel density estimate (kde).
        
        Optional
        --------
        hist : bool
            If true, a density histogram is displayed (histtype='stepfilled')
        kwargs : any valid matplotlib.pyplot.plot or .hist kwarg
        
        """
        vals = self._mcpts
        low = min(vals)
        high = max(vals)

        p = ss.kde.gaussian_kde(vals)
        xp = np.linspace(low,high,100)

#        plt.figure()
        if hist:
            h = plt.hist(vals, bins=np.round(np.sqrt(len(vals))), 
                     histtype='stepfilled', normed=True, **kwargs)
            if self.tag is not None:
                plt.suptitle('Histogram of ('+self.tag+')')
                plt.title(str(self), fontsize=12)
            else:
                plt.suptitle('Histogram of')
                plt.title(str(self), fontsize=12)
            plt.ylim(0, 1.1*h[0].max())
        else:
            plt.plot(xp,p.evaluate(xp))
            if self.tag is not None:
                plt.suptitle('KDE of ('+self.tag+')')
                plt.title(str(self), fontsize=12)
            else:
                plt.suptitle('KDE of')
                plt.title(str(self), fontsize=12)

        plt.xlim(low - (high - low)*0.1, high + (high - low)*0.1)

    def show(self):
        plt.show()

    def __add__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts + uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __radd__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts + uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __mul__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts * uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __rmul__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts * uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __sub__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts - uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __rsub__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts - uf[0]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __div__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts/uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __rdiv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts/uf[0]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __truediv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts/uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __rtruediv__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts/uf[0]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __pow__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[0]._mcpts**uf[1]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)

    def __rpow__(self, val):
        uf = map(to_uncertain_func, [self, val])
        mcpts = uf[1]._mcpts**uf[0]._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
    
    def __neg__(self):
        mcpts = -self._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
        
    def __pos__(self):
        mcpts = self._mcpts
        return UncertainFunction(np.mean(mcpts), mcpts)
    
    def __abs__(self):
        mcpts = np.abs(self._mcpts)
        return UncertainFunction(np.mean(mcpts), mcpts)
    
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

################################################################################

class UncertainVariable(UncertainFunction):
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
    
    +---------------------------+-------------+-------------------+-----+---------+
    | Distribution              | scipy.stats |  args             | loc | scale   |
    |                           | class name  | (shape params)    |     |         |
    +===========================+=============+===================+=====+=========+
    | Normal(mu, sigma)         | norm        |                   | mu  | sigma   | 
    +---------------------------+-------------+-------------------+-----+---------+
    | Uniform(a, b)             | uniform     |                   | a   | b-a     |
    +---------------------------+-------------+-------------------+-----+---------+
    | Exponential(lamda)        | expon       |                   |     | 1/lamda |
    +---------------------------+-------------+-------------------+-----+---------+
    | Gamma(k, theta)           | gamma       | k                 |     | theta   |
    +---------------------------+-------------+-------------------+-----+---------+
    | Beta(alpha, beta, [a, b]) | beta        | alpha, beta       | a   | b-a     |
    +---------------------------+-------------+-------------------+-----+---------+
    | Log-Normal(mu, sigma)     | lognorm     | sigma             | mu  |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Chi-Square(k)             | chi2        | k                 |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | F(d1, d2)                 | f           | d1, d2            |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Triangular(a, b, c)       | triang      | c                 | a   | b-a     |
    +---------------------------+-------------+-------------------+-----+---------+
    | Student-T(v)              | t           | v                 |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Weibull(lamda, k)         | exponweib   | lamda, k          |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Bernoulli(p)              | bernoulli   | p                 |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Binomial(n, p)            | binomial    | n, p              |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Geometric(p)              | geom        | p                 |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Hypergeometric(M, n, N)   | hypergeom   | M, n, N           |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    | Poisson(lamda)            | poisson     | lamda             |     |         |
    +---------------------------+-------------+-------------------+-----+---------+
    
    Thus, each distribution above would have the same call signature::
        
        >>> import scipy.stats as ss
        >>> ss.your_dist_here(args, loc=loc, scale=scale)
        
    ANY SCIPY.STATS.DISTRIBUTION SHOULD WORK! IF ONE DOESN'T, PLEASE LET ME
    KNOW!
    
    Convenient constructors have been created to make assigning these 
    distributions easier. They follow the parameter notation found in the
    respective Wikipedia articles:
    
    +---------------------------+---------------------------------------------------------------+
    | MCERP Distibution         | Wikipedia page                                                |
    +===========================+===============================================================+
    | N(mu, sigma)              | http://en.wikipedia.org/wiki/Normal_distribution              |
    +---------------------------+---------------------------------------------------------------+
    | U(a, b)                   | http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)|
    +---------------------------+---------------------------------------------------------------+
    | Exp(lamda, [mu])          | http://en.wikipedia.org/wiki/Exponential_distribution         |
    +---------------------------+---------------------------------------------------------------+
    | Gamma(k, theta)           | http://en.wikipedia.org/wiki/Gamma_distribution               |
    +---------------------------+---------------------------------------------------------------+
    | Beta(alpha, beta, [a, b]) | http://en.wikipedia.org/wiki/Beta_distribution                |
    +---------------------------+---------------------------------------------------------------+
    | LogN(mu, sigma)           | http://en.wikipedia.org/wiki/Log-normal_distribution          |
    +---------------------------+---------------------------------------------------------------+
    | X2(df)                    | http://en.wikipedia.org/wiki/Chi-squared_distribution         |
    +---------------------------+---------------------------------------------------------------+
    | F(dfn, dfd)               | http://en.wikipedia.org/wiki/F-distribution                   |
    +---------------------------+---------------------------------------------------------------+
    | Tri(a, b, c)              | http://en.wikipedia.org/wiki/Triangular_distribution          |
    +---------------------------+---------------------------------------------------------------+
    | T(df)                     | http://en.wikipedia.org/wiki/Student's_t-distribution         |
    +---------------------------+---------------------------------------------------------------+
    | Weib(lamda, k)            | http://en.wikipedia.org/wiki/Weibull_distribution             |
    +---------------------------+---------------------------------------------------------------+
    | Bern(p)                   | http://en.wikipedia.org/wiki/Bernoulli_distribution           |
    +---------------------------+---------------------------------------------------------------+
    | B(n, p)                   | http://en.wikipedia.org/wiki/Binomial_distribution            |
    +---------------------------+---------------------------------------------------------------+
    | G(p)                      | http://en.wikipedia.org/wiki/Geometric_distribution           |
    +---------------------------+---------------------------------------------------------------+
    | H(M, n, N)                | http://en.wikipedia.org/wiki/Hypergeometric_distribution      |
    +---------------------------+---------------------------------------------------------------+
    | Pois(lamda)               | http://en.wikipedia.org/wiki/Poisson_distribution             |
    +---------------------------+---------------------------------------------------------------+

    Thus, the following are equivalent::

        >>> x = N(10, 1)
        >>> x = uv(ss.norm(loc=10, scale=1))

    Examples
    --------
    A three-part assembly
        
        >>> x1 = N(24, 1)
        >>> x2 = N(37, 4)
        >>> x3 = Exp(2)  # Exp(mu=0.5) works too
        
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
    calculation to the ne        if hist:
            plt.hist(vals, bins=np.round(np.sqrt(vals)), histtype='stepfilled',
                **kwargs)
        else:
xt, this parameter should be changed before any 
    calculations are performed to ensure parameter compatibility (this may 
    change to be more dynamic in the future, but for now this is how it is).
    
    Also, to see the underlying distribution of the variable, and if matplotlib
    is installed, simply call its plot method::
        
        >>> x1.plot()
    
    Optional kwargs can be any valid kwarg used by matplotlib.pyplot.plot
    
    See Also
    --------
    N, U, Exp, Gamma, Beta, LogN, X2, F, Tri, T, Weib, Bern, B, G, H, Pois
        
    """
    
    def __init__(self, rv, tag=None):
        
        assert hasattr(rv, 'dist'), 'Input must be a  distribution from ' + \
            'the scipy.stats module.'
        self.rv = rv
        
        # generate the latin-hypercube points
        self._mcpts = lhd(dist=self.rv, size=npts).flatten()
        self.tag = tag
        
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
    
    def plot(self, hist=False, **kwargs):
        """
        Plot the distribution of the UncertainVariable. Continuous 
        distributions are plotted with a line plot and discrete distributions
        are plotted with discrete circles.
        
        Optional
        --------
        hist : bool
            If true, a histogram is displayed
        kwargs : any valid matplotlib.pyplot.plot kwarg
        
        """
        
        if hist:
            vals = self._mcpts
            low = vals.min()
            high = vals.max()
            h = plt.hist(vals, bins=np.round(np.sqrt(len(vals))), 
                     histtype='stepfilled', normed=True, **kwargs)

            if self.tag is not None:
                plt.suptitle('Histogram of (' + self.tag + ')')
                plt.title(str(self), fontsize=12)
            else:
                plt.suptitle('Histogram of')
                plt.title(str(self), fontsize=12)

            plt.ylim(0, 1.1*h[0].max())
        else:
            bound = 1e-6
            low = self.rv.ppf(bound)
            high = self.rv.ppf(1 - bound)
            if hasattr(self.rv.dist, 'pmf'):
                low = int(low)
                high = int(high)
                vals = range(low, high + 1)
                plt.plot(vals, self.rv.pmf(vals), 'o', **kwargs)

                if self.tag is not None:
                    plt.suptitle('PMF of (' + self.tag + ')')
                    plt.title(str(self), fontsize=12)
                else:
                    plt.suptitle('PMF of')
                    plt.title(str(self), fontsize=12)

            else:
                vals = np.linspace(low, high, 500)
                plt.plot(vals, self.rv.pdf(vals), **kwargs)

                if self.tag is not None:
                    plt.suptitle('PDF of ('+self.tag+')')
                    plt.title(str(self), fontsize=12)
                else:
                    plt.suptitle('PDF of')
                    plt.title(str(self), fontsize=12)

        plt.xlim(low - (high - low)*0.1, high + (high - low)*0.1)

                
        
uv = UncertainVariable # a nicer form for the user

###############################################################################
# Define some convenience constructors for common statistical distributions.
# Hopefully these are a little easier/more intuitive to use than the 
# scipy.stats.distributions.
###############################################################################

def N(mu, sigma, tag=None):
    """
    A Normal (or Gaussian) random variate
    
    Parameters
    ----------
    mu : scalar
        The mean value of the distribution
    sigma : scalar
        The standard deviation (must be positive and non-zero)
    """
    assert sigma>0, 'Sigma must be positive'
    return uv(ss.norm(loc=mu, scale=sigma), tag=tag)

###############################################################################

def U(a, b, tag=None):
    """
    A Uniform random variate
    
    Parameters
    ----------
    low : scalar
        Lower bound of the distribution support.
    high : scalar
        Upper bound of the distribution support.
    """
    assert a<b, 'Lower bound must be less than the upper bound'
    return uv(ss.uniform(loc=a, scale=b-a), tag=tag)

###############################################################################

def Exp(lamda, mu=None, tag=None):
    """
    An Exponential random variate
    
    Parameters
    ----------
    lamda : scalar
        The inverse scale (as shown on Wikipedia). 
    
    Optional
    --------
    mu : scalar
        The mean value of the distribution (must be positive and non-zero). If 
        this is given, ``lamda`` is ignored. (FYI: mu = 1/lamda.)
    """
    if mu is not None:
        assert mu>0, 'Mean must be positive and not zero'
        return uv(ss.expon(scale=mu), tag=tag)
    else:
        return uv(ss.expon(scale=1./lamda), tag=tag)

###############################################################################

def Gamma(k, theta, tag=None):
    """
    A Gamma random variate
    
    Parameters
    ----------
    k : scalar
        The shape parameter (must be positive and non-zero)
    theta : scalar
        The scale parameter (must be positive and non-zero)
    """
    assert k>0 and theta>0, 'Gamma parameters must be greater than zero'
    return uv(ss.gamma(k, scale=theta), tag=tag)

###############################################################################

def Beta(alpha, beta, a=0, b=1, tag=None):
    """
    A Beta random variate
    
    Parameters
    ----------
    alpha : scalar
        The first shape parameter
    beta : scalar
        The second shape parameter
    
    Optional
    --------
    a : scalar
        Lower bound of the distribution support (default=0)
    b : scalar
        Upper bound of the distribution support (default=1)
    """
    assert alpha>0 and beta>0, 'Shape parameters must be greater than zero'
    return uv(ss.beta(alpha, beta, loc=a, scale=b-a), tag=tag)

###############################################################################

def LogN(mu, sigma, tag=None):
    """
    A Log-Normal random variate
    
    Parameters
    ----------
    mu : scalar
        The location parameter
    sigma : scalar
        The scale parameter (must be positive and non-zero)
    """
    assert sigma>0, 'Sigma must be positive'
    return uv(ss.lognorm(sigma, loc=mu), tag=tag)

###############################################################################

def X2(k, tag=None):
    """
    A Chi-Squared random variate
    
    Parameters
    ----------
    k : int
        The degrees of freedom of the distribution (must be greater than one)
    """
    assert isinstance(df, int) and df>1, 'DF must be an int greater than 1'
    return uv(ss.chi2(df), tag=tag)

###############################################################################

def F(d1, d2, tag=None):
    """
    An F (fisher) random variate
    
    Parameters
    ----------
    d1 : int
        Numerator degrees of freedom
    d2 : int
        Denominator degrees of freedom
    """
    assert isinstance(d1, int) and d1>1, 'DFN must be an int greater than 1'
    assert isinstance(d2, int) and d2>1, 'DFD must be an int greater than 1'
    return uv(ss.f(d1, d2), tag=tag)

###############################################################################

def Tri(a, b, c, tag=None):
    """
    A triangular random variate
    
    Parameters
    ----------
    a : scalar
        Lower bound of the distribution support (default=0)
    b : scalar
        Upper bound of the distribution support (default=1)
    c : scalar
        The location of the triangle's peak (a <= c <= b)
    """
    assert a<=c<=b, 'peak must lie in between low and high'
    return uv(ss.triang(c, loc=a, scale=b-a), tag=tag)

###############################################################################

def T(v, tag=None):
    """
    A Student-T random variate
    
    Parameters
    ----------
    v : int
        The degrees of freedom of the distribution (must be greater than one)
    """
    assert isinstance(v, int) and v>1, 'DF must be an int greater than 1'
    return uv(ss.t(v), tag=tag)

###############################################################################

def Weib(lamda, k, tag=None):
    """
    A Weibull random variate
    
    Parameters
    ----------
    lamda : scalar
        The scale parameter
    k : scalar
        The shape parameter
    """
    assert lamda>0 and k>0, 'Weibull scale and shape parameters must be greater than zero'
    return uv(ss.exponweib(lamda, k), tag=tag)

###############################################################################

def Bern(p, tag=None):
    """
    A Bernoulli random variate
    
    Parameters
    ----------
    p : scalar
        The probability of success
    """
    assert 0<p<1, 'Bernoulli probability must be between zero and one'
    return uv(ss.bernoulli(p), tag=tag)

###############################################################################

def B(n, p, tag=None):
    """
    A Binomial random variate
    
    Parameters
    ----------
    n : int
        The number of trials
    p : scalar
        The probability of success
    """
    assert int(n)==n and n>0, 'Binomial number of trials must be an integer greater than zero'
    assert 0<p<1, 'Binomial probability must be between zero and one'
    return uv(ss.binom(n, p), tag=tag)

###############################################################################

def G(p, tag=None):
    """
    A Geometric random variate
    
    Parameters
    ----------
    p : scalar
        The probability of success
    """
    assert 0<p<1, 'Geometric probability must be between zero and one'
    return uv(ss.geom(p), tag=tag)

###############################################################################

def H(p, tag=None):
    """
    A Hypergeometric random variate
    
    Parameters
    ----------
    M : int
        The total population size
    n : int
        The number of individuals of interest in the population
    N : int
        The number of individuals that will be chosen from the population
    """
    assert int(M)==M and M>0, 'Hypergeometric total population size must be an integer greater than zero.'
    assert int(n)==n and 0<n<=M, 'Hypergeometric interest population size must be an integer greater than zero and no more than the total population size.'
    assert int(N)==N and 0<N<=M, 'Hypergeometric chosen population size must be an integer greater than zero and no more than the total population size.'
    return uv(ss.hypergeom(M, n, N), tag=tag)

###############################################################################

def Pois(lamda, tag=None):
    """
    A Poisson random variate
    
    Parameters
    ----------
    lamda : scalar
        The rate of an occurance within a specified interval of time or space.
    """
    assert lamda>0, 'Poisson rate must be greater than zero.'
    return uv(ss.poisson(lamda), tag=tag)

###############################################################################

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
    
        >>> x = N(1, 0.1)
        >>> y = N(10, 0.1)
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
    
        >>> x = N(1, 0.1)
        >>> y = N(10, 0.1)
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
    x1 = N(24, 1)
    x2 = N(37, 4)
    x3 = Exp(2)  # Exp(mu=0.5) is the same
    Z = (x1*x2**2)/(15*(1.5 + x3))
    Z.describe()
    
    print '*'*80
    print 'Example of volumetric gas flow through orifice meter'
    H = N(64, 0.5)
    M = N(16, 0.1)
    P = N(361, 2)
    t = N(165, 0.5)
    C = 38.4
    Q = C*umath.sqrt((520*H*P)/(M*(t + 460)))
    Q.describe()

    print '*'*80
    print 'Example of manufacturing tolerance stackup'
    # for a gamma distribution we need the following conversions:
    # shape = mean**2/var
    # scale = var/mean
    mn = 1.5
    vr = 0.25
    k = mn**2/vr
    theta = vr/mn
    x = Gamma(k, theta)
    y = Gamma(k, theta)
    z = Gamma(k, theta)
    w = x + y + z
    w.describe()

    print '*'*80
    print 'Example of scheduling facilities (six stations)'
    s1 = N(10, 1)
    s2 = N(20, 2**0.5)
    mn1 = 1.5
    vr1 = 0.25
    k1 = mn1**2/vr1
    theta1 = vr1/mn1
    s3 = Gamma(k1, theta1)
    mn2 = 10
    vr2 = 10
    k2 = mn2**2/vr2
    theta2 = vr2/mn2
    s4 = Gamma(k2, theta2)
    s5 = Exp(5)  # Exp(mu=0.2) is the same
    s6 = X2(10)
    T = s1 + s2 + s3 + s4 + s5 + s6
    T.describe()

    print '*'*80
    print 'Example of two-bar truss stress/deflection analysis'
    H = N(30, 5/3., tag='H')
    B = N(60, 0.5/3., tag='B')
    d = N(3, 0.1/3, tag='d')
    t = N(0.15, 0.01/3, tag='t')
    E = N(30000, 1500/3., tag='E')
    rho = N(0.3, 0.01/3., tag='rho')
    P = N(66, 3/3., tag='P')
    pi = np.pi
    wght = 2*pi*rho*d*t*umath.sqrt((B/2)**2 + H**2)
    strs = (P*umath.sqrt((B/2)**2 + H**2))/(2*pi*d*t*H)
    buck = (pi**2*E*(d**2 + t**2))/(8*((B/2)**2 + H**2))
    defl = (P*((B/2)**2 + H**2)**(1.5))/(2*pi*d*t*H**2*E)
    print wght.describe('wght')
    print strs.describe('strs')
    print buck.describe('buck')
    print defl.describe('defl')

    print '** TESTS COMPLETE **'
