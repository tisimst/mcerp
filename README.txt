Overview
=========

``mcerp`` is a combination of on-the-fly `Monte Carlo methods`_, `latin-hypercube sampling`_, and second-order error propagation (see `soerp`_ for the Python implementation of the original Fortran code `SOERP` by N. D. Cox) to perform non-order specific `error propagation`_ (or uncertainty analysis). The ``soerp`` package allows you to **easily** and **transparently** track the effects of uncertainty through mathematical calculations. Advanced mathematical functions, similar to those in the standard `math`_ module can also be evaluated directly. 

Due to the nature of sampling techniques, calculation results will vary from session to session (but consistent within the session) since new latin-hypercube samples are only generated when variables are newly defined or redifined. By default, each variable uses 10,000 latin-hypercube samples that are sufficiently random. This can be changed by assigning an integer value to the ``mcerp.npts`` object (typically, values between 1,000 and 1,000,000 are sufficiently large to insure small errors in the resulting statistics). This should only be changed prior to performing calculations since all subsequent calculations assume that each input has the same number of sampled points (this could be changed through resampling, I suppose...).

In order to correctly use ``mcerp``, knowledge of the distributions from the ``scipy.stats`` sub-module is required. See the documentation of the ``uv`` constructor for how to use the distributions most effectively. The result of all calculations generates a *mean*, *variance*, and *standardized skewness and kurtosis* coefficients. Running the source code by itself in a terminal with ``python mcerp.py`` will demonstrate some more advanced examples.


Required Packages
=================

- `NumPy`_ : Numeric Python
- `SciPy`_ : Scientific Python
- `ad`_ : For automatic differentiation

Suggested Packages
==================

- `uncertainties`_ : First-order uncertainty propagation
- `Matplotlib`_ : Python plotting library (required for the ``mcerp.plot`` command)

Basic examples
==============
::

    >>> from mcerp import uv  # the constructor for uncertain variables
    >>> from mcerp import ss  # for using scipy.stats distributions
    >>> x1 = uv(ss.norm(loc=24, scale=1)) # normally distributed
    >>> x2 = uv(ss.norm(loc=37, scale=4)) # normally distributed
    >>> x3 = uv(ss.expon(scale=0.5))  # exponentially distributed

    >>> Z = (x1*x2**2)/(15*(1.5 + x3))
    >>> print Z
    uv(1161.14296871, 116093.134064, 0.361152281239, 3.05247793644)

    >>> Z.describe()  # explains what the moments are
    MCERP Uncertain Value:
     > Mean...................  1161.14296871
     > Variance...............  116093.134064
     > Skewness Coefficient...  0.361152281239
     > Kurtosis Coefficient...  3.05247793644

    >>> x1.stats  # the eight moments can be accessed at any time
    [24.0, 1.0, 0.0, 3.0]
    
    >>> x1.plot()  # if matplotlib is installed it shows the underlying distribution

    >>> x1-x1  # correlations are correctly handled
    0.0
    
    # convenient access to derivatives
    >>> Z.d(x1)  # first derivative wrt x1 (returns all if no input, 0 if derivative doesn't exist)
    45.63333333333333
    >>> Z.d2(x2)  # second derivative wrt x2
    1.6
    >>> Z.d2c(x1,x3) # second cross-derivative wrt x1 and x3 (order doesn't matter)
    -22.816666666666666
    
    >>> Z.gradient([x1, x2, x3])  # convenience function, useful for optimization
    [  45.63333333   59.2        -547.6       ]

    >>> Z.hessian([x1, x2, x3])   # another convenience function
    [[   0.            2.46666667  -22.81666667]
     [   2.46666667    1.6         -29.6       ]
     [ -22.81666667  -29.6         547.6       ]]

    # Example of volumetric gas flow through orifice meter
    >>> import mcerp.umath as umath  # sin, exp, sqrt, etc.
    >>> H = uv(ss.norm(loc=64, scale=0.5))
    >>> M = uv(ss.norm(loc=16, scale=0.1))
    >>> P = uv(ss.norm(loc=361, scale=2))
    >>> t = uv(ss.norm(loc=165, scale=0.5))
    >>> C = 38.4
    >>> Q = C*umath.sqrt((520*H*P)/(M*(t + 460)))
    >>> Q.describe()
    MCERP Uncertain Value:
     > Mean...................  1330.9997362
     > Variance...............  57.5497899824
     > Skewness Coefficient...  0.0229295468388
     > Kurtosis Coefficient...  2.99662898689
    
    >>> Q.plot()  # currently shows a histogram of the resultant calculations

Distribution Creation
=====================

Since variables are created using the ``scipy.stats`` distributions, here is a table that conveniently describes how to construct many of the most common statistical continuous distributions::

=======================  =============  ===================  ====  =========
Distribution             scipy.stats    args                 loc   scale
                         distribution   (shape parameters)
                         name 
=======================  =============  ===================  ====  =========
Normal(mu,sigma)         norm                                mu    sigma
Uniform(a,b)             uniform                             a     b-a
Exponential(lambda)      expon                                     1/lambda
Gamma(k,theta)           gamma          k                          theta
Beta(alpha,beta,[a,b])   beta           alpha,beta           a     b-a
Log-Normal(mu,sigma)     lognorm        sigma                mu
Chi-Square(dv)           chi2           dv
F(df_numer,df_denom)     f              df_numer,df_denom
Triangular(a,b,peak)     triang         peak                 a     b-a
Student-T(df)            t              df
=======================  =============  ===================  ====  =========

Thus, each distribution above would have the same call signature::
    
    >>> import scipy.stats as ss
    >>> ss.your_dist_here(args, loc=loc, scale=scale)
    

Main Features
=============

1. **Transparent calculations** with derivatives automatically calculated. **No or little modification** to existing code required.
2. Basic `NumPy` support without modification. Vectorized calculations built-in to the ``ad`` package.
3. Nearly all standard `math`_ module functions supported through the ``soerp.umath`` sub-module. If you think a function is in there, it probably is.
4. Nearly all derivatives calculated analytically using ``ad`` functionality.

Installation
============

**Make sure you have the** `ad`_ **,** `SciPy`_ **,** `NumPy`_ **packages installed!**

You have several easy, convenient options to install the ``mcerp`` package 
(administrative privileges may be required)

1. Download the package files below, unzip to any directory, and run ``python setup.py install`` from the command-line
2. Simply copy the unzipped ``mcerp-XYZ`` directory to any other location that python can find it and rename it ``mcerp``
3. If ``setuptools`` is installed, run ``easy_install --upgrade mcerp`` from the command-line
4. If ``pip`` is installed, run ``pip --upgrade mcerp`` from the command-line

Contact
=======

Please send **feature requests, bug reports, or feedback** to 
`Abraham Lee`_.


  
.. _Monte Carlo methods: http://en.wikipedia.org/wiki/Monte_Carlo_method
.. _latin-hypercube sampling: http://en.wikipedia.org/wiki/Latin_hypercube_sampling
.. _soerp: http://pypi.python.org/pypi/soerp
.. _error propagation: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _math: http://docs.python.org/library/math.html
.. _ad: http://pypi.python.org/pypi/ad
.. _NumPy: http://www.numpy.org/
.. _SciPy: http://scipy.org
.. _Matplotlib: http://matplotlib.org/
.. _uncertainties: http://pypi.python.org/pypi/uncertainties
.. _Abraham Lee: mailto: tisimst@gmail.com
.. _PEP8: http://www.python.org/dev/peps/pep-0008
