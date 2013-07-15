``mcerp`` Package Documentation
===============================

Overview
--------

``mcerp`` is an on-the-fly calculator for `Monte Carlo methods`_ that uses 
`latin-hypercube sampling`_ (see soerp_ for the Python implementation of the
analytical second-order error propagation original Fortran code `SOERP` by 
N. D. Cox) to perform non-order specific `error propagation`_ (or uncertainty
analysis). The ``mcerp`` package allows you to **easily** and 
**transparently** track the effects of uncertainty through mathematical 
calculations. Advanced mathematical functions, similar to those in the 
standard math_ module can also be evaluated directly. 

Due to the nature of random sampling techniques, calculation results will vary 
from session to session (but be consistent within the session) since new 
latin-hypercube samples are only generated when variables are newly defined or 
re-defined. By default, each variable uses 10,000 latin-hypercube samples that 
are sufficiently random. The number of samples can be changed by assigning an 
integer value to the ``mcerp.npts`` object (typically, values between 1,000 and 
1,000,000 are sufficiently large to insure small errors in the resulting 
statistics). This should only be changed prior to performing calculations since 
all subsequent calculations assume that each input has the same number of 
sampled points (this could be changed through resampling, I suppose...).

In order to correctly use ``mcerp``, knowledge of the distributions from the 
``scipy.stats`` sub-module is required. The result of all calculations 
generates a *mean*, *variance*, and *standardized skewness and kurtosis* 
coefficients. 


Required Packages
-----------------

- NumPy_ : Numeric Python

- SciPy_ : Scientific Python

Suggested Packages
------------------

- Matplotlib_ : Python plotting library (required for ``mcerp.plot``)

Basic examples
--------------
::

    >>> from mcerp import *  # N, U, Gamma, Beta, etc.
    >>> x1 = N(24, 1) # normally distributed
    >>> x2 = N(37, 4) # normally distributed
    >>> x3 = Exp(2)  # exponentially distributed

    >>> Z = (x1*x2**2)/(15*(1.5 + x3))
    >>> Z
    uv(1161.14296871, 116093.134064, 0.361152281239, 3.05247793644)

    >>> Z.describe()  # explains what the moments are
    MCERP Uncertain Value:
     > Mean...................  1161.14296871
     > Variance...............  116093.134064
     > Skewness Coefficient...  0.361152281239
     > Kurtosis Coefficient...  3.05247793644

    >>> x1.stats  # the eight moments can be accessed at any time
    [24.0, 1.0, 0.0, 3.0]
    
    >>> x1.plot()  # if matplotlib is installed it shows the distribution

    >>> x1-x1  # correlations are correctly handled
    0.0
    
    # Example of volumetric gas flow through orifice meter
    >>> import mcerp.umath as umath  # sin, exp, sqrt, etc.
    >>> H = N(64, 0.5)
    >>> M = N(16, 0.1)
    >>> P = N(361, 2)
    >>> t = N(165, 0.5)
    >>> C = 38.4
    >>> Q = C*umath.sqrt((520*H*P)/(M*(t + 460)))
    >>> Q.describe()
    MCERP Uncertain Value:
     > Mean...................  1330.9997362
     > Variance...............  57.5497899824
     > Skewness Coefficient...  0.0229295468388
     > Kurtosis Coefficient...  2.99662898689
    
    >>> Q.plot()  # shows a kde of the resultant calculations

Using Distributions
-------------------

Since all of the variables in MCERP are statistical distributions, they are 
created using the ``scipy.stats`` distributions. There are also some 
convenience constructors that should make defining a distribution easier, 
though it's not necessary to use them. See the source code for info that 
describes how to construct many of the most common statistical continuous 
and discrete distributions using the ``scipy.stats`` distributions.

Though its not entirely discouraged to use the ``scipy.stats`` distributions
directly, here are the **equivalent constructors** that I've found to be 
**easier to use** (the location, scale, and shape parameters are described in 
the respective Wikipedia pages):

Continuous Distributions

- ``N(mu, sigma)`` : `Normal distribution`_

- ``U(a, b)`` : `Uniform distribution`_

- ``Exp(lamda, [mu])`` : `Exponential distribution`_

- ``Gamma(k, theta)`` : `Gamma distribution`_

- ``Beta(alpha, beta, [a, b])`` : `Beta distribution`_

- ``LogN(mu, sigma)`` : `Log-normal distribution`_

- ``X2(k)`` : `Chi-squared distribution`_

- ``F(d1, d2)`` : `F-distribution`_

- ``Tri(a, b, c)`` : `Triangular distribution`_

- ``T(v)`` : `T-distribution`_

- ``Weib(lamda, k)`` : `Weibull distribution`_

Discrete Distributions

- ``Bern(p)`` : `Bernoulli distribution`_

- ``B(n, p)`` : `Binomial distribution`_

- ``G(p)`` : `Geometric distribution`_

- ``H(M, n, N)`` : `Hypergeometric distribution`_

- ``Pois(lamda)`` : `Poisson distribution`_

For example, the following constructions are equivalent::

    >>> x = uv(ss.norm(loc=10, scale=1))  # scipy.stats distribution
    >>> x = N(10, 1)  # nicer constructor IMHO :)

Main Features
-------------

1. **Transparent calculations**. **No or little modification** to existing 
   code required.
    
2. Basic NumPy_ support without modification. (I haven't done extensive 
   testing, so please let me know if you encounter bugs.)

3. Advanced mathematical functions supported through the ``mcerp.umath`` 
   sub-module. If you think a function is in there, it probably is. If it 
   isn't, please request it!

4. **Easy continuous distribution constructors**. The location, scale, 
   and shape parameters follow the notation in the respective Wikipedia 
   articles.

Installation
------------

**Make sure you have the**  SciPy_ **and** NumPy_ **packages installed!**
This package won't work without them.

You have several easy, convenient options to install the ``mcerp`` package 
(administrative privileges may be required)

1. Manually download the package files below, unzip to any directory, and run 
   ``python setup.py install`` from the command-line.

2. Simply copy the unzipped ``mcerp-XYZ`` directory to any other location that
   python can find it and rename it ``mcerp``.
    
3. If ``setuptools`` is installed, run ``easy_install --upgrade mcerp`` from 
   the command-line.
    
4. If ``pip`` is installed, run ``pip --upgrade mcerp`` from the command-line.

Python 3
--------

To use this package with Python 3.x, you will need to run the ``2to3`` 
conversion tool at the command-line using the following syntax while in the 
unzipped ``mcerp`` directory::

    $ 2to3 -w -f all *.py
    
This should take care of the main changes required. If bugs continue to pop up,
please email the author.
    
See also
--------

- uncertainties_ : First-order uncertainty propagation

- soerp_ : Second Order ERror Propagation

Contact
-------

Please send **feature requests, bug reports, or feedback** to 
`Abraham Lee`_.


    
.. _Monte Carlo methods: http://en.wikipedia.org/wiki/Monte_Carlo_method
.. _latin-hypercube sampling: http://en.wikipedia.org/wiki/Latin_hypercube_sampling
.. _soerp: http://pypi.python.org/pypi/soerp
.. _error propagation: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _math: http://docs.python.org/library/math.html
.. _NumPy: http://www.numpy.org/
.. _SciPy: http://scipy.org
.. _Matplotlib: http://matplotlib.org/
.. _uncertainties: http://pypi.python.org/pypi/uncertainties
.. _Abraham Lee: mailto: tisimst@gmail.com
.. _PEP8: http://www.python.org/dev/peps/pep-0008
.. _Normal distribution: http://en.wikipedia.org/wiki/Normal_distribution
.. _Uniform distribution: http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
.. _Exponential distribution: http://en.wikipedia.org/wiki/Exponential_distribution
.. _Gamma distribution: http://en.wikipedia.org/wiki/Gamma_distribution
.. _Beta distribution: http://en.wikipedia.org/wiki/Beta_distribution
.. _Log-normal distribution: http://en.wikipedia.org/wiki/Log-normal_distribution
.. _Chi-squared distribution: http://en.wikipedia.org/wiki/Chi-squared_distribution
.. _F-distribution: http://en.wikipedia.org/wiki/F-distribution
.. _Triangular distribution: http://en.wikipedia.org/wiki/Triangular_distribution
.. _T-distribution: http://en.wikipedia.org/wiki/Student's_t-distribution
.. _Weibull distribution: http://en.wikipedia.org/wiki/Weibull_distribution
.. _Bernoulli distribution: http://en.wikipedia.org/wiki/Bernoulli_distribution
.. _Binomial distribution: http://en.wikipedia.org/wiki/Binomial_distribution
.. _Geometric distribution: http://en.wikipedia.org/wiki/Geometric_distribution
.. _Hypergeometric distribution: http://en.wikipedia.org/wiki/Hypergeometric_distribution
.. _Poisson distribution: http://en.wikipedia.org/wiki/Poisson_distribution