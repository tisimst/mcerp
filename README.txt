Overview
========

``mcerp`` is an on-the-fly calculator for `Monte Carlo methods`_ that uses 
`latin-hypercube sampling`_ (see `soerp`_ for the Python implementation of the
analytical second-order error propagation original Fortran code `SOERP` by 
N. D. Cox) to perform non-order specific `error propagation`_ (or uncertainty
 analysis). The ``mcerp`` package allows you to **easily** and 
 **transparently** track the effects of uncertainty through mathematical 
 calculations. Advanced mathematical functions, similar to those in the 
 standard `math`_ module can also be evaluated directly. 

Due to the nature of random sampling techniques, calculation results will vary 
from session to session (but be consistent within the session) since new 
latin-hypercube samples are only generated when variables are newly defined or 
re-defined. By default, each variable uses 10,000 latin-hypercube samples that 
are sufficiently random. This can be changed by assigning an integer value to 
the ``mcerp.npts`` object (typically, values between 1,000 and 1,000,000 are 
sufficiently large to insure small errors in the resulting statistics). This 
should only be changed prior to performing calculations since all subsequent 
calculations assume that each input has the same number of sampled points 
(this could be changed through resampling, I suppose...).

In order to correctly use ``mcerp``, knowledge of the distributions from the 
``scipy.stats`` sub-module is required. The result of all calculations 
generates a *mean*, *variance*, and *standardized skewness and kurtosis* 
coefficients. Running the source code by itself in a terminal with ``python 
mcerp.py`` will demonstrate some more advanced examples.


Required Packages
=================

- `NumPy`_ : Numeric Python
- `SciPy`_ : Scientific Python

Suggested Packages
==================

- `Matplotlib`_ : Python plotting library (required for ``mcerp.plot``)

Other Recommendations
=====================

- `uncertainties`_ : First-order uncertainty propagation
- `soerp`_ : Second Order ERror Propagation

Basic examples
==============
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

Distribution Creation
=====================

Since all of the variables in MCERP are statistical distributions, they are 
created using the ``scipy.stats`` distributions. There are also some 
convenience constructors that should make defining a distribution easier, 
though it's not necessary to use them. Here is a table that describes how to 
construct many of the most common statistical continuous and discrete 
distributions using the ``scipy.stats`` distributions. See the source code for
help using these distributions.

Though its not entirely discouraged to use the ``scipy.stats`` distributions
directly, here are the **equivalent constructors** that I've found to be 
**easier to use** (the location, scale, and shape parameters are described in 
the respective Wikipedia pages):

- Continuous Distributions
    - `N(mu, sigma)`_ : Normal distribution
    - `U(a, b)`_ : Uniform distribution (continuous)
    - `Exp(lamda, [mu])`_ : Exponential distribution
    - `Gamma(k, theta)`_ : Gamma distribution
    - `Beta(alpha, beta, [a, b])`_ : Beta distribution
    - `LogN(mu, sigma)`_ : Log-normal distribution
    - `X2(df)`_ : Chi-squared distribution
    - `F(dfn, dfd)`_ : F (or, Fisher) distribution
    - `Tri(a, b, c)`_ : Triangular distribution
    - `T(df)`_ : Student's t-distribution
    - `Weib(lamda, k)`_ : Weibull distribution
- Discrete Distributions
    - `Bern(p)`_ : Bernoulli distribution
    - `B(n, p)`_ : Binomial distribution
    - `G(p)`_ : Geometric distribution
    - `H(M, n, N)`_ : Hypergeometric distribution
    - `Pois(lamda)`_ : Poisson distribution

Thus, the following are equivalent::

    >>> x = uv(ss.norm(loc=10, scale=1))  # scipy.stats distribution
    >>> x = N(10, 1)  # nice constructor :)

Main Features
=============

1. **Transparent calculations**. **No or little modification** to existing 
   code required.
   
2. Basic `NumPy` support without modification.

3. Nearly all standard `math`_ module functions supported through the 
   ``mcerp.umath`` sub-module. If you think a function is in there, it 
   probably is. If it isn't, please request it!

4. **Easy continuous distribution constructors**. The location, scale, 
   and shape parameters follow the notation in the respective Wikipedia 
   articles.

Installation
============

**Make sure you have the**  `SciPy`_ **,** `NumPy`_ **packages installed!**

You have several easy, convenient options to install the ``mcerp`` package 
(administrative privileges may be required)

1. Download the package files below, unzip to any directory, and run 
   ``python setup.py install`` from the command-line.

2. Simply copy the unzipped ``mcerp-XYZ`` directory to any other location that
   python can find it and rename it ``mcerp``.
   
3. If ``setuptools`` is installed, run ``easy_install --upgrade mcerp`` from 
   the command-line.
   
4. If ``pip`` is installed, run ``pip --upgrade mcerp`` from the command-line.

Python 3
========

To use this package with Python 3.x, you will need to run the ``2to3`` tool at
the command-line using the following syntax while in the unzipped ``ad`` 
directory::

    $ 2to3 -w -f all *.py
    
This should take care of the main changes required. If bugs continue to pop up,
please email the author.
    
Contact
=======

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
.. _N(mu, sigma): http://en.wikipedia.org/wiki/Normal_distribution
.. _U(a, b): http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
.. _Exp(lamda, [mu]): http://en.wikipedia.org/wiki/Exponential_distribution
.. _Gamma(k, theta): http://en.wikipedia.org/wiki/Gamma_distribution
.. _Beta(alpha, beta, [a, b]): http://en.wikipedia.org/wiki/Beta_distribution
.. _LogN(mu, sigma): http://en.wikipedia.org/wiki/Log-normal_distribution
.. _X2(df): http://en.wikipedia.org/wiki/Chi-squared_distribution
.. _F(dfn, dfd): http://en.wikipedia.org/wiki/F-distribution
.. _Tri(a, b, c): http://en.wikipedia.org/wiki/Triangular_distribution
.. _T(df): http://en.wikipedia.org/wiki/Student's_t-distribution
.. _Weib(lamda, k): http://en.wikipedia.org/wiki/Weibull_distribution
.. _Bern(p): http://en.wikipedia.org/wiki/Bernoulli_distribution
.. _B(n, p): http://en.wikipedia.org/wiki/Binomial_distribution
.. _G(p): http://en.wikipedia.org/wiki/Geometric_distribution
.. _H(M, n, N): http://en.wikipedia.org/wiki/Hypergeometric_distribution
.. _Pois(lamda): http://en.wikipedia.org/wiki/Poisson_distribution
