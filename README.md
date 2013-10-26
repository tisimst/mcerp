===============================
``mcerp`` Package Documentation
===============================

Overview
========

``mcerp`` is a stochastic calculator for `Monte Carlo methods`_ that uses 
`latin-hypercube sampling`_ to perform non-order specific 
`error propagation`_ (or uncertainty analysis). With this package you can 
**easily** and **transparently** track the effects of uncertainty 
through mathematical calculations. Advanced mathematical functions, similar 
to those in the standard `math`_ module can also be evaluated directly. 

Due to the nature of random sampling techniques, calculation results will vary 
from session to session (but be consistent within the session) since new 
samples are only generated when variables are newly defined or re-defined. 
By default, each variable uses 10,000 samples that are sufficiently random. 
The number of samples can be changed by assigning an integer value to the 
``mcerp.npts`` object (typically, values between 1,000 and 1,000,000 are 
sufficiently large to ensure small errors in the resulting statistics). 
This should only be changed prior to performing calculations since 
all subsequent calculations assume that each input has the same number of 
sampled points (this could be changed through resampling, I suppose...).

In order to correctly use ``mcerp``, knowledge of the distributions from 
the ``scipy.stats`` sub-module is helpful, but not required since many
of the most common distributions can be created with convenient functions
(see table below). 
The result of all calculations generates a *mean*, *variance*, and 
*standardized skewness and kurtosis* coefficients (this means that a 
*Normal* distribution has a kurtosis of 3, **NOT** 0).

What's New In This Release
==========================

- New ways that *equality* and *inequality* comparison operators function
  (see `the examples`_ below).
  
  - **Two uncertain objects** created through MCERP (like 
    ``N(4, 1)<=Exp(3)``): Statistical tests are used to decide whether 
    ``<=``, ``>``, etc. hold for the means of the two objects.
  - **One uncertain object and a scalar** (like ``N(3, 0.2)<2``): An
    estimated probability that the distribution lies above/below that value
    is calculated.

- Some important bug fixes in the `distribution constructors`_. These
  should all be working correctly now.

Required Packages
=================

The following packages should be installed automatically (if using ``pip``
or ``easy_install``), otherwise they will need to be installed manually:

- `NumPy`_ : Numeric Python
- `SciPy`_ : Scientific Python
- `Matplotlib`_ : Python plotting library

These packages come standard in *Python(x,y)*, *Spyder*, and other 
scientific computing python bundles.

Basic examples
==============

Let's start with the main import::

    >>> from mcerp import *  # N, U, Gamma, Beta, correlate, etc.

Now, we can construct many kinds of statistical distributions (both 
continuous and discrete). Here's a basic example that involves a 
three-part assembly::

    >>> x1 = N(24, 1) # normally distributed
    >>> x2 = N(37, 4)
    >>> x3 = Exp(2)  # exponentially distributed

The first four moments can be accessed at any time::

    >>> x1.mean
    24.000015983319528
    >>> x1.var
    1.0001064993938098
    >>> x1.skew
    0.00053245280109193898
    >>> x1.kurt
    3.0017697746273813
    >>> x1.stats  # this returns all four as a list
    
Now we'll compute the actual stack-up using normal mathematics and see what 
happens::

    >>> Z = (x1*x2**2)/(15*(1.5 + x3))

We can see how the statistics of each of these distributions propagated 
through the calculations in two basic ways:

#. Telling python to print the object::

    >>> Z  # Explicit "print" command not necessary at the command-line
    uv(1161.35518507, 116688.945979, 0.353867228823, 3.00238273799)

#. Using the ``describe`` class method (provides a more detailed explanation)::

    >>> Z.describe()
    MCERP Uncertain Value:
     > Mean...................  1161.35518507
     > Variance...............  116688.945979
     > Skewness Coefficient...  0.353867228823
     > Kurtosis Coefficient...  3.00238273799
 
Viewing the distribution
------------------------

We can also plot the distributions using the ``plot`` class method (NOTE:
the main import above brings in ``matplotlib.pyplot`` as ``plt`` for
further plotting use)::

    >>> x1.plot()  # No inputs shows the distribution's PDF
    >>> plt.show()  # explicit 'show()' required to display to screen

.. image:: https://raw.github.com/tisimst/mcerp/master/x1.png
   :scale: 60%
   
and for the outputs::

    >>> Z.plot()  # shows the Kernel Density Estimate (KDE) of the data
    >>> plt.show()

.. image:: https://raw.github.com/tisimst/mcerp/master/Z_kde.png
   :scale: 60%

::

    >>> Z.plot(hist=True)  # shows a histogram instead of a KDE
    >>> plt.show()

.. image:: https://raw.github.com/tisimst/mcerp/master/Z_hist.png
   :scale: 60%

or both::

    >>> Z.plot()
    >>> Z.plot(hist=True)
    >>> plt.show()

.. image:: https://raw.github.com/tisimst/mcerp/master/Z_kde_hist.png
   :scale: 60%

.. _the examples:

Probabilities
-------------

To estimate the probability that a value lies above or below some point
in a distribution, we can use the standard inequality comparison 
operators (<, <=, >, >=)::

    >>> x1<21
    0.0014
    
    >>> Z>=1000
    0.6622
    
    >>> x2<x1
    False
    
On the otherhand, if we are comparing distributions to see if one is
"less" or "more" than the other, we actually perform a T-test of the two
objects to compare the two sample means. If the p-value is greater than
0.05 AND the t-statistic has the correct sign, then the comparison will
return ``True``. Let's first create some new samples (the actual values
are contained in the ``_mcpts`` member of the ``UncertainFunction`` class::

    >>> rvs1 = N(5, 10)
    >>> rvs2 = N(5, 10) + N(0, 0.2)
    >>> rvs3 = N(8, 10) + N(0, 0.2)
    
Now, let's compare ``rvs1`` and ``rvs2``. They are similar, but with slightly
different variances, so we would expect the p-value to be large::

    >>> from scipy.stats import ttest_rel
    >>> tstat, pval = ttest_rel(rvs1._mcpts, rvs2._mcpts)
    >>> pval
    0.99888340212679583
    
As expected, because the p-value is essentially, 1.0, the test couldn't tell
them apart, so our comparison returns::

    >>> rvs1<rvs2
    False

However, let's try distributions that are a more different, ``rvs1`` and
``rvs3``. This test should return a smaller p-value and a t-statistic that
we will get the sign from to check the orientation of the comparison::

    >>> tstat, pval = ttest_rel(rvs1._mcpts, rvs3._mcpts)
    >>> pval
    3.0480674044727307e-97

That's a very small p-value, indicating that the distributions are
separated from each other distinctly enough that the test could tell them
apart. Now we need to check the sign of the t-statistic to see if 
``rvs1`` is on the "left" of ``rvs3`` for the comparison::

    >>> float(tstat)
    -21.158661004433682

Because we are using the *less than* comparison and the sign of the 
t-statistic is negative, then we say that this is "oriented" correctly
and, no surprise, we get::

    >>> rvs1<rvs3
    True

If we had tried *greater than*, then we would have gotten the wrong sign
on the t-statistic and the comparison evaluates to ``False``.

One interesting thing about this way of testing two distributions is that
it's possible to get the following::

    >>> x = N(0, 1)
    >>> y = N(0, 10)
    >>> x<y
    False
    >>> x>y
    False
    >>> x==y
    False
    
The equality comparison operators (== and !=) actually test to see if 
the distributions are identical, thus::
    
    >>> x1==x1
    True

    >>> n1 = N(0, 1)
    >>> n2 = N(0, 1)
    >>> n1==n2  # n1 and n2 are independently sampled, so they are not equal
    False
    
    >>> Z*Z==Z**2  # Both sides have the same root samples, so they are equal
    True

Correlations
------------

By default, the samples try to be as uncorrelated and independent as
possible from any other inputs. However, sometimes inputs to have some
degree of correlation between them. If it is desired to force a set of 
variables to have certain correlations, which is not uncommon in 
real-life situations, this can be done with the ``correlate`` function 
(NOTE: this should be done BEFORE any calculations have taken place in 
order to work properly).

For example, let's look at our original example with inputs ``x1``, ``x2``,
and ``x3``::

    # The correlation coefficients before adjusting
    >>> print correlation_matrix([x1, x2, x3])
    [[ 1.          0.00558381  0.01268168]
     [ 0.00558381  1.          0.00250815]
     [ 0.01268168  0.00250815  1.        ]]

You'll notice a few things about the correlation matrix. First, the 
diagonals are all 1.0 (they always are). Second, the matrix is symmetric.
Third, the correlation coefficients in the upper and lower triangular
parts are relatively small. This is how ``mcerp`` is designed. Here is 
what the actual samples looks like in a matrix plot form:

.. image:: https://raw.github.com/tisimst/mcerp/master/before_correlation_matrixplot.png
    :scale: 60%

Now, let's say we desire to impose a -0.75 correlation between ``x1``
and ``x2``. Let's create the desired correlation matrix (note that all 
diagonal elements should be 1.0)::

    # The desired correlation coefficient matrix
    >>> c = np.array([[  1.0, -0.75, 0.0],
    ...               [-0.75,   1.0, 0.0],
    ...               [  0.0,   0.0, 1.0]])

Using the ``mcerp.correlate`` function, we can now apply the desired
correlation coefficients to our samples to try and force the inputs
to be correlated::
    
    # Apply the correlations into the samples (works in-place)
    >>> correlate([x1, x2, x3], c)
    
    # Show the new correlation coefficients
    >>> print correlation_matrix([x1, x2, x3])
    [[  1.00000000e+00  -7.50010477e-01   1.87057576e-03]
     [ -7.50010477e-01   1.00000000e+00   8.53061774e-04]
     [  1.87057576e-03   8.53061774e-04   1.00000000e+00]]
 
The correlation matrix is roughly what we expected within a few percent.
Even the other correlation coefficients are closer to zero than before. If 
any exceptions appear during this operation, it is likely because the
correlation matrix is either not **symmetric**, not **positive-definite**, 
or both.
    
The newly correlated samples will now look something like:

.. image:: https://raw.github.com/tisimst/mcerp/master/after_correlation_matrixplot.png
    :scale: 60%

.. note: This correlation operation doesn't change any of the original sampled
   values, it simply re-organizes them in such a way that they closely
   match the desired correlations.

Now that the inputs' relations have been modified, let's check how 
the output of our stack-up has changed (sometimes the correlations won't
change the output much, but others can change a lot!)::

    # Z should now be a little different
    >>> Z = (x1*x2**2)/(15*(1.5 + x3))
    >>> Z.describe
    MCERP Uncertain Value:
     > Mean...................  1153.710442
     > Variance...............  97123.3417748
     > Skewness Coefficient...  0.211835225063
     > Kurtosis Coefficient...  2.87618465139
 
We can also see what adding that correlation did: reduced the mean,
reduced the variance, increased the skewness, increased the kurtosis.

Advanced Example
================

Here's a *slightly* more advanced engineering example, showing how the 
random effects of input parameters propagates through the calculation of 
the volumetric gas flow through an orifice meter::

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

Interestingly enough, even though the calculations involve multiplication, 
division, and a square-root, the result appears to be very close to a Normal 
distribution (Q ~ N(1331, 7.6)).

.. _distribution constructors:

Using Distributions
===================

Since all of the variables in ``mcerp`` are statistical distributions, they 
are created using the `scipy.stats`_ distributions. There are also some 
convenience constructors that should make defining a distribution easier, 
though it's not necessary to use them. See the `source code`_ of the
``UncertainVariable`` class for info that describes how to construct many 
of the most common statistical continuous and discrete distributions using 
the `scipy.stats`_ distributions.

Here are a set of **equivalent constructors** that I've found to be 
**easier to use** for the most common kinds of distributions (the location, 
scale, and shape parameters are described in the respective Wikipedia pages):

+----------------------------------------------------------------------------------+
| **Continuous Distributions**                                                     |
+-----------------------------------------+----------------------------------------+
| ``N(mu, sigma)``                        | `Normal distribution`_                 |
+-----------------------------------------+----------------------------------------+
| ``U(a, b)``                             | `Uniform distribution`_                |
+-----------------------------------------+----------------------------------------+
| ``Exp(lamda)``                          | `Exponential distribution`_            |
+-----------------------------------------+----------------------------------------+
| ``Gamma(k, theta)``                     | `Gamma distribution`_                  |
+-----------------------------------------+----------------------------------------+
| ``Beta(alpha, beta, [a, b])``           | `Beta distribution`_                   |
+-----------------------------------------+----------------------------------------+
| ``LogN(mu, sigma)``                     | `Log-normal distribution`_             |
+-----------------------------------------+----------------------------------------+
| ``Chi2(k)``                             | `Chi-squared distribution`_            |
+-----------------------------------------+----------------------------------------+
| ``F(d1, d2)``                           | `F-distribution`_                      |
+-----------------------------------------+----------------------------------------+
| ``Tri(a, b, c)``                        | `Triangular distribution`_             |
+-----------------------------------------+----------------------------------------+
| ``T(v)``                                | `T-distribution`_                      |
+-----------------------------------------+----------------------------------------+
| ``Weib(lamda, k)``                      | `Weibull distribution`_                |
+-----------------------------------------+----------------------------------------+
| **Discrete Distributions**                                                       |
+-----------------------------------------+----------------------------------------+
| ``Bern(p)``                             | `Bernoulli distribution`_              |
+-----------------------------------------+----------------------------------------+
| ``B(n, p)``                             | `Binomial distribution`_               |
+-----------------------------------------+----------------------------------------+
| ``G(p)``                                | `Geometric distribution`_              |
+-----------------------------------------+----------------------------------------+
| ``H(N, n, K)``                          | `Hypergeometric distribution`_         |
+-----------------------------------------+----------------------------------------+
| ``Pois(lamda)``                         | `Poisson distribution`_                |
+-----------------------------------------+----------------------------------------+

For example, the following constructions are equivalent::

    # explicitly calling out the scipy.stats distribution
    >>> import scipy.stats as ss
    >>> x = uv(ss.norm(loc=10, scale=1))

    # using a built-in constructor
    >>> x = N(10, 1)

From my experience, the first option can be tedious and difficult to work 
with, but it does allow you to input any distribution defined in the 
`scipy.stats`_ sub-module, both continuous and discrete, if you know how. 
For the most common distributions, the MCERP constructors are hard to beat.
If you feel like another distribution should be included in the "common"
list, let me know!

Main Features
=============

1. **Transparent calculations**. **No or little modification** to existing 
   code required.
    
2. Basic `NumPy`_ support without modification. (I haven't done extensive 
   testing, so please let me know if you encounter bugs.)

3. Advanced mathematical functions supported through the ``mcerp.umath`` 
   sub-module. If you think a function is in there, it probably is. If it 
   isn't, please request it!

4. **Easy statistical distribution constructors**. The location, scale, 
   and shape parameters follow the notation in the respective Wikipedia 
   articles.

5. **Correlation enforcement** and variable sample visualization capabilities.

6. **Probability calculations** using conventional comparison operators.

Installation
============

**Make sure you have the**  `SciPy`_ **and** `NumPy`_ **and** `Matplotlib`_ **packages installed!**
This package won't work without them.

You have **several easy, convenient options** to install the ``mcerp`` 
package (administrative privileges may be required)

#. Simply copy the unzipped ``mcerp-XYZ`` directory to any other location that
   python can find it and rename it ``mcerp``.
    
#. From the command-line, do one of the following:
   
   a. Manually download the package files below, unzip to any directory, and run::
   
       $ [sudo] python setup.py install

   b. If ``setuptools`` is installed, run::

       $ [sudo] easy_install --upgrade mcerp
    
   c. If ``pip`` is installed, run::

       $ [sudo] pip install --upgrade mcerp

Python 3
--------

To use this package with Python 3.x, you will need to run the ``2to3`` 
conversion tool at the command-line using the following syntax while in the 
unzipped ``mcerp`` directory::

    $ 2to3 -w .
    
This should take care of the main changes required. Then, run::

    $ python3 setup.py install

If bugs continue to pop up, please email the author.
    
See also
========

- `uncertainties`_ : First-order error propagation
- `soerp`_ : Second-order error propagation

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
.. _scipy.stats: http://docs.scipy.org/doc/scipy/reference/stats.html
.. _uncertainties: http://pypi.python.org/pypi/uncertainties
.. _source code: https://github.com/tisimst/mcerp
.. _Abraham Lee: mailto:tisimst@gmail.com
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