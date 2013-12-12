
.. index:: Constructing Distributions

.. _using distributions:

Using Distributions
===================

Since all of the variables in ``mcerp`` are statistical distributions, they 
are created internally using the `scipy.stats`_ distributions. There are also 
some convenience constructors that should make defining a distribution easier, 
though it's not necessary to use them. See the `source code`_ of the
``UncertainVariable`` class for info that describes how to construct many 
of the most common statistical continuous and discrete distributions using 
the `scipy.stats`_ distributions (and some others not currently part of
`scipy.stats`).

Tracking
--------

To track variables that may share the same distribution, each of the
constructors below accepts an optional ``tag=...`` kwarg that can be
accessed at anytime using ``x.tag``.

Constructors
------------

Now, if you are like me, the easier the syntax--the better, thus, here are a 
set of **equivalent constructors** that I've found to be **easier to use** 
for the most common kinds of distributions (the location, scale, and shape 
parameters are described in their respective web pages and the `source code`_):

+--------------------------------------------------------------------------------------------------------+
| **Continuous Distributions**                                                                           |
+---------------------------------------------------------------+----------------------------------------+
| ``Beta(alpha, beta, [low, high])``                            | `Beta distribution`_                   |
+---------------------------------------------------------------+----------------------------------------+
| ``Bradford(q, [low, high])``                                  | `Bradford distribution`_               |
+---------------------------------------------------------------+----------------------------------------+
| ``Burr(c, k)``                                                | `Burr distribution`_                   |
+---------------------------------------------------------------+----------------------------------------+
| ``ChiSquared(k)`` or ``Chi2(k)``                              | `Chi-squared distribution`_            |
+---------------------------------------------------------------+----------------------------------------+
| ``Erf(h)``                                                    | `Error Function distribution`_         |
+---------------------------------------------------------------+----------------------------------------+
| ``Erlang(m, b)``                                              | `Erlang distribution`_                 |
+---------------------------------------------------------------+----------------------------------------+
| ``Exponential(lamda)`` or ``Exp(lamda)``                      | `Exponential distribution`_            |
+---------------------------------------------------------------+----------------------------------------+
| ``ExtValueMax(mu, sigma)`` or ``EVMax(mu, sigma)``            | `Extreme Value Maximum distribution`_  |
+---------------------------------------------------------------+----------------------------------------+
| ``ExtValueMin(mu, sigma)`` or ``EVMin(mu, sigma)``            | `Extreme Value Minimum distribution`_  |
+---------------------------------------------------------------+----------------------------------------+
| ``Fisher(d1, d2)`` or ``F(d1, d2)``                           | `F-distribution`_                      |
+---------------------------------------------------------------+----------------------------------------+
| ``Gamma(k, theta)``                                           | `Gamma distribution`_                  |
+---------------------------------------------------------------+----------------------------------------+
| ``LogNormal(mu, sigma)`` or ``LogN(mu, sigma)``               | `Log-normal distribution`_             |
+---------------------------------------------------------------+----------------------------------------+
| ``Normal(mu, sigma)`` or ``N(mu, sigma)``                     | `Normal distribution`_                 |
+---------------------------------------------------------------+----------------------------------------+
| ``Pareto(q, a)`` (first kind)                                 | `Pareto distribution`_                 |
+---------------------------------------------------------------+----------------------------------------+
| ``Pareto2(q, b)`` (second kind)                               | `Pareto2 distribution`_                |
+---------------------------------------------------------------+----------------------------------------+
| ``PERT(low, peak, high)``                                     | `PERT distribution`_                   |
+---------------------------------------------------------------+----------------------------------------+
| ``StudentT(v)`` or ``T(v)``                                   | `T-distribution`_                      |
+---------------------------------------------------------------+----------------------------------------+
| ``Triangular(low, peak, high)`` or ``Tri(low, peak, high)``   | `Triangular distribution`_             |
+---------------------------------------------------------------+----------------------------------------+
| ``Uniform(low, high)`` or ``U(low, high)``                    | `Uniform distribution`_                |
+---------------------------------------------------------------+----------------------------------------+
| ``Weibull(lamda, k)`` or ``Weib(lamda, k)``                   | `Weibull distribution`_                |
+---------------------------------------------------------------+----------------------------------------+
| **Discrete Distributions**                                                                             |
+---------------------------------------------------------------+----------------------------------------+
| ``Bernoulli(p)`` or ``Bern(p)``                               | `Bernoulli distribution`_              |
+---------------------------------------------------------------+----------------------------------------+
| ``Binomial(n, p)`` or ``B(n, p)``                             | `Binomial distribution`_               |
+---------------------------------------------------------------+----------------------------------------+
| ``Geometric(p)`` or ``G(p)``                                  | `Geometric distribution`_              |
+---------------------------------------------------------------+----------------------------------------+
| ``Hypergeometric(N, n, K)`` or ``H(N, n, K)``                 | `Hypergeometric distribution`_         |
+---------------------------------------------------------------+----------------------------------------+
| ``Poisson(lamda)`` or ``Pois(lamda)``                         | `Poisson distribution`_                |
+---------------------------------------------------------------+----------------------------------------+

For example, the following constructions are equivalent::

    # explicitly calling out the scipy.stats distribution
    >>> import scipy.stats as ss
    >>> x = uv(ss.norm(loc=10, scale=1))

    # using a built-in constructor
    >>> x = Normal(10, 1)
    
    # and if there's a short-name alias available
    >>> x = N(10, 1)

From my experience, the first option can be tedious and difficult to work 
with, but it does allow you to input any distribution defined in the 
`scipy.stats`_ sub-module, both continuous and discrete, if you know how. 
For the most common distributions, the MCERP constructors are hard to beat.
If you feel like another distribution should be included in the "common"
list, let me know!


.. _scipy.stats: http://docs.scipy.org/doc/scipy/reference/stats.html
.. _source code: https://github.com/tisimst/mcerp/blob/master/mcerp/__init__.py
.. _Beta distribution: http://en.wikipedia.org/wiki/Beta_distribution
.. _Bradford distribution: http://www.vosesoftware.com/ModelRiskHelp/index.htm#Distributions/Continuous_distributions/Bradford_distribution.htm
.. _Burr distribution: http://en.wikipedia.org/wiki/Burr_distribution
.. _Chi-squared distribution: http://en.wikipedia.org/wiki/Chi-squared_distribution
.. _Error Function distribution: http://www.mathwave.com/articles/error_function_distribution.html
.. _Erlang distribution: http://en.wikipedia.org/wiki/Erlang_distribution
.. _Exponential distribution: http://en.wikipedia.org/wiki/Exponential_distribution
.. _Extreme Value Maximum distribution: http://www.math.uah.edu/stat/special/ExtremeValue.html
.. _Extreme Value Minimum distribution: http://www.math.uah.edu/stat/special/ExtremeValue.html
.. _F-distribution: http://en.wikipedia.org/wiki/F-distribution
.. _Gamma distribution: http://en.wikipedia.org/wiki/Gamma_distribution
.. _Log-normal distribution: http://en.wikipedia.org/wiki/Log-normal_distribution
.. _Normal distribution: http://en.wikipedia.org/wiki/Normal_distribution
.. _Pareto distribution: http://www.vosesoftware.com/ModelRiskHelp/Distributions/Continuous_distributions/Pareto_(first_kind)_distribution.htm
.. _Pareto2 distribution: http://www.vosesoftware.com/ModelRiskHelp/Distributions/Continuous_distributions/Pareto_(second_kind)_distribution.htm
.. _PERT distribution: http://www.vosesoftware.com/ModelRiskHelp/index.htm#Distributions/Continuous_distributions/PERT_distribution.htm
.. _T-distribution: http://en.wikipedia.org/wiki/Student's_t-distribution
.. _Triangular distribution: http://en.wikipedia.org/wiki/Triangular_distribution
.. _Uniform distribution: http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
.. _Weibull distribution: http://en.wikipedia.org/wiki/Weibull_distribution
.. _Bernoulli distribution: http://en.wikipedia.org/wiki/Bernoulli_distribution
.. _Binomial distribution: http://en.wikipedia.org/wiki/Binomial_distribution
.. _Geometric distribution: http://en.wikipedia.org/wiki/Geometric_distribution
.. _Hypergeometric distribution: http://en.wikipedia.org/wiki/Hypergeometric_distribution
.. _Poisson distribution: http://en.wikipedia.org/wiki/Poisson_distribution
