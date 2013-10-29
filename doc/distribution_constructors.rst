
.. _using distributions:

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


.. _scipy.stats: http://docs.scipy.org/doc/scipy/reference/stats.html
.. _source code: https://github.com/tisimst/mcerp
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