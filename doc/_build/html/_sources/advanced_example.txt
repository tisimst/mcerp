
.. index:: Advanced Examples

.. _advanced example:

Advanced Examples
=================

In this section, some more advanced/complex usages will be presented.

Volumetric Gas Flow Through an Orifice Meter
--------------------------------------------

Here's an interesting engineering example, showing how the 
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

Manufacturing Tolerance Stackup
-------------------------------

In a certain manufacturing plant, a record of process data is taken and
just happens to fit a Gamma distribution with a mean of 1.5 and a
variance of 0.25. With an assembly of three of these parts, our analysis
for tolerance stackup would proceed as follows. We first need to convert
the distribution statistics to the shape parameters k and theta, which 
are found by the relations:

    * ``k = mean**2/variance = (1.5)**2/(0.25) = 9``
    * ``theta = variance/mean = (0.25)/(1.5) = 1/6``

We can now proceed with the analysis::

    
    >>> k = 9.0
    >>> theta = 1/6.0
    >>> x = Gamma(k, theta)
    >>> y = Gamma(k, theta)
    >>> z = Gamma(k, theta)
    >>> w = x + y + z
    >>> w.describe()
    MCERP Uncertain Value:
     > Mean...................  4.50000470462
     > Variance...............  0.76376726781
     > Skewness Coefficient...  0.368543723948
     > Kurtosis Coefficient...  3.18692837067

Due to the skewed nature of the inputs, the output is also slightly skewed.

Scheduling Facilities
---------------------

At a manufacturing plant, that has six test-and-repair stations, data 
is collected about the amount of time that product is present at each 
station:

    1. Station 1: Normally distributed, mean of 10 hours, variance of
       1 hour.
    2. Station 2: Normally distributed, mean of 20 hours, variance of
       2 hours.
    3. Station 3: Gamma distributed, mean of 1.5 hours, variance of
       0.25 hours.
    4. Station 4: Gamma distributed, mean of 10 hours, variance of 10
       hours.
    5. Station 5: Exponentially distributed, mean of 0.2 hours, variance
       of 0.04 hours (decay constant lamda=5).
    6. Station 6: Chi-squared distributed, mean of 10 hours, variance of
       20 hours (degree of freedom v=10).
    
Management wants to understand the uncertainty associated with the whole
process::

    # Station 1
    >>> s1 = N(10, 1)
    
    # Station 2
    >>> s2 = N(20, 2**0.5)
    
    # Station 3
    >>> mn1 = 1.5
    >>> vr1 = 0.25
    >>> k1 = mn1**2/vr1
    >>> theta1 = vr1/mn1
    >>> s3 = Gamma(k1, theta1)
    
    # Station 4
    >>> mn2 = 10
    >>> vr2 = 10
    >>> k2 = mn2**2/vr2
    >>> theta2 = vr2/mn2
    >>> s4 = Gamma(k2, theta2)
    
    # Station 5
    >>> s5 = Exp(5)
    
    # Station 6
    >>> s6 = Chi2(10)
    
    >>> T = s1 + s2 + s3 + s4 + s5 + s6
    >>> T.describe()
    MCERP Uncertain Value:
     > Mean...................  51.6999259156
     > Variance...............  33.6983675299
     > Skewness Coefficient...  0.520212339449
     > Kurtosis Coefficient...  3.52754453865

From this analysis, we see that the average process time is about 51.7 hours,
but the variance is quite large (standard deviation = 5.8 hours). This
gives management the desire to understand which is the greatest contributors,
so we can analyze the standard deviations of each process step::

    >>> for i, si in enumerate([s1, s2, s3, s4, s5, s6]):
    ...     print 'Station', i + 1, ':', si.std
    ...
    Station 1 : 0.9998880644
    Station 2 : 1.41409415266
    Station 3 : 0.499878358909
    Station 4 : 3.16243741632
    Station 5 : 0.199970343107
    Station 6 : 4.47143708522    

This would seem to indicate that management could focus their efforts on
making the cycle times of stations 4 and 6 more consistent.

It may also be useful to understand the probability that a complete cycle
will exceed a certain amount, say at 59, 62 and 68 hours::

    >>> prob = [T>hr for hr in [59, 62, 68]]
    >>> prob
    [0.1091, 0.0497, 0.0083]
    
That is to say that it is expected that the entire process will take 59
hours approximately 11% of the time, 62 hours 5% of the time, and 68 hours
about 1% of the time.

                                                                                                                                                                                                                                                                
    
