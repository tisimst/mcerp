
.. _basic examples:

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
 