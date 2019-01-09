========================
``mcerp`` Python Package
========================

.. image:: https://dev.azure.com/tisimst/tisimst/_apis/build/status/tisimst.mcerp
    :target: https://dev.azure.com/tisimst/tisimst/_build/latest?definitionId=1
    :alt: Build Status

- Code: https://github.com/tisimst/mcerp
- Documentation: (not online yet, for now, see the doc folder on Github)
- License: BSD-3-Clause

Overview
========

``mcerp`` is a stochastic calculator for `Monte Carlo methods`_ that uses 
`latin-hypercube sampling`_ to perform non-order specific 
`error propagation`_ (or uncertainty analysis). 

With this package you can **easily** and **transparently** track the effects
of uncertainty through mathematical calculations. Advanced mathematical 
functions, similar to those in the standard `math`_ module, and statistical
functions like those in the `scipy.stats`_ module, can also be evaluated 
directly.

If you are familiar with Excel-based risk analysis programs like *@Risk*, 
*Crystal Ball*, *ModelRisk*, etc., this package **will work wonders** for you
(and probably even be faster!) and give you more modelling flexibility with 
the powerful Python language. This package also *doesn't cost a penny*, 
compared to those commercial packages which cost *thousands of dollars* for a 
single-seat license. Feel free to copy and redistribute this package as much 
as you desire!

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
   articles and other relevant web pages.

5. **Correlation enforcement** and variable sample visualization capabilities.

6. **Probability calculations** using conventional comparison operators.

7. Advanced Scipy **statistical function compatibility** with package 
   functions. Depending on your version of Scipy, some functions might not
   work.

Installation
============

``mcerp`` works on Linux, MacOS and Windows, with Python 2.7 or with Python 3.5 or later.

To install it, use ``pip``::

    pip install mcerp

The ``mcerp`` dependencies should be installed automatically if using ``pip``,
otherwise they will need to be installed manually:

- `NumPy`_ : Numeric Python
- `SciPy`_ : Scientific Python
- `Matplotlib`_ : Python plotting library

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
.. _package documentation: http://pythonhosted.org/mcerp
.. _GitHub: http://github.com/tisimst/mcerp
