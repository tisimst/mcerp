"""
================================================================================
mcerp: Real-time latin-hypercube-sampling-based Monte Carlo Error Propagation
================================================================================

Generalizes mathematical operators that work on numeric objects (from the math
module or numpy) compatible with objects with uncertainty distributions

Author: Abraham Lee
Copyright: 2013
"""
from mcerp import UncertainFunction, to_uncertain_func
import numpy as np
import umath as um

__author__ = 'Tim Galvin'

#!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Conflicts with the abs keyword
#!!!!!!!!!!!!!!!!!!!!!!!!!!!

# abs = np.vectorize(um.abs)




"""
Inverse cosine
"""
acos = np.vectorize(um.acos)

"""
Inverse hyperbolic cosine
"""
acosh = np.vectorize(um.acosh)

"""
Inverse sine
"""
asin = np.vectorize(um.asin)

"""
Inverse hyperbolic sine
"""
asinh = np.vectorize(um.asinh)

"""
Inverse tangent
"""
atan = np.vectorize(um.atan)

"""
Inverse hyperbolic tangent
"""
atanh = np.vectorize(um.atanh)

"""
Ceiling function (round towards positive infinity)
"""
ceil = np.vectorize(um.ceil)

"""
Cosine
"""
cos = np.vectorize(um.cos)

"""
Hyperbolic cosine
"""
cosh = np.vectorize(um.cosh)

"""
Convert radians to degrees
"""
degrees = np.vectorize(um.degrees)

"""
Exponential function
"""
exp = np.vectorize(um.exp)

"""
Calculate exp(x) - 1
"""
expm1 = np.vectorize(um.expm1)

"""
Absolute value function
"""
fabs = np.vectorize(um.fabs)

"""
Floor function (round towards negative infinity)
"""
floor = np.vectorize(um.floor)

"""
Calculate the hypotenuse given two "legs" of a right triangle
"""
hypot = np.vectorize(um.hypot)

"""
Natural logarithm (same as "log(x)")
"""
ln = np.vectorize(um.ln)

"""
Natural logarithm
"""
log = np.vectorize(um.log)

"""
Base-10 logarithm
"""
log10 = np.vectorize(um.log10)

"""
Natural logarithm of (1 + x)
"""
log1p = np.vectorize(um.log1p)

"""
Convert degrees to radians
"""
radians = np.vectorize(um.radians)

"""
Sine
"""
sin = np.vectorize(um.sin)

"""
Hyperbolic sine
"""
sinh = np.vectorize(um.sinh)

"""
Square-root function
"""
sqrt = np.vectorize(um.sqrt)

"""
Tangent
"""
tan = np.vectorize(um.tan)

"""
Hyperbolic tangent
"""
tanh = np.vectorize(um.tanh)

"""
Truncate the values to the integer value without rounding
"""
trunc = np.vectorize(um.trunc)
