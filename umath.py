"""
Generalizes mathematical operators that work on numeric objects (from the math
module) compatible with objects with uncertainty distributions
"""
from mcerp import to_uncertain_func,_make_UF_compatible_object
import ad.admath as admath
import numpy as np
#import sys

__author__ = 'Abraham Lee'

__all__ = admath.__all__

def abs(x):
    return fabs(x)
    
def acos(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acos(uf._mcpts))
    return _make_UF_compatible_object(admath.acos(x),mcpts)

def acosh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acosh(uf._mcpts))
    return _make_UF_compatible_object(admath.acosh(x),mcpts)

def acot(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acot(uf._mcpts))
    return _make_UF_compatible_object(admath.acot(x),mcpts)

def acoth(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acoth(uf._mcpts))
    return _make_UF_compatible_object(admath.acoth(x),mcpts)

def acsc(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acsc(uf._mcpts))
    return _make_UF_compatible_object(admath.acsc(x),mcpts)

def acsch(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.acsch(uf._mcpts))
    return _make_UF_compatible_object(admath.acsch(x),mcpts)

def asec(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.asec(uf._mcpts))
    return _make_UF_compatible_object(admath.asec(x),mcpts)

def asech(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.asech(uf._mcpts))
    return _make_UF_compatible_object(admath.asech(x),mcpts)

def asin(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.asin(uf._mcpts))
    return _make_UF_compatible_object(admath.asin(x),mcpts)

def asinh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.asinh(uf._mcpts))
    return _make_UF_compatible_object(admath.asinh(x),mcpts)

def atan(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.atan(uf._mcpts))
    return _make_UF_compatible_object(admath.atan(x),mcpts)

def atanh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.atanh(uf._mcpts))
    return _make_UF_compatible_object(admath.atanh(x),mcpts)

def ceil(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.ceil(uf._mcpts))
    return _make_UF_compatible_object(admath.ceil(x),mcpts)

def cos(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.cos(uf._mcpts))
    return _make_UF_compatible_object(admath.cos(x),mcpts)

def cosh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.cosh(uf._mcpts))
    return _make_UF_compatible_object(admath.cosh(x),mcpts)

def cot(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.cot(uf._mcpts))
    return _make_UF_compatible_object(admath.cot(x),mcpts)

def coth(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.coth(uf._mcpts))
    return _make_UF_compatible_object(admath.coth(x),mcpts)

def csc(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.csc(uf._mcpts))
    return _make_UF_compatible_object(admath.csc(x),mcpts)

def csch(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.csch(uf._mcpts))
    return _make_UF_compatible_object(admath.csch(x),mcpts)

def degrees(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.degrees(uf._mcpts))
    return _make_UF_compatible_object(admath.degrees(x),mcpts)

def erf(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.erf(uf._mcpts))
    return _make_UF_compatible_object(admath.erf(x),mcpts)

def erfc(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.erfc(uf._mcpts))
    return _make_UF_compatible_object(admath.erfc(x),mcpts)

def exp(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.exp(uf._mcpts))
    return _make_UF_compatible_object(admath.exp(x),mcpts)

def expm1(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.expm1(uf._mcpts))
    return _make_UF_compatible_object(admath.expm1(x),mcpts)

def fabs(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.fabs(uf._mcpts))
    return _make_UF_compatible_object(admath.fabs(x),mcpts)

def factorial(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.factorial(uf._mcpts))
    return _make_UF_compatible_object(admath.factorial(x),mcpts)

def floor(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.floor(uf._mcpts))
    return _make_UF_compatible_object(admath.floor(x),mcpts)

def gamma(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.gamma(uf._mcpts))
    return _make_UF_compatible_object(admath.gamma(x))

def lgamma(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.lgamma(uf._mcpts))
    return _make_UF_compatible_object(admath.lgamma(x))

def hypot(x,y):
    ufx = to_uncertain_func(x)
    ufy = to_uncertain_func(y)
    mcpts = np.array(admath.hypot(ufx._mcpts,ufy._mcpts))
    return _make_UF_compatible_object(admath.hypot(x,y),mcpts)

def ln(x):
    return log(x)

def log(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.log(uf._mcpts))
    return _make_UF_compatible_object(admath.log(x),mcpts)

def log10(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.log10(uf._mcpts))
    return _make_UF_compatible_object(admath.log10(x),mcpts)

def log1p(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.log1p(uf._mcpts))
    return _make_UF_compatible_object(admath.log1p(x),mcpts)

def pow(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.pow(uf._mcpts))
    return _make_UF_compatible_object(admath.pow(x),mcpts)

def radians(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.radians(uf._mcpts))
    return _make_UF_compatible_object(admath.radians(x),mcpts)

def sec(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.sec(uf._mcpts))
    return _make_UF_compatible_object(admath.sec(x),mcpts)

def sech(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.sech(uf._mcpts))
    return _make_UF_compatible_object(admath.sech(x),mcpts)

def sin(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.sin(uf._mcpts))
    return _make_UF_compatible_object(admath.sin(x),mcpts)

def sinh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.sinh(uf._mcpts))
    return _make_UF_compatible_object(admath.sinh(x),mcpts)

def sqrt(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.sqrt(uf._mcpts))
    return _make_UF_compatible_object(admath.sqrt(x),mcpts)

def tan(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.tan(uf._mcpts))
    return _make_UF_compatible_object(admath.tan(x),mcpts)

def tanh(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.tanh(uf._mcpts))
    return _make_UF_compatible_object(admath.tanh(x),mcpts)

def trunc(x):
    uf = to_uncertain_func(x)
    mcpts = np.array(admath.trunc(uf._mcpts))
    return _make_UF_compatible_object(admath.trunc(x),mcpts)


    
