"""
===============================================================================
mcerp: Real-time latin-hypercube-sampling-based Monte Carlo Error Propagation
===============================================================================

Generalizes many statistical functions that work on numeric objects (from the 
scipy.stats module) to be compatible with objects defined by statistical
distributions. 

NOTE: Although all of these functions can be used without this import, this 
package was created for convenience and transparent operation. For usage,
see the respective documentation at 

http://docs.scipy.org/doc/scipy/reference/stats.html#statistical-functions

Author: Abraham Lee
Copyright: 2013
"""
from mcerp import UncertainFunction
import numpy as np
import scipy.stats as ss

__author__ = 'Abraham Lee'

def wrap(func):
    def wrappedfunc(*args, **kwargs):
        """
        Wraps a Scipy.Stats (or any) function, checking for MCERP objects 
        as non-keyword arguments
        """
        tmpargs = []
        for arg in args:
            if isinstance(arg, UncertainFunction):
                tmpargs += [arg._mcpts]
            else:
                tmpargs += [arg]
        args = tuple(tmpargs)

        return func(*args, **kwargs)
    
    wrappedfunc.__name__ = func.__name__
    wrappedfunc.__doc__ = func.__doc__
    
    return wrappedfunc
    
describe = wrap(ss.describe)
gmean = wrap(ss.gmean)
hmean = wrap(ss.hmean)
kurtosis = wrap(ss.kurtosis)
kurtosistest = wrap(ss.kurtosistest)
mode = wrap(ss.mode)
moment = wrap(ss.moment)
normaltest = wrap(ss.normaltest)
skew = wrap(ss.skew)
skewtest = wrap(ss.skewtest)
tmean = wrap(ss.tmean)
tvar = wrap(ss.tvar)
tmin = wrap(ss.tmin)
tmax = wrap(ss.tmax)
tstd = wrap(ss.tstd)
tsem = wrap(ss.tsem)
variation = wrap(ss.variation)
percentileofscore = wrap(ss.percentileofscore)
scoreatpercentile = wrap(ss.scoreatpercentile)
bayes_mvs = wrap(ss.bayes_mvs)
sem = wrap(ss.sem)
zmap = wrap(ss.zmap)
zscore = wrap(ss.zscore)
f_oneway = wrap(ss.f_oneway)
pearsonr = wrap(ss.pearsonr)
spearmanr = wrap(ss.spearmanr)
pointbiserialr = wrap(ss.pointbiserialr)
kendalltau = wrap(ss.kendalltau)
linregress = wrap(ss.linregress)
ttest_1samp = wrap(ss.ttest_1samp)
ttest_ind = wrap(ss.ttest_ind)
ttest_rel = wrap(ss.ttest_rel)
kstest = wrap(ss.kstest)
chisquare = wrap(ss.chisquare)
# power_divergence = wrap(ss.power_divergence)
ks_2samp = wrap(ss.ks_2samp)
mannwhitneyu = wrap(ss.mannwhitneyu)
rankdata = wrap(ss.rankdata)
ranksums = wrap(ss.ranksums)
wilcoxon = wrap(ss.wilcoxon)
kruskal = wrap(ss.kruskal)
friedmanchisquare = wrap(ss.friedmanchisquare)
ansari = wrap(ss.ansari)
bartlett = wrap(ss.bartlett)
levene = wrap(ss.levene)
shapiro = wrap(ss.shapiro)
anderson = wrap(ss.anderson)
binom_test = wrap(ss.binom_test)
fligner = wrap(ss.fligner)
mood = wrap(ss.mood)


