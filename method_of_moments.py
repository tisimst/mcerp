import numpy as np
from copy import copy

assume_linear = False # if True, the sqc and scp parts are ignored

################################################################################
def linear_coef(func,vars):
    lc = np.empty(len(vars),dtype=object)
    for i,v in enumerate(vars):
        lc[i] = func.diff(v)
    return lc
################################################################################
def quad_coef(func,vars,lc=False):
    qc = np.empty(len(vars),dtype=object)
    for i,v in enumerate(vars):
        if lc:
            qc[i] = func[i].diff(v)
        else:
            qc[i] = func.diff(v,v)
    return qc/2
################################################################################
def cross_coef(func,vars,lc=False):
    cp = np.empty((len(vars),len(vars)),dtype=object)
    if len(vars) >= 2:
        for i,vi in enumerate(vars):
            for j,vj in enumerate(vars):
                if j>i:
                    if lc:
                        cp[i,j] = func[i].diff(vj)
                    else:
                        cp[i,j] = func.diff(vi,vj)
                    cp[j,i] = cp[i,j]
    return cp
################################################################################
def standard_vars(vars,means,stdevs):
    return np.array([(v-m)/s for v,m,s in zip(vars,means,stdevs)])
################################################################################
def standard_lc(lc,stdevs):
    return np.array([coef*stdev for coef,stdev in zip(lc,stdevs)])
################################################################################
def standard_qc(qc,stdevs):
    return np.array([coef*stdev**2 for coef,stdev in zip(qc,stdevs)])
################################################################################
def standard_cp(cp,stdevs):
    nvars = cp.shape[1]
    scp = np.empty_like(cp)
    if nvars >= 2:
        for i in range(nvars):
            for j in range(i+1,nvars):
                scp[i,j] = cp[i,j]*stdevs[i]*stdevs[j]
                scp[j,i] = scp[i,j]
    return scp
################################################################################
def standardize(lc,qc,cp,stdevs):#,means,vars):
    #w = standard_vars(vars,means,stdevs)
    slc = standard_lc(lc,stdevs)
    sqc = standard_qc(qc,stdevs)
    scp = standard_cp(cp,stdevs)
    #return w,slc,sqc,scp
    return slc,sqc,scp
###############################################################################
def makenumeric(lc,qc,cp,vars,nominals):
    n = len(vars)
    noms = dict(zip(vars,nominals))
    lc_n = np.array([lc[i].subs(noms) for i in range(n)])
    qc_n = np.array([qc[i].subs(noms) for i in range(n)])
    cp_n = np.empty_like(cp)
    for i in range(n):
        for j in range(i+1,n):
            cp_n[i,j] = cp[i,j].subs(noms)
            cp_n[j,i] = cp_n[i,j]
    return lc_n,qc_n,cp_n
###############################################################################
def rawmoment(slc,sqc,scp,vm,k):
    lc = copy(slc)
    qc = copy(sqc)
    cp = copy(scp)
    n = len(lc)
    
    if assume_linear:
        qc[:] = 0.0
        cp[:,:] = 0.0
        
    ans = 0.0
    if k==0:
        ans = 1
    ############################
    elif k==1:
        for i in range(n):
            ans += qc[i]*vm[i,2]
    ############################
    elif k==2:
        for i in range(n):
            ans += lc[i]**2*vm[i,2] + 2*lc[i]*qc[i]*vm[i,3] + qc[i]**2*vm[i,4]
        
        if n>=2:
            for i in range(n-1):
                for j in range(i+1,n):
                    ans += (2*qc[i]*qc[j] + cp[i,j]**2)*vm[i,2]*vm[j,2]
    ############################
    elif k==3:
        for i in range(n):
            ans += lc[i]**3*vm[i,3] + qc[i]**3*vm[i,6] + \
                   3*lc[i]**2*qc[i]*vm[i,4] + 3*lc[i]*qc[i]**2*vm[i,5]
        
        if n>=2:
            for i in range(n-1):
                for j in range(i+1,n):
                    ans += cp[i,j]**3*vm[i,3]*vm[j,3] + \
                           6*lc[i]*lc[j]*cp[i,j]*vm[i,2]*vm[j,2] + \
                           6*qc[i]*qc[j]*cp[i,j]*vm[i,3]*vm[j,3]
            for i in range(n):
                for j in range(n):
                    if j!=i:
                        ans += 3*qc[i]**2*vm[i,4]*qc[j]*vm[j,2] + \
                               6*lc[i]*qc[j]*cp[i,j]*vm[i,2]*vm[j,3] +\
                               3*qc[i]*lc[j]**2*vm[i,2]*vm[j,2] + \
                               6*lc[i]*qc[i]*qc[j]*vm[i,3]*vm[j,2] + \
                               3*lc[i]*cp[i,j]**2*vm[i,3]*vm[j,2] + \
                               3*qc[i]*cp[i,j]**2*vm[i,4]*vm[j,2]
        
        if n>=3:
            for i in range(n-2):
                for j in range(i+1,n-1):
                    for k in range(j+1,n):
                        ans += (6*qc[i]*qc[j]*qc[k] + \
                                6*cp[i,j]*cp[i,k]*cp[j,k] + 
                                3*(qc[i]*cp[j,k]**2 + \
                                qc[j]*cp[i,k]**2 + \
                                qc[k]*cp[i,j]**2))*vm[i,2]*vm[j,2]*vm[k,1]
    ############################
    elif k==4:
        for i in range(n):
            ans += lc[i]**4*vm[i,4] + qc[i]**4*vm[i,8] + \
                   4*lc[i]**3*qc[i]*vm[i,5] + 4*lc[i]*qc[i]**3*vm[i,7] + \
                   6*lc[i]**2*qc[i]**2*vm[i,6]
        
        if n>=2:
            for i in range(n-1):
                for j in range(i+1,n):
                    ans += 6*lc[i]**2*lc[j]**2*vm[i,2]*vm[j,2] + \
                           6*qc[i]**2*qc[j]**2*vm[i,4]*vm[j,4] + \
                           cp[i,j]**4*vm[i,4]*vm[j,4] + \
                           12*cp[i,j]*(lc[i]**2*lc[j]*vm[i,3]*vm[j,2] + \
                                       lc[i]*lc[j]**2*vm[i,2]*vm[j,3]) + \
                           12*cp[i,j]*qc[i]*qc[j]*(qc[i]*vm[i,5]*vm[j,3] + \
                                                   qc[j]*vm[i,3]*vm[j,6]) + \
                           12*qc[i]*qc[j]*(lc[i]**2*vm[i,4]*vm[j,2] + \
                                           lc[j]**2*vm[i,2]*vm[j,4] + \
                                           2*lc[i]*lc[j]*vm[i,3]*vm[j,3]) + \
                           6*cp[i,j]**2*(lc[i]**2*vm[i,4]*vm[j,2] + \
                                         lc[j]**2*vm[i,2]*vm[j,4] + \
                                         2*lc[i]*lc[j]*vm[i,3]*vm[j,3]) + \
                           6*cp[i,j]**2*(qc[i]**2*vm[i,6]*vm[j,2] + \
                                         qc[j]**2*vm[i,2]*vm[j,6] + \
                                         2*qc[i]*qc[j]*vm[i,4]*vm[j,4]) + \
                           12*cp[i,j]*(lc[j]*qc[i]*(lc[j]*vm[i,3]*vm[j,3] + \
                                                    2*lc[i]*vm[i,4]*vm[j,2]) + \
                                       lc[i]*qc[j]*(lc[i]*vm[i,3]*vm[j,3] + \
                                                    2*lc[j]*vm[i,2]*vm[j,4])) +\
                           12*cp[i,j]*(lc[i]*qc[j]*(qc[j]*vm[i,2]*vm[j,5] + \
                                                    2*qc[i]*vm[i,4]*vm[j,3]) + \
                                       lc[j]*qc[i]*(qc[i]*vm[i,5]*vm[j,2] + \
                                                    2*qc[j]*vm[i,3]*vm[j,4])) +\
                           12*cp[i,j]**2*(qc[i]*(lc[i]*vm[i,5]*vm[i,2] + \
                                                 lc[j]*vm[i,4]*vm[j,3]) + \
                                          qc[j]*(lc[i]*vm[i,3]*vm[j,4] + \
                                                 lc[j]*vm[i,2]*vm[j,5]))
            for i in range(n):
                for j in range(n):
                    if i!=j:
                        ans += 4*qc[i]**3*qc[j]*vm[i,6]*vm[j,2] + \
                               4*qc[i]*lc[j]**3*vm[i,2]*vm[j,3] + \
                               12*lc[i]*qc[i]*lc[j]**2*vm[i,3]*vm[i,2] + \
                               12*lc[i]*qc[i]**2*qc[j]*vm[i,5]*vm[i,2] + \
                               12*lc[i]*qc[i]*qc[j]**2*vm[i,3]*vm[j,4] + \
                               4*lc[i]*cp[i,j]**3*vm[i,4]*vm[j,3] + \
                               4*qc[i]*cp[i,j]**3*vm[i,5]*vm[j,3] + \
                               6*qc[i]**2*lc[j]**2*vm[i,4]*vm[j,2]
            
        if n>=3:
            for i in range(n-2):
                for j in range(i+1,n-1):
                    for k in range(j+1,n):
                        ans += (12*qc[i]**2*qc[j]*qc[k] + \
                               6*cp[i,j]**2*cp[i,k]**2 + \
                               12*qc[i]*(qc[k]*cp[i,j]**2 + qc[j]*cp[i,k]**2) + \
                               6*qc[i]**2*cp[j,k])*vm[i,4]*vm[j,2]*vm[k,2]
                        ans += (12*qc[i]*qc[j]**2*qc[k] + \
                               6*cp[i,j]**2*cp[j,k]**2 + \
                               12*qc[j]*(qc[k]*cp[i,j]**2 + qc[i]*cp[j,k]**2) + \
                               6*qc[j]**2*cp[i,k]**2)*vm[i,2]*vm[j,4]*vm[k,2]
                        ans += (12*qc[i]*qc[j]*qc[k]**2 + \
                               6*cp[i,k]**2*cp[j,k]**2 + \
                               12*qc[k]*(qc[i]*cp[j,k]**2 + qc[j]*cp[i,k]**2) + \
                               6*qc[k]**2*cp[i,j]**2)*vm[i,2]*vm[j,2]*vm[k,4]
                        ans += (12*cp[i,j]**2*cp[i,k]*cp[j,k] + \
                               24*qc[i]*qc[j]*qc[k]*cp[i,j] + \
                               4*qc[k]*cp[i,j]**3 + \
                               24*qc[i]*qc[k]*cp[i,k]*cp[j,k])*vm[i,3]*vm[j,3]*vm[k,2]
                        ans += (12*cp[i,j]*cp[i,k]**2*cp[j,k] + \
                               24*qc[i]*qc[j]*qc[k]*cp[i,k] + \
                               4*qc[j]*cp[i,k]**3 + \
                               24*qc[i]*qc[k]*cp[i,j]*cp[j,k])*vm[i,3]*vm[j,2]*vm[k,3]
                        ans += (12*cp[i,j]*cp[i,k]*cp[j,k]**2 + \
                               24*qc[i]*qc[j]*qc[k]*cp[j,k] + \
                               4*qc[i]*cp[j,k]**3 + \
                               24*qc[j]*qc[j]*cp[i,j]*cp[i,k])*vm[i,2]*vm[j,3]*vm[k,3]
                        ans += (12*cp[i,j]*cp[i,k]*cp[j,k] + \
                               24*qc[i]*qc[j]*qc[k]*cp[j,k] + \
                               4*qc[i]*cp[j,k]**3 + \
                               24*qc[j]*qc[k]*cp[i,j]*cp[i,k])*vm[i,2]*vm[j,3]*vm[k,3]
                        ans += 24*(qc[i]*qc[j]*qc[k] + \
                               cp[i,j]*cp[i,k]*cp[j,k])*(lc[i]*vm[i,3]*vm[j,2]*vm[k,2] + \
                                                         lc[j]*vm[i,2]*vm[j,3]*vm[k,2] + \
                                                         lc[k]*vm[i,2]*vm[j,2]*vm[k,3])
                        ans += 12*(lc[i]*cp[j,k]**2*vm[i,2]*(cp[i,j]*vm[j,3]*vm[k,2] + \
                               cp[i,k]*vm[j,2]*vm[k,3]) + lc[j]*cp[i,k]**2*vm[j,2]*(cp[i,j]*vm[i,3]*vm[k,2] + \
                               cp[j,k]*vm[i,2]*vm[k,3]) + \
                               lc[k]*cp[i,j]**2*vm[k,2]*(cp[i,k]*vm[i,3]*vm[j,2] + \
                               cp[j,k]*vm[i,2]*vm[j,3]))
                        ans += 12*(qc[i]*cp[j,k]**2*vm[i,3]*(cp[i,j]*vm[j,3]*vm[k,2] + \
                               cp[i,k]*vm[j,2]*vm[k,3]) + qc[j]*cp[i,k]**2*vm[j,3]*(cp[i,j]*vm[i,3]*vm[k,2] + \
                               cp[j,k]*vm[i,2]*vm[k,3]) + \
                               qc[k]*cp[i,j]**2*vm[k,3]*(cp[i,k]*vm[i,3]*vm[j,2] + \
                               cp[j,k]*vm[i,2]*vm[j,3]))
                        ans += 24*cp[i,j]*cp[i,k]*cp[j,k]*(qc[i]*vm[i,4]*vm[j,2]*vm[k,2] + \
                               qc[j]*vm[i,2]*vm[j,4]*vm[k,2] + qc[k]*vm[i,2]*vm[j,2]*vm[k,4])
                        ans += vm[i,2]*vm[j,2]*vm[k,2]*(12*(qc[i]*qc[j]*lc[k]**2 + \
                               qc[i]*qc[k]*lc[j]**2 + qc[j]*qc[k]*lc[i]**2) + \
                               6*(lc[i]**2*cp[j,k]**2 + lc[j]**2*cp[i,k]**2 + \
                               lc[k]**2*cp[i,j]**2) + \
                               24*(cp[i,j]*cp[i,k]*lc[j]*lc[k] + \
                                   cp[i,j]*cp[j,k]*lc[i]*lc[k] + \
                                   cp[i,k]*cp[j,k]*lc[i]*lc[j]) + \
                               24*(lc[i]*lc[j]*qc[k]*cp[i,j] + \
                                   lc[i]*lc[k]*qc[j]*cp[i,k] + \
                                   lc[j]*lc[k]*qc[i]*cp[j,k]))
                        ans += vm[i,3]*vm[j,2]*vm[k,2]*(24*lc[j]*cp[i,j]*qc[i]*qc[k] + \
                               24*lc[k]*cp[i,k]*qc[i]*qc[j] + \
                               12*lc[i]*cp[j,k]**2*qc[i] + \
                               24*lc[j]*cp[i,k]*cp[j,k]*qc[i] + \
                               24*lc[k]*cp[i,j]*cp[j,k]*qc[i] + \
                               12*lc[i]*cp[i,k]**2*qc[j] + \
                               12*lc[i]*cp[i,j]**2*qc[k])
                        ans += vm[i,2]*vm[j,3]*vm[k,2]*(24*lc[i]*cp[i,j]*qc[j]*qc[k] + \
                               24*lc[k]*cp[j,k]*qc[i]*qc[j] + \
                               12*lc[j]*cp[i,k]**2*qc[j] + \
                               24*lc[i]*cp[i,k]*cp[j,k]*qc[j] + \
                               24*lc[k]*cp[i,j]*cp[i,k]*qc[j] + \
                               12*lc[j]*cp[j,k]**2*qc[i] + \
                               12*lc[j]*cp[i,j]**2*qc[k])
                        ans += vm[i,2]*vm[j,2]*vm[k,3]*(24*lc[i]*cp[i,k]*qc[j]*qc[k] + \
                               24*lc[j]*cp[j,k]*qc[i]*qc[k] + \
                               12*lc[k]*cp[i,j]**2*qc[k] + \
                               24*lc[i]*cp[i,j]*cp[j,k]*qc[k] + \
                               24*lc[j]*cp[i,j]*cp[i,k]*qc[k] + \
                               12*lc[k]*cp[j,k]**2*qc[i] + \
                               12*lc[k]*cp[i,k]**2*qc[j])
        
        if n>=4:
            for i in range(n-3):
                for j in range(i+1,n-2):
                    for k in range(j+1,n-1):
                        for m in range(k+1,n):
                            ans += vm[i,2]*vm[j,2]*vm[k,2]*vm[m,2]*(24*(qc[i]*qc[j]*qc[k]*qc[m] +\
                                   cp[i,j]*cp[i,k]*cp[j,m]*cp[k,m] + \
                                   cp[i,j]*cp[i,m]*cp[j,k]*cp[k,m] + \
                                   cp[i,k]*cp[i,m]*cp[j,k]*cp[j,m] + \
                                   qc[i]*cp[j,k]*cp[j,m]*cp[k,m] + \
                                   qc[j]*cp[i,k]*cp[i,m]*cp[i,m] + \
                                   qc[k]*cp[i,j]*cp[i,m]*cp[j,m] + \
                                   qc[m]*cp[i,j]*cp[i,k]*cp[j,k]) + \
                                   12*(qc[i]*qc[j]*cp[k,m]**2 + \
                                       qc[i]*qc[k]*cp[j,m]**2 + \
                                       qc[i]*qc[m]*cp[j,k]**2 + \
                                       qc[j]*qc[k]*cp[i,m]**2 + \
                                       qc[j]*qc[m]*cp[i,k]**2 + \
                                       qc[k]*qc[m]*cp[i,j]**2) + \
                                   6*(cp[i,j]**2*cp[k,m]**2 + \
                                      cp[i,k]**2*cp[j,m]**2 + \
                                      cp[i,m]**2*cp[j,k]**2))
    else:
        print 'Can only calculate raw moments k = 0 to 4. Sorry.'
        ans = None
    
    return ans
    
################################################################################
def centralmoment(vi,k):
    if k==0:
        ans = 1
    elif k==1:
        ans = 0
    elif k==2:
        ans = vi[2]-vi[1]**2
    elif k==3:
        ans = vi[3]-3*vi[2]*vi[1] + 2*vi[1]**3
    elif k==4:
        ans = vi[4] - 4*vi[3]*vi[1] + 6*vi[2]*vi[1]**2 - 3*vi[1]**4
    else:
        print 'Can only calculate central moments k = 0 to 4. Sorry.'
        ans = None
    return ans
################################################################################
def variance_components(slc,sqc,scp,var_moments,vz):
    n = len(slc)
    var_comp_lc = np.empty(n)
    var_comp_qc = np.empty(n)
    var_comp_cp = np.empty((n,n))

    for i in range(n):
        slc_tmp = copy(slc)
        slc_tmp[i] = 0.0
        vy_tmp = [rawmoment(slc_tmp,sqc,scp,var_moments,k) for k in range(3)]
        var_comp_lc[i] = vz[2]-centralmoment(vy_tmp,2)
#        if not silent:
#            print 'Variance component of lc[{:s}]: {:7.5%}'.format(x[i],var_comp_lc[i])

    for i in range(n):
        sqc_tmp = copy(sqc)
        sqc_tmp[i] = 0.0
        vy_tmp = [rawmoment(slc,sqc_tmp,scp,var_moments,k) for k in range(3)]
        var_comp_qc[i] = vz[2]-centralmoment(vy_tmp,2)
#        if not silent:
#            print 'Variance component of qc[{:s}]: {:7.5%}'.format(x[i],var_comp_qc[i])

    for i in range(n-1):
        for j in range(i+1,n):
            scp_tmp = copy(scp)
            scp_tmp[i,j] = 0.0
            scp_tmp[j,i] = 0.0
            vy_tmp = [rawmoment(slc,sqc,scp_tmp,var_moments,k) for k in range(3)]
            var_comp_cp[i,j] = vz[2]-centralmoment(vy_tmp,2)
    #        if not silent:
    #            print 'Variance component of cp[{:s},{:s}]: {:7.5%}'.format(x[i],x[j],var_comp_cp[i,j])

    return var_comp_lc,var_comp_qc,var_comp_cp
################################################################################
def variance_contrib(var_comp_lc,var_comp_qc,var_comp_cp,vz):
    n = len(var_comp_lc)
    var_contrib_lc = np.empty_like(var_comp_lc)
    var_contrib_qc = np.empty_like(var_comp_qc)
    var_contrib_cp = np.empty_like(var_comp_cp)

    for i in range(n):
        if vz[2]:
            var_contrib_lc[i] = np.abs(var_comp_lc[i]/vz[2])
        else:
            var_contrib_lc[i] = 0.0
#        if not silent:
#            print 'Variance contribution of lc[{:s}]: {:7.5%}'.format(x[i],var_contrib_lc[i])

    for i in range(n):
        if vz[2]:
            var_contrib_qc[i] = np.abs(var_comp_qc[i]/vz[2])
        else:
            var_contrib_qc[i] = 0.0
#        if not silent:
#            print 'Variance contribution of qc[{:s}]: {:7.5%}'.format(x[i],var_contrib[i])

    for i in range(n-1):
        for j in range(i+1,n):
            if vz[2]:
                var_contrib_cp[i,j] = np.abs(var_comp_cp[i,j]/vz[2])
            else:
                var_contrib_cp[i,j] = 0.0
#            if not silent:
#                print 'Variance contribution of cp[{:s},{:s}]: {:7.5%}'.format(x[i],x[j],var_contrib[i])
        
    return var_contrib_lc,var_contrib_qc,var_contrib_cp
################################################################################
def soerp_symbolic(func,x,means,stdevs,var_moments,title=None,debug=False,
                   silent=False):
    means = np.array(means)
    stdevs = np.array(stdevs)
    if not silent:
        print '\n','*'*80
        if title:
            print '{:*^80}'.format(' SOERP: '+title+' ')
        else:
            print '*'*80
        print '*'*80

    func0 = func.subs(dict(zip(x,means)))
    if debug and not silent:
        print '*'*80
        print 'System Function :',func
        print 'System Variables:',x
        print 'Variable Means  :',means
        print 'Variable Stdevs :',stdevs
        print 'Variable Moments:\n',var_moments
        print 'System Nominal  :',func0
    ############################
    lc = linear_coef(func,x)
    if debug and not silent:
        print '*'*80
        print 'Linear coefficients:\n',lc
    ############################
    qc = quad_coef(func,x)
    if debug and not silent:
        print '*'*80
        print '1/2*Quadratic coefficients:\n',qc
    ############################
    cp = cross_coef(func,x)
    if debug and not silent:
        print '*'*80
        print 'Cross-product coefficients:\n',cp
    ############################
    lc,qc,cp = makenumeric(lc,qc,cp,x,means)
    if debug and not silent:
        print '*'*80
        print 'Numerical linear coefficients:\n',lc
        print 'Numerical quadratic coefficients:\n',qc
        print 'Numerical cross-product coefficients:\n',cp
    ############################
    slc,sqc,scp = standardize(lc,qc,cp,stdevs)
    
    if debug and not silent:
        print '*'*80
        print 'Standardized linear coefficients:\n',slc
        print 'Standardized quadratic coefficients:\n',sqc
        print 'Standardized cross-product coefficients:\n',scp
    ############################
    vy = np.empty(5)
    if debug and not silent:
        print '*'*80
    for k in range(5):
        vy[k] = rawmoment(slc,sqc,scp,var_moments,k)
        if debug and not silent:
            print 'Raw Moment',k,':',vy[k]
    ############################
    vz = np.empty(5)
    if debug and not silent:
        print '*'*80
    for k in range(5):
        vz[k] = centralmoment(vy,k)
        if debug and not silent:
            print 'Central Moment',k,':',vz[k]
    sysmean = float(vy[1]+func0)
    ############################
    # Calculate variance contributions
    vc_lc,vc_qc,vc_cp = variance_components(slc,sqc,scp,var_moments,vz)
    vlc,vqc,vcp = variance_contrib(vc_lc,vc_qc,vc_cp,vz)
    n = len(x)
#    var_contrib = np.empty(n)
    if not silent:
        print '*'*80
        for i in range(n):
#        slc_tmp = copy.copy(slc)
#        slc_tmp[i] = 0.0
#        vy_tmp = [rawmoment(slc_tmp,sqc,scp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        var_contrib[i] = np.abs(1-(vz_tmp)/vz[2])
            print 'Variance Contribution of lc[{:s}]: {:7.5%}'.format(x[i],vlc[i])

    if not silent:
        print '*'*80
        for i in range(n):
#        sqc_tmp = copy.copy(sqc)
#        sqc_tmp[i] = 0.0
#        vy_tmp = [rawmoment(slc,sqc_tmp,scp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        var_contrib[i] = np.abs(1-(vz_tmp)/vz[2])
            print 'Variance Contribution of qc[{:s}]: {:7.5%}'.format(x[i],vqc[i])

    if not silent:
        print '*'*80
        for i in range(n-1):
            for j in range(i+1,n):
#            scp_tmp = copy.copy(scp)
#            scp_tmp[i,j] = 0.0
#            scp_tmp[j,i] = 0.0
#            vy_tmp = [rawmoment(slc,sqc,scp_tmp,var_moments,k) for k in range(3)]
#            vz_tmp = centralmoment(vy_tmp,2)
#            var_contrib[i] = np.abs(1-(vz_tmp)/vz[2])
                print 'Variance Contribution of cp[{:s},{:s}]: {:7.5%}'.format(x[i],x[j],vcp[i,j])
        
#    if not silent:
#        print '*'*80
#    for i in range(n):
#        slc_tmp = copy.copy(slc)
#        sqc_tmp = copy.copy(sqc)
#        scp_tmp = copy.copy(scp)
#        slc_tmp[i] = 0.0
#        sqc_tmp[i] = 0.0
#        scp_tmp[i,:] = 0.0
#        scp_tmp[:,i] = 0.0
#        vy_tmp = [rawmoment(slc_tmp,sqc_tmp,scp_tmp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        var_contrib[i] = np.abs(1-(vz_tmp)/vz[2])
#        if not silent:
#            print 'Approx. Variance Contribution of {:s}: {:7.5%}'.format(x[i],var_contrib[i])
        
        
    ############################
    stdev = vz[2]**(0.5)
    rtbt1 = vz[3]/vz[2]**(1.5)
    beta1 = rtbt1**2
    beta2 = vz[4]/vz[2]**2
    if not silent:
        print '*'*80
        print 'MEAN-INTERCEPT (EDEL1).................... =','{: 8.7E}'.format(vy[1])
        print 'MEAN...................................... =','{: 8.7E}'.format(sysmean)
        print 'SECOND MOMENT (EDEL2)..................... =','{: 8.7E}'.format(vy[2])
        print 'VARIANCE (VARDL).......................... =','{: 8.7E}'.format(vz[2])
        print 'STANDARD DEVIATION.(RTVAR)................ =','{: 8.7E}'.format(stdev)
        print 'THIRD MOMENT (EDEL3)...................... =','{: 8.7E}'.format(vy[3])
        print 'THIRD CENTRAL MOMENT (MU3DL).............. =','{: 8.7E}'.format(vz[3])
        print 'COEFFICIENT OF SKEWNESS SQUARED (BETA1)... =','{: 8.7E}'.format(beta1)
        print 'COEFFICIENT OF SKEWNESS (RTBT1)........... =','{: 8.7E}'.format(rtbt1)
        print 'FOURTH MOMENT (EDEL4)..................... =','{: 8.7E}'.format(vy[4])
        print 'FOURTH CENTRAL MOMENT (MU4DL)............. =','{: 8.7E}'.format(vz[4])
        print 'COEFFICIENT OF KURTOSIS (BETA2)........... =','{: 8.7E}'.format(beta2)
        print '*'*80
    
    return [sysmean,vz[2],rtbt1,beta2]
    
################################################################################
def soerp_numeric(slc,sqc,scp,var_moments,func0,
                  title=None,debug=False,silent=False):
    """
    This performs the same moment calculations, but expects that all input
    derivatives and moments have been put in standardized form.
    
    Parameters
    ----------
    slc : array
        1st-order standardized derivatives (i.e., multiplied by the standard 
        deviation of the related input)
    sqc : array
        2nd-order derivatives (i.e., multiplied by the standard 
        deviation squared, or variance, of the related input)
    scp : 2d-array
        2nd-order cross-derivatives (i.e., multiplied by the two standard 
        deviations of the related inputs)
    var_moments : 2-d array
        Standardized moments where row[i] contains the first 9 moments of 
        variable x[i]. FYI: the first 3 values should always be [1, 0, 1]
    func0 : scalar
        System mean (i.e. value of the system evaluated at the means of all
        the input variables)
        
    Optional
    --------
    title : str
        Identifier for results that get printed to the screen
    debug : bool, false by default
        If true, all intermediate calculation results get printed to the screen
    silent : bool, false by default
        If true, nothing gets printed to the screen (overrides debug).
    
    Returns
    -------
    moments : list
        The first four standard moments (mean, variance, skewness and kurtosis
        coefficients)
    
    Example
    -------
    Example taken from the original SOERP user guide by N. D. Cox:
        >>> lc = [-802.65,-430.5]
        >>> qc = [205.54, 78.66]
        >>> cp = np.array([[0,-216.5],[-216.5,0]])
        >>> vm = np.array([norm_moments,norm_moments])
        >>> f0 = 4152
        >>> soerp_numeric(lc,qc,cp,vm,f0,
        ...               title='EXAMPLE FROM ORIGINAL SOERP USER GUIDE')
        ********************************************************************************
        **************** SOERP: EXAMPLE FROM ORIGINAL SOERP USER GUIDE *****************
        ********************************************************************************
        ********************************************************************************
        Variance Contribution of lc[x0]: 66.19083%
        Variance Contribution of lc[x1]: 19.04109%
        Variance Contribution of qc[x0]: 8.68097%
        Variance Contribution of qc[x1]: 1.27140%
        Variance Contribution of cp[x0,x1]: 4.81572%
        ********************************************************************************
        Approx. Variance Contribution of x0: 79.68751%
        Approx. Variance Contribution of x1: 25.12821%
        ********************************************************************************
        MEAN-INTERCEPT (EDEL1).................... =  2.8420000E+02
        MEAN...................................... =  4.4362000E+03
        SECOND MOMENT (EDEL2)..................... =  1.0540873E+06
        VARIANCE (VARDL).......................... =  9.7331770E+05
        STANDARD DEVIATION (RTVAR)................ =  9.8656865E+02
        THIRD MOMENT (EDEL3)...................... =  1.4392148E+09
        THIRD CENTRAL MOMENT (MU3DL).............. =  5.8640938E+08
        COEFFICIENT OF SKEWNESS SQUARED (BETA1)... =  3.7293913E-01
        COEFFICIENT OF SKEWNESS (RTBT1)........... =  6.1068742E-01
        FOURTH MOMENT (EDEL4)..................... =  5.0404781E+12
        FOURTH CENTRAL MOMENT (MU4DL)............. =  3.8956371E+12
        COEFFICIENT OF KURTOSIS (BETA2)........... =  4.1121529E+00
        ********************************************************************************

    """
    if not silent:
        print '\n','*'*80
        if title:
            print '{:*^80}'.format(' SOERP: '+title+' ')
        else:
            print '*'*80
        print '*'*80
    ############################
    vy = np.empty(5)
    if debug and not silent:
        print '*'*80
    for k in range(5):
        vy[k] = rawmoment(slc,sqc,scp,var_moments,k)
        if debug and not silent:
            print 'Raw Moment',k,':',vy[k]
    ############################
    vz = np.empty(5)
    if debug and not silent:
        print '*'*80
    for k in range(5):
        vz[k] = centralmoment(vy,k)
        if debug and not silent:
            print 'Central Moment',k,':',vz[k]
    sysmean = float(vy[1]+func0)
    ############################
    # Calculate variance contributions
    vc_lc,vc_qc,vc_cp = variance_components(slc,sqc,scp,var_moments,vz)
    vlc,vqc,vcp = variance_contrib(vc_lc,vc_qc,vc_cp,vz)
    n = len(slc)

    if not silent:
        print '*'*80
        for i in range(n):
#        slc_tmp = copy.copy(slc)
#        slc_tmp[i] = 0.0
#        vy_tmp = [rawmoment(slc_tmp,sqc,scp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        if vz[2]:
#            var_contrib = np.abs(1-(vz_tmp)/vz[2])
#        else:
#            var_contrib = 0.0
#        if not silent:
            print 'Variance Contribution of lc[x{:d}]: {:7.5%}'.format(i,vlc[i])

        for i in range(n):
#        sqc_tmp = copy.copy(sqc)
#        sqc_tmp[i] = 0.0
#        vy_tmp = [rawmoment(slc,sqc_tmp,scp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        if vz[2]:
#            var_contrib = np.abs(1-(vz_tmp)/vz[2])
#        else:
#            var_contrib = 0.0
#        if not silent:
            print 'Variance Contribution of qc[x{:d}]: {:7.5%}'.format(i,vqc[i])

        for i in range(n-1):
            for j in range(i+1,n):
#            scp_tmp = copy.copy(scp)
#            scp_tmp[i,j] = 0.0
#            scp_tmp[j,i] = 0.0
#            vy_tmp = [rawmoment(slc,sqc,scp_tmp,var_moments,k) for k in range(3)]
#            vz_tmp = centralmoment(vy_tmp,2)
#            if vz[2]:
#                var_contrib = np.abs(1-(vz_tmp)/vz[2])
#            else:
#                var_contrib = 0.0
#            if not silent:
                print 'Variance Contribution of cp[x{:d},x{:d}]: {:7.5%}'.format(i,j,vcp[i,j])
        
#    if not silent:
#        print '*'*80
#    for i in range(n):
#        slc_tmp = copy.copy(slc)
#        sqc_tmp = copy.copy(sqc)
#        scp_tmp = copy.copy(scp)
#        slc_tmp[i] = 0.0
#        sqc_tmp[i] = 0.0
#        scp_tmp[i,:] = 0.0
#        scp_tmp[:,i] = 0.0
#        vy_tmp = [rawmoment(slc_tmp,sqc_tmp,scp_tmp,var_moments,k) for k in range(3)]
#        vz_tmp = centralmoment(vy_tmp,2)
#        if vz[2]:
#            var_contrib[i] = np.abs(1-(vz_tmp)/vz[2])
#        else:
#            var_contrib[i] = 0.0
#        if not silent:
#            print 'Approx. Variance Contribution of x{:d}: {:7.5%}'.format(i,var_contrib[i])
        
    ############################
    stdev = vz[2]**(0.5)
    if stdev:
        rtbt1 = vz[3]/vz[2]**(1.5)
        beta2 = vz[4]/vz[2]**2
    else:
        rtbt1 = 0.0
        beta2 = 0.0
    beta1 = rtbt1**2
    if not silent:
        print '*'*80
        print 'MEAN-INTERCEPT (EDEL1).................... =','{: 8.7E}'.format(vy[1])
        print 'MEAN...................................... =','{: 8.7E}'.format(sysmean)
        print 'SECOND MOMENT (EDEL2)..................... =','{: 8.7E}'.format(vy[2])
        print 'VARIANCE (VARDL).......................... =','{: 8.7E}'.format(vz[2])
        print 'STANDARD DEVIATION (RTVAR)................ =','{: 8.7E}'.format(stdev)
        print 'THIRD MOMENT (EDEL3)...................... =','{: 8.7E}'.format(vy[3])
        print 'THIRD CENTRAL MOMENT (MU3DL).............. =','{: 8.7E}'.format(vz[3])
        print 'COEFFICIENT OF SKEWNESS SQUARED (BETA1)... =','{: 8.7E}'.format(beta1)
        print 'COEFFICIENT OF SKEWNESS (RTBT1)........... =','{: 8.7E}'.format(rtbt1)
        print 'FOURTH MOMENT (EDEL4)..................... =','{: 8.7E}'.format(vy[4])
        print 'FOURTH CENTRAL MOMENT (MU4DL)............. =','{: 8.7E}'.format(vz[4])
        print 'COEFFICIENT OF KURTOSIS (BETA2)........... =','{: 8.7E}'.format(beta2)
        print '*'*80
    
    return [sysmean,vz[2],rtbt1,beta2]

###############################################################################

if __name__=='__main__':
    import sympy as sym
    
    x1,x2,x3 = sym.symbols('x1,x2,x3')
    Z = (x1*x2**2)/(15*(1.5+x3))
    
    # standardized moments of a normal distribution: 
    norm_moments = [1,0,1,0,3,0,15,0,105]
    # standardized moments of a uniform distribution: 
    unif_moments = [1,0,1,0,1.8,0,3.857,0,9]
    # standardized moments of an exponential distribution: 
    expn_moments = [1,0,1,2,9,44,265,1854,14833]
    
    x_nom = [24, 37, 0.5]
    x_std = [1, 4, 0.5]
    x_vm = np.array([norm_moments,
                     norm_moments,
                     expn_moments]) 
                     
    
    soerp_symbolic(Z,[x1,x2,x3],x_nom,x_std,x_vm,title='THREE PART ASSEMBLY')
    
    ########################
    C,H,M,P,t = sym.symbols('C,H,M,P,t')
    Q = C*sym.sqrt((520.0*H*P)/(M*(t+460.0)))
    Q = Q.subs({C:38.4})
    
    x_nom = [64.0,16.0,361.0,165.0]
    x_std = [0.5,0.1,2.0,0.5]
    x_vm = np.array([norm_moments,
                     norm_moments,
                     norm_moments,
                     norm_moments])
    
    soerp_symbolic(Q,[H,M,P,t],x_nom,x_std,x_vm,title='VOLUMETRIC GAS FLOW THROUGH ORIFICE METER')
    
    ########################
    x = sym.Symbol('x')
    
    soerp_symbolic(x**2,[x],[0],[1],np.array([norm_moments]),title='SIMPLE QUADRATIC: X^2',debug=True)
    ########################
    x = sym.Symbol('x')
    
    soerp_symbolic(x*sym.sin(x)+1,[x],[0],[1],np.array([norm_moments]),title='X*SIN(X)+1')
    
    ########################
    V,Ra,Rb,Rc = sym.symbols('V,Ra,Rb,Rc')
    I = V*(1/Ra+1/Rb+1/Rc)
    
    x_nom = [120.0,2.98,15.0,20.0]
    x_std = np.sqrt([15.0,1.0,1.0,2.0])
    x_vm = np.array([[1,0,1,0,1.8,0,0,0,0],
                     [1,0,1,1.044,5,0,0,0,0],
                     [1,0,1,-1,3,0,0,0,0],
                     [1,0,1,2,9,0,0,0,0]])
    
    soerp_symbolic(I,[V,Ra,Rb,Rc],x_nom,x_std,x_vm,title='ELECTRICAL CIRCUIT')
    
    ########################
    S0,S1,t1,f2,I,t2,t3,f4,f5,D,r = sym.symbols('S0,S1,t1,f2,I,t2,t3,f4,f5,D,r')
    W0 = S0
    W1 = S1*(1-sym.exp(-r*t1))/(r*t1)
    S = (f2*D*I/t2*r**2)*(-1-r*t2+sym.exp(r*t2))
    W2 = S*sym.exp(-r*(t2+t1))
    W3 = (f2*D*I/r)*(-1+sym.exp(r*t3))*(sym.exp(-r*(t3+t2+t1)))
    W4 = (f4*S1/r)*(-1+sym.exp(r*(t3+t2)))*(sym.exp(-r*(t3+t2+t1)))
    W5 = (1-f5)**(t3+t2)*S1*sym.exp(-r*(t3+t2+t1))
    W = -W0 - W1 + W2 + W3 - W4 + W5
    W = W.subs({D:182.5e9,r:0.2})
    
    x_nom = [100000,1e9,5,0.9,4/1000.,5,15,0.05,0.10]
    x_tol = np.array([100000*0.2,1e9*0.25,1,0.10,1/1000.,1,10,0.01,0.05])
    x_var = (2*x_tol)**2/12
    x_std = np.sqrt(x_var)
    x_vm = np.array([unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments,
                     unif_moments])
    noms = dict(zip([S0,S1,t1,f2,I,t2,t3,f4,f5],x_nom))
    #print '-W0 =',-W0.subs(noms).subs({D:182.5e9,r:0.2})
    #print '-W1 =',-W1.subs(noms).subs({D:182.5e9,r:0.2})
    #print ' W2 =',W2.subs(noms).subs({D:182.5e9,r:0.2})
    #print ' W3 =',W3.subs(noms).subs({D:182.5e9,r:0.2})
    #print '-W4 =',-W4.subs(noms).subs({D:182.5e9,r:0.2})
    #print ' W5 =',W5.subs(noms).subs({D:182.5e9,r:0.2})
    
    #soerp_symbolic(W,[S0,S1,t1,f2,I,t2,t3,f4,f5],x_nom,x_std,x_vm,title='PROFITABILITY ANALYSIS')
    
    
    ########################
    lc = [-802.65,-430.5]
    qc = [205.54, 78.66]
    cp = np.array([[0,-216.5],[-216.5,0]])
    vm = np.array([norm_moments,norm_moments])
    f0 = 4152
    soerp_numeric(lc,qc,cp,vm,f0,title='EXAMPLE FROM ORIGINAL SOERP USER GUIDE')

    ########################
    H,B,d,t,E,rho,P = sym.symbols('H,B,d,t,E,rho,P')
    wght = 2*sym.pi*rho*d*t*sym.sqrt((B/2)**2+H**2)
    strs = (P*sym.sqrt((B/2)**2+H**2))/(2*sym.pi*d*t*H)
    buck = (sym.pi**2*E*(d**2+t**2))/(8*((B/2)**2+H**2))
    defl = (P*((B/2)**2+H**2)**(1.5))/(2*sym.pi*d*t*H**2*E)
    
    x_nom = [30,60,3,0.15,30000,0.3,66]
    x_std = np.array([5,0.5,0.1,0.01,1500,0.01,3])/3
    x_vm = np.array([norm_moments,
                     norm_moments,
                     norm_moments,
                     norm_moments,
                     norm_moments,
                     norm_moments,
                     norm_moments])
    
    soerp_symbolic(wght,[H,B,d,t,E,rho,P],x_nom,x_std,x_vm,title='TWO-BAR TRUSS (WEIGHT)')
    soerp_symbolic(strs,[H,B,d,t,E,rho,P],x_nom,x_std,x_vm,title='TWO-BAR TRUSS (STRESS)')
    soerp_symbolic(buck,[H,B,d,t,E,rho,P],x_nom,x_std,x_vm,title='TWO-BAR TRUSS (BUCKLING)')
    soerp_symbolic(defl,[H,B,d,t,E,rho,P],x_nom,x_std,x_vm,title='TWO-BAR TRUSS (DEFLECTION)')
    











