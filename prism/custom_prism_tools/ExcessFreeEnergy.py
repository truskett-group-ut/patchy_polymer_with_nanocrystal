import pyPRISM
import numpy as np
from itertools import combinations_with_replacement
from copy import deepcopy
import warnings

def ExcessFreeEnergy(PRISM, default_method='HNC'): 
    
    #check on the closures and raise a warning of incompatible
    closures = set()
    for type_1, type_2 in combinations_with_replacement(PRISM.sys.types, 2):
        closures.add(type(PRISM.sys.closure[type_1, type_2]))
    
    if not (closures - set([type(KovalenkoHirata())])):
        method = 'KH'
    elif not (closures - set([type(pyPRISM.closure.HyperNettedChain())])):
        method = 'HNC'
    else:
        warnings.warn('This calculation is only theoretically valid if either \\
                      the KovalenkoHirata (KH) or HyperNettedChain (HNC) \\
                      closures are used for all pairs. Defaulting to {} \\
                      functional form for use with mixed closures.'.format(default_method), UserWarning)
        method = default_method
    
    #extra quantities for k-space calculation
    k = PRISM.sys.domain.k
    I = pyPRISM.IdentityMatrixArray(length=len(PRISM.sys.domain.k), rank=3, 
                                    types=PRISM.sys.types, space=pyPRISM.Space.Fourier)
    
    #extra quantities for r-space calculation
    r = PRISM.sys.domain.r
    p = np.diag(np.diag(PRISM.sys.density.site.data[0]))  #this is a diagonal density matrix
    
    #direct correlation function in r- and k-space
    if PRISM.directCorr.space == pyPRISM.Space.Real:
        Cr = deepcopy(PRISM.directCorr)
        Ck = deepcopy(Cr)
        PRISM.sys.domain.MatrixArray_to_fourier(Ck)
    elif PRISM.directCorr.space == pyPRISM.Space.Fourier:
        Ck = deepcopy(PRISM.directCorr)
        Cr = deepcopy(Ck)
        PRISM.sys.domain.MatrixArray_to_real(Cr)
    
    #intramolecular structure in k-space
    Wk = deepcopy(PRISM.omega)
    if Wk.space == pyPRISM.Space.Real:
        PRISM.sys.domain.MatrixArray_to_fourier(Wk)
    
    #total correlation function in r-space
    Hr = deepcopy(PRISM.totalCorr)
    if Hr.space == pyPRISM.Space.Fourier:
        PRISM.sys.domain.MatrixArray_to_real(Hr)
    
    #calculate the k-space contribution
    tr_WkCk = np.trace(Wk.dot(Ck).data, axis1=1, axis2=2)
    lndet_WkCk = np.log( np.linalg.det( (I-Wk.dot(Ck)).data ) )
    F_ex_k = 4.0*np.pi/(2.0*(2.0*np.pi)**3)*np.trapz(k*k*(tr_WkCk + lndet_WkCk), k)
    
    #calculate the r-space contribution
    if method == 'KH':
        structure = ((1.0/2.0)*(Hr*Hr).data*np.heaviside(-Hr.data, 0) - Cr.data)
    elif method == 'HNC':
        structure = ((1.0/2.0)*(Hr*Hr).data - Cr.data)
    kernel_r = np.sum(np.matmul(np.matmul(p, structure), p), axis=(1,2))
    F_ex_r = (4.0*np.pi/2.0)*np.trapz(r*r*kernel_r, r)
    
    #total free energy per unit volume and kBT
    F_ex = F_ex_r + F_ex_k
    
    return F_ex