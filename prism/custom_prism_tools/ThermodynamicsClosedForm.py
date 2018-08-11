import pyPRISM
import numpy as np
from itertools import combinations_with_replacement
from copy import deepcopy
import inspect
import warnings

def ThermodynamicsClosedForm(PRISM, default_method='HNC'): 
    
    r'''Calculate the Real-space *inter*-molecular pair correlation function 
    Parameters
    ----------
    PRISM: pyPRISM.core.PRISM
        A **solved** PRISM object.
    
    Returns
    -------
    pairCorr: pyPRISM.core.MatrixArray
        The full MatrixArray of pair correlation functions.
    
    **Mathematical Definition**
    .. math::
         
        g_{\alpha,\beta}(r) = h_{\alpha,\beta}(r) + 1.0
    
    **Variable Definitions**
        - :math:`g_{\alpha,\beta}(r)`
            Pair correlation function between site types :math:`\alpha` and
            :math:`\beta` at a distance :math:`r`
        - :math:`h_{\alpha,\beta}(r)`
            Total correlation function between site types :math:`\alpha` and
            :math:`\beta` at a distance :math:`r`
    **Description**
        The pair correlation function describes the spatial correlations
        between pairs of sites in Real-space. Also known as the *radial
        distribution function* (rdf), the :math:`g(r)` function is
        related to the underlying spatial probability distributions of a given
        system. In a PRISM calculation, :math:`g(r)` is strictly an
        *inter*-molecular quantity.
        After convergence of a PRISM object, the stored total correlation
        attribute function can simply be shifted to obtain the :math:`g(r)` 
    .. warning::
        Passing an unsolved PRISM object to this function will still produce
        output based on the default values of the attributes of the PRISM
        object.
    
    Example
    -------
    .. code-block:: python
        import pyPRISM
        sys = pyPRISM.System(['A','B'])
        
        # ** populate system variables **
        
        PRISM = sys.createPRISM()
        PRISM.solve()
        rdf = pyPRISM.calculate.pair_correlation(PRISM)
        rdf_AA = rdf['A','A']
        rdf_AB = rdf['A','B']
        rdf_BB = rdf['B','B']
    
    '''
    
    #Check on the closures and raise a warning of incompatible.
    #These formulas can be used anyways but they are not theoretically 
    #consistent with the input structure anymore.
    closures = set()
    for type_1, type_2 in combinations_with_replacement(PRISM.sys.types, 2):
        closures.add(dict(inspect.getmembers(PRISM.sys.closure[type_1, type_2]))['__module__'])
    
    if not (closures - set(['custom_prism_tools.KovalenkoHirata'])):
        method = 'KH'
    elif not (closures - set(['pyPRISM.closure.HyperNettedChain'])):
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
    Rho = np.diag(np.diag(PRISM.sys.density.site.data[0]))  #this is a diagonal density matrix
    
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
    #two quantites of relevance to the thermodynamics in k-space
    K1 = Wk.dot(Ck)
    K2 = (I-K1)
    #work with numpy arrays from here on out
    K1 = K1.data                 #WkCk
    K2 = K2.data                 #I-WkCk
    K2inv = np.linalg.inv(K2)    #(I-WkCk)^-1
    
    #specific quantities used in the k-space integrals
    tr_K1 = np.trace(K1, axis1=1, axis2=2)
    tr_K2invK1 = np.trace(np.matmul(K2inv, K1), axis1=1, axis2=2)
    lndet_K2 = np.log(np.linalg.det(K2))
    
    #free energy and pressure k-space contributions
    prefactor_k = 4.0*np.pi/(2.0*(2.0*np.pi)**3)
    F_ex_k = prefactor_k*np.trapz(k*k*(tr_K1 + lndet_K2), k)
    P_ex_k = -prefactor_k*np.trapz(k*k*(tr_K2invK1 + lndet_K2), k)
    
    #calculate the r-space contribution
    if method == 'HNC':
        Structure = ((1.0/2.0)*(Hr*Hr).data*np.heaviside(-Hr.data, 0) - Cr.data)
    elif method == 'KH':
        Structure = ((1.0/2.0)*(Hr*Hr).data - Cr.data)
    kernel_r = np.sum(np.matmul(np.matmul(Rho, Structure), Rho), axis=(1,2))
    F_and_P_ex_r = (4.0*np.pi/2.0)*np.trapz(r*r*kernel_r, r)
    
    #total free energy per volume and kBT
    F_ex = F_and_P_ex_r + F_ex_k
    
    #pressure per kBT (NOT divided by density)
    P_ex = F_and_P_ex_r + P_ex_k
    
    return (F_ex, P_ex)