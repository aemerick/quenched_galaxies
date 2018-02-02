import numpy as _np
from letstalkaboutquench import fstarforms as _fstarforms # fitting routine for SFMS

def _bradford_fgas_limit(mstar):
    """
    Limits pulled from digitized plot
    """
    data = _np.genfromtxt('./Data/bradford_fgas_limit.dat',names=True)

    return _np.interp(mstar, data['log_Mstar'], data['fgas'])


# generalize this a bit to be able to hand
# a data set and just cut based on keyword arguements
# for now, don't do this
def fgas_limits(log_Mstar, fgas, obs_data = 'Bradford'):
    """
    Given stellar mass and gas fraction of simulated dataset
    returns the array that would cut out galaxies according
    to the desired observational limits.
    """

    if obs_data  == 'Bradford':
        fgas_cut = lambda x : _bradford_fgas_limit(x)

    selection    = fgas >= fgas_cut(log_Mstar)

    return selection

def fit_SFMS(log_mstar, log_SFR, *args, **kwargs):
    """
    Convenience wrapper on Chang's routine
    """

    fit = _fstarforms.fstarforms()

    # set default behavior unless overridden
    if not 'fit_range' in kwargs.keys():
        min_bin = _np.floor(_np.min(log_mstar))
        max_bin = _np.ceil(_np.max(log_mstar))  
        kwargs['fit_range'] = [min_bin, max_bin]

    if not 'method' in kwargs.keys():
        kwargs['method'] = 'negbinomfit'
        kwargs['method'] = 'gaussmix'

    if not 'Nbin_thresh' in kwargs.keys():
        kwargs['Nbin_thresh'] = 20

    if not 'dlogm' in kwargs.keys():
        kwargs['dlogm'] = 0.1

    m_fit, sfr_fit = fit.fit(log_mstar, log_SFR, *args, **kwargs)

    # return an array of distance to SFMS
    D_SFMS = _np.ones(_np.size(log_mstar)) * -99 # init to flag

    D_SFMS[ log_SFR > -99] = log_SFR[ log_SFR > -99] - _np.interp(log_mstar[ log_SFR > -99], m_fit, sfr_fit)


    return D_SFMS, m_fit, sfr_fit
    
        
        



    
    
    
