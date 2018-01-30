import numpy as _np


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
