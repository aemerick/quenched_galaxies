"""
Pre-load file paths and data
"""
import numpy as _np

import matplotlib as _mpl
_mpl.use('Agg')

from matplotlib import rc as _rc

fsize = 17
_rc('text', usetex=False)
_rc('font', size=fsize)#, ftype=42)
line_width = 3.5
point_size = 30

_home = "/home/emerick/work/quenched_galaxies/"
_data_path = _home + "data/"
#eagle_ms="EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt"
#eagle_

colors = {'illustris' : 'C0', 'SAM' : 'C1', 'MUFASA' : 'C2', 'Brooks' : 'C3', 'EAGLE' : 'C4'}

data_files = {'illustris_extended': "Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv",
              'illustris_all'     : "Illustris1_extended_individual_galaxy_values_z0.csv",
              'illustris'         : "Illustris1_individual_galaxy_values_z0.csv",
              'SAM'               : "SAM_group_catalog.dat",
              'SAM_cen_and_sat'   : "SAM_group_catalog_cenANDsat.dat",
              'Brooks'            : "brooks_updated_catalog.dat",
              'MUFASA'            : "MUFASA_GALAXY.txt"}

#
# Yes, this is gross
#
_data = {} # to hold actual data
for k in data_files.keys():
    data_files[k] = _data_path + data_files[k] # append full path

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH','ic_subfind'],
                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8','u1']})
_data['illustris_extended'] = _np.genfromtxt(data_files['illustris_extended'], delimiter=',', skip_header = 1, dtype =_dt)

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH'],
                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8']})
_data['illustris_extended']          = _np.genfromtxt(data_files['illustris_extended'],delimiter=',',skip_header=1,dtype=_dt)

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_MBH'],
                'formats': ['f8','f8','f8','f8','f8','f8']})
_data['illustris']          = _np.genfromtxt(data_files['illustris'],delimiter=',',skip_header=1,dtype=_dt)

for x in [k for k in data_files.keys() if 'SAM' in k]:
    _data[x] = _np.genfromtxt(data_files[x], names = True)

# load MUFASA and remove satellites
_data['MUFASA'] = _np.genfromtxt(data_files['MUFASA'],names=True)
_data['MUFASA'] = _data['MUFASA'][ _data['MUFASA']['cen_sat'] == 1]

# brooks
_data['Brooks'] = _np.genfromtxt(data_files['Brooks'], names = True)

#
# Maybe pre-compute gas fractions and things here for all the data
#

#
data = {} # now convert ndarrays of data into dictionaries for each data set
for k in _data.keys():
    data[k] = {}
    for l in _data[k].dtype.names:
        data[k][l] = _data[k][l]

#
# Pre-compute some things for each data set
#
def _compute_fgas(mstar, mgas, log = True):
    if log:
        fgas = 10**(mgas) / (10**(mstar) + 10**(mgas))
    else:
        fgas = mgas / (mstar + mgas)
    return fgas
data['illustris']['fgas'] = _compute_fgas(data['illustris']['log_Mstar'],data['illustris']['log_MHI'])
data['SAM']['fgas'] = _compute_fgas(data['SAM']['mstar'],data['SAM']['mcold'])
data['MUFASA']['fgas'] = _compute_fgas(data['MUFASA']['log_Mstar'],data['MUFASA']['log_Mcold'])
data['Brooks']['fgas'] = _compute_fgas(data['Brooks']['Mstar'],data['Brooks']['HI_Mass'],log=False)

#
#
for k in ['_10Myr','_1Gyr']:
    data['MUFASA']['SFR' + k] = 10.0**(data['MUFASA']['log_SFR' + k])
    data['MUFASA']['SFR' + k][ data['MUFASA']['log_SFR' + k] == -99.0 ] = 0.0 # flag for zero
    data['MUFASA']['sSFR' + k] = data['MUFASA']['SFR' + k] / 10.0**(data['MUFASA']['log_Mstar'])

# del(_data) - take off memory since we don't need this anymore

