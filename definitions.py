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
_data = _home + "data/"
#eagle_ms="EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt"
#eagle_

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
data = {} # to hold actual data
for k in data_files.keys():
    data_files[k] = _data + data_files[k] # append full path

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH','ic_subfind'],
                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8','u1']})
data['illustris_extended'] = _np.genfromtxt(data_files['illustris_extended'], delimiter=',', skip_header = 1, dtype =_dt)

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH'],
                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8']})
data['illustris_extended']          = _np.genfromtxt(data_files['illustris_extended'],delimiter=',',skip_header=1,dtype=_dt)

_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_MBH'],
                'formats': ['f8','f8','f8','f8','f8','f8']})
data['illustris']          = _np.genfromtxt(data_files['illustris'],delimiter=',',skip_header=1,dtype=_dt)

for x in [k for k in data_files.keys() if 'SAM' in k]:
    data[x] = _np.genfromtxt(data_files[x], names = True)

# load MUFASA and remove satellites
data['MUFASA'] = _np.genfromtxt(data_files['MUFASA'],names=True)
data['MUFASA'] = data['MUFASA'][ data['MUFASA']['cen_sat'] == 1]

# brooks
data['Brooks'] = _np.genfromtxt(data_files['Brooks'], names = True)

#
# Maybe pre-compute gas fractions and things here for all the data
#
