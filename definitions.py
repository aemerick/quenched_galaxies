"""
Pre-load file paths and data
"""
import numpy as _np

import matplotlib as _mpl
_mpl.use('Agg')

from matplotlib import rc as _rc
from astropy.io import fits

fsize = 17
_rc('text', usetex=False)
_rc('font', size=fsize)#, ftype=42)
line_width = 3.5
point_size = 30

_home = "/home/aemerick/code/quenched_galaxies/"
_data_path = _home + "Data/"


#eagle_ms="EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt"
#eagle_

colors = {'illustris' : 'C0', 'SAM' : 'C1', 'MUFASA' : 'C2', 'Brooks' : 'C3', 'EAGLE' : 'C4', 'Bradford2015' : 'C5'}

data_files = { # 'illustris_extended': "Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv",
               # 'illustris_all'     : "Illustris1_extended_individual_galaxy_values_z0.csv",
              'illustris'         : "Illustris1_individual_galaxy_values_z0.csv",
              'SAM'               : 'sc_sam_cat_day2.txt', #"SAM_group_catalog.dat",
              #'SAM_cen_and_sat'   : "SAM_group_catalog_cenANDsat.dat",
              #'Brooks'            : "brooks_data_day2.dat",
              'MUFASA'            : "MUFASA_GALAXY.txt",
              'EAGLE'             : ['EAGLE_RefL0100_Mcoldgas_allabove1.8e8Msun.txt', 
                                     'EAGLE_RefL0100_Mstarin12HalfMassRadStars_allabove1.8e8Msun.txt',
                                     'EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt']}

data_dtypes = {'EAGLE'     : [ _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mcold'], 'formats' : ['f8','f8','f8']}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'StellarRhalf', 'log_Mstar', 'log_Mstar1Rhalf', 'log_Mstar2Rhalf'], 'formats':['f8']*6}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mstar', 'SFR_10Myr', 'SFR_1Gyr'], 'formats' : ['f8','f8','f8','f8','f8']})],
               'illustris' :   _np.dtype( {'names' : ['log_Mstar','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_MBH'], 'formats': ['f8','f8','f8','f8','f8','f8']}),
               'SAM'       :   _np.dtype( {'names' : ['log_Mstar','log_sSFR_100kyr','log_sSFR_100Myr','log_Mcold','sigma_8','log_MBH'], 'formats': ['f8']*6}),
               'MUFASA'    :   _np.dtype( {'names' : ['x','y','z','vx','vy','vz','log_Mstar','log_SFR_10Myr','log_SFR_1Gyr','log_Mcold','log_Z_SFR','cen_sat'], 'formats': ['f8']*12 + ['u1']})
              }

delimiters = {} # assign delimiters and exceptions
for k in data_dtypes.keys():
    delimiters[k] = None
delimiters['illustris']  = ','

skip_headers = {} # assign skip headers and exceptions
for k in data_dtypes.keys():
    skip_headers[k] = 1
skip_headers['MUFASA'] = 13

#
# This is maybe gross, but performs things consistently for all data 
#
_data = {} # to hold actual data
for k in data_files.keys():

    if not isinstance(data_files[k], basestring): # if multiple files
        _data[k] = [None]*len(data_files[k])
        for i in _np.arange(len(data_files[k])):
            data_path = _data_path + data_files[k][i]
            _data[k][i]  = _np.genfromtxt( data_path, delimiter = delimiters[k], skip_header = skip_headers[k], dtype = data_dtypes[k][i])
    else: # single data file to load
        data_path = _data_path + data_files[k]
        _data[k]  = _np.genfromtxt( data_path, delimiter = delimiters[k], skip_header = skip_headers[k], dtype = data_dtypes[k])

# Filter out some things for specific data sets before making into a dictionary
_data['MUFASA'] = _data['MUFASA'][ _data['MUFASA']['cen_sat'] == 1]


#
# now for each case, stitch together the files using
#
data = {}
for k in data_files.keys():
    data[k] = {}
    if not isinstance(data_files[k], basestring): # if multiple files
        for i in _np.arange(len(data_files[k])):
            for l in _data[k][i].dtype.names:
                data[k][l] = _data[k][i][l]        
    else: # single data file to load
        for l in _data[k].dtype.names:
            data[k][l] = _data[k][l]


#_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH','ic_subfind'],
#                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8','u1']})
#_data['illustris_extended'] = _np.genfromtxt(data_files['illustris_extended'], delimiter=',', skip_header = 1, dtype =_dt)

#_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI', 'Z_SF','log_MBH'],
#                'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8']})
#_data['illustris_extended']          = _np.genfromtxt(data_files['illustris_extended'],delimiter=',',skip_header=1,dtype=_dt)

#_dt = _np.dtype( {'names' : ['log_Mstar','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_MBH'],
#                'formats': ['f8','f8','f8','f8','f8','f8']})
#_data['illustris']          = _np.genfromtxt(data_files['illustris'],delimiter=',',skip_header=1,dtype=_dt)

#for x in [k for k in data_files.keys() if 'SAM' in k]:
#    _data[x] = _np.genfromtxt(data_files[x], names = True)

# load MUFASA and remove satellites
#_data['MUFASA'] = _np.genfromtxt(data_files['MUFASA'],names=True)
#_data['MUFASA'] = _data['MUFASA'][ _data['MUFASA']['cen_sat'] == 1]

# brooks
# _data['Brooks'] = _np.genfromtxt(data_files['Brooks'], names = True)

#
# Maybe pre-compute gas fractions and things here for all the data
#

#
#data = {} # now convert ndarrays of data into dictionaries for each data set
#for k in _data.keys():
#    data[k] = {}
#    for l in _data[k].dtype.names:
#        data[k][l] = _data[k][l]


# Load up the observational sample:
_Bradford_2015 = fits.open(_data_path + 'table_1_bradford_2015.fits')

data['Bradford2015'] = {}
Bradford_keys = {'R_eff' : 'R_EFF', 'R_eff_err' : 'R_eff_err',
                 'MHI'   : 'M_HI' , 'Mstar'     : 'M_STAR',
                 'Mstar_err' : 'M_STAR_ERR'}
for k in ['R_eff', 'R_eff_err', 'MHI', 'Mstar', 'Mstar_err']:
    data['Bradford2015'][k] =  _Bradford_2015[1].data[Bradford_keys[k]]
data['Bradford2015']['log_Mstar'] = _np.log10(data['Bradford2015']['Mstar'])
data['Bradford2015']['log_MHI']   = _np.log10(data['Bradford2015']['MHI'])

#
# Pre-compute some things for each data set
#
def _compute_fgas(mstar, mgas, log = True):
    if log:
        fgas = 10**(mgas) / (10**(mstar) + 10**(mgas))
    else:
        fgas = mgas / (mstar + mgas)
    return fgas

data['illustris']['fgas'] = _compute_fgas(data['illustris']['log_Mstar'], data['illustris']['log_MHI'])
data['SAM']['fgas'] = _compute_fgas(data['SAM']['log_Mstar'], data['SAM']['log_Mcold'])
data['MUFASA']['fgas'] = _compute_fgas(data['MUFASA']['log_Mstar'], data['MUFASA']['log_Mcold'])
data['EAGLE']['fgas']  = _compute_fgas(data['EAGLE']['log_Mstar'], data['EAGLE']['log_Mcold'])
data['Bradford2015']['fgas'] = _compute_fgas(data['Bradford2015']['log_Mstar'], data['Bradford2015']['log_MHI'])

# data['Brooks']['fgas'] = _compute_fgas(data['Brooks']['Mstar'], data['Brooks']['HI_Mass'],log=False)

#
#
for k in ['_10Myr','_1Gyr']:
    data['MUFASA']['SFR' + k] = 10.0**(data['MUFASA']['log_SFR' + k])
    data['MUFASA']['SFR' + k][ data['MUFASA']['log_SFR' + k] == -99.0 ] = 0.0 # flag for zero
    data['MUFASA']['sSFR' + k] = data['MUFASA']['SFR' + k] / 10.0**(data['MUFASA']['log_Mstar'])

    data['EAGLE']['sSFR' + k] = data['EAGLE']['SFR' + k] / 10.0**(data['EAGLE']['log_Mstar'])
    data['EAGLE']['log_sSFR' + k] = _np.log10(data['EAGLE']['sSFR' + k])

for k in ['_100kyr', '_100Myr']:
    data['SAM']['sSFR' + k] = 10.0**(data['SAM']['log_sSFR' + k])

#    data['SAM']['sSFR' + k][ data['SAM']['log_sSFR' + k] == -

data['illustris']['sSFR_10Myr'] = data['illustris']['sSFR_20Myr']
data['SAM']['sSFR_10Myr'] = data['SAM']['sSFR_100kyr']
data['SAM']['sSFR_1Gyr'] = data['SAM']['sSFR_100Myr']

# del(_data) - take off memory since we don't need this anymore
#data['Brooks']['sSFR_20Myr'] = data['Brooks']['SFR_20Myr'] / data['Brooks']['Mstar']
#data['Brooks']['sSFR_1Gyr']  = data['Brooks']['SFR_1Gyr'] / data['Brooks']['Mstar']
#data['Brooks']['log_Mstar']  = _np.log10(data['Brooks']['Mstar'])
#data['Brooks']['log_MHI']    = _np.ones(_np.size(data['Brooks']['HI_Mass'])) * -99
#data['Brooks']['log_MHI'][ data['Brooks']['HI_Mass'] > 0] = _np.log10(data['Brooks']['HI_Mass'][data['Brooks']['HI_Mass']>0])
#data['Brooks']['log_MHI'][ data['Brooks']['HI_Mass'] == 0.0] = -99.0


#
# SANDBOX: play with some of the data here to see what happens
#
f_HI = 0.735
HI_APPROXIMATION = True
if HI_APPROXIMATION: # compute M_HI from M_cold
    print "CAUTION: HI Approximation on - cold gas mass converted to HI with constant factor of ", f_HI
    for k in data.keys():
        if not 'log_MHI' in data[k].keys():
            data[k]['log_MHI'] = data[k]['log_Mcold'] + _np.log10(f_HI)
