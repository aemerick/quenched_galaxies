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

colors = {'Illustris' : 'C0', 'SCSAM' : 'C1', 'MUFASA' : 'C2', 'Brooks' : 'C3', 'EAGLE' : 'C4', 'Bradford2015' : 'C5'}

data_files = { # 'Illustris_extended': "Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv",
               # 'Illustris_all'     : "Illustris1_extended_individual_galaxy_values_z0.csv",
#              'Illustris'         :   "Illustris1_individual_galaxy_values_z0.csv",
              'Illustris'         : ['Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv',
                                     'project3_Illustris1_individual_galaxy_values_all1e8Msunh_z0.csv'],
              #'SAM'               : 'sc_sam_cat_day2.txt', #"SAM_group_catalog.dat",
              'SCSAM'                : 'newSCSAMgalprop.dat',
              #'SAM_cen_and_sat'   : "SAM_group_catalog_cenANDsat.dat",
              'Brooks'            : ["brooks_updated_catalog.dat","brooks.cca.centrals_info.dat"],
              'MUFASA'            : "MUFASA_GALAXY.txt",
              'EAGLE'             : ['EAGLE_RefL0100_Mcoldgas_allabove1.8e8Msun.txt', 
                                     'EAGLE_RefL0100_Mstarin12HalfMassRadStars_allabove1.8e8Msun.txt',
                                     'EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt']}

data_dtypes = {'EAGLE'     : [ _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mcold'], 'formats' : ['f8','f8','f8']}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'StellarRhalf', 'log_Mstar', 'log_Mstar1Rhalf', 'log_Mstar2Rhalf'], 'formats':['f8']*6}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mstar', 'SFR_10Myr', 'SFR_1Gyr'], 'formats' : ['f8','f8','f8','f8','f8']})],
               'Illustris' : [ _np.dtype( {'names' : ['log_Mstar','sSFR_10Myr','sSFR_20Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI','Z_SF_gas','log_MBH','cen_sat'], 'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8','u1']}),
                               _np.dtype( {'names' : ['r_half', 'log_Mstar_1Rh', 'log_Mstar_2Rh', 'log_MHI_1Rh', 'log_MHI_2Rh', 'log_Mcold_1Rh', 'log_Mcold_2Rh', 'SFR_0_1Rh',' SFR_0_2Rh', 'SFR_100Myr_1Rh','SFR_100Myr_2Rh','SFR_1Gyr_1Rh','SFR_2Gyr_2Rh'], 'formats' : ['f8']*13})],
               #'SAM'       :   _np.dtype( {'names' : ['log_Mstar','log_sSFR_100kyr','log_sSFR_100Myr','log_Mcold','sigma_8','log_MBH'], 'formats': ['f8']*6}),
               'SCSAM' : None, # load elsewhere

               'MUFASA'    :   _np.dtype( {'names' : ['x','y','z','vx','vy','vz','log_Mstar','log_SFR_10Myr','log_SFR_1Gyr','log_Mcold','log_Z_SFR','cen_sat'], 'formats': ['f8']*12 + ['u1']}),
               'Brooks'    : [ _np.dtype( {'names' : ['Sim','Grp','Mstar','Mvir','SFR_10Myr','SFR_1Gyr','time_last_SF','MHI', 'X','Y','Z','VX','VY','VZ','parent'], 'formats' : ['U5','U3'] + ['f8']*12 + ['u1']}),
                               _np.dtype( {'names' : ['Sim','Grp','r_half','Mstar_1Rh','Mstar_2Rh','MHI_1Rh','MHI_2Rh','MHI','Mcold_1Rh','Mcold_2Rh','Mcold','Mstar_100Myr_1Rh','Mstar_1Gyr_1Rh','Mstar_100Myr_2Rh','Mstar_1Gyr_2Rh','Mstar_100Myr','Mstar_1Gyr'], 'formats' : ['U5','U3'] + ['f8']*15})]
              }

delimiters = {} # assign delimiters and exceptions
for k in data_dtypes.keys():
    delimiters[k] = None
delimiters['Illustris']  = ','

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
_data['MUFASA']    = _data['MUFASA'][ _data['MUFASA']['cen_sat'] == 1]             # remove satellites

_data['Illustris'][1] = _data['Illustris'][1][ _data['Illustris'][0]['cen_sat'] == 1]       # remove satellites
_data['Illustris'][0] = _data['Illustris'][0][ _data['Illustris'][0]['cen_sat'] == 1]

#  Collate the two Brooks data sets together grabbing only centrals that 
#  exist in both files
#  - kinda gross, but only option really
#

# generate list of unique names for all galaxies (Sim name + group number)
full_names = [None,None]
for i in [0,1]:
    full_names[i] = [str(_data['Brooks'][i]['Sim'][j])  + str(_data['Brooks'][i]['Grp'][j]) for j in _np.arange(_np.size(_data['Brooks'][i]['Sim']))]

#  Cross-match: second data set is centrals only --- remove entries from first data set if they DNE in second
select = _np.zeros(_np.size(full_names[0]))
for i,k in enumerate(full_names[0]):
    select[i] = (k in full_names[1])
select = select.astype(bool)
_data['Brooks'][0] = _data['Brooks'][0][ select ]


# load scdata
scdata = _np.genfromtxt('./Data/newSCSAMgalprop.dat', names = True, skip_header = 46)
scdata = scdata[scdata['sat_type'] == 0]
scdata = scdata[scdata['mstar'] > 0]
_data['SCSAM'] = scdata

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

# screw with SCSAM data
data['SCSAM']['Mhalo']       = data['SCSAM']['mhalo']        * 1.0E9
data['SCSAM']['Mstar']       = data['SCSAM']['mstar']        * 1.0E9
data['SCSAM']['Mcold']       = data['SCSAM']['mcold']        * 1.0E9
data['SCSAM']['SFR_20Myr']   = data['SCSAM']['sfr_ave20M']
data['SCSAM']['SFR_100Myr']  = data['SCSAM']['sfr_ave100M']
data['SCSAM']['SFR_1Gyr']    = data['SCSAM']['sfr_ave1G']
data['SCSAM']['sSFR_100Myr'] = data['SCSAM']['SFR_100Myr']   / data['SCSAM']['Mstar']
data['SCSAM']['sSFR_1Gyr']   = data['SCSAM']['SFR_1Gyr']     / data['SCSAM']['Mstar']
data['SCSAM']['sSFR_20Myr']  = data['SCSAM']['SFR_20Myr']    / data['SCSAM']['Mstar']


for k in ['SFR_100Myr','sSFR_100Myr','SFR_1Gyr','sSFR_1Gyr', 'Mhalo','Mstar','Mcold']:
    data['SCSAM']['log_' + k] = _np.log10(data['SCSAM'][k])

# compute some things for Brooks dataset:
data['Brooks']['SFR_100Myr']     = data['Brooks']['Mstar_100Myr']     / ( 100.0E6)  #  SFR in Msun / yr
data['Brooks']['SFR_100Myr_1Rh'] = data['Brooks']['Mstar_100Myr_1Rh'] / ( 100.0E6)
data['Brooks']['SFR_100Myr_2Rh'] = data['Brooks']['Mstar_100Myr_2Rh'] / ( 100.0E6)
data['Brooks']['SFR_1Gyr_1Rh']   = data['Brooks']['Mstar_1Gyr_1Rh']   / (1000.0E6)
data['Brooks']['SFR_1Gyr_2Rh']   = data['Brooks']['Mstar_1Gyr_2Rh']   / (1000.0E6)
data['Brooks']['sSFR_100Myr']     = data['Brooks']['SFR_100Myr']     / data['Brooks']['Mstar']
data['Brooks']['sSFR_1Gyr']     = data['Brooks']['SFR_1Gyr']     / data['Brooks']['Mstar']
data['Brooks']['sSFR_100Myr_1Rh'] = data['Brooks']['SFR_100Myr_1Rh'] / data['Brooks']['Mstar_1Rh']
data['Brooks']['sSFR_100Myr_2Rh'] = data['Brooks']['SFR_100Myr_2Rh'] / data['Brooks']['Mstar_2Rh']
data['Brooks']['sSFR_1Gyr_1Rh']   = data['Brooks']['SFR_1Gyr_1Rh']   / data['Brooks']['Mstar_1Rh']
data['Brooks']['sSFR_1Gyr_2Rh']   = data['Brooks']['SFR_1Gyr_2Rh']   / data['Brooks']['Mstar_2Rh']

for k in ['sSFR_100Myr','sSFR_1Gyr','sSFR_100Myr_1Rh','sSFR_1Gyr_2Rh', 'sSFR_100Myr_1Rh','sSFR_1Gyr_2Rh']:
    data['Brooks']['log_' + k] = data['Brooks'][k]

for k in data['Brooks'].keys():
    if any( [x in k for x in ['MHI','Mcold','Mstar']]):
        data['Brooks']['log_' + k] = _np.log10(data['Brooks'][k])


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

data['Illustris']['fgas'] = _compute_fgas(data['Illustris']['log_Mstar'], data['Illustris']['log_MHI'])
#data['SAM']['fgas'] = _compute_fgas(data['SAM']['log_Mstar'], data['SAM']['log_Mcold'])
data['MUFASA']['fgas'] = _compute_fgas(data['MUFASA']['log_Mstar'], data['MUFASA']['log_Mcold'])
data['EAGLE']['fgas']  = _compute_fgas(data['EAGLE']['log_Mstar'], data['EAGLE']['log_Mcold'])
data['Bradford2015']['fgas'] = _compute_fgas(data['Bradford2015']['log_Mstar'], data['Bradford2015']['log_MHI'])
data['Brooks']['fgas'] = _compute_fgas(data['Brooks']['log_Mstar'], data['Brooks']['log_MHI'])
data['SCSAM']['fgas'] = _compute_fgas(data['SCSAM']['log_Mstar'], data['SCSAM']['log_Mcold'])

#
#
for k in ['_10Myr','_1Gyr']:
    data['MUFASA']['SFR' + k] = 10.0**(data['MUFASA']['log_SFR' + k])
    data['MUFASA']['SFR' + k][ data['MUFASA']['log_SFR' + k] == -99.0 ] = 0.0 # flag for zero
    data['MUFASA']['sSFR' + k] = data['MUFASA']['SFR' + k] / 10.0**(data['MUFASA']['log_Mstar'])

    data['EAGLE']['sSFR' + k] = data['EAGLE']['SFR' + k] / 10.0**(data['EAGLE']['log_Mstar'])
    data['EAGLE']['log_sSFR' + k] = _np.log10(data['EAGLE']['sSFR' + k])

#for k in ['_100kyr', '_100Myr']:
#    data['SAM']['sSFR' + k] = 10.0**(data['SAM']['log_sSFR' + k])


data['Illustris']['sSFR_10Myr'] = data['Illustris']['sSFR_20Myr']

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


# Points are approximate only, and should be used with care.
#
# From Cantinella et. al. 2013 and Brown et. al. 2015
#   as taken from a plot digitizer from points in Fig. 4
#   of Crain et. al. 2016 (EAGLE HI in galaxies paper)
#

cantinella_13 = {'log_Mstar': _np.array([10.1251, 10.3728, 10.6289, 10.86818, 11.12782]),
                 'log_MHI': _np.array([9.5548, 9.5548, 9.67491, 9.64664, 9.752650])}
cantinella_13['MHI']   = 10**(cantinella_13['log_MHI'])
cantinella_13['Mstar'] = 10**(cantinella_13['log_Mstar'])
cantinella_13['fgas']      = _compute_fgas(cantinella_13['log_Mstar'], cantinella_13['log_MHI'])

brown_15 = {'log_Mstar' : _np.array([9.1984, 9.6298, 10.1451, 10.6245, 11.08389]),
            'log_MHI'   : _np.array([9.3639, 9.45583, 9.51590, 9.61484, 9.67378])}
brown_15['MHI']   = 10.0**(brown_15['log_MHI'])
brown_15['Mstar'] = 10.0**(brown_15['log_Mstar'])
brown_15['fgas']  = _compute_fgas(brown_15['log_Mstar'], brown_15['log_MHI'])

