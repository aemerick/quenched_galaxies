"""
   This file defines the filepaths and data structures for each of the
   data files, loads the data into a single dictionary with uniform
   keyword names across simulations for *most* parameters, and calculates
   additional dataset information to either generalize existing information
   (e.x. SFR's to sSFR's and vice-versa) or compute new properties (e.x. gas
   fraction).

   Warning: This is is a very ugly script. I did some hard coding for
   specific data files in order to minimize time spent coding this up.
   Generalizing this to readily read in any file (provided one defined
   column names for that file first) and auto-compute all releveant information
   would be really nice.

   To Do: Re-write this to use a data class / functions...
"""
import matplotlib as _mpl
_mpl.use('Agg', warn = False)

import numpy as _np
from matplotlib import rc as _rc
from astropy.io import fits

from quenched_galaxies.function_definitions import fit_SFMS, compute_fgas, LOGZERO
from quenched_galaxies.tools import NoStdStreams


# define general plot styles:
fsize = 17
_rc('text', usetex=False)
_rc('font', size=fsize)#, ftype=42)
line_width = 3.5
point_size = 30

INCLUDE_SATELLITES = True

#
# set data paths - CHANGE THIS
#
_home = "/home/aemerick/code/quenched_galaxies/"
_data_path = _home + "Data/"

# set colors associated with datasets
colors = {'Illustris' : 'C9', 'SCSAM' : 'C1', 'MUFASA' : 'C2', 'Brooks' : 'C3', 'EAGLE' : 'C0', 'Bradford2015' : 'C5',
          'catinella13' : 'C4', 'brown15' : 'black', 'MUFASA_ari' : 'navy', 'xGASS' : 'black'}

#
# Data file paths. If multiple for a single dataset, give as list. This is handled later to combine to single
#                  dictionary
#
data_files = { # 'Illustris_extended': "Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv",
               # 'Illustris_all'     : "Illustris1_extended_individual_galaxy_values_z0.csv",
               # 'Illustris'         :   "Illustris1_individual_galaxy_values_z0.csv",
              'Illustris'         : ['Illustris1_extended_individual_galaxy_values_all1e8Msunh_z0.csv',  #
                                     'project3_Illustris1_individual_galaxy_values_all1e8Msunh_z0.csv'], #
              #'SAM'               : 'sc_sam_cat_day2.txt', #"SAM_group_catalog.dat",
              'SCSAM'                : 'SCSAMgalprop_updatedVersion.dat',
              'xGASS'             : 'xGASS_representative_sample.ascii.txt',
              #'SAM_cen_and_sat'   : "SAM_group_catalog_cenANDsat.dat",
              'Brooks'            : ["brooks_updated_catalog.dat","brooks.cca.centrals_info.dat"],
              # MUFASA dataset obtained from Romeel / Mika
              'MUFASA'            : "MUFASA_GALAXY_extra.txt",
              'EAGLE'             : ['0-29810EAGLE_RefL0100Hash_MHIH2Rhalfmass_allabove1.8e8Msun.txt',
                                     '0-29810EAGLE_galIDs_MvirsMdm_RefL0100.txt',
                                     'EAGLE_RefL0100_Mcoldgas_allabove1.8e8Msun.txt',
                                     'EAGLE_RefL0100_Mstarin12HalfMassRadStars_allabove1.8e8Msun.txt',
                                     'EAGLE_RefL0100_MstarSFR_allabove1.8e8Msun.txt',
                                     '0-29810EAGLE_RefL0100Hash_XYZMHIH2RhalfmassOnEoST1e4_allabove1.8e8Msun.txt',
                                     '0-29810EAGLE_galIDs_Aperture30kpc70kpcMstarMgas.txt',
                                     'EAGLE_RefL0100_MstarSFR100Myr_allabove1.8e8Msun.txt'],
              'NSA_catalog'       : 'dickey_NSA_iso_lowmass_gals.txt',
              # MUFASA_ari: mufasa dataset obtained from Ari
              'MUFASA_ari'        : 'halos_m50n512_z0.0.txt'}

#
#
# Gross bit. Set dtype and format for every column in every data file
#            this is done to have a somewhat consistent naming scheme
#            for the data dictionary, even though names are different among
#            files.
#
data_dtypes = {'EAGLE'     : [ _np.dtype( {'names' : ['GroupNum','SubGroupNum','r_half_old','log_MHI_old','log_MHI_1Rh_old','log_MHI_2Rh_old','log_MH2_old','log_MH2_1Rh_old','log_MH2_2Rh_old'], 'formats' : ['f8']*9}),
                               _np.dtype( {'names' : ['galaxyID', 'log_Mvir', 'log_MvirMean', 'log_MDM_total'], 'formats' : ['f8']*4}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mcold'], 'formats' : ['f8','f8','f8']}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'r_half', 'log_Mstar', 'log_Mstar_1Rh', 'log_Mstar_2Rh'], 'formats':['f8']*6}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum', 'log_Mstar', 'SFR_10Myr', 'SFR_1Gyr', 'cen_sat'], 'formats' : ['f8','f8','f8','f8','f8','f8']}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum','r_half','log_Mtot','log_Mcold','log_MH_p','log_MHe_p','log_MMet_p','log_MHI','log_MHI_1Rh','log_MHI_2Rh','log_MHI_70','log_MH2','log_MH2_1Rh','log_MH2_2Rh','log_MH2_70'], 'formats' : ['f8']*16  }),
                               _np.dtype( {'names' : ['GalaxyID','log_Mstar_30','log_Mstar_70','log_Mgas_30','log_Mgas_70'], 'formats' : ['f8']*5}),
                               _np.dtype( {'names' : ['GroupNum','SubGroupNum','log_Mstar','SFR_10Myr','SFR_100Myr','cen_sat'], 'formats' : ['f8']*6})],
               'Illustris' : [ _np.dtype( {'names' : ['log_Mstar','sSFR_0Myr','sSFR_10Myr','sSFR_20Myr','sSFR_100Myr','sSFR_1Gyr','log_MHI','sigma_8','log_SF_MHI','Z_SF_gas','log_MBH','cen_sat'], 'formats': ['f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','u1']}),
                               _np.dtype( {'names' : ['r_half', 'log_Mstar_1Rh', 'log_Mstar_2Rh', 'log_MHI_1Rh', 'log_MHI_2Rh', 'log_Mcold_1Rh', 'log_Mcold_2Rh', 'SFR_0_1Rh',' SFR_0_2Rh', 'SFR_100Myr_1Rh','SFR_100Myr_2Rh','SFR_1Gyr_1Rh','SFR_2Gyr_2Rh','log_Mvir'], 'formats' : ['f8']*14})],
               #'SAM'       :   _np.dtype( {'names' : ['log_Mstar','log_sSFR_100kyr','log_sSFR_100Myr','log_Mcold','sigma_8','log_MBH'], 'formats': ['f8']*6}),
               'SCSAM' : None, # load elsewhere
               'xGASS' : None,
               'MUFASA'    :   _np.dtype( {'names' : ['x','y','z','vx','vy','vz','log_Mstar','log_SFR_10Myr','log_SFR_100Myr','log_SFR_1Gyr','log_Mcold','log_Z_SFR','cen_sat','log_MHI', 'log_MH2', 'log_Mhalo'], 'formats': ['f8']*12 + ['u1'] + ['f8']*3}),
               'Brooks'    : [ _np.dtype( {'names' : ['Sim','Grp','Mstar','Mvir','SFR_10Myr','SFR_1Gyr','time_last_SF','MHI', 'X','Y','Z','VX','VY','VZ','parent'], 'formats' : ['U5','U3'] + ['f8']*12 + ['u1']}),
                               _np.dtype( {'names' : ['Sim','Grp','r_half','Mstar_1Rh','Mstar_2Rh','MHI_1Rh','MHI_2Rh','MHI','Mcold_1Rh','Mcold_2Rh','Mcold','Mstar_100Myr_1Rh','Mstar_1Gyr_1Rh','Mstar_100Myr_2Rh','Mstar_1Gyr_2Rh','Mstar_100Myr','Mstar_1Gyr'], 'formats' : ['U5','U3'] + ['f8']*15})],
               'NSA_catalog' : _np.dtype({'names' : ['NSAID', 'log_Mstar', 'DHOST', 'D4000', 'HAEW', 'HALPHA_SFR', 'HALPHA_SSFR'] , 'formats' : ['int'] + ['f8']*6}),
               'MUFASA_ari' : _np.dtype({'names' : ['Mhalo','Mstar','MHI','MH2','SFR'], 'formats' : ['f8']*5})
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
# This is maybe gross, but performs things consistently for all data.
#   Loop through all of the datasets above and load them into a dictionary
#   the temporary dictionary is then used to construct the final dictionary
#
_data = {} # to hold actual data
for k in data_files.keys():

    if not isinstance(data_files[k], basestring): # check if multiple files
        _data[k] = [None]*len(data_files[k])
        for i in _np.arange(len(data_files[k])):
            data_path = _data_path + data_files[k][i]
            _data[k][i]  = _np.genfromtxt( data_path, delimiter = delimiters[k],
                                           skip_header = skip_headers[k], dtype = data_dtypes[k][i],
                                           encoding = None)
    else: # single data file to load
        data_path = _data_path + data_files[k]
        _data[k]  = _np.genfromtxt( data_path, delimiter = delimiters[k],
                                    skip_header = skip_headers[k], dtype = data_dtypes[k],
                                    encoding = None)

#
# Filter out galaxies for specific data sets.
#
if not INCLUDE_SATELLITES:
    _data['MUFASA']    = _data['MUFASA'][ _data['MUFASA']['cen_sat'] == 1] # remove satellites
    _data['Illustris'][1] = _data['Illustris'][1][ _data['Illustris'][0]['cen_sat'] == 1] # remove satellites
    _data['Illustris'][0] = _data['Illustris'][0][ _data['Illustris'][0]['cen_sat'] == 1]

#  Collate the two Brooks data sets together grabbing only centrals that
#  exist in both files
#  - kinda gross, but only option really
# generate list of unique names for all galaxies (Sim name + group number)
full_names = [None,None]
for i in [0,1]:
    full_names[i] = [str(_data['Brooks'][i]['Sim'][j])  + str(_data['Brooks'][i]['Grp'][j]) for j in _np.arange(_np.size(_data['Brooks'][i]['Sim']))]

#  Cross-match: second Brooks data set is centrals only --- remove entries from first data set if they DNE in second
select = _np.zeros(_np.size(full_names[0]))
for i,k in enumerate(full_names[0]):
    select[i] = (k in full_names[1])
select = select.astype(bool)
_data['Brooks'][0] = _data['Brooks'][0][ select ]


# load SAM data separately
scdata = _np.genfromtxt(_data_path + data_files['SCSAM'] , names = True, skip_header = 39 ) # 46)
if not INCLUDE_SATELLITES:
    scdata = scdata[scdata['sat_type'] == 0]  # filter out satellites
scdata = scdata[scdata['Mstar'] > 0]      # filter out galaxies with no stars
_data['SCSAM'] = scdata

#    remove xGASS detections that are marginal (2), confused (5), or both (3)
#    take only isolated centrals
#    including upper limits for now
xgass = _np.genfromtxt(_data_path + data_files['xGASS'], names = True)
xgass = xgass[ (xgass['HI_FLAG'] != 2) *\
               (xgass['HI_FLAG'] != 3) *\
               (xgass['HI_FLAG'] != 5)]
xgass = xgass[ xgass['env_code_B'] == 1 ] # isolated centrals
#  data['XGASS'] = data['xGASS'][ data['xGASS']['HIsrc'] != 4] # remove non-detections
_data['xGASS'] = xgass

#
# now for each case, stitch together all data files and datasets to
# a single, cohesive dictionary. This will be structured so a given PROPERTY
# of all galaxies for a dataset of SIMULATION_NAME is given as:
#
#    data[SIMULATION_NAME][PROPERTY]
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

#
# Some more filtering. select centrals only from EAGLE dataset
#
if not INCLUDE_SATELLITES:
    select = data['EAGLE']['cen_sat'] == 1
    for k in data['EAGLE'].keys():
        data['EAGLE'][k] = data['EAGLE'][k][select]

#
#
# The below contains additional field definitions for each data set. Relables
# some fields and adds new fields of similar names for making life easy later on
# converts all un-logged fields to logged fields and vice-versa to have stored
# copies of both in each dataset
#
#


#
# WARNING: Hard coding fix to current eagle data. HI and H2 masses are
#          below what they should be by the below factor. This is only approximate
#
for k in ['log_MHI','log_MHI_70']:
    data['EAGLE'][k] += 0.675

# make sure logged values that are zero are non-zero (set to LOGZERO as flag)
#
for k in ['MHI','MH2','MHI_70','MH2_70','Mstar_30','Mstar_70','Mstar','Mgas_30','Mgas_70','Mcold','MHI_old','MH2_old','MH_p','MHe_p']:
    data['EAGLE'][k] = 10.0**(data['EAGLE']['log_' + k])
    data['EAGLE']['log_' + k] = LOGZERO * _np.ones(_np.size(data['EAGLE'][k]))
    select = data['EAGLE'][k] > 0
    data['EAGLE']['log_' + k][select] = _np.log10(data['EAGLE'][k][select])

# make a combined HI and H2 for all types of HI and H2
for k,k1,k2 in [('MHI_MH2','MHI','MH2'),('MHI_MH2_70','MHI_70','MH2_70')]:
    data['EAGLE'][k] = data['EAGLE'][k1] + data['EAGLE'][k2]
    data['EAGLE']['log_' + k] = LOGZERO * _np.ones(_np.size(data['EAGLE'][k]))
    select = data['EAGLE'][k] > 0
    data['EAGLE']['log_' + k][select] = _np.log10(data['EAGLE'][k][select])

#
# Do some field defines for the MUFASA datasets
#
data['MUFASA_ari']['Mcold'] = data['MUFASA_ari']['MHI'] + data['MUFASA_ari']['MH2']

for k in ['log_MHI','log_MH2','log_Mstar']:
    data['MUFASA'][k][  data['MUFASA'][k] == -99.0 ] = LOGZERO

data['MUFASA']['MHI'] = 10.0**(data['MUFASA']['log_MHI'])
data['MUFASA']['MH2'] = 10.0**(data['MUFASA']['log_MH2'])
data['MUFASA']['Mcold'] = data['MUFASA']['MHI'] + data['MUFASA']['MH2']
data['MUFASA']['log_Mcold'] = _np.log10(data['MUFASA']['Mcold'])
for k in data['MUFASA_ari'].keys():
    data['MUFASA_ari']['log_' + k] = LOGZERO * _np.ones(_np.size(data['MUFASA_ari'][k]))
    select = data['MUFASA_ari'][k] > 0
    data['MUFASA_ari']['log_' + k][select] = _np.log10(data['MUFASA_ari'][k][select])


#
# Hack the mufasa ari temporary data for now - just use same SFR for everything
#
for k in ['10Myr','100Myr','1Gyr']:
    x = data['MUFASA_ari']['SFR']

    data['MUFASA_ari']['SFR_' + k] = 1.0 * x
    data['MUFASA_ari']['sSFR_' + k] = x / data['MUFASA_ari']['Mstar']

    data['MUFASA_ari']['log_SFR_' + k] = LOGZERO * _np.ones(_np.size(x))
    data['MUFASA_ari']['log_sSFR_' + k] = LOGZERO * _np.ones(_np.size(x))

    data['MUFASA_ari']['log_SFR_' + k][x>0] = _np.log10(x[x>0])
    data['MUFASA_ari']['log_sSFR_' + k][x>0] = _np.log10(x[x>0]) / data['MUFASA_ari']['Mstar'][x>0]






# adjust SCSAM data naming to fit conventions above
data['SCSAM']['Mhalo']       = data['SCSAM']['Mhalo']        * 1.0E9 # convert to Msun
data['SCSAM']['Mstar']       = data['SCSAM']['Mstar']        * 1.0E9 # convert to Msun
data['SCSAM']['Mcold']       = data['SCSAM']['Mcold']        * 1.0E9 # convert to Msun
data['SCSAM']['MHI']         = data['SCSAM']['MHI']          * 1.0E9
data['SCSAM']['MH2']         = data['SCSAM']['MH2']
data['SCSAM']['Mcold']       = data['SCSAM']['MHI'] + data['SCSAM']['MH2']
#data['SCSAM']['SFR_20Myr']   = data['SCSAM']['sfr_ave20M']
#data['SCSAM']['SFR_100Myr']  = data['SCSAM']['sfr_ave100M']
#data['SCSAM']['SFR_1Gyr']    = data['SCSAM']['sfr_ave1G']
data['SCSAM']['sSFR_100Myr'] = data['SCSAM']['SFR_100Myr']   / data['SCSAM']['Mstar']
data['SCSAM']['sSFR_1Gyr']   = data['SCSAM']['SFR_1Gyr']     / data['SCSAM']['Mstar']
data['SCSAM']['sSFR_20Myr']  = data['SCSAM']['SFR_20Myr']    / data['SCSAM']['Mstar']
data['SCSAM']['log_SFR_20Myr']   = _np.log10(data['SCSAM']['SFR_20Myr'])
data['SCSAM']['log_SFR_100Myr']  = _np.log10(data['SCSAM']['SFR_100Myr'])
data['SCSAM']['log_SFR_1Gyr']    = _np.log10(data['SCSAM']['SFR_1Gyr'])
for k in ['20Myr','100Myr','1Gyr']:
    data['SCSAM']['log_SFR_' + k][ data['SCSAM']['log_SFR_' + k] == -_np.inf] = LOGZERO

#
# WARNING: Using 20 Myr data as the 10 Myr!!!!
#
data['SCSAM']['sSFR_10Myr'] = data['SCSAM']['sSFR_20Myr']
data['SCSAM']['SFR_10Myr'] = data['SCSAM']['SFR_20Myr']
data['SCSAM']['log_SFR_10Myr'] = data['SCSAM']['log_SFR_20Myr']

for k in ['SFR_100Myr','sSFR_100Myr','SFR_1Gyr','sSFR_1Gyr', 'Mhalo','Mstar','Mcold']:
    data['SCSAM']['log_' + k] = LOGZERO * _np.ones(_np.size(data["SCSAM"][k]))
    data['SCSAM']['log_' + k][ data['SCSAM'][k] > 0.0 ]= _np.log10(data['SCSAM'][k][ data['SCSAM'][k] > 0.0])

#
# rename some things in the xGASS dataset
#
data['xGASS']['log_Mstar'] = data['xGASS']['lgMstar']
data['xGASS']['SFR']       = data['xGASS']['SFR_best']
data['xGASS']['log_MHI']   = data['xGASS']['lgMHI']
data['xGASS']['log_Mhalo'] = data['xGASS']['logMh_Mst_B'] # using group catalog and stellar mass function

for k in ['Mstar','MHI', 'Mhalo']:
    data['xGASS'][k] = 10.0**(data['xGASS']['log_' + k])

for k in ['SFR']:
    data['xGASS']['log_' + k] = LOGZERO * _np.ones(_np.size(data['xGASS'][k]))
    select = data['xGASS'][k] > 0
    data['xGASS']['log_' + k][select] = _np.log10(data['xGASS'][k][select])
data['xGASS']['log_SFR_1Gyr']       = data['xGASS']['log_SFR']
data['xGASS']['log_SFR_100Myr']       = data['xGASS']['log_SFR']


# compute some things for Brooks dataset:
data['Brooks']['SFR_100Myr']     = data['Brooks']['Mstar_100Myr']     / ( 100.0E6)  #  SFR in Msun / yr
data['Brooks']['SFR_100Myr_1Rh'] = data['Brooks']['Mstar_100Myr_1Rh'] / ( 100.0E6)
data['Brooks']['SFR_100Myr_2Rh'] = data['Brooks']['Mstar_100Myr_2Rh'] / ( 100.0E6)
data['Brooks']['SFR_1Gyr_1Rh']   = data['Brooks']['Mstar_1Gyr_1Rh']   / (1000.0E6)
data['Brooks']['SFR_1Gyr_2Rh']   = data['Brooks']['Mstar_1Gyr_2Rh']   / (1000.0E6)
data['Brooks']['sSFR_10Myr']     = data['Brooks']['SFR_10Myr']     / data['Brooks']['Mstar']
data['Brooks']['sSFR_100Myr']     = data['Brooks']['SFR_100Myr']     / data['Brooks']['Mstar']
data['Brooks']['sSFR_1Gyr']     = data['Brooks']['SFR_1Gyr']     / data['Brooks']['Mstar']
data['Brooks']['sSFR_100Myr_1Rh'] = data['Brooks']['SFR_100Myr_1Rh'] / data['Brooks']['Mstar_1Rh']
data['Brooks']['sSFR_100Myr_2Rh'] = data['Brooks']['SFR_100Myr_2Rh'] / data['Brooks']['Mstar_2Rh']
data['Brooks']['sSFR_1Gyr_1Rh']   = data['Brooks']['SFR_1Gyr_1Rh']   / data['Brooks']['Mstar_1Rh']
data['Brooks']['sSFR_1Gyr_2Rh']   = data['Brooks']['SFR_1Gyr_2Rh']   / data['Brooks']['Mstar_2Rh']

for k in ['sSFR_100Myr','sSFR_1Gyr','sSFR_100Myr_1Rh','sSFR_1Gyr_2Rh', 'sSFR_100Myr_1Rh','sSFR_1Gyr_2Rh', 'SFR_100Myr','SFR_1Gyr','SFR_10Myr']:
    data['Brooks']['log_' + k] = _np.log10(data['Brooks'][k])
    data['Brooks']['log_' + k][ data['Brooks']['log_' + k] == -_np.inf ] = LOGZERO

for k in data['Brooks'].keys():
    if any( [x in k for x in ['MHI','Mcold','Mstar']]):
        data['Brooks']['log_' + k] = _np.log10(data['Brooks'][k])

for k in ['0Myr','10Myr','20Myr','100Myr','1Gyr']:
    data['Illustris']['SFR_' + k] = data['Illustris']['sSFR_' + k] * 10.0**(data['Illustris']['log_Mstar'])
    data['Illustris']['log_SFR_' + k] = _np.log10(data['Illustris']['SFR_' + k])
    data['Illustris']['log_SFR_' + k][ data['Illustris']['log_SFR_' + k] == -_np.inf] = LOGZERO # flag


# Load up the observational sample:
_Bradford_2015 = fits.open(_data_path + 'table_1_bradford_2015.fits')

data['Bradford2015'] = {}
Bradford_keys = {'R_eff' : 'R_EFF', 'R_eff_err' : 'R_eff_err',
                 'MHI'   : 'M_HI' , 'Mstar'     : 'M_STAR',
                 'Mstar_err' : 'M_STAR_ERR', 'NSAID' : 'NSAID'}

for k in ['R_eff', 'R_eff_err', 'MHI', 'Mstar', 'Mstar_err', 'NSAID']:
    data['Bradford2015'][k] =  _Bradford_2015[1].data[Bradford_keys[k]]
data['Bradford2015']['log_Mstar'] = _np.log10(data['Bradford2015']['Mstar'])
data['Bradford2015']['log_MHI']   = _np.log10(data['Bradford2015']['MHI'])

data['Bradford2015']['log_SFR_Halpha'] = _np.ones(_np.size( data['Bradford2015']['Mstar'])) * LOGZERO
data['Bradford2015']['log_sSFR_Halpha'] = _np.ones(_np.size( data['Bradford2015']['Mstar'])) * LOGZERO
for select, i in enumerate(data['Bradford2015']['NSAID']):
    if (i in data['NSA_catalog']['NSAID']):
        select2 = _np.where(data['NSA_catalog']['NSAID'] == i)[0]

#        print select, select2, i, data['Bradford2015']['NSAID'][select], data['NSA_catalog']['NSAID'][select2]

        data['Bradford2015']['log_SFR_Halpha'][select]  = data['NSA_catalog']['HALPHA_SFR'][select2]
        data['Bradford2015']['log_sSFR_Halpha'][select] = data['NSA_catalog']['HALPHA_SSFR'][select2]

data['Bradford2015']['SFR_Halpha']  = 10.0**(data['Bradford2015']['log_SFR_Halpha'])
data['Bradford2015']['sSFR_Halpha'] = 10.0**(data['Bradford2015']['log_sSFR_Halpha'])

data['Bradford2015']['log_SFR_1Gyr'] = data['Bradford2015']['log_SFR_Halpha']
data['Bradford2015']['log_SFR_100Myr'] = data['Bradford2015']['log_SFR_Halpha']


#print "BRADFORD: ", _np.size(data['Bradford2015']['Mstar'][data['Bradford2015']['log_SFR_Halpha'] < -90]), _np.size(data['Bradford2015']['Mstar'][data['Bradford2015']['log_SFR_Halpha'] > -90])






#
#
for k in ['_10Myr','_100Myr','_1Gyr']:
    data['MUFASA']['SFR' + k] = 10.0**(data['MUFASA']['log_SFR' + k])
    data['MUFASA']['SFR' + k][ data['MUFASA']['log_SFR' + k] == LOGZERO ] = 0.0 # flag for zero
    data['MUFASA']['sSFR' + k] = data['MUFASA']['SFR' + k] / 10.0**(data['MUFASA']['log_Mstar'])

    #if not k == '_100Myr':
    data['EAGLE']['log_SFR' + k] = _np.log10(data['EAGLE']['SFR' + k])
    data['EAGLE']['log_SFR' + k][ data['EAGLE']['log_SFR' + k] == -_np.inf] = LOGZERO

    data['EAGLE']['sSFR' + k] = data['EAGLE']['SFR' + k] / 10.0**(data['EAGLE']['log_Mstar'])
    data['EAGLE']['log_sSFR' + k] = _np.log10(data['EAGLE']['sSFR' + k])

#for k in ['_100kyr', '_100Myr']:
#    data['SAM']['sSFR' + k] = 10.0**(data['SAM']['log_sSFR' + k])

# generate logged sSFR fields for all data sets and flag those with zero SFR
for simname in data.keys():
    sSFR_types = [k for k in data[simname].keys() if 'sSFR' in k]
    for k in sSFR_types:
        if not ( ('log_' + k) in data[simname].keys() ):
            data[simname]['log_' + k] = _np.log10(data[simname][k])
            data[simname]['log_' + k][ data[simname][k] == 0.0 ] = LOGZERO

#
#
# Fit the SFMS for each data set:
#
#     Be careful with this
#
#

def _fit_sfms(sim_name, years):
    """
    Function to fit the SFMS for each data set and store the result in
    the dataset dictionary
    """
    for k in years:
        log_mstar = data[sim_name]['log_Mstar']
        log_sfr   = data[sim_name]['log_SFR_' + k]

        # print sim_name, k, _np.min(log_sfr), _np.max(log_sfr), _np.median(log_sfr), _np.min(log_sfr[log_sfr>LOGZERO])
        # print sim_name, k, _np.min(log_mstar), _np.max(log_mstar), _np.median(log_mstar)
        with NoStdStreams(): # this function call can be loud... keep it quiet
            D, m_fit, sfr_fit = fit_SFMS(log_mstar, log_sfr)
        data[sim_name]['D_SFMS_' + k] = D
        data[sim_name]['SFMS_fit_' + k] = [m_fit, sfr_fit]

        #print sim_name, k, _np.min(D), _np.min(D[D>LOGZERO]), _np.max(D), _np.median(D)

    return


_fit_sfms('Illustris', ['0Myr','10Myr','20Myr','100Myr','1Gyr'])
# _fit_sfms('Brooks'   , ['100Myr','1Gyr'])
_fit_sfms('EAGLE'    , ['10Myr','100Myr','1Gyr'])
_fit_sfms('MUFASA'   , ['10Myr','100Myr','1Gyr'])
_fit_sfms('MUFASA_ari'   , ['10Myr','100Myr','1Gyr'])
_fit_sfms('SCSAM'    , ['10Myr','20Myr','100Myr','1Gyr'])



#
# SANDBOX: play with some of the data here to see what happens
#
f_HI = 0.735
HI_APPROXIMATION = True
if HI_APPROXIMATION: # compute M_HI from M_cold
    for k in data.keys():
        if k == 'NSA_catalog':
            continue
        if k == 'MUFASA_ari':
            continue

        if not 'log_MHI' in data[k].keys():
            print "CAUTION: HI Approximation on - For data with no explicit HI mass, HI mass is taken as the cold gas mass times the constant factor ", f_HI

            data[k]['log_MHI'] = data[k]['log_Mcold'] + _np.log10(f_HI)

#
# FOR ALL DATASETS, WE SHOULD BE LOOKING AT THE GAS PROPERTIES
# AND COLD GAS FRACTIONS UNIFORMLY. This means adopting:
#    1) M_HI + M_H2  if H2 is computed
#    2) M_HI if only M_HI is provided
#    3) M_cold * 0.735 if neither is computed
#
#   call this: M_cold_gas

for simname in data.keys():
    if simname == 'NSA_catalog':
        continue

    for k in ['MHI','MH2']:
        if (not ('log_' + k) in data[simname].keys()) and k in data[simname].keys():
            data[simname]['log_' + k] = LOGZERO * _np.ones(_np.size(data[simname][k]))
            data[simname]['log_' + k][data[simname][k]>0.0] = _np.log10(data[simname][k][data[simname][k]>0.0])

    if 'log_MHI' in data[simname].keys():
        mass = 10.0**(data[simname]['log_MHI'])
        mass[ data[simname]['log_MHI'] == LOGZERO] = 0.0

        if 'log_MH2' in data[simname].keys():
            temp = 10.0**(data[simname]['log_MH2'])
            temp[ data[simname]['log_MH2'] == LOGZERO ] = 0.0
            mass = mass + temp

    else:
        mass = 10.0**(data[simname]['log_Mcold'])
        mass[ data[simname]['log_Mcold'] == LOGZERO ] = 0.0


    data[simname]['Mcold_gas'] = mass * 1.0
    data[simname]['log_Mcold_gas'] = LOGZERO * _np.ones(_np.size(data[simname]['Mcold_gas']))
    data[simname]['log_Mcold_gas'][ data[simname]['Mcold_gas'] > 0.0 ] = _np.log10( data[simname]['Mcold_gas'][data[simname]['Mcold_gas']>0.0])


for simname in data.keys():
    if simname == 'NSA_catalog':
        continue

    # print simname, _np.size(data[simname]['log_Mstar']), _np.size(data[simname]['log_Mcold_gas'])
    data[simname]['fgas'] = compute_fgas(data[simname]['log_Mstar'], data[simname]['log_Mcold_gas'])
    data[simname]['fgas_obs'] = compute_fgas(data[simname]['log_Mstar'], data[simname]['log_MHI'])

    data[simname]['log_fgas'] = LOGZERO * _np.ones(_np.size(data[simname]['fgas']))
    data[simname]['log_fgas'][ data[simname]['fgas']  > 0.0] =_np.log10(data[simname]['fgas'][ data[simname]['fgas'] > 0.0])
    data[simname]['log_fgas_obs'] = LOGZERO * _np.ones(_np.size(data[simname]['fgas_obs']))
    data[simname]['log_fgas_obs'][ data[simname]['fgas_obs']  > 0.0] =_np.log10(data[simname]['fgas_obs'][ data[simname]['fgas_obs'] > 0.0])

for k in data.keys():
    if simname == 'NSA_catalog':
        continue

    if ('log_Mstar_1Rh' in data[k].keys()) and (('log_MHI_1Rh' in data[k].keys()) or ('log_Mcold_1Rh' in data[k].keys())):
        #gas_key = 'log_MHI'
        #if not (gas_key in data[k].keys()):
        #    gas_key = 'log_Mcold'

        if 'log_MHI_1Rh' in data[k].keys():
            mass = 10.0**(data[k]['log_MHI_1Rh'])
            mass[ data[k]['log_MHI_1Rh'] == LOGZERO] = 0.0

            if 'log_MH2_1Rh' in data[k].keys():
                temp = 10.0**(data[k]['log_MH2_1Rh'])
                temp[ data[k]['log_MH2_1Rh'] == LOGZERO ] = 0.0
                mass = mass + temp
        data[k]['Mcold_gas_1Rh']     = mass * 1.0
        data[k]['log_Mcold_gas_1Rh'] = LOGZERO * _np.ones(_np.size(data[k]['Mcold_gas_1Rh']))
        data[k]['log_Mcold_gas_1Rh'] [ data[k]['Mcold_gas_1Rh'] > 0.0]= _np.log10(mass[ data[k]['Mcold_gas_1Rh']>0.0 ] * 1.0)
#########
        if 'log_MHI_2Rh' in data[k].keys():
            mass = 10.0**(data[k]['log_MHI_2Rh'])
            mass[ data[k]['log_MHI_2Rh'] == LOGZERO] = 0.0

            if 'log_MH2_2Rh' in data[k].keys():
                temp = 10.0**(data[k]['log_MH2_2Rh'])
                temp[ data[k]['log_MH2_2Rh'] == LOGZERO ] = 0.0
                mass = mass + temp
        data[k]['Mcold_gas_2Rh']     = mass * 1.0
        data[k]['log_Mcold_gas_2Rh'] = LOGZERO * _np.ones(_np.size(data[k]['Mcold_gas_2Rh']))
        data[k]['log_Mcold_gas_2Rh'] [ data[k]['Mcold_gas_2Rh'] > 0.0]= _np.log10(mass[ data[k]['Mcold_gas_2Rh']>0.0 ] * 1.0)

        gas_key = 'log_Mcold_gas'
        data[k]['fgas_1Rh'] = compute_fgas(data[k]['log_Mstar_1Rh'], data[k][gas_key + '_1Rh'])
        data[k]['fgas_2Rh'] = compute_fgas(data[k]['log_Mstar_2Rh'], data[k][gas_key + '_2Rh'])
        data[k]['log_fgas_1Rh'] = LOGZERO * _np.ones(_np.size(data[k]['fgas_1Rh']))
        data[k]['log_fgas_1Rh'][ data[k]['fgas_1Rh']  > 0.0] =_np.log10(data[k]['fgas_1Rh'][ data[k]['fgas_1Rh'] > 0.0])
        data[k]['log_fgas_2Rh'] = LOGZERO * _np.ones(_np.size(data[k]['fgas_2Rh']))
        data[k]['log_fgas_2Rh'][ data[k]['fgas_2Rh']  > 0.0] =_np.log10(data[k]['fgas_2Rh'][ data[k]['fgas_2Rh'] > 0.0])

        gas_key = 'log_MHI'
        data[k]['fgas_obs_1Rh'] = compute_fgas(data[k]['log_Mstar_1Rh'], data[k][gas_key + '_1Rh'])
        data[k]['fgas_obs_2Rh'] = compute_fgas(data[k]['log_Mstar_2Rh'], data[k][gas_key + '_2Rh'])
        data[k]['log_fgas_obs_1Rh'] = LOGZERO * _np.ones(_np.size(data[k]['fgas_obs_1Rh']))
        data[k]['log_fgas_obs_1Rh'][ data[k]['fgas_obs_1Rh']  > 0.0] = _np.log10(data[k]['fgas_obs_1Rh'][ data[k]['fgas_obs_1Rh'] > 0.0])
        data[k]['log_fgas_obs_2Rh'] = LOGZERO * _np.ones(_np.size(data[k]['fgas_obs_2Rh']))
        data[k]['log_fgas_obs_2Rh'][ data[k]['fgas_obs_2Rh']  > 0.0] = _np.log10(data[k]['fgas_obs_2Rh'][ data[k]['fgas_obs_2Rh'] > 0.0])

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
cantinella_13['fgas']      = compute_fgas(cantinella_13['log_Mstar'], cantinella_13['log_MHI'])

brown_15 = {'log_Mstar' : _np.array([9.1984, 9.6298, 10.1451, 10.6245, 11.08389]),
            'log_MHI'   : _np.array([9.3639, 9.45583, 9.51590, 9.61484, 9.67378])}
brown_15['MHI']   = 10.0**(brown_15['log_MHI'])
brown_15['Mstar'] = 10.0**(brown_15['log_Mstar'])
brown_15['fgas']  = compute_fgas(brown_15['log_Mstar'], brown_15['log_MHI'])
#
#
#
#
#
#
data['EAGLE']['volume']     = (100.0)**3
data['MUFASA']['volume']    = (50.0)**3
data['Illustris']['volume'] = (100.0)**3
data['SCSAM']['volume']     = (100.0)**3
