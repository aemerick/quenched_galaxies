from definitions import * # -- local - predefines

# -- need to port this to local file so no one needs to download this code:
from galaxy_analysis.plot.plot_styles import plot_histogram


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.stats import binned_statistic_2d

from quenched_galaxies.function_definitions import fgas_limits, LOGZERO

_output_dir = './plots/'


MSTAR_BINS = np.arange(7.0,13.1,0.2)
DEFAULT_SFR = '100Myr' # default sfr to use (10Myr, 100Myr, 1Gyr, etc.)
#
# Add simulation names to these global dicts for plotting
#
ALL_DATA = ['Illustris','SCSAM' ,'EAGLE','MUFASA','Romulus25', 'Bradford2015']
SIM_DATA = ['Illustris','SCSAM' ,'EAGLE','MUFASA','Romulus25'] # Brooks
LSIM_DATA = ['Illustris','EAGLE','MUFASA','SCSAM','Romulus25'] # Brooks
OBS_DATA = ['Bradford2015', 'xGASS']

point_size = 40        # default point size for the above


#def compute_fgas(mstar, mgas, log = True):
#    """
#    Compute gas fraction as gas_mass / baryon_mass
#    """
#
#    if log:
#        fgas = 10**(mgas) / (10**(mstar) + 10**(mgas))
#    else:
#        fgas = mgas / (mstar + mgas)
#    return fgas

def _check_bins(x,xbins):
    """
    Return number of points in each bin
    """
    N = np.zeros(np.size(xbins) -1)
    for i in np.arange(np.size(xbins)-1):
        N[i] = np.size( x[ (x >= xbins[i]) * (x < xbins[i+1])])
    return N

def _select_scatter_points(x, xbins, threshold = 10):
    """
    Provide cut-arrays that can be used to selectively plot points
    that sit in bins with low counts. Also provides selection for
    points that do sit in bins with enough counts to do statistics.
    """
    N_per_bin = _check_bins(x, xbins)

    scatter_select = np.zeros(np.size(x))
    for i, N in enumerate(N_per_bin):
        if N < threshold:
            scatter_select += ( x >= xbins[i] ) * ( x < xbins[i+1] )
    scatter_select = scatter_select.astype(bool)

    return scatter_select, scatter_select == 0



def _compute_statistics(x, y, xbins, return_dict = False):
    """
    Function to compute statistics on 'y' in each 'xbin'. Original
    version of this returns:

         bin_centers, median, std, Q1, Q3, average, Num_points_per_bin

    and removes bins that are empty. Does not have any catches for low
    number statistics in bins (i.e. will compute these values even for
    N = 1). return_dict = False operates this way. However, return_dict = True
    is new api functionality which returns a dictionary of all the above
    (and more) which makes adding new stats easy in the future. At the moment
    this includes the above and the percentiles for 10% and 90% (poorly named
    'Q10' and 'Q90' respectively).
    """
    flag = -99999999
    # generic function to compute median and statistics
    median = np.ones(np.size(xbins)-1) * flag # leave neg as flag
    average = np.zeros(np.size(median))
    Q1 = np.zeros(np.size(median)); Q3 = np.zeros(np.size(median))
    Q10 = np.zeros(np.size(median)); Q90 = np.zeros(np.size(median))
    std = np.zeros(np.size(median)); N = np.zeros(np.size(median))
    for i in np.arange(np.size(xbins)-1):
        y_select  = y[(x>=xbins[i])*(x<xbins[i+1])]

        if np.size(y_select) >= 1:
            median[i] = np.median( y_select )
            average[i] = np.average( y_select)
            Q10[i]    = np.percentile( y_select, 10)
            Q1[i]     = np.percentile( y_select, 25)
            Q3[i]     = np.percentile( y_select, 75)
            Q90[i]    = np.percentile( y_select, 90)
            std[i]    = np.std(y_select)
            N[i]      = np.size( y_select )

    select = median > flag
    centers = (xbins[1:] + xbins[:-1])*0.5

    if return_dict:

        rdict = {'x' : centers[select], 'median' : median[select], 'std' : std[select],
                 'Q1' : Q1[select], 'Q3' : Q3[select], 'average' : average[select], 'N' : N[select],
                 'Q10' : Q10[select], 'Q90' : Q90[select] }
        return rdict

    else:
        return centers[select], median[select], std[select], Q1[select], Q3[select], average[select], N[select]


def plot_SFMS(mstar_bins = MSTAR_BINS, remove_zero = True, SFR_type = DEFAULT_SFR,
              datasets = SIM_DATA, figdim = (2,2)):
    """
    Plot the SFMS fits for each simulation individually as panels. For simplicity,
    will just be plotting median and IQR for now, but should move to doing
    contours / points.
    """

    if (figdim == (2,2)) and (len(datasets) >  (figdim[0] * figdim[1])):
        figdim = (1, len(datasets))

    fig, ax = plt.subplots(figdim[0],figdim[1]) # hard coded for now - not good
    fig.set_size_inches(12,12)

    axi, axj = 0, 0
    for k in datasets:
        if figdim[0] > 1:
            axindex = (axi,axj)
        else:
            axindex = axi

        SFR   = data[k]['log_SFR_' + SFR_type]
        Mstar = data[k]['log_Mstar']

        # filter
        select = SFR > LOGZERO
        SFR    = SFR[select]
        Mstar  = Mstar[select]

        scatter_select, line_select = _select_scatter_points(Mstar, mstar_bins)
        ax[axindex].scatter(Mstar[scatter_select], SFR[scatter_select], s = point_size, color = colors[k])
        computed_stats = _compute_statistics(Mstar[line_select], SFR[line_select], mstar_bins, return_dict = True)

        ax[axindex].plot(computed_stats['x'], computed_stats['median'], lw = line_width, color = colors[k], ls = '-', label = k)

        ax[axindex].fill_between(computed_stats['x'], computed_stats['Q1'], computed_stats['Q3'], facecolor = colors[k],
                             interpolate = True, lw = line_width, alpha = 0.4)
        ax[axindex].fill_between(computed_stats['x'], computed_stats['Q10'], computed_stats['Q90'], facecolor = colors[k],
                             interpolate = True, lw = line_width, alpha = 0.1)

        try:
            ax[axindex].plot(data[k]['SFMS_fit_' + SFR_type][0], data[k]['SFMS_fit_' + SFR_type][1],
                         lw = line_width, ls = '--', color = 'black')
        except:
            print "I'm a lazy programmer"

        ax[axindex].set_xlim(np.min(mstar_bins), np.max(mstar_bins))
        ax[axindex].set_ylim(-4, 2)
        ax[axindex].legend(loc='upper left')

        ax[axindex].set_xlabel(r'log( M$_{*}$ [M$_{\odot}$])')
        ax[axindex].set_ylabel(r'log( SFR (M$_{\odot}$ yr$^{-1}$)')

        if figdim[0] > 1:
            axj = axj + 1
            if axj >= figdim[1]:
                axj = 0
                axi = axi + 1
        else:
            axi = axi + 1

    plt.tight_layout()
    plt.minorticks_on()
    fig.savefig('SFMS_fits.png')

    return


def plot_halo_stellar_2D_hist(mstar_bins = np.arange(8.0, 12.6, 0.05), halo_bins = np.arange(8.0, 15.0, 0.05),
                              statistic = 'fraction',
                              log_fgas = False, datasets = SIM_DATA, cmap = 'viridis',figdim=(1,4)):
    """
    Plot 2D histograms for each simulation indvidually as panels of gas fraction vs. stellar
    mass. By default, shading is by fraction of galaxies in a given bin.
    """

    if (figdim == (1,4)) and (len(datasets) > figdim[1]):
        figdim = (1,len(datasets))

    ylabel = r'log(M$_{*}$ [M$_{\odot}$])'

    fig, ax = plt.subplots(figdim[0],figdim[1], sharey=True, sharex=True)  # change hard coding
    fig.set_size_inches(figdim[1]*6,figdim[0]*6)

    if figdim == (2,2):
        ax_indexes = [(0,0),(0,1),(1,0),(1,1)]
    else:
        ax_indexes = np.arange(figdim[1])

    stat_name = statistic
    if statistic == 'count':
        cbar_label = r'log(Number of Galaxies)'
        vmin = 0
        vmax = 3
        cmap = 'viridis'
    elif statistic == 'fraction':
        cbar_label = r'log(Fraction of Galaxies)'
        # may need to log the fraction
        vmin =  -4.0
        vmax =  -2.0
        stat_name = 'count'
        cmap = 'magma'
    elif statistic == 'median':
        cbar_label = r'f$_{\rm gas}$'
        vmin = 0.0
        vmax = 1.0
        cmap = 'PRGn'


    xpos_label = 13.75
    ypos_label = 8.25


    for axi, k in zip( ax_indexes, datasets):

        field = None
        if ('log_Mhalo' in data[k].keys()):
            field = 'log_Mhalo'
        elif ('log_Mvir') in data[k].keys():
            field = 'log_Mvir'
        else:
            if ('Mhalo' in data[k].keys()):
                data[k]['log_Mhalo'] = np.log10(data[k]['Mhalo'])
                field = 'log_Mhalo'
            elif ('Mvir' in data[k].keys()):
                data[k]['log_Mvir'] = np.log10(data[k]['Mvir'])
                field = 'log_Mvir'
            else:
                ax[0].set_ylabel(ylabel)
                continue

        Mhalo  = data[k][field]
        Mstar  = data[k]['log_Mstar']


        ngal   = np.size(Mstar)

        if statistic == 'count' or statistic == 'fraction':
            statistic_data = np.ones(ngal)
        else:
            statistic_data = data[k]['fgas']

        N, x_edge, y_edge, binnum = binned_statistic_2d(Mhalo, Mstar, statistic_data,
                                                        statistic = stat_name, bins = (halo_bins, mstar_bins))

        if statistic == 'fraction':
            fraction = N / (1.0 * ngal)
            fraction[fraction <= 0] = LOGZERO
            fraction[fraction >  0] = np.log10(fraction[fraction > 0])
            fraction = np.log10(N / (1.0 * ngal))
            plot_val = fraction
        else:
            plot_val = N


        xmesh, ymesh = np.meshgrid(x_edge, y_edge)



        img1 = ax[axi].pcolormesh(xmesh, ymesh, plot_val.T,
                                  cmap = cmap, vmin = vmin, vmax = vmax)

        if axi == ax_indexes[0]:
            ax[axi].set_ylabel(ylabel)
        elif axi == ax_indexes[-1]:
            divider = make_axes_locatable(ax[axi])
            cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(img1, cax=cax1, label = cbar_label)


        dx = mstar_bins[1] - mstar_bins[0]
        dx = 0
        ax[axi].set_ylim(8.0,12.0) #np.min(mstar_bins) - dx, np.max(mstar_bins) + dx)
        dx = halo_bins[1] - halo_bins[0]
        dx = 0
        ax[axi].set_xlim(9.0, 15.0) #np.min(halo_bins)  - dx, np.max(halo_bins) + dx)


        ax[axi].set_xlabel(r'log( M$_{\rm vir}$ [M$_{\odot}$])')

        ax[axi].text(xpos_label, ypos_label, k)
        plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    plt.minorticks_on()

    outname = 'halo_stellar_2D_' + statistic

    if log_fgas:
        outname += '_log_fgas'

    fig.savefig(outname + '.png')

    plt.close()

    return



def plot_fgas_mstar_2D_hist(mstar_bins = np.arange(8.0, 12.6, 0.1), fgas_bins = None, halo_ratio = False,
                            log_fgas = False, datasets = SIM_DATA, cmap = 'viridis',figdim=(1,4),
                            sSFR_cut = None, outstr = ''):
    """
    Plot 2D histograms for each simulation indvidually as panels of gas fraction vs. stellar
    mass. By default, shading is by fraction of galaxies in a given bin.
    """

    if (figdim == (1,4)) and (len(datasets) > figdim[1]):
        figdim = (1,len(datasets))

    if fgas_bins is None:

        if log_fgas:
            fgas_bins = np.linspace(-4, 0, np.size(mstar_bins))
        else:
            fgas_bins = np.linspace(0.0, 1.0, np.size(mstar_bins))
    if log_fgas:
        ylabel = r'log(f$_{\rm gas}$)'
        xpos_label = 8.5
        ypos_label = -3.75
    else:
        ylabel = r'f$_{\rm gas}$'
        xpos_label = 11.
        ypos_label = 0.925

    fig, ax = plt.subplots(figdim[0],figdim[1], sharey=True, sharex=True)  # change hard coding
    fig.set_size_inches(figdim[1]*6,figdim[0]*6)

    if figdim == (2,2):
        ax_indexes = [(0,0),(0,1),(1,0),(1,1)]
    else:
        ax_indexes = np.arange(figdim[1])

    if halo_ratio:
        xlabel = r'log(M$_{*}$/M$_{\rm vir}$)'
        mstar_bins = np.arange(-4,0,0.1)
    else:
        xlabel = r'log(M$_{*}$ [M$_{\odot}$])'


    for axi, k in zip( ax_indexes, datasets):
        fgas   = data[k]['fgas']
        Mstar  = data[k]['log_Mstar']

        if halo_ratio:
            if ('log_Mhalo' in data[k].keys()):
                field = 'log_Mhalo'
            elif ('log_Mvir') in data[k].keys():
                field = 'log_Mvir'
            else:
                if ('Mhalo' in data[k].keys()):
                    data[k]['log_Mhalo'] = np.log10(data[k]['Mhalo'])
                    field = 'log_Mhalo'
                elif ('Mvir' in data[k].keys()):
                    data[k]['log_Mvir'] = np.log10(data[k]['Mvir'])
                    field = 'log_Mvir'
                else:
                    continue

            Mhalo  = data[k][field]

            x_data = Mstar - Mhalo
        else:
            x_data = Mstar


        if log_fgas:
            select = fgas > 0
            fgas  = np.log10(fgas[select])
            x_data = x_data[select]

        ngal   = np.size(fgas)

        N, x_edge, y_edge, binnum = binned_statistic_2d(x_data, fgas, np.ones(np.shape(fgas)),
                                                        statistic = 'count', bins = (mstar_bins, fgas_bins))

        fraction = N / (1.0 * ngal)
        fraction[fraction <= 0] = LOGZERO
        fraction[fraction >  0] = np.log10(fraction[fraction > 0])
        fraction = np.log10(N / (1.0 * ngal))

        xmesh, ymesh = np.meshgrid(x_edge, y_edge)

        # may need to log the fraction
        vmin =  -4.0
        vmax =  -1

        img1 = ax[axi].pcolormesh(xmesh, ymesh, fraction.T,
                                  cmap = cmap, vmin = vmin, vmax = vmax)

        if axi == ax_indexes[0]:
            ax[axi].set_ylabel(ylabel)
        elif axi == ax_indexes[-1]:
            divider = make_axes_locatable(ax[axi])
            cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
            cbar_label = r'log(Fraction of Galaxies)'
            fig.colorbar(img1, cax=cax1, label = cbar_label)


        dx = mstar_bins[1] - mstar_bins[0]
        dx = 0
        ax[axi].set_xlim( np.min(mstar_bins) - dx, np.max(mstar_bins) + dx)
        dx = fgas_bins[1] - fgas_bins[0]
        dx = 0
        ax[axi].set_ylim( np.min(fgas_bins)  - dx, np.max(fgas_bins) + dx)

        ax[axi].set_xlabel(xlabel)


        ax[axi].text(xpos_label, ypos_label, k)
        plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    plt.minorticks_on()

    if halo_ratio:
        outname = 'fgas_msmh_2D'
    else:
        outname = 'fgas_mstar_2D'

    if log_fgas:
        outname += '_log_fgas'

    outname += outstr

    fig.savefig(outname + '.png')
    plt.close()

    return

def plot_fgas_DSFMS_2d_hist(DSFMS_bins = np.arange(-5,3,0.1), fgas_bins = None,
                            SFR_type = DEFAULT_SFR, log_fgas = True, statistic = 'fraction',
                            datasets = SIM_DATA, figdim = (1,4), outstr = ''):
    """
    Plot the gas fraction of galaxies as a function of the distance to the SFMS
    using the SFR_type as given ('1Gyr' defualt) in a 2D histogram showing the
    fraction of galaxies in  a given pixel.

    statistic sets what statistic to compute for gas fraction. Default is fraction.
    """
    if fgas_bins is None:

        if log_fgas:
            fgas_bins = np.linspace(-4, 0, np.size(DSFMS_bins))
        else:
            fgas_bins = np.linspace(0.0, 1.0, np.size(DSFMS_bins))

    if (figdim == (1,4)) and (len(datasets) > figdim[1]):
        figdim = (1, len(datasets))

    fig, ax = plt.subplots(figdim[0],figdim[1], sharey=True, sharex=True)  # change hard coding
    fig.set_size_inches(figdim[1]*6,figdim[0]*6)

    if figdim == (2,2):
        ax_indexes = [(0,0),(0,1),(1,0),(1,1)]
    else:
        ax_indexes = np.arange(figdim[1])

    stat_name = statistic
    if statistic == 'count':
        cbar_label = r'log(Number of Galaxies)'
        vmin = 0
        vmax = 3
    elif statistic == 'fraction':
        cbar_label = r'log(Fraction of Galaxies)'
        # may need to log the fraction
        vmin =  -4.0
        vmax =  -1
        stat_name = 'count'
        cmap = 'viridis'
    elif statistic == 'median':
        if log_fgas:
            cbar_label = r'log(f$_{\rm gas}$)'
            vmin   = -2.5
            vmax   = 0
            cmap   = 'plasma'
        else:
            cbar_label = r'f$_{\rm gas}$'
            vmin = 0.0
            vmax = 1.0
            cmap = "PRGn"

    stat_name == 'statistic'
    if statistic == 'fraction':
        stat_name = 'count'

    for axi, k in zip( ax_indexes, datasets):
        DSFMS = data[k]['D_SFMS_' + SFR_type]
        fgas  = data[k]['fgas']

        if log_fgas:
            fgas   = np.log10(fgas)

        ngal = np.size(fgas)
        if statistic == 'count' or statistic == 'fraction':
            stat_field = np.ones(ngal)
        else:
            stat_field = fgas

        stat, x_edge, y_edge, binnum = binned_statistic_2d(DSFMS,fgas, stat_field, statistic = stat_name,
                                                             bins = (DSFMS_bins,fgas_bins))

        if statistic == 'fraction':
            fraction = stat / (1.0 * ngal)
            fraction[fraction <= 0] = LOGZERO
            fraction[fraction >  0] = np.log10(fraction[fraction > 0])
            fraction = np.log10(stat / (1.0 * ngal))
            plot_val = fraction
        else:
            plot_val = stat


        xmesh, ymesh = np.meshgrid(x_edge, y_edge)

        img1 = ax[axi].pcolormesh(xmesh, ymesh, plot_val.T,
                                  cmap = cmap, vmin = vmin, vmax = vmax)

        #ax[axi].set_xlabel(r'log( f$_{\rm gas}$)')
        ax[axi].set_xlabel(r'log( DSFMS)')
        if axi == ax_indexes[0]:
            if log_fgas:
                #ax[axi].set_ylabel(' log(DSFMS)')
                ax[axi].set_ylabel(r'log( f$_{\rm gas}$)')
            else:
                ax[axi].set_ylabel(r'f$_{\rm gas}$')
                #ax[axi].set_ylabel(' log(DSFMS)')
        elif axi == ax_indexes[-1]:
            divider = make_axes_locatable(ax[axi])
            cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(img1, cax=cax1, label = cbar_label)

        ax[axi].set_xlim(DSFMS_bins[0], DSFMS_bins[-1])
        ax[axi].set_ylim(fgas_bins[0], fgas_bins[-1])
        #ax[axi].text(8.2, 1.5, k)

        plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    plt.minorticks_on()

    outname = 'fgas_DSFMS_2dhist_' + stat_name

    if log_fgas:
        outname += '_log_fgas'

    outname += outstr

    fig.savefig(outname + '.png')
    plt.close()

    return

def plot_SFMS_2D_hist(mstar_bins = np.arange(8.0,12.6,0.1), sfr_bins = np.arange(-4, 2, 0.1),
                      remove_zero = True, SFR_type = DEFAULT_SFR, log_fgas = False, statistic = 'median',
                      show_fit = True, datasets = SIM_DATA, figdim = (1,4), outstr = ''):
    """
    Plot the SFMS fits for each simulation individually as panels. For simplicity,
    will just be plotting median and IQR for now, but should move to doing
    contours / points.

    statistic sets what statistic to compute for gas fraction. Default is median. This can be set
    to "fraction", for example, to plot density of points in a bin
    """

    if (figdim == (1,4)) and (len(datasets) > figdim[1]):
        figdim = (1, len(datasets))

    fig, ax = plt.subplots(figdim[0],figdim[1], sharey=True, sharex=True)  # change hard coding
    fig.set_size_inches(figdim[1]*6,figdim[0]*6)

    if figdim == (2,2):
        ax_indexes = [(0,0),(0,1),(1,0),(1,1)]
    else:
        ax_indexes = np.arange(figdim[1])

    stat_name = statistic
    if statistic == 'count':
        cbar_label = r'log(Number of Galaxies)'
        vmin = 0
        vmax = 3
    elif statistic == 'fraction':
        cbar_label = r'log(Fraction of Galaxies)'
        # may need to log the fraction
        vmin =  -4.0
        vmax =  -1
        stat_name = 'count'
        cmap = 'viridis'
    elif statistic == 'median':
        if log_fgas:
            cbar_label = r'log(f$_{\rm gas}$)'
            vmin   = -2.5
            vmax   = 0
            cmap   = 'plasma'
        else:
            cbar_label = r'f$_{\rm gas}$'
            vmin = 0.0
            vmax = 1.0
            cmap = "PRGn"

    stat_name == 'statistic'
    if statistic == 'fraction':
        stat_name = 'count'


    for axi, k in zip( ax_indexes, datasets):
        SFR   = data[k]['log_SFR_' + SFR_type]
        Mstar = data[k]['log_Mstar']
        fgas  = data[k]['fgas']

        # filter
        select = SFR > LOGZERO
        SFR    = SFR[select]
        Mstar  = Mstar[select]

        if log_fgas:
            fgas   = np.log10(fgas[select])

        else:
            fgas   = fgas[select]

        ngal = np.size(SFR)
        if statistic == 'count' or statistic == 'fraction':
            stat_field = np.ones(ngal)
        else:
            stat_field = fgas



        stat, x_edge, y_edge, binnum = binned_statistic_2d(Mstar, SFR, stat_field, statistic = stat_name,
                                                             bins = (mstar_bins, sfr_bins))

        if statistic == 'fraction':
            fraction = stat / (1.0 * ngal)
            fraction[fraction <= 0] = LOGZERO
            fraction[fraction >  0] = np.log10(fraction[fraction > 0])
            fraction = np.log10(stat / (1.0 * ngal))
            plot_val = fraction
        else:
            plot_val = stat


        xmesh, ymesh = np.meshgrid(x_edge, y_edge)

        img1 = ax[axi].pcolormesh(xmesh, ymesh, plot_val.T,
                                  cmap = cmap, vmin = vmin, vmax = vmax)

        if axi == ax_indexes[0]:
            ax[axi].set_ylabel(r'log( SFR (M$_{\odot}$ yr$^{-1}$)')
        elif axi == ax_indexes[-1]:
            divider = make_axes_locatable(ax[axi])
            cax1 = divider.append_axes('right', size = '5%', pad = 0.05)

            fig.colorbar(img1, cax=cax1, label = cbar_label)


        if show_fit:
            try:
                ax[axi].plot(data[k]['SFMS_fit_' + SFR_type][0], data[k]['SFMS_fit_' + SFR_type][1],
                             lw = line_width, ls = '--', color = 'black')
            except:
                print "I'm a lazy programmer"

        ax[axi].set_xlim(7.8, 12.5)
        ax[axi].set_ylim(-3.5, 2)
        ax[axi].text(8.2, 1.5, k)

        ax[axi].set_xlabel(r'log( M$_{*}$ [M$_{\odot}$])')
        plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    plt.minorticks_on()

    outname = 'SFMS_fits_2dhist_' + stat_name

    if log_fgas:
        outname += '_log_fgas'

    outname += outstr

    fig.savefig(outname + '.png')
    plt.close()

    return


def plot_gas_mass_stellar_mass(method = 'scatter', include_range = None,
                               mstar_bins = MSTAR_BINS, remove_zero = False,
                               rhalf = None, observational_limits = None,
                               datasets = ALL_DATA):
    """
    Plot the M_HI vs. M_star relationship for all data. When M_HI is not available,
    M_cold is used. The default behavior is to plot a scatter plot (which looks terrible),
    but can also plot the median / average with shaded regions showing the IQR / standard
    deviation.



    method   :  string. 'scatter' or 'binned', latter bins data with the actual values
                determined by include_range. Default: 'scatter'
    include_range : When method = 'binned', select either 'IQR' or 'std' to plot the
                    median with IQR shading or average with standard deviation shading.
                    The default behavior, if include_range is left as None, is to
                    do median with IQR shading. Default : None.
    mstar_bins  : bins in logspace to use to bin by stellar mass
    remove_zero : filter out galaxies with NO gas mass first
    rhalf       : Show quantities within "1" or "2" half light radii (must be 1 or 2).
                  Default shows within virial radius. Default: None

    observational_limits : apply observational limits on f_gas for better comparison
                           to observational work. Options are 'Bradford' or None.

    """

    if method == 'binned' and include_range is None:
        include_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '_%iRh'%(rhalf)
    else:
        rhalf_str = ''

    # assert that datasets contain the right information
    #     if they do not, do not plot them
#    for k in datasets:
#        key_check = all([ ('log_Mcold' + rhalf_str in datasets[k].keys()) or ('log_MHI' + rhalf_str in datasets[k].keys()),
#                          ('fgas' in datasets[k].keys())])
#        if not key_check:


    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        for k in datasets:


            x = data[k]['log_Mstar' + rhalf_str]
            if 'log_MHI' in data[k].keys():
                y = data[k]['log_MHI' + rhalf_str]
            else:
                y = data[k]['log_Mcold' + rhalf_str]

            if remove_zero:
                x = x[y>-2]
                y = y[y>-2]

            if observational_limits == 'Bradford':
                cut = fgas_limits(x, data[k]['fgas'])
                x   = x[cut]
                y   = y[cut]

            ax.scatter(x, y, s = point_size, lw = 0, label = k, color = colors[k], alpha = 0.5)
    elif method == 'binned':
        def _compute_and_plot(axis, xdata, ydata, fgas, xbins, label, *args, **kwargs):

            # assume ydata is logged - compute stats on un-logged data
            ydata = 10.0**(ydata)
            if remove_zero:
                xdata = xdata[ydata>0.01] # ONLY those with gas
                ydata = ydata[ydata>0.01]

            if observational_limits == 'Bradford':
                cut = fgas_limits(xdata, fgas)
                xdata   = xdata[cut]
                ydata   = ydata[cut]

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)
            axis.scatter(xdata[scatter_select], np.log10(ydata[scatter_select]), s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            # scatter plot points that don't have proper statistics


            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                axis.plot(x, np.log10(average), lw = line_width, label = label, *args, **kwargs)
                #print label, x, average
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0
                axis.plot(x, np.log10(median), lw = line_width, label = label, *args, **kwargs)
                #print label, x, median

            if not (fill_low is None):
                if 'color' in kwargs.keys():
                    facecolor = kwargs['color']
                else:
                    facecolor = 'black'
                axis.fill_between(x, np.log10(fill_low), np.log10(fill_up), facecolor = facecolor,
                                interpolate = True, lw = line_width, alpha = 0.25, *args, **kwargs)

            return

        # plot each data set
        for k in datasets:

            if not (rhalf is None):

                if ( not data[k].has_key('log_Mstar' + rhalf_str)):
                    print 'missing Rhalf info for ', k, ' continuing'
                    continue

            x = data[k]['log_Mstar' + rhalf_str]

            #if k == 'EAGLE' and rhalf_str == '':
            #    y = data[k]['log_MHI_MH2' + rhalf_str]
            #    y = data[k]['log_Mcold']

            if 'log_MHI' + rhalf_str in data[k].keys():
                y = data[k]['log_MHI' + rhalf_str]
            else:

                if ( not data[k].has_key('log_Mcold' + rhalf_str)):
                    print 'missing Rhalf info for ', k,' continuing'
                    continue

                y = data[k]['log_Mcold' + rhalf_str]

            # function assumes that y is logged
            _compute_and_plot(ax, x, y, data[k]['fgas'], mstar_bins, k, color = colors[k])

        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'log(Average M$_{\rm gas}$[M$_{\odot}$])' + rhalf_label
            ylabel += r'with STD'
        else:
            ylabel = r'log(Median M$_{\rm gas}$ [M$_{\odot}$])' + rhalf_label
            ylabel += r'with IQR'
        ax.set_ylabel(ylabel)


    # add aditional observational lines:
    if not (observational_limits is None):
        ax.legend(loc = 'upper left', ncol = 2)
    elif False: # old method --- ignore these for now
    #not (observational_limits is None):

        ax.plot(brown_15['log_Mstar'], brown_15['log_MHI'], color = colors['brown15'], lw = 3)
        ax.scatter(brown_15['log_Mstar'], brown_15['log_MHI'], color = colors['brown15'],
                marker = 's', s = 3*point_size, label = 'Brown-15')

        ax.plot(cantinella_13['log_Mstar'], cantinella_13['log_MHI'], color = colors['catinella13'], lw = 3)
        ax.scatter(cantinella_13['log_Mstar'], cantinella_13['log_MHI'], color = colors['catinella13'],
                   marker = 'o', s = 3*point_size, label = 'Catinella-13')

        ax.legend(loc = 'upper left', ncol=2)
    else:
    # axis labels and legend
        ax.legend(loc = 'best')
    ax.set_xlabel(r'log( M$_{\rm *}$ / M$_{\odot}$)')
    plt.minorticks_on()
    ax.set_ylim(7, 11.4)
    ax.set_xlim(7, 12.75)
    if method == 'scatter':
        outname = 'MHI_Mstar_scatter'
    elif method == 'binned':
        outname = 'MHI_Mstar_binned'


    outname += rhalf_str
    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'IQR':
        outname += '_IQR'
    elif include_range == 'std':
        outname += '_std'

    if remove_zero:
        outname += '_nonzero_gas'
        ax.set_title(r"Galaxies with (M$_{\rm gas} > 0$)")

    if observational_limits == 'Bradford':
        outname += '_bradford_fgas_cut'
        ax.set_title(r'With Bradford 2015 f$_{\rm gas}$ Cut')
    else:
        ax.set_title(r'Including M$_{\rm gas} = 0$')

    ax.plot([4,14],[4,14], lw = 0.5*line_width, ls = '--', color = 'black')
    plt.tight_layout()

    fig.savefig(_output_dir + outname + '.png')
    plt.close()
    return


def plot_fgas_histograms(mass_bins = np.array([8, 9, 10, 11]),
                         fgas_bins = np.arange(0,1.1,0.05) , norm = 'fraction',
                         log_fgas = False, datasets = ALL_DATA):
    """
    Constructs a panel plot of histograms showing f_gas distributions, with panels
    selecting galaxies over a given stellar mass range.

    mass_bins   :    stellar mass bins to break data into; sets number of panels.
                     last number is used to make a final, open-ended bin for all
                     galaxies above that mass. This is not true for lower mass bin
    fgas_bins   :    plots give histograms as function of f_gas. The bins for the
                     histograms.
    norm        :    Histogram normalization. Options are 'fraction' or None. 'fraction'
                     gives number in bin / total number in panel. None gives number in bin

    log_fgas    :    plot linear or log f_gas bins. User must supply own f_gas bins if
                     True.
    """

    _coldict = {4 : (2,3), 5 : (2,3)}
    nrow, ncol = _coldict[np.size(mass_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _mass_bins = np.zeros(np.size(mass_bins)+1)
    _mass_bins[:-1] = mass_bins
    _mass_bins[-1] = np.inf
    mass_bins = 1.0 * _mass_bins
    axi,axj = 0,0

    def _compute_and_plot(x, y, label, axis, *args, **kwargs):

        select = (x > mass_bins[ibin]) * (x<mass_bins[ibin+1])

        hist, bins  = np.histogram( y, bins = fgas_bins)
        A = 1.0
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])

        plot_histogram(axis, bins, hist * A, label = label, lw = line_width, *args, **kwargs)
        return

    def _plot_all(fgas, label, axis, *args, **kwargs):

        hist, bins  = np.histogram(fgas, bins = fgas_bins)
        A = 1.0
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])

        plot_histogram(axis, bins, hist * A, label = label, lw = line_width, *args, **kwargs)
        return


    for ibin in np.arange(np.size(mass_bins) - 1):
        axind = (axi, axj)

        for k in datasets:
            x = data[k]['log_Mstar']
            #if 'log_MHI' in data[k].keys():
            #    y = data[k]['log_MHI']
            #else:
            #    y = data[k]['log_Mcold']
            if log_fgas:
                y = data[k]['log_fgas']
            else:
                y = data[k]['fgas']

            _compute_and_plot(x, y, k, ax[axind], color = colors[k])

        if log_fgas:
            ax[axind].set_xlabel(r'log(f$_{\rm gas}$)')
            ax[axind].set_xlim(np.min(fgas_bins),np.max(fgas_bins))

        else:
            ax[axind].set_xlabel(r'f$_{\rm gas}$')
            ax[axind].set_xlim(0.0,1.0)
            ylim = ax[axind].get_ylim()
            ax[axind].set_ylim(0.0, ylim[1])
        ax[axind].set_ylim(0.0,0.9)

        if norm == 'fraction':
            ax[axind].set_ylabel(r'Fraction')
        else:
            ax[axind].set_ylabel(r'Number')
        ax[axind].minorticks_on()

        if ibin < np.size(mass_bins) - 2:
            ax[axind].set_title(r'%i < log(M$_{\rm *}$ [M$_{\odot}]$) < %i'%(int(mass_bins[ibin]),int(mass_bins[ibin+1])))
        else:
            ax[axind].set_title(r'log(M$_{\rm *}$ [M$_{\odot}$]) > %i'%(int(mass_bins[ibin])))


        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    for k in datasets:
        if log_fgas:
            x = data[k]['log_fgas'] # np.log10(data[k]['fgas'][data[k]['fgas'] > 0])
        else:
            x = data[k]['fgas']
        _plot_all(x, k, ax[(axi,axj)], color = colors[k])

    ax[(axi,axj)].set_title(r'All Galaxies')
    ax[(axi,axj)].set_ylabel(ax[(0,0)].get_ylabel())
    ax[(axi,axj)].set_xlabel(ax[(0,0)].get_xlabel())
    ax[(axi,axj)].set_xlim(ax[(0,0)].get_xlim())
    ax[(axi,axj)].set_ylim(ax[(0,0)].get_ylim())

    ax[(0,0)].legend(loc='best')
    fig.set_size_inches(ncol*5, nrow*5)
    plt.tight_layout()
    plt.minorticks_on()

    outname = 'fgas_mstar_histogram_panel.png'
    if log_fgas:
        outname = 'log_' + outname
    fig.savefig(_output_dir + outname)
    plt.close()
    return

def plot_fgas_DSFMS(method   = 'binned', include_range = 'IQR', DSFMS_bins = np.arange(-5,3,0.2),
                    log_fgas = False,    rhalf = None, remove_zero = False,
                    datasets = ALL_DATA, SFR_type = DEFAULT_SFR, single_SFMS = False,
                    recalculate_SFMS = False, *args, **kwargs):
    """
    Plot gas fraction as a function of the distance to the star forming main sequence
    of galaxies. This is difficult to do in reality, as it requires first defining the SFMS
    which may or may not be uniform across simulations.

    single_sfms : boolean, optional
        If True, use the distance to SFMS provided for all data

    """

    if recalculate_SFMS:
        print "This will recalculate the pre-computed SFMS and distances to each, "
        print "but is not yet implemented. When it is, args and kwargs will be passed "
        print "to the fitting routine."
        raise NotImplementedError
#
#
#

    if log_fgas:
        remove_zero = True
    if method == 'binned' and include_range is None:
        include_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '_%iRh'%(rhalf)
    else:
        rhalf_str = ''

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'binned':
        def _compute_and_plot(axis, xdata, ydata, fgas, xbins, label, *args, **kwargs):

            # assume ydata is logged - compute stats on un-logged data
            # ydata = 10.0**(ydata)
            if remove_zero:
                xdata = xdata[ydata > 1.0E-8] # ONLY those with gas
                ydata = ydata[ydata > 1.0E-8]

            #if observational_limits == 'Bradford':
            #    cut = fgas_limits(xdata, fgas)
            #    xdata   = xdata[cut]
            #    ydata   = ydata[cut]

            if log_fgas:
                yscatter = np.log10(ydata)
                yscatter[ydata == 0.0] = LOGZERO
            else:
                yscatter = ydata

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)
            axis.scatter(xdata[scatter_select], yscatter[scatter_select], s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            # scatter plot points that don't have proper statistics

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                axis.plot(x, np.log10(average), lw = line_width, label = label, *args, **kwargs)
                #print label, x, average
            elif include_range == 'IQR':

                fill_up = Q3
                fill_low = Q1

                if log_fgas: # and False:
                    fill_low[fill_low <= 0] = 1.0E-8
                    fill_up = np.log10(fill_up)
                    fill_low = np.log10(fill_low)
                    median   = np.log10(median)
                else:
                    fill_low[fill_low < 0] = 0

                axis.plot(x, median, lw = line_width, label = label, *args, **kwargs)

            if not (fill_low is None):
                if 'color' in kwargs.keys():
                    facecolor = kwargs['color']
                else:
                    facecolor = 'black'

                axis.fill_between(x, fill_low, fill_up, facecolor = facecolor,
                                interpolate = True, lw = line_width, alpha = 0.25, *args, **kwargs)

            return

        # plot each data set
        for k in datasets:

            if not (rhalf is None):

                if ( not data[k].has_key('log_Mstar' + rhalf_str)):
                    print 'missing Rhalf info for ', k, ' continuing'
                    continue

            x = data[k]['D_SFMS_' + SFR_type]                # distance is in dex
            if 'fgas' + rhalf_str in data[k].keys():
                y = data[k]['fgas' + rhalf_str]

            # remove those with distance_SFMS < LOGZERO9
            select = x > LOGZERO                        # hard coded flag

            _compute_and_plot(ax, x[select], y[select], y[select], DSFMS_bins, k, color = colors[k])


        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'log(Average M$_{\rm gas}$ [M$_{\odot}$])' + rhalf_label
            ylabel += r'with STD'
        else:
            if log_fgas:
                ylabel = r'log(Median f$_{\rm gas}$)' + rhalf_label
            else:
                ylabel = r'Median f$_{\rm gas}$' + rhalf_label
            ylabel += r'with IQR'
        ax.set_ylabel(ylabel)


    # axis labels and legend
    ax.legend(loc = 'best')
    ax.set_xlabel(r'Distance to SFMS (dex)')
    plt.minorticks_on()

    if log_fgas:
        ax.set_ylim(-4, 0)
    else:
        ax.set_ylim(0.0, 1.0)

    ax.set_xlim(np.min(DSFMS_bins), np.max(DSFMS_bins))

#    if method == 'scatter':
#        outname = 'MHI_Mstar_scatter'
    if method == 'binned':
        outname = 'fgas_DSFMS_binned'

    outname += rhalf_str
    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'IQR':
        outname += '_IQR'
    elif include_range == 'std':
        outname += '_std'

    if log_fgas:
        outname += '_logged'

    ax.plot([0,0], ax.get_ylim(), lw = 0.5*line_width, ls = '--', color = 'black')
#   shade in width of SFMS here

    plt.tight_layout()

    fig.savefig(_output_dir + outname + '.png')

    plt.close()
    return

def plot_fgas_DSFMS_panel(method   = 'binned', include_range = 'IQR', DSFMS_bins = np.arange(-5,3,0.2),
                          mass_bins = [8,9,10,11],
                          log_fgas = False,    rhalf = None, remove_zero = False,
                          datasets = ALL_DATA, SFR_type = DEFAULT_SFR, single_SFMS = False,
                          recalculate_SFMS = False, *args, **kwargs):
    """
    Plot gas fraction as a function of the distance to the star forming main sequence
    of galaxies. This is difficult to do in reality, as it requires first defining the SFMS
    which may or may not be uniform across simulations.

    single_sfms : boolean, optional
        If True, use the distance to SFMS provided for all data

    """

    if recalculate_SFMS:
        print "This will recalculate the pre-computed SFMS and distances to each, "
        print "but is not yet implemented. When it is, args and kwargs will be passed "
        print "to the fitting routine."
        raise NotImplementedError
#
#
#

    if log_fgas:
        remove_zero = True
    if method == 'binned' and include_range is None:
        include_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '_%iRh'%(rhalf)
    else:
        rhalf_str = ''

    _coldict = {4 : (2,3), 5 : (2,3)}
    nrow, ncol = _coldict[np.size(mass_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _mass_bins = np.zeros(np.size(mass_bins)+1)
    _mass_bins[:-1] = mass_bins
    _mass_bins[-1] = np.inf
    mass_bins = 1.0 * _mass_bins
    axi,axj = 0,0

    fig, ax = plt.subplots(nrow,ncol)
    fig.set_size_inches(5*ncol,5*nrow)

    if method == 'binned':
        def _compute_and_plot(axis, xdata, ydata, fgas, xbins, label, *args, **kwargs):

            # assume ydata is logged - compute stats on un-logged data
            # ydata = 10.0**(ydata)
            if remove_zero:
                xdata = xdata[ydata > 1.0E-8] # ONLY those with gas
                ydata = ydata[ydata > 1.0E-8]

            #if observational_limits == 'Bradford':
            #    cut = fgas_limits(xdata, fgas)
            #    xdata   = xdata[cut]
            #    ydata   = ydata[cut]

            if log_fgas:
                yscatter = np.log10(ydata)
            else:
                yscatter = ydata

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)
            axis.scatter(xdata[scatter_select], yscatter[scatter_select], s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            # scatter plot points that don't have proper statistics

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                axis.plot(x, np.log10(average), lw = line_width, label = label, *args, **kwargs)
                #print label, x, average
            elif include_range == 'IQR':

                fill_up = Q3
                fill_low = Q1

                if log_fgas: # and False:
                    fill_low[fill_low <= 0] = 1.0E-8
                    fill_up = np.log10(fill_up)
                    fill_low = np.log10(fill_low)
                    median   = np.log10(median)
                else:
                    fill_low[fill_low < 0] = 0

                axis.plot(x, median, lw = line_width, label = label, *args, **kwargs)

            if not (fill_low is None):
                if 'color' in kwargs.keys():
                    facecolor = kwargs['color']
                else:
                    facecolor = 'black'

                axis.fill_between(x, fill_low, fill_up, facecolor = facecolor,
                                interpolate = True, lw = line_width, alpha = 0.25, *args, **kwargs)

            return

        # plot each data set
        for k in datasets:

            if not (rhalf is None):

                if ( not data[k].has_key('log_Mstar' + rhalf_str)):
                    print 'missing Rhalf info for ', k, ' continuing'
                    continue

            x = data[k]['D_SFMS_' + SFR_type]                # distance is in dex
            mstar = data[k]['log_Mstar']
            if 'fgas' + rhalf_str in data[k].keys():
                y = data[k]['fgas' + rhalf_str]

            # remove those with distance_SFMS < LOGZERO9
            select = x > LOGZERO                        # hard coded flag

            axi, axj = 0, 0
            for i in np.arange(np.size(mass_bins)-1):
                x_select = x[ select * (mstar >= mass_bins[i]) * (mstar < mass_bins[i+1])]
                y_select = y[ select * (mstar >= mass_bins[i]) * (mstar < mass_bins[i+1])]

                _compute_and_plot(ax[(axi,axj)], x_select, y_select, y_select, DSFMS_bins, k, color = colors[k])
                ax[(axi,axj)].set_title(r'%.1f < log(M$_{*}$) < %.1f'%(mass_bins[i], mass_bins[i+1]))
                # iterate axi axj
                axj = axj + 1
                if axj >= ncol:
                    axj = 0
                    axi = axi + 1

            # plot all data
            _compute_and_plot(ax[(axi,axj)], x[select], y[select],y[select], DSFMS_bins, k, color = colors[k])
            ax[(axi,axj)].set_title('All Galaxies')

        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'log(Average M$_{\rm gas}$ [M$_{\odot}$])' + rhalf_label
            ylabel += r'with STD'
        else:
            if log_fgas:
                ylabel = r'log(Median f$_{\rm gas}$)' + rhalf_label
            else:
                ylabel = r'Median f$_{\rm gas}$' + rhalf_label
            ylabel += r'with IQR'



    # axis labels and legend
    ax[(0,2)].legend(loc = 'best')

    plt.minorticks_on()

    for axind in [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]:
        ax[axind].set_xlabel(r'Distance to SFMS (dex)')
        if log_fgas:
            ax[axind].set_ylim(-4, 0)
        else:
            ax[axind].set_ylim(0.0, 1.0)

        ax[axind].set_xlim(np.min(DSFMS_bins), np.max(DSFMS_bins))
        ax[axind].plot([0,0], ax[axind].get_ylim(), lw = 0.5*line_width, ls = '--', color = 'black')
        ax[axind].set_ylabel(ylabel)

#    if method == 'scatter':
#        outname = 'MHI_Mstar_scatter'
    if method == 'binned':
        outname = 'fgas_DSFMS_panel_plot'

    outname += rhalf_str
    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'IQR':
        outname += '_IQR'
    elif include_range == 'std':
        outname += '_std'

    if log_fgas:
        outname += '_logged'

    plt.tight_layout()

    fig.savefig(_output_dir + outname + '.png')

    plt.close()
    return

def plot_no_gas(mstar_bins = MSTAR_BINS, datasets = ALL_DATA):

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)


    def _compute_and_plot(ax, x, all_x, bins, label, *args, **kwargs):

        hist, bins = np.histogram(x, bins = bins)
        full_hist, bins = np.histogram(all_x, bins = bins)

        hist = hist / ( 1.0 * full_hist)

#        hist = hist
        plot_histogram(ax, bins, hist, label = label, **kwargs)


        return

    for k in datasets:
        all_Mstar = data[k]['log_Mstar']
        fgas  = data[k]['fgas']

        gasless = fgas <= 1.0E-10

        Mstar     = all_Mstar[gasless]
        N_gasless = np.size(Mstar)

        print '---------------', k, N_gasless

        label = k + r'  (%0.1f %%)'%(100.0 * N_gasless / (1.0 * np.size(all_Mstar)))
        _compute_and_plot(ax, Mstar, all_Mstar, mstar_bins, label,
                          lw = line_width, color = colors[k], ls = '-')

    ax.set_xlabel(r'log(M$_{*}$) [M$_{\odot}$])')
    ax.set_ylabel(r'Fraction Per Bin')
    ax.semilogy()

    ax.set_xlim(8.0, 12.0)
    ax.set_ylim(1.0E-4, 0.6)

    plt.minorticks_on()
    ax.legend(loc='best')
    plt.tight_layout()

    fig.savefig('gasless_fraction_mstar.png')
    plt.close()
    return



def plot_fgas_mstar(method = 'scatter', include_range = None,
                    mstar_bins = MSTAR_BINS, log_fgas = False, rhalf = None,
                    remove_zero = None,
                    observational_limits = None,
                    datasets = ALL_DATA,
                    sSFR_cut = None):
    """
    Plot fgas as function of stellar mass for all data. Default behavior is to
    make a scatter plot (which looks terrible). Alternatively, can plot median
    in bins of mstar with options to include shaded regions showing either the
    IQR in that bin or the standard deviation in that bin. This functionality can
    be expanded in the future.

    method : 'scatter' or 'median'. Latter plots a line where data is binned using
             mstar_bins
    mstar_bins : stellar mass (horizontal axis) bins if using 'median' method

    include_range : If None, 'binned' method shows only median or average line. If "IQR" or
                    "std", also shades in the IQR or std (respectively) of the
                    data in the given stellar mass bin
    """

    if method == 'binned' and include_range is None:
        inlcude_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '_%iRh'%(rhalf)
    else:
        rhalf_str = ''

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        for k in datasets:
            x = data[k]['log_Mstar']
            y = data[k]['fgas']

            if remove_zero:
                x = x[data[k]['log_MHI'] >-2]
                y = y[data[k]['log_MHI'] >-2]

            if observational_limits == 'Bradford':
                cut = fgas_limits(x, data[k]['fgas'])
                x   = x[cut]
                y   = y[cut]

            if log_fgas:
                y = np.log10(y)

            ax.scatter(x, y, s = point_size, lw = 0, label = k, color = colors[k], alpha = 0.5)

        ax.set_ylabel(r'f$_{\rm gas}$')
    elif method == 'binned':

        def _compute_and_plot(axis,xdata,ydata, mgas, SFR,
                              xbins, label, *args, **kwargs):
            # helper generic function to compute and then plot the data with fills

            log_mgas = 1.0 * mgas
            log_mgas[mgas==0.0] = LOGZERO

            cut = (xdata == xdata)

            if remove_zero:
                cut *= (log_mgas > -2)

            if observational_limits == 'Bradford':
                cut *= (fgas_limits(xdata, ydata))

            if not (sSFR_cut is None):
                cut *= ((sSFR >= sSFR_cut[0]) * (sSFR < sSFR_cut[1]))

            xdata = xdata[cut] * 1.0
            ydata = ydata[cut] * 1.0

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)

            if log_fgas:
                yscatter = np.log10(ydata[scatter_select])
            else:
                yscatter = ydata[scatter_select]

            axis.scatter(xdata[scatter_select], yscatter, s = point_size, **kwargs)

            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)


            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                if True:
                    fill_low[fill_low < 0] = 0
                if log_fgas: # and False:
                    fill_up  = np.log10(fill_up)
                    fill_low[fill_low == 0.0] = 1.0E-10 * fill_up[fill_low == 0.0] # need to do somehing here
                    fill_low = np.log10(fill_low)
                    average  = np.log10(average)

                axis.plot(x, average, lw = line_width, label = label, *args, **kwargs)
                #print label, x, average
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                if True:
                    fill_low[fill_low < 0] = 0

                # print label, np.min(ydata), np.max(ydata), np.min(Q1), np.max(Q3), np.size(ydata[ydata<1.0E-4])
                # print label, Q1
                # print label, Q3

                if log_fgas: # and False:
                    fill_up = np.log10(fill_up)
                    fill_low = np.log10(fill_low)
                    median   = np.log10(median)

                axis.plot(x, median, lw = line_width, label = label, *args, **kwargs)
                #print label, x, median

            if not (fill_low is None):
                if 'color' in kwargs.keys():
                    facecolor = kwargs['color']
                else:
                    facecolor = 'black'

                axis.fill_between(x, fill_low, fill_up, facecolor = facecolor,
                                interpolate = True, lw = line_width, alpha = 0.25, *args, **kwargs)

            return

        # plot each data source
        for k in datasets:

            if not (rhalf is None):
                if ( not data[k].has_key('fgas' + rhalf_str)):
                    print 'plot_fgas_mstar: missing Rhalf info for ', k, ' continuing'
                    continue

            x    = data[k]['log_Mstar' + rhalf_str]
            y    = data[k]['fgas' + rhalf_str]
            if 'log_MHI' + rhalf_str in data[k].keys():
                mgas = data[k]['log_MHI' + rhalf_str]
            else:

                if ( not data[k].has_key('log_Mcold' + rhalf_str)):
                    print 'plot_fgas_mstar: missing Rhalf info for ', k,' continuing'
                    continue

                mgas = data[k]['log_Mcold' + rhalf_str]

            if 'log_sSFR_1Gyr' in data[k].keys():
                sSFR = data[k]['log_sSFR_1Gyr']
            elif not (sSFR_cut is None):
                continue # skip if no data available on SFR and we make SFR cut

            # send off to plotter routine - mgas is for making cuts
            _compute_and_plot(ax, x, y, mgas, sSFR, mstar_bins, k, color = colors[k])


        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'Average f$_{\rm gas}$ ' + rhalf_label
            ylabel += r' with STD'
        else:
            ylabel = r'Median f$_{\rm gas}$ ' + rhalf_label
            ylabel += r' with IQR'
        ax.set_ylabel(ylabel)
        # end plot median and fill between


    legend_col = 1
    if not (observational_limits is None):

        fgas = brown_15['fgas']
        if log_fgas:
            fgas = np.log10(fgas)

        ax.plot(brown_15['log_Mstar'], fgas, color = colors['brown15'], lw = 3)
        ax.scatter(brown_15['log_Mstar'], fgas, color = colors['brown15'],
                marker = 's', s = 3*point_size, label = 'Brown-15')

        fgas = cantinella_13['fgas']
        if log_fgas:
             fgas = np.log10(fgas)

        ax.plot(cantinella_13['log_Mstar'], fgas, color = colors['catinella13'], lw = 3)
        ax.scatter(cantinella_13['log_Mstar'], fgas, color = colors['catinella13'],
                   marker = 'o', s = 3*point_size, label = 'Catinella-13')
        legend_col = 2

    # axis labels
    if log_fgas:
        ax.legend(loc='lower left', ncol = legend_col)
    else:
        ax.legend(loc='upper right', ncol = legend_col)

    ax.set_xlabel(r'log( M$_{\rm *}$ / M$_{\odot}$)')
    plt.minorticks_on()
    if log_fgas:
        if not (sSFR_cut is None):
            ax.set_ylim(-8,0)
        else:
            ax.set_ylim(-3, 0)
    else:
        ax.set_ylim(0,1.0)
    ax.set_xlim(7,12.75)


    # set output filename and save
    if method == 'scatter':
        outname = 'fgas_mstar_scatter'
    elif method == 'binned':
        outname = 'fgas_mstar_binned'

    outname += rhalf_str
    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'std':
        outname = outname + '_std'
    elif include_range == 'IQR':
        outname = outname + '_IQR'

    if log_fgas:
        outname += '_logged'

    if remove_zero:
        outname += '_nonzero_gas'
        ax.set_title(r"Galaxies with (M$_{\rm gas} > 0$)")

    if observational_limits == 'Bradford':
        outname += '_bradford_fgas_cut'
        ax.set_title(r'With Bradford 2015 f$_{\rm gas}$ Cut')
    else:
        ax.set_title(r'Including M$_{\rm gas} = 0$')

    if not (sSFR_cut is None):
        if (sSFR_cut[0] == -np.inf) and (sSFR_cut[1] < -15):
            outname += '_NO_SF_'
        else:
            outname += '_sSFR_cut_'

    plt.tight_layout()
    fig.savefig(_output_dir + outname + '.png')
    plt.close()
    return

def plot_fgas_ssfr(method = 'scatter', include_range = None,
                   ssfr_bins = np.arange(-15,-8.9,0.2),
                   ssfr_type = DEFAULT_SFR, plot_zeros = False, remove_zeros = False, log_fgas = False,
                   datasets = ALL_DATA, rhalf = None, observational_limits = None, extra_label = '',
                   add_observations = False):
    """
    Plot fgas as function of sSFR

    method : 'scatter' or 'median'. Latter plots a line where data is binned using
             mstar_bins
    mstar_bins : stellar mass (horizontal axis) bins if using 'median' method

    include_range : If None, 'median' method shows only the median or averageline. If "IQR" or
                    "std", also shades in the IQR or std (respectively) of the
                    data in the given stellar mass bin.
    """

    if (include_range is None) and method == 'binned':
        include_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '%iRh'%(rhalf)
    else:
        rhalf_str = ''

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        for k in datasets:
            x = data[k]['sSFR_' + ssfr_type]
            y = data[k]['fgas']

            ax.scatter( np.log10(x[x>0]), y[x>0], s = point_size, lw = 0, label = k, alpha = 0.5, color = colors[k])

        ax.set_ylabel(r'f$_{\rm gas}$')

    elif method == 'binned':

        def _compute_and_plot(axis, input_xdata, ydata, log_Mgas, xbins, label, *args, **kwargs):

            selection = input_xdata > 0
            if log_fgas:
                selection = selection * (log_Mgas > -2)

            xdata = np.log10(input_xdata[selection]) # log the ssfr's
            ydata = ydata[selection]

            scatter_select, line_select = _select_scatter_points(xdata, xbins)

            if log_fgas:
                yscatter = np.log10( ydata[scatter_select])
            else:
                yscatter = ydata[scatter_select]

            axis.scatter(xdata[scatter_select], yscatter, s = point_size, **kwargs)

            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0

                if log_fgas:
                    fill_up  = np.log10(fill_up)
                    fill_low[fill_low == 0.0] = 1.0E-10 * fill_up[fill_low == 0.0] # need to do so$
                    fill_low = np.log10(fill_low)
                    average  = np.log10(average)

                axis.plot(x, average, lw = line_width, label = label, *args, **kwargs)
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0

                if log_fgas: # and False:
                    fill_up = np.log10(fill_up)
                    fill_low = np.log10(fill_low)
                    median   = np.log10(median)

                axis.plot(x, median, lw = line_width, label = label, *args, **kwargs)

            if not (fill_low is None):
                axis.fill_between(x, fill_low, fill_up, facecolor = colors[label],
                                interpolate = True, lw = line_width, alpha = 0.25)

            # plot median and relevant error bars for galaxies with zero ssfr
            zero_ssfr = input_xdata[input_xdata == 0.0]
            if plot_zeros and np.size(zero_ssfr) >= 3:
                zero_ssfr_fgas = ydata[input_xdata == 0.0]
                median  = np.median(zero_ssfr_fgas)
                average = np.average(zero_ssfr_fgas)

                ye1 = np.percentile(zero_ssfr_fgas,75) - median
                ye2 = median - np.percentile(zero_ssfr_fgas,25)

                yerr = np.array([[ ye2, ye1]]).T

                axis.scatter(np.min(ssfr_bins) - 0.5, median, s = point_size*2, *args, **kwargs)

                axis.errorbar(np.min(ssfr_bins) - 0.5, median, yerr = yerr, markersize = point_size*4,
                            color = colors[label], elinewidth = 0.75 * line_width, capsize = 4)

            return

        # plot each data source
        for k in datasets:
            #print k
            x = data[k]['sSFR_' + ssfr_type]
            y = data[k]['fgas']

            if 'log_MHI' in data[k].keys():
                log_Mgas = data[k]['log_MHI']
            else:
                log_Mgas = data[k]['log_Mcold']

            _compute_and_plot(ax, x, y, log_Mgas, ssfr_bins, k, color = colors[k])

        if add_observations:
            x = data['Bradford2015']['sSFR_Halpha']
            y = data['Bradford2015']['fgas']
            log_Mgas = data['Bradford2015']['log_MHI']

            _compute_and_plot(ax, x, y, log_Mgas, ssfr_bins, 'Bradford2015', color = colors['Bradford2015'])

        if include_range == 'std':
            ylabel = r'Average f$_{\rm gas}$'
            ylabel += r' with STD'
        else:
            ylabel = r'Median f$_{\rm gas}$'
            ylabel += r' with IQR'
        ax.set_ylabel(ylabel)
        # end plot median and fill between

    # axis labels
    ax.legend(loc='best')
    ax.set_xlabel(r'log( sSFR yr$^{-1}$)')
    plt.minorticks_on()

    if log_fgas:
        ax.set_ylim(-4, 0)
    else:
        ax.set_ylim(0,1.0)

    if plot_zeros:
        ax.set_xlim(np.min(ssfr_bins)-1, np.max(ssfr_bins))
    else:
        ax.set_xlim(np.min(ssfr_bins), np.max(ssfr_bins))

    # set output filename and save
    if method == 'scatter':
        outname = 'fgas_ssfr_' + ssfr_type + '_scatter'
    elif method == 'binned':
        outname = 'fgas_ssfr_' + ssfr_type + '_binned'
    outname += '_' + rhalf_str
    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'IQR':
        outname += 'IQR'
    elif include_range == 'std':
        outname += 'std'

    if remove_zeros:
        outname += '_nonzero_gas'
    ax.set_title(r"M$_{\rm gas} > 0$ and sSFR > 0")
    #else:
    #    ax.set_title(r"Including Galaxies with M$_{\rm gas} = 0$")

    if log_fgas:
        outname += '_logged'

    if observational_limits == 'Bradford':
        outname += '_bradford_fgas_cut'
        ax.set_title(r'With Observational f$_{\rm gas}$ Cut')
    #else:
    #    ax.set_title(r'Including M$_{\rm gas} = 0$')

    if add_observations:
        outname += '_obs'

    plt.tight_layout()
    fig.savefig(_output_dir + outname + extra_label + '.png')
    plt.close()
    return

def plot_fgas_ssfr_histograms(ssfr_bins = np.array([-20,-13,-12,-11,-10,-9]),
                              fgas_bins = None, norm = 'fraction',
                              sSFR_type = DEFAULT_SFR, sSFR_alternate = None, log_fgas = False, datasets = ALL_DATA):
    """
    Make a panel plot showing distributions of f_gas for each simulation, grouped
    by sSFR.

    ssfr_bins  :  sets bins to panel. First panel will show distribution where NO SF is
                  measured in the time range. Second and on will show all galaxies in the
                  range. For example, with default values, panels are: No SF, -20 to -13,
                  -13 to -12, -12 to -11, -11 to -10, and -10 to -9.

    fgas_bins  : We are plotting distributions of f_gas. This sets the horizontal axis bins.

    norm       : Currently only 'fraction' is supported:
                    'fraction' : fraction of galaxies in bin relative to total number in panel

    sSFR_type  : time range over which sSFR is measured. Currently does nothing, but will be
                 used to select between 1Gyr and 10Myr sSFR's
    """
    _coldict = {5 : (2,3), 6 : (2,3)}  # use this to set key for number of ssfr bins and number of panels

    if fgas_bins is None:
        if log_fgas:
            fgas_bins = np.arange(-4, 0.1, 0.2)
        else:
            fgas_bins = np.arange(0, 1.10, 0.1)

    nrow, ncol = _coldict[np.size(ssfr_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _ssfr_bins = np.zeros(np.size(ssfr_bins)+1)
    _ssfr_bins[1:] = ssfr_bins

    _ssfr_bins[0] = -np.inf
    ssfr_bins = 1.0 * _ssfr_bins
    axi,axj = 0,0

    def _compute_and_plot(axis, fgas, ssfr, label, ibin, *args, **kwargs):

        if log_fgas:
            #ssfr = ssfr[fgas > 0]
            temp = np.log10(fgas)
            temp[fgas == 0.0] = LOGZERO
            fgas = 1.0*temp

        # generic helper function to compute data and plot
        if ibin > 1:
            logssfr  = np.log10( ssfr[ssfr>0] )
            select = (logssfr > ssfr_bins[ibin-1]) * (logssfr<=ssfr_bins[ibin])
            fgas_data = fgas[ssfr>0]
            fgas_data = fgas_data[select]
        else:
            select    = ssfr == 0.0
            fgas_data = fgas[select] # fgas
            if not log_fgas:
                print "NO SF in last    type    label    size_full_data   size fgas < 0.1     fraction f_gas > 0    median f_gas   max f_gas"                
                print 'NO SF in last ', sSFR_type, label, np.size(fgas_data), np.size(fgas_data[fgas_data < 0.1]), np.size(fgas_data)/(1.0*np.size(fgas)), np.median(fgas_data), np.max(fgas_data)

        hist, bins  = np.histogram( fgas_data, bins = fgas_bins)
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, 1.0*np.sum(hist)])
        plot_histogram(axis, bins, hist * A, label = label, lw = line_width, *args, **kwargs)

        return

    for ibin in np.arange(1, np.size(ssfr_bins)): # loop over each bin / panel
        axind = (axi, axj)

        for k in datasets:
            if not (('sSFR_' + sSFR_type) in data[k].keys()):
                y = data[k]['sSFR_' + sSFR_alternate]
            else:
                y = data[k]['sSFR_' + sSFR_type]

            x = data[k]['fgas']

            _compute_and_plot(ax[axind], x, y, k, ibin, color = colors[k])

        # set appropriate axis labels and plot styles
        if ibin == 1:
            ax[axind].set_title(r'No SF in last ' + sSFR_type)
        elif ibin == 2:
            ax[axind].set_title(r'log(sSFR yr$^{-1}$) < -13')
        else:
            ax[axind].set_title(r'%i < log(sSFR yr$^{-1}$) < %i'%( int(ssfr_bins[ibin-1]), int(ssfr_bins[ibin])))

        ax[axind].set_xlabel(r'f$_{\rm gas}$')
        if norm == 'fraction':
            ax[axind].set_ylabel(r'Fraction of Total')
        else:
            ax[axind].set_ylabel(r'Number')

        if log_fgas:
            ax[axind].set_xlim(-4,0)
        else:
            ax[axind].set_xlim(0.0, 1.0)
        ax[axind].set_ylim(0, 0.5)

        ax[axind].minorticks_on()

        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1
        # end loop over panels

    ax[(0,0)].legend(loc='best')
    fig.set_size_inches(ncol*5, nrow*5)
    plt.tight_layout()
    plt.minorticks_on()
    if log_fgas:
        outname = 'log_'
    else:
        outname = ''

    fig.savefig(_output_dir + outname + 'fgas_ssfr_histogram_panel_' + sSFR_type + '.png')
    for a1 in ax:
        for a in a1:
            a.set_ylim(0.0,0.2)
    fig.savefig(_output_dir + outname + 'fgas_ssfr_histogram_panel_' + sSFR_type + '_zoom.png')
    return


def plot_sfr_mgas(method = 'binned', include_range = 'IQR',
                  mgas_bins = np.arange(3,12,0.2),
                  sfr_type = DEFAULT_SFR, plot_zeros = False, remove_zeros = False, log_fgas = False,
                  datasets = ALL_DATA, rhalf = None, observational_limits = None, extra_label = '' ):
    """
    Plot SFR vs Mgas
    """

    if (include_range is None) and method == 'binned':
        include_range = 'IQR'

    if not (rhalf is None):
        rhalf_str = '%iRh'%(rhalf)
    else:
        rhalf_str = ''

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'binned':

        def _compute_and_plot(axis, xdata, ydata, log_Mgas, xbins, label, *args, **kwargs):

            selection = ydata > LOGZERO

            xdata = xdata[selection]
            ydata = ydata[selection]

            scatter_select, line_select = _select_scatter_points(xdata, xbins)

            axis.scatter(xdata[scatter_select], ydata[scatter_select], s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , 10.0**(ydata[line_select]), xbins)

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0

                axis.plot(x, average, lw = line_width, label = label, *args, **kwargs)
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0

                axis.plot(x, np.log10(median), lw = line_width, label = label, *args, **kwargs)

            if not (fill_low is None):
                axis.fill_between(x, np.log10(fill_low), np.log10(fill_up), facecolor = colors[label],
                                interpolate = True, lw = line_width, alpha = 0.25)


            return

        # plot each data source
        for k in datasets:
            #print k
            if 'log_MHI' in data[k].keys():
                x = data[k]['log_MHI']
            else:
                x = data[k]['log_Mcold']

            try:
                y = data[k]['log_SFR_' + sfr_type]
            except:
                y = data[k]['log_SFR_' + sfr_type]

            mstar = data[k]['log_Mstar']

            x = x[mstar < 9.5]
            y = y[mstar < 9.5]

            _compute_and_plot(ax, x, y, x, mgas_bins, k, color = colors[k])

        if include_range == 'std':
            ylabel = r'Average f$_{\rm gas}$'
            ylabel += r' with STD'
        else:
            ylabel = r'log(Median SFR [M$_{\odot}$ yr$^{-1}$])'
            ylabel += r' with IQR'
        ax.set_ylabel(ylabel)
        # end plot median and fill between

    # axis labels
    ax.legend(loc='best')
    ax.set_xlabel(r'log(M$_{\rm gas}$ [M$_{\odot}$])')
    plt.minorticks_on()

    if log_fgas:
        ax.set_ylim(4, 12)
    else:
        ax.set_ylim(-6, 2)


    # set output filename and save

    if method == 'binned':
        outname = 'sfr_mgas_' + sfr_type + '_binned'

    outname += '_' + rhalf_str

    if not (rhalf is None):
        outname = rhalf_str.strip('_') + '/' + outname

    if include_range == 'IQR':
        outname += 'IQR'
    elif include_range == 'std':
        outname += 'std'


    # ax.set_title(r"M$_{\rm gas} > 0$ and sSFR > 0")
    # else:
    #    ax.set_title(r"Including Galaxies with M$_{\rm gas} = 0$")

#    if observational_limits == 'Bradford':
#        outname += '_bradford_fgas_cut'
#        ax.set_title(r'With Observational f$_{\rm gas}$ Cut')
    #else:
    #    ax.set_title(r'Including M$_{\rm gas} = 0$')

    plt.tight_layout()
    fig.savefig(_output_dir + outname + extra_label + '.png')
    plt.close()
    return


if __name__ == "__main__":

    plot_fgas_DSFMS_2d_hist(datasets = LSIM_DATA, log_fgas = True)
    plot_fgas_DSFMS_2d_hist(datasets = LSIM_DATA, log_fgas = False)


    plot_no_gas(datasets = LSIM_DATA)

    plot_halo_stellar_2D_hist(datasets=LSIM_DATA)
    plot_halo_stellar_2D_hist(datasets=LSIM_DATA, statistic = 'median')

    plot_sfr_mgas(datasets=['Illustris','SCSAM','EAGLE', 'MUFASA', 'Romulus25']) #Brooks

    plot_fgas_mstar_2D_hist(datasets = LSIM_DATA)
    plot_fgas_mstar_2D_hist(datasets = LSIM_DATA, log_fgas = True)
    plot_fgas_mstar_2D_hist(datasets = LSIM_DATA, log_fgas = True, sSFR_cut = (-np.inf,-20))
    plot_fgas_mstar_2D_hist(datasets = LSIM_DATA, log_fgas = True, halo_ratio = True)



    #
    # Plot various versions of the SFMS fit. Include 2D hist with gas fraction
    #
    plot_SFMS(datasets = LSIM_DATA)
    plot_SFMS_2D_hist(datasets = LSIM_DATA)
    plot_SFMS_2D_hist(datasets = LSIM_DATA, statistic = 'fraction', show_fit = False)
    plot_SFMS_2D_hist(datasets = LSIM_DATA, log_fgas = True)

    #plot_fgas_mstar_2D_hist(datasets = ['Illustris','SCSAM','EAGLE', 'MUFASA'],figdim=(1,4), log_fgas = True)
    plot_fgas_mstar_2D_hist(datasets = ['Bradford2015','xGASS'],figdim=(1,2), log_fgas = True, outstr = '_obs')

    #plot_SFMS_2D_hist(datasets = ['Illustris','SCSAM','EAGLE', 'MUFASA'],figdim=(1,4))
    plot_SFMS_2D_hist(datasets = ['Bradford2015','xGASS'],figdim=(1,2), outstr = '_obs')


    #
    # plot gas fraction as function of distance to SFMS
    #
    plot_fgas_DSFMS_panel(method = 'binned', include_range = 'IQR', datasets = LSIM_DATA, remove_zeros = False)

    plot_fgas_DSFMS_panel(method = 'binned', include_range = 'IQR', datasets = LSIM_DATA, remove_zeros = False, log_fgas = True)

    plot_fgas_DSFMS(method = 'binned', include_range = 'IQR', datasets = LSIM_DATA, remove_zeros = False)
    plot_fgas_DSFMS(method = 'binned', include_range = 'IQR', datasets = LSIM_DATA, log_fgas = True, remove_zeros=True)

    #
    # gas mass vs stellar mass
    #
    GM_SM_data = SIM_DATA #+ ['MUFASA']
    GM_SM_data_obs = ALL_DATA #+ ['MUFASA']
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = GM_SM_data)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = GM_SM_data, rhalf = 1)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = GM_SM_data, rhalf = 2)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', observational_limits = 'Bradford', datasets = GM_SM_data_obs + ['xGASS'])

    #
    # gas fraction vs sSFR and variations
    #
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, add_observations=True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, remove_zeros = True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True, remove_zeros=True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True, remove_zeros=True, add_observations=True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True, remove_zeros=True, ssfr_bins = np.arange(-20.5,-8.9,0.2), extra_label = '_extended')
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True)

    #
    # make all gas fraction vs stellar mass plots and variants
    #
    fgas_MStar_data     = SIM_DATA #+ ['MUFASA']
    fgas_MStar_data_obs = ALL_DATA #+ ['MUFASA']
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = LSIM_DATA,
                    log_fgas = True, sSFR_cut = (-np.inf, -90), remove_zero = True)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = fgas_MStar_data)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', observational_limits = 'Bradford', datasets = fgas_MStar_data_obs + ['xGASS'])
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = fgas_MStar_data)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, observational_limits = 'Bradford', datasets = fgas_MStar_data_obs + ['xGASS'])
        # with spatial information
    # plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = SIM_DATA, rhalf = 1)
    # plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 1)
    # plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = SIM_DATA, rhalf = 2)
    # plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 2)

    # do above as histograms
    #
    plot_fgas_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE', 'MUFASA' , 'Romulus25', 'Bradford2015'])
    plot_fgas_histograms(fgas_bins = np.arange(-3, 0.01, 0.1) , log_fgas = True, datasets = ['Illustris', 'SCSAM', 'EAGLE',  'MUFASA' , 'Romulus25', 'Bradford2015'])

    #
    # gas fraction as a function of ssfr histograms
    #
    plot_fgas_ssfr_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE', 'MUFASA' , 'Romulus25'], log_fgas = True)
    plot_fgas_ssfr_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE', 'MUFASA' , 'Romulus25' ], log_fgas = False)
    plot_fgas_ssfr_histograms(sSFR_type = '100Myr', sSFR_alternate = '100Myr', datasets =['Illustris', 'SCSAM', 'EAGLE', 'MUFASA' , 'Romulus25' ])
    plot_fgas_ssfr_histograms(sSFR_type = '100Myr', sSFR_alternate = '100Myr', datasets =['Illustris', 'SCSAM', 'EAGLE', 'MUFASA', 'Romulus25'  ], log_fgas = True)
