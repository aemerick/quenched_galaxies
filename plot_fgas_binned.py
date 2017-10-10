from definitions import * # -- local - predefines
from galaxy_analysis.plot.plot_styles import plot_histogram

import numpy as np
import matplotlib.pyplot as plt


def compute_fgas(mstar, mgas, log = True):
    if log:
        fgas = 10**(mgas) / (10**(mstar) + 10**(mgas))
    else:
        fgas = mgas / (mstar + mgas)
    return fgas

def plot_fgas_histograms(mass_bins = np.array([8, 9, 10, 11, 12]),
                         fgas_bins = np.arange(0,1.05,0.025) , norm = 'fraction',
                         log_fgas = False):
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

    _coldict = {5 : (2,3)}
    nrow, ncol = _coldict[np.size(mass_bins)]

    fig, ax = plt.subplots(nrow,ncol)
    
    # now loop through and plot data
    _mass_bins = np.zeros(np.size(mass_bins)+1)
    _mass_bins[:-1] = mass_bins
    _mass_bins[-1] = np.inf
    mass_bins = 1.0 * _mass_bins
    axi,axj = 0,0

    def _compute_and_plot(x, y, label):

        select = (x > mass_bins[ibin]) * (x<mass_bins[ibin+1])

        fgas = compute_fgas(x[select],y[select])
        if log_fgas:
            fgas   = np.log10(fgas[fgas>0])

        hist, bins  = np.histogram( fgas, bins = fgas_bins)
        A = 1.0
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])

        plot_histogram(ax[axind], bins, hist * A, label = label, lw = line_width)


    for ibin in np.arange(np.size(mass_bins) - 1):
        axind = (axi, axj)

        _compute_and_plot(data['illustris']['log_Mstar'],
                             data['illustris']['log_MHI'], 'illustris')
        _compute_and_plot(data['SAM']['mstar'], data['SAM']['mcold'], 'SAM')
        _compute_and_plot(data['MUFASA']['log_Mstar'], data['MUFASA']['log_Mcold'],'MUFASA')

        if log_fgas:
            ax[axind].set_xlabel(r'log(f$_{\rm gas}$)')
            ax[axind].set_xlim(np.min(fgas_bins),np.max(fgas_bins))
        else:
            ax[axind].set_xlabel(r'f$_{\rm gas}$')
            ax[axind].set_xlim(0.0,1.0)

        if norm == 'fraction':
            ax[axind].set_ylabel(r'Fraction of Total')
        else:
            ax[axind].set_ylabel(r'Number')
        ax[axind].minorticks_on()

        if ibin < np.size(mass_bins) - 2:
            ax[axind].set_title(r'%i < log(M$_{\rm *}$ [M$_{\odot}]$) < %i'%(int(mass_bins[ibin]),int(mass_bins[ibin+1])))
        else:
            ax[axind].set_title(r'log(M$_{\rm *}$ [M$_{\odot}$]) > %i'%(int(mass_bins[ibin])))

        ylim = ax[axind].get_ylim()
        ax[axind].set_ylim(0.0, ylim[1])

        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    ax[(0,0)].legend(loc='best')
    fig.set_size_inches(ncol*5, nrow*5)
    plt.tight_layout()
    plt.minorticks_on()
 
    outname = 'fgas_mstar_histogram_panel.png'
    if log_fgas:
        outname = 'log_' + outname
    fig.savefig(outname)
    return

def plot_fgas_mstar(method = 'scatter', include_range = None,
                    mstar_bins = np.arange(8.0,13.1,0.1) ):
    """
    Plot fgas as function of stellar mass for all data. Default behavior is to
    make a scatter plot (which looks terrible). Alternatively, can plot median
    in bins of mstar with options to include shaded regions showing either the 
    IQR in that bin or the standard deviation in that bin. This functionality can
    be expanded in the future. 

    method : 'scatter' or 'median'. Latter plots a line where data is binned using
             mstar_bins
    mstar_bins : stellar mass (horizontal axis) bins if using 'median' method

    include_range : If None, 'median' method shows only the median line. If "IQR" or
                    "std", also shades in the IQR or std (respectively) of the
                    data in the given stellar mass bin  
    """

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        ax.scatter( data['illustris']['log_Mstar'],
                    compute_fgas(data['illustris']['log_Mstar'],data['illustris']['log_MHI']),
                    s = point_size, lw = 0, label = 'illustris', color=colors['illustris'], alpha=0.5)
        ax.scatter( data['SAM']['mcold'],
                    compute_fgas(data['SAM']['mstar'],data['SAM']['mcold']),
                    s = point_size, lw = 0, label = 'SAM', color=colors['SAM'], alpha = 0.5)
        ax.scatter( data['MUFASA']['log_Mstar'],
                    compute_fgas(data['MUFASA']['log_Mstar'],data['MUFASA']['log_Mcold']),
                    s = point_size, lw = 0, label = 'MUFASA', color=colors['MUFASA'], alpha = 0.5)

    elif method == 'median':
        # use bins to bin data and compute average in bin
        def _compute_median(x, y, xbins):
            median = np.ones(np.size(xbins)-1) * -1 # leave neg as flag
            Q1 = np.zeros(np.size(median)); Q3 = np.zeros(np.size(median))
            std = np.zeros(np.size(median))
            for i in np.arange(np.size(xbins)-1):
                y_select  = y[(x>=xbins[i])*(x<xbins[i+1])] 
                if np.size(y_select) >= 3:
                    median[i] = np.median( y_select )
                    Q1[i]     = np.percentile( y_select, 25)
                    Q3[i]     = np.percentile( y_select, 75)
                    std[i]    = np.std(y_select)

            select = median > 0
            centers = (xbins[1:] + xbins[:-1])*0.5
            return centers[select], median[select], std[select], Q1[select], Q3[select]

        def _compute_and_plot(axis,xdata,ydata,xbins, label):

            x,median,std,Q1,Q3 = _compute_median(xdata,ydata,xbins)
            ax.plot(x,median, lw = line_width, color = colors[label], label = label)
            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0

            if not (fill_low is None):
                axis.fill_between(x, fill_low, fill_up, facecolor = colors[label],
                                interpolate = True, lw = line_width, alpha = 0.25)
                 
            return

        _compute_and_plot(ax, data['illustris']['log_Mstar'], data['illustris']['fgas'],
                              mstar_bins, 'illustris')
        _compute_and_plot(ax, data['SAM']['mstar'], data['SAM']['fgas'], mstar_bins, "SAM")
        _compute_and_plot(ax, data['MUFASA']['log_Mstar'], data['MUFASA']['fgas'],
                              mstar_bins, 'MUFASA')

        # end plot median and fill between

    # axis labels
    ax.legend(loc='best')
    ax.set_ylabel(r'f$_{\rm gas}$')
    ax.set_xlabel(r'log( M$_{\rm *}$ / M$_{\odot}$)')
    plt.minorticks_on()
    ax.set_ylim(0,1.0)
    ax.set_xlim(8,12.75)

    # set output filename and save
    if method == 'scatter':
        outname = 'fgas_mstar_scatter'
    elif method == 'median':
        outname = 'fgas_mstar_median'
    if include_range == 'std':
        outname = outname + '_std'
    elif include_range == 'IQR':
        outname = outname + '_IQR'
    fig.savefig(outname + '.png')

    return
 
def plot_fgas_ssfr_histograms(ssfr_bins = np.array([-20,-13,-12,-11,-10,-9]),
                              fgas_bins = np.arange(0,1.05,0.05) , norm = 'fraction',
                              sSFR_type = '1Gyr'):
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

    nrow, ncol = _coldict[np.size(ssfr_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _ssfr_bins = np.zeros(np.size(ssfr_bins)+1)
    _ssfr_bins[1:] = ssfr_bins

    _ssfr_bins[0] = -np.inf
    ssfr_bins = 1.0 * _ssfr_bins
    axi,axj = 0,0

    def _compute_and_plot(fgas, ssfr, label):
        
        if ibin > 1:
            logssfr  = np.log10( ssfr[ssfr>0] )
            select = (logssfr > ssfr_bins[ibin-1]) * (logssfr<=ssfr_bins[ibin])
            fgas_data = fgas[ssfr>0]
            fgas_data = fgas_data[select]
        else:
            select    = ssfr == 0.0
            fgas_data = fgas[select]

        hist, bins  = np.histogram( fgas_data, bins = fgas_bins)
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])
        plot_histogram(ax[axind], bins, hist * A, label = label, lw = line_width, color = colors[label])    

        return

    for ibin in np.arange(1, np.size(ssfr_bins)):
        axind = (axi, axj)

        _compute_and_plot(data['illustris']['fgas'], data['illustris']['sSFR_1Gyr'], 'illustris')
#        _compute_and_plot(data['SAM']['fgas'], data['SAM']['sSFR_1Gyr'], 'SAM')
        _compute_and_plot(data['MUFASA']['fgas'], data['MUFASA']['sSFR_1Gyr'], 'MUFASA')

        if ibin == 1:
            ax[axind].set_title(r'No SF in last Gyr')
        elif ibin == 2:
            ax[axind].set_title(r'log(sSFR yr$^{-1}$) < -13')
        else:
            ax[axind].set_title(r'%i < log(sSFR yr$^{-1}$) < %i'%( int(ssfr_bins[ibin-1]), int(ssfr_bins[ibin])))

        ax[axind].set_xlabel(r'f$_{\rm gas}$')
        if norm == 'fraction':
            ax[axind].set_ylabel(r'Fraction of Total')
        else:
            ax[axind].set_ylabel(r'Number')
        ax[axind].set_xlim(0.0, 1.0)
        ax[axind].minorticks_on()

        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    ax[(0,0)].legend(loc='best')
    fig.set_size_inches(ncol*5, nrow*5)
    plt.tight_layout()
    plt.minorticks_on()
    fig.savefig('fgas_ssfr_histogram_panel.png')
    return


if __name__ == "__main__":
    plot_fgas_mstar(method = 'scatter')
    plot_fgas_mstar(method = 'median', include_range = 'std')
    plot_fgas_mstar(method = 'median', include_range = 'IQR')
    plot_fgas_histograms()
    plot_fgas_histograms(fgas_bins = np.arange(-4, 0.1, 0.1) , log_fgas = True)
    plot_fgas_ssfr_histograms()

