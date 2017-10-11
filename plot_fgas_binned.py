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

    def _compute_and_plot(x, y, label, axis):

        select = (x > mass_bins[ibin]) * (x<mass_bins[ibin+1])

        fgas = compute_fgas(x[select],y[select])
        if log_fgas:
            fgas   = np.log10(fgas[fgas>0])

        hist, bins  = np.histogram( fgas, bins = fgas_bins)
        A = 1.0
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])

        plot_histogram(axis, bins, hist * A, label = label, lw = line_width)
        return

    def _plot_all(fgas, label, axis):

        hist, bins  = np.histogram(fgas, bins = fgas_bins)
        A = 1.0
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])

        plot_histogram(axis, bins, hist * A, label = label, lw = line_width)     
        return


    for ibin in np.arange(np.size(mass_bins) - 1):
        axind = (axi, axj)

        _compute_and_plot(data['illustris']['log_Mstar'],
                             data['illustris']['log_MHI'], 'illustris',ax[axind])
        _compute_and_plot(data['SAM']['mstar'], data['SAM']['mcold'], 'SAM', ax[axind])
        _compute_and_plot(data['MUFASA']['log_Mstar'], data['MUFASA']['log_Mcold'],'MUFASA', ax[axind])
        _compute_and_plot(data['Brooks']['log_Mstar'], data['Brooks']['log_MHI'], 'Brooks', ax[axind])

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

    _plot_all(data['illustris']['fgas'], 'illustris', ax[(axi,axj)])
    _plot_all(data['SAM']['fgas'], 'SAM', ax[(axi,axj)])
    _plot_all(data['MUFASA']['fgas'], "MUFASA", ax[(axi,axj)])
    ax[(axi,axj)].set_title(r'All Galaxies')
    ax[(axi,axj)].set_ylabel(ax[(0,0)].get_ylabel())
    ax[(axi,axj)].set_xlabel(ax[(0,0)].get_xlabel())

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

    include_range : If None, 'binned' method shows only median or average line. If "IQR" or
                    "std", also shades in the IQR or std (respectively) of the
                    data in the given stellar mass bin  
    """

    if method == 'binned' and include_range is None:
        inlcude_range = 'IQR'

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        ax.scatter( data['illustris']['log_Mstar'], data['illustris']['fgas'],
                    s = point_size, lw = 0, label = 'illustris', color=colors['illustris'], alpha=0.5)
        ax.scatter( data['SAM']['mcold'], data['SAM']['fgas'],
                    s = point_size, lw = 0, label = 'SAM', color=colors['SAM'], alpha = 0.5)
        ax.scatter( data['MUFASA']['log_Mstar'], data['MUFASA']['fgas'],
                    s = point_size, lw = 0, label = 'MUFASA', color=colors['MUFASA'], alpha = 0.5)
        ax.scatter( data['Brooks']['log_Mstar'], data['Brooks']['fgas'],
                    s = point_size, lw = 0, label = 'Brooks', color=colors['Brooks'], alpha = 0.5)

        ax.set_ylabel(r'f$_{\rm gas}$')
    elif method == 'binned':
        def _compute_median(x, y, xbins):
            # generic function to compute median and statistics
            median = np.ones(np.size(xbins)-1) * -1 # leave neg as flag
            average = np.zeros(np.size(median))
            Q1 = np.zeros(np.size(median)); Q3 = np.zeros(np.size(median))
            std = np.zeros(np.size(median))
            for i in np.arange(np.size(xbins)-1):
                y_select  = y[(x>=xbins[i])*(x<xbins[i+1])] 
                if np.size(y_select) >= 1:
                    median[i] = np.median( y_select )
                    average[i] = np.average( y_select)
                    Q1[i]     = np.percentile( y_select, 25)
                    Q3[i]     = np.percentile( y_select, 75)
                    std[i]    = np.std(y_select)

            select = median > -1
            centers = (xbins[1:] + xbins[:-1])*0.5
            return centers[select], median[select], std[select], Q1[select], Q3[select], average[select]


        def _compute_and_plot(axis,xdata,ydata,xbins, label):
            # helper generic function to compute and then plot the data with fills
            x,median,std,Q1,Q3,average = _compute_median(xdata,ydata,xbins)

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                ax.plot(x, average, lw = line_width, color = colors[label], label = label)
                #print label, x, average
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0
                ax.plot(x, median, lw = line_width, color = colors[label], label = label)
                #print label, x, median

            if not (fill_low is None):
                axis.fill_between(x, fill_low, fill_up, facecolor = colors[label],
                                interpolate = True, lw = line_width, alpha = 0.25)
                 
            return

        # plot each data source
        _compute_and_plot(ax, data['illustris']['log_Mstar'], data['illustris']['fgas'],
                              mstar_bins, 'illustris')
        _compute_and_plot(ax, data['SAM']['mstar'], data['SAM']['fgas'], mstar_bins, "SAM")
        _compute_and_plot(ax, data['MUFASA']['log_Mstar'], data['MUFASA']['fgas'],
                              mstar_bins, 'MUFASA')
        _compute_and_plot(ax, data['Brooks']['log_Mstar'], data['Brooks']['fgas'], mstar_bins, 'Brooks')
       # print data['Brooks']['log_Mstar']
        #print data['Brooks']['fgas']
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
    ax.set_xlabel(r'log( M$_{\rm *}$ / M$_{\odot}$)')
    plt.minorticks_on()
    ax.set_ylim(0,1.0)
    ax.set_xlim(8,12.75)

    # set output filename and save
    if method == 'scatter':
        outname = 'fgas_mstar_scatter'
    elif method == 'binned':
        outname = 'fgas_mstar_binned'
    if include_range == 'std':
        outname = outname + '_std'
    elif include_range == 'IQR':
        outname = outname + '_IQR'
    fig.savefig(outname + '.png')

    return

def plot_fgas_ssfr(method = 'scatter', include_range = None,
                   ssfr_bins = np.arange(-14,-8.9,0.1),
                   ssfr_type = '1Gyr', plot_zeros = False):
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

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if method == 'scatter':
        x = data['illustris']['sSFR_' + ssfr_type]
        ax.scatter( np.log10(x[x>0]), data['illustris']['fgas'][x>0],
                    s = point_size, lw = 0, label = 'illustris', color=colors['illustris'], alpha=0.5)
        x = data['SAM']['sSFR_' + ssfr_type]
        ax.scatter( np.log10(x[x>0]), data['SAM']['fgas'][x>0],
                    s = point_size, lw = 0, label = 'SAM', color=colors['SAM'], alpha = 0.5)
        x = data['MUFASA']['sSFR_' + ssfr_type]
        ax.scatter( np.log10(x[x>0]), data['MUFASA']['fgas'][x>0],
                    s = point_size, lw = 0, label = 'MUFASA', color=colors['MUFASA'], alpha = 0.5)
        x = data['Brooks']['sSFR_' + ssfr_type]
        ax.scatter( np.log10(x[x>0]), data['Brooks']['fgas'][x>0], s = point_size, lw = 0,
                                label = 'Brooks', color = colors['Brooks'], alpha = 0.5)

        ax.set_ylabel(r'f$_{\rm gas}$')

    elif method == 'binned':
        def _compute_median(x, y, xbins):
            # generic function to compute median and statistics
            median = np.ones(np.size(xbins)-1) * -1 # leave neg as flag
            average = np.zeros(np.size(median))
            Q1 = np.zeros(np.size(median)); Q3 = np.zeros(np.size(median))
            std = np.zeros(np.size(median))
            for i in np.arange(np.size(xbins)-1):
                y_select  = y[(x>=xbins[i])*(x<xbins[i+1])] 
                if np.size(y_select) >= 1:
                    median[i] = np.median( y_select )
                    Q1[i]     = np.percentile( y_select, 25)
                    Q3[i]     = np.percentile( y_select, 75)
                    std[i]    = np.std(y_select)
                    average[i]   = np.average(y_select)

            select = median > -1
            centers = (xbins[1:] + xbins[:-1])*0.5
            return centers[select], median[select], std[select], Q1[select], Q3[select], average[select]


        def _compute_and_plot(axis, input_xdata, ydata,xbins, label):
            xdata = np.log10(input_xdata[input_xdata>0]) # log the ssfr's

            # helper generic function to compute and then plot the data with fills
            x,median,std,Q1,Q3,average = _compute_median(xdata,ydata[input_xdata>0],xbins)

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                ax.plot(x, average, lw = line_width, color = colors[label], label = label)
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0
                ax.plot(x, median, lw = line_width, color = colors[label], label = label)

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
#                print yerr, np.max(zero_ssfr_fgas), np.min(zero_ssfr_fgas)
                ax.scatter(np.min(ssfr_bins) - 0.5, median, s = point_size*2)

                ax.errorbar(np.min(ssfr_bins) - 0.5, median, yerr = yerr, markersize = point_size*4,
                            color = colors[label], elinewidth = 0.75 * line_width, capsize = 4)

#                low = None; up = None
#                if include_range == 'std':
#                    low = np.max([0.0,median - np.std(zero_ssfr)])
#                    up  = np.min([1.0,median + np.std(zero_ssfr)])
#                elif include_range == 'IQR':
#                    low = median - np.percentile( zero_ssfr, 25)
#                    up  = np.percentile( zero_ssfr, 75) - median
#
#                if not (low is None):
#
#                    ax.errorbar(np.min(ssfr_bins) - 0.5, median, yerr=np.array([[low,up]]).T, color = colors[label], 
#                                markersize = point_size, elinewidth = 0.5*line_width, capsize=4, capthick=line_width)
#
#                    ax.errorbar(np.min(ssfr_bins) - 0.5, average, np.array([[low,up]]).T, color = colors[label], 
#                                markersize = point_size, elinewidth = 0.5*line_width, capsize=4, capthick=line_width)
#                    print median, average, low,up
#                else:
#                    ax.scatter(np.min(ssfr_bins) - 0.5, median, color = colors[label], s = point_size)
#                    ax.scatter(np.min(ssfr_bins) - 0.5, average, color = colors[label], s = point_size)
                 
            return

        # plot each data source
        _compute_and_plot(ax, data['illustris']['sSFR_' + ssfr_type], data['illustris']['fgas'],
                              ssfr_bins, 'illustris')
        _compute_and_plot(ax, data['SAM']['sSFR_' + ssfr_type], data['SAM']['fgas'], ssfr_bins, "SAM")
        _compute_and_plot(ax, data['MUFASA']['sSFR_' + ssfr_type], data['MUFASA']['fgas'],
                              ssfr_bins, 'MUFASA')
        _compute_and_plot(ax, data['Brooks']['sSFR_' + ssfr_type], data['Brooks']['fgas'], ssfr_bins, 'Brooks')


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

    def _compute_and_plot(fgas, ssfr, label, ibin):
        # generic helper function to compute data and plot
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

    for ibin in np.arange(1, np.size(ssfr_bins)): # loop over each bin / panel
        axind = (axi, axj)

        if sSFR_type == '10Myr':
            il_ssfr_type = '20Myr'
        else:
            il_ssfr_type = sSFR_type
        _compute_and_plot(data['illustris']['fgas'], data['illustris']['sSFR_' + il_ssfr_type], 'illustris',ibin)
        _compute_and_plot(data['SAM']['fgas'], data['SAM']['sSFR_' + sSFR_type], 'SAM', ibin)
        _compute_and_plot(data['MUFASA']['fgas'], data['MUFASA']['sSFR_' + sSFR_type], 'MUFASA',ibin)
        _compute_and_plot(data['Brooks']['fgas'], data['Brooks']['sSFR_' + il_ssfr_type], 'Brooks', ibin)

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
        ax[axind].set_xlim(0.0, 1.0)
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
    fig.savefig('fgas_ssfr_histogram_panel_' + sSFR_type + '.png')
    for a1 in ax:
        for a in a1:
            a.set_ylim(0.0,0.2)
    fig.savefig('fgas_ssfr_histogram_panel_' + sSFR_type + '_zoom.png')
    return


if __name__ == "__main__":
    plot_fgas_ssfr(method = 'scatter')
    plot_fgas_ssfr(method = 'binned', include_range = 'std')
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR')


    plot_fgas_mstar(method = 'scatter')
    plot_fgas_mstar(method = 'binned', include_range = 'std')
    plot_fgas_mstar(method = 'binned', include_range = 'IQR')


    plot_fgas_histograms()
    plot_fgas_histograms(fgas_bins = np.arange(-4, 0.1, 0.1) , log_fgas = True)


    plot_fgas_ssfr_histograms()
    plot_fgas_ssfr_histograms(sSFR_type = '10Myr')

