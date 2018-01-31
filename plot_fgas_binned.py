from definitions import * # -- local - predefines
from galaxy_analysis.plot.plot_styles import plot_histogram

import numpy as np
import matplotlib.pyplot as plt

from function_definitions import fgas_limits

_output_dir = './plots/'

MSTAR_BINS = np.arange(7.0,13.1,0.2)
ALL_DATA = ['Illustris','SCSAM','MUFASA','EAGLE','Brooks', 'Bradford2015']
SIM_DATA = ['Illustris','SCSAM','MUFASA','EAGLE','Brooks']
OBS_DATA = ['Bradford2015']

SCATTER_THRESHOLD = 10 # threshold to plot scatter poitns instead of median curves

point_size = 40

#
# - really just need to handle this as general, only making median 
#   curves when number in bin is > 10 or something like that
#
scatter_only_dataset = {}
for k in ALL_DATA:
    scatter_only_dataset[k] = False
scatter_only_dataset['Brooks'] = True

def compute_fgas(mstar, mgas, log = True):
    if log:
        fgas = 10**(mgas) / (10**(mstar) + 10**(mgas))
    else:
        fgas = mgas / (mstar + mgas)
    return fgas

def _check_bins(x,xbins):
    N = np.zeros(np.size(xbins) -1)
    for i in np.arange(np.size(xbins)-1):
        N[i] = np.size( x[ (x >= xbins[i]) * (x < xbins[i+1])])
    return N

def _select_scatter_points(x, xbins):
    N_per_bin = _check_bins(x, xbins)

    scatter_select = np.zeros(np.size(x))
    for i, N in enumerate(N_per_bin):
        if N < SCATTER_THRESHOLD:
            scatter_select += ( x >= xbins[i] ) * ( x < xbins[i+1] )
    scatter_select = scatter_select.astype(bool)

    return scatter_select, scatter_select == 0

def _compute_statistics(x, y, xbins):
    flag = -99999999
    # generic function to compute median and statistics
    median = np.ones(np.size(xbins)-1) * flag # leave neg as flag
    average = np.zeros(np.size(median))
    Q1 = np.zeros(np.size(median)); Q3 = np.zeros(np.size(median))
    std = np.zeros(np.size(median)); N = np.zeros(np.size(median))
    for i in np.arange(np.size(xbins)-1):
        y_select  = y[(x>=xbins[i])*(x<xbins[i+1])] 

        if np.size(y_select) >= 1:
            median[i] = np.median( y_select )
            average[i] = np.average( y_select)
            Q1[i]     = np.percentile( y_select, 25)
            Q3[i]     = np.percentile( y_select, 75)
            std[i]    = np.std(y_select)
            N[i]      = np.size( y_select )

    select = median > flag
    centers = (xbins[1:] + xbins[:-1])*0.5
    return centers[select], median[select], std[select], Q1[select], Q3[select], average[select], N[select]



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
            ax.scatter(xdata[scatter_select], np.log10(ydata[scatter_select]), s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            # scatter plot points that don't have proper statistics
            

            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                ax.plot(x, np.log10(average), lw = line_width, label = label, *args, **kwargs)
                #print label, x, average
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0
                ax.plot(x, np.log10(median), lw = line_width, label = label, *args, **kwargs)
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
            if 'log_MHI' + rhalf_str in data[k].keys():
                y = data[k]['log_MHI' + rhalf_str]
            else:

                if ( not data[k].has_key('log_Mcold' + rhalf_str)):
                    print 'missing Rhalf info for ', k,' continuing'
                    continue

                y = data[k]['log_Mcold' + rhalf_str]
 
            #if scatter_only_dataset[k]:
            #    if observational_limits == 'Bradford':
            #        cut = fgas_limits(x, data[k]['fgas'])
            #        x   = x[cut]
            #        y   = y[cut]    

#                ax.scatter(x, y, color = colors[k], lw = 0, s = point_size, label = k, alpha = 1.0)
            #else:
            # function assumes that y is logged
            _compute_and_plot(ax, x, y, data[k]['fgas'], mstar_bins, k, color = colors[k])

        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'log(Average M$_{\rm gas}$/M$_{\odot}$)' + rhalf_label
            ylabel += r'with STD'
        else:
            ylabel = r'log(Median M$_{\rm gas}$/M$_{\odot}$)' + rhalf_label
            ylabel += r'with IQR'
        ax.set_ylabel(ylabel)


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

    if include_range == 'IQR':
        outname += '_IQR'
    elif include_range == 'std':
        outname += '_std'

    if remove_zero:
        outname += '_nonzero_gas'
        ax.set_title(r"Galaxies with (M$_{\rm gas} > 0$)")

    if observational_limits == 'Bradford':
        outname += '_bradford_fgas_cut'
        ax.set_title(r'With Observational f$_{\rm gas}$ Cut')
    else:
        ax.set_title(r'Including M$_{\rm gas} = 0$')

    ax.plot(ax.get_xlim(), ax.get_ylim(), lw = 0.5*line_width, ls = '--', color = 'black')


    fig.savefig(_output_dir + outname + '.png')

    return
    

def plot_fgas_histograms(mass_bins = np.array([8, 9, 10, 11, 12]),
                         fgas_bins = np.arange(0,1.05,0.025) , norm = 'fraction',
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

    _coldict = {5 : (2,3)}
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

        fgas = compute_fgas(x[select],y[select])
        if log_fgas:
            fgas   = np.log10(fgas[fgas>0])

        hist, bins  = np.histogram( fgas, bins = fgas_bins)
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
            if 'log_MHI' in data[k].keys():
                y = data[k]['log_MHI']
            else:
                y = data[k]['log_Mcold']

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
            x = np.log10(data[k]['fgas'][data[k]['fgas'] > 0])
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
                    datasets = ALL_DATA, SFR_type = '1Gyr', single_SFMS = False,
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
            else:
                yscatter = ydata

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)
            ax.scatter(xdata[scatter_select], yscatter[scatter_select], s = point_size, **kwargs)
            x,median,std,Q1,Q3,average, N = _compute_statistics(xdata[line_select] , ydata[line_select], xbins)

            # scatter plot points that don't have proper statistics
            
            fill_low = None ; fill_up = None
            if include_range == 'std':
                fill_up = median + std
                fill_low = median - std
                fill_low[fill_low < 0] = 0
                ax.plot(x, np.log10(average), lw = line_width, label = label, *args, **kwargs)
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

                ax.plot(x, median, lw = line_width, label = label, *args, **kwargs)

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

            # remove those with distance_SFMS < -999
            select = x > -99                        # hard coded flag
         
            _compute_and_plot(ax, x[select], y[select], y[select], DSFMS_bins, k, color = colors[k])


        if not (rhalf is None):
            rhalf_label = r' Interior to %i R$_{*,\rm{half}}$ '%(rhalf)
        else:
            rhalf_label = r' '

        if include_range == 'std':
            ylabel = r'log(Average M$_{\rm gas}$/M$_{\odot}$)' + rhalf_label
            ylabel += r'with STD'
        else:
            ylabel = r'log(Median M$_{\rm gas}$/M$_{\odot}$)' + rhalf_label
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

    if include_range == 'IQR':
        outname += '_IQR'
    elif include_range == 'std':
        outname += '_std'

    if log_fgas:
        outname += '_logged'
    


    ax.plot([0,0], ax.get_ylim(), lw = 0.5*line_width, ls = '--', color = 'black')
#   shade in width of SFMS here

    fig.savefig(_output_dir + outname + '.png')

    plt.close()
    return
    
    

def plot_fgas_mstar(method = 'scatter', include_range = None,
                    mstar_bins = MSTAR_BINS, log_fgas = False, rhalf = None,
                    remove_zero = None,
                    observational_limits = None, datasets = ALL_DATA):
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

        def _compute_and_plot(axis,xdata,ydata, mgas, xbins, label, *args, **kwargs):
            # helper generic function to compute and then plot the data with fills

            if remove_zero:
                xdata = xdata[log_mgas > -2] # ONLY those with gas
                ydata = ydata[log_mgas > -2]

            if observational_limits == 'Bradford':
                cut = fgas_limits(xdata, ydata)
                xdata   = xdata[cut]
                ydata   = ydata[cut]

            # check number of points in bins - plot low number counts as scatter
            scatter_select, line_select = _select_scatter_points(xdata, xbins)

            if log_fgas:
                yscatter = np.log10(ydata[scatter_select])
            else:
                yscatter = ydata[scatter_select]

            ax.scatter(xdata[scatter_select], yscatter, s = point_size, **kwargs)

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

                ax.plot(x, average, lw = line_width, label = label, *args, **kwargs)
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

                ax.plot(x, median, lw = line_width, label = label, *args, **kwargs)
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
            
            # send off to plotter routine - mgas is for making cuts
            _compute_and_plot(ax, x, y, mgas, mstar_bins, k, color = colors[k])


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

    # axis labels
    if log_fgas:
        ax.legend(loc='lower left')
    else:
        ax.legend(loc='upper right')
    ax.set_xlabel(r'log( M$_{\rm *}$ / M$_{\odot}$)')
    plt.minorticks_on()
    if log_fgas:
        ax.set_ylim(-3, 0)
    else:
        ax.set_ylim(0,1.0)
    ax.set_xlim(7,12.75)
    plt.tight_layout()

    # set output filename and save
    if method == 'scatter':
        outname = 'fgas_mstar_scatter'
    elif method == 'binned':
        outname = 'fgas_mstar_binned'

    outname += rhalf_str

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
        ax.set_title(r'With Observational f$_{\rm gas}$ Cut')
    else:
        ax.set_title(r'Including M$_{\rm gas} = 0$')

    fig.savefig(_output_dir + outname + '.png')
    plt.close()
    return

def plot_fgas_ssfr(method = 'scatter', include_range = None,
                   ssfr_bins = np.arange(-15,-8.9,0.2),
                   ssfr_type = '1Gyr', plot_zeros = False, remove_zeros = False, log_fgas = False,
                   datasets = ALL_DATA, rhalf = None, observational_limits = None, extra_label = '' ):
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

            ax.scatter(xdata[scatter_select], yscatter, s = point_size, **kwargs)

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

                ax.plot(x, average, lw = line_width, label = label, *args, **kwargs)
            elif include_range == 'IQR':
                fill_up = Q3
                fill_low = Q1
                fill_low[fill_low < 0] = 0

                if log_fgas: # and False:
                    fill_up = np.log10(fill_up)
                    fill_low = np.log10(fill_low)
                    median   = np.log10(median)

                ax.plot(x, median, lw = line_width, label = label, *args, **kwargs)

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

                ax.scatter(np.min(ssfr_bins) - 0.5, median, s = point_size*2, *args, **kwargs)

                ax.errorbar(np.min(ssfr_bins) - 0.5, median, yerr = yerr, markersize = point_size*4,
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

    fig.savefig(_output_dir + outname + extra_label + '.png')
    plt.close()
    return

def plot_fgas_ssfr_histograms(ssfr_bins = np.array([-20,-13,-12,-11,-10,-9]),
                              fgas_bins = None, norm = 'fraction',
                              sSFR_type = '1Gyr', sSFR_alternate = None, log_fgas = False, datasets = ALL_DATA):
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
            fgas_bins = np.arange(-4, 0.1, 0.1)
        else:
            fgas_bins = np.arange(0, 1.05, 0.05)

    nrow, ncol = _coldict[np.size(ssfr_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _ssfr_bins = np.zeros(np.size(ssfr_bins)+1)
    _ssfr_bins[1:] = ssfr_bins

    _ssfr_bins[0] = -np.inf
    ssfr_bins = 1.0 * _ssfr_bins
    axi,axj = 0,0

    def _compute_and_plot(fgas, ssfr, label, ibin, *args, **kwargs):

        if log_fgas:
            ssfr = ssfr[fgas > 0]
            fgas = np.log10(fgas[fgas>0])
        # generic helper function to compute data and plot
        if ibin > 1:
            logssfr  = np.log10( ssfr[ssfr>0] )
            select = (logssfr > ssfr_bins[ibin-1]) * (logssfr<=ssfr_bins[ibin])
            fgas_data = fgas[ssfr>0]
            fgas_data = fgas_data[select]
        else:
            select    = ssfr == 0.0
            fgas_data = fgas[select]
            print label, np.size(fgas_data)

        hist, bins  = np.histogram( fgas_data, bins = fgas_bins)
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])
        plot_histogram(ax[axind], bins, hist * A, label = label, lw = line_width, *args, **kwargs)

        return 

    for ibin in np.arange(1, np.size(ssfr_bins)): # loop over each bin / panel
        axind = (axi, axj)

        for k in datasets:
            if not (('sSFR_' + sSFR_type) in data[k].keys()):
                y = data[k]['sSFR_' + sSFR_alternate]
            else:
                y = data[k]['sSFR_' + sSFR_type]

            x = data[k]['fgas']

            _compute_and_plot(x, y, k, ibin, color = colors[k])

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


if __name__ == "__main__":

    plot_fgas_DSFMS(method = 'binned', include_range = 'IQR', datasets = ['Illustris','SCSAM','EAGLE','MUFASA'], remove_zeros = False)
    plot_fgas_DSFMS(method = 'binned', include_range = 'IQR', datasets = ['Illustris','SCSAM','EAGLE','MUFASA'], log_fgas = True, remove_zeros=True)

    #plot_gas_mass_stellar_mass(method = 'scatter')
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = SIM_DATA)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 1)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 2)
    plot_gas_mass_stellar_mass(method = 'binned', include_range = 'IQR', observational_limits = 'Bradford')

    # plot_fgas_ssfr(method = 'scatter')
    # plot_fgas_ssfr(method = 'binned', include_range = 'std')
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, remove_zeros = True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True, remove_zeros=True)
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True, remove_zeros=True, ssfr_bins = np.arange(-20.5,-8.9,0.2), extra_label = '_extended')
    plot_fgas_ssfr(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, log_fgas = True)

    #plot_fgas_mstar(method = 'scatter')
#    plot_fgas_mstar(method = 'binned', include_range = 'std')
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = SIM_DATA)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', observational_limits = 'Bradford')
#    plot_fgas_mstar(method = 'binned', include_range = 'std', log_fgas = True)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = SIM_DATA)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, observational_limits = 'Bradford')

    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = SIM_DATA, rhalf = 1)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 1)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', log_fgas = True, datasets = SIM_DATA, rhalf = 2)
    plot_fgas_mstar(method = 'binned', include_range = 'IQR', datasets = SIM_DATA, rhalf = 2)

    plot_fgas_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE','MUFASA', 'Bradford2015'])
    plot_fgas_histograms(fgas_bins = np.arange(-3, 0.01, 0.1) , log_fgas = True, datasets = ['Illustris', 'SCSAM', 'EAGLE', 'MUFASA','Bradford2015'])


    plot_fgas_ssfr_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE','MUFASA'], log_fgas = True)
    plot_fgas_ssfr_histograms(datasets = ['Illustris', 'SCSAM', 'EAGLE','MUFASA'], log_fgas = False)
    plot_fgas_ssfr_histograms(sSFR_type = '10Myr', sSFR_alternate = '20Myr', datasets =['Illustris', 'SCSAM', 'EAGLE','MUFASA'])
    plot_fgas_ssfr_histograms(sSFR_type = '10Myr', sSFR_alternate = '20Myr', datasets =['Illustris', 'SCSAM', 'EAGLE','MUFASA'], log_fgas = True)



