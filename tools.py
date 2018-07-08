import numpy as np
import sys, os

from quenched_galaxies import plot_definitions

class NoStdStreams(object):
    def __init__(self,stdout = None, stderr = None):
        self.devnull = open(os.devnull,'w')
        self._stdout = stdout or self.devnull or sys.stdout
        self._stderr = stderr or self.devnull or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()

def check_bins(x,xbins):
    """
    Return number of points in each bin
    """
    N = np.zeros(np.size(xbins) -1)
    for i in np.arange(np.size(xbins)-1):
        N[i] = np.size( x[ (x >= xbins[i]) * (x < xbins[i+1])])
    return N

def select_scatter_points(x, xbins, threshold = 10):
    """
    Provide cut-arrays that can be used to selectively plot points
    that sit in bins with low counts. Also provides selection for
    points that do sit in bins with enough counts to do statistics.
    """
    N_per_bin = check_bins(x, xbins)

    scatter_select = np.zeros(np.size(x))
    for i, N in enumerate(N_per_bin):
        if N < threshold:
            scatter_select += ( x >= xbins[i] ) * ( x < xbins[i+1] )
    scatter_select = scatter_select.astype(bool)

    return scatter_select, scatter_select == 0

def compute_statistics(x, y, xbins, return_dict = False):
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

def bin_and_plot(axis, xdata, ydata, xbins, label = '', n_thresh = 10,
                     remove_zero = False, include_range = 'IQR',
                     observational_limits = None, fgas = None,
                     *args, **kwargs):
    """
    Helper function designed to facilitate plotting median / average curves
    with standard devation / IQR shading in a general way with the option to
    add in a filter to `remove_zero` values, for example, or to inlcude
    observational limits on the gas fraction (`observational_limits` and
    `fgas`).

    xdata and ydata are assumed to be logged and figure is assumed
    to have log axis. This can be adopted and messed with a bit if undesired.

    Input:
    --------

    axis : matplotlib axis object
        The axis to plot the curve on

    xdata : numpy array
        The horizontal axis data for all points. Must be the same size
        as `ydata`.

    ydata : numpy array
       The vertial axis data for all points. Must be the same size as `xdata`.
       THIS MUST BE THE LOGGED (base 10) DATA

    xbins : numpy array
       The bins to bin along in the x-axis. The statistics (median, average,
       standard deviation, IQR, etc.) will be computed in the horizontal axis
       according to these bins. If there are less than 10 galaxies in a bin,
       they are plotted instead as points.

    label : string, optional
        What should the line be called (for legend purposes). Default : ''

    n_thresh : int, optional
        Threshold number of data points below which points will be placed on
        graph instead of the median line or shading. Default : 10

    remove_zero : bool, optional
        Filter out zero values (or log(x) = -99) values. Optional. Default False

    observational_limits: string, optional
        Place an observational limit on the gas fraction if desired. Gas fraction
        must be provided via `fgas`. Current options are "Bradford". Default: None

    fgas : numpy array, optional
        Gas fraction data associated with xdata and ydata. Only used when placing
        observational limits. Default: None

    Additional args and kwargs are passed to the ax.plot calls

    Return:
    --------
    Void
    """

            # assume ydata is logged - compute stats on un-logged data
    ydata = 10.0**(ydata)
    if remove_zero:
        xdata = xdata[ydata>0.0000001]
        ydata = ydata[ydata>0.0000001]

    if observational_limits == 'Bradford':
        if fgas is None:
            print "fgas must be provided if applying observational limit cuts"
            raise ValueError
        cut = fgas_limits(xdata, fgas)
        xdata   = xdata[cut]
        ydata   = ydata[cut]

    # check number of points in bins - plot low number counts as scatter
    scatter_select, line_select = select_scatter_points(xdata, xbins, threshold = n_thresh)
    axis.scatter(xdata[scatter_select], np.log10(ydata[scatter_select]), s = plot_definitions.point_size, **kwargs)
    x,median,std,Q1,Q3,average, N = compute_statistics(xdata[line_select] , ydata[line_select], xbins)

    # scatter plot points that don't have proper statistics

    fill_low = None ; fill_up = None
    if include_range == 'std':
        fill_up = median + std
        fill_low = median - std
        fill_low[fill_low < 0] = 0
        axis.plot(x, np.log10(average), lw = plot_definitions.line_width, label = label, *args, **kwargs)
        #print label, x, average
    elif include_range == 'IQR':
        fill_up = Q3
        fill_low = Q1
        fill_low[fill_low < 0] = 0
        axis.plot(x, np.log10(median), lw = plot_definitions.line_width, label = label, *args, **kwargs)
        #print label, x, median

    #
    if not (fill_low is None):
        if 'color' in kwargs.keys():
            facecolor = kwargs['color']
        else:
            facecolor = 'black'
        axis.fill_between(x, np.log10(fill_low), np.log10(fill_up), facecolor = facecolor,
                        interpolate = True, lw = plot_definitions.line_width, alpha = 0.25, *args, **kwargs)

    return
