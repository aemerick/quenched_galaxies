import numpy as np

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
