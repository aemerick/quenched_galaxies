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


#        x      = data['illustris']['log_Mstar']
#        y      = data['illustris']['log_MHI']
#        select = (x > mass_bins[ibin]) * (x<mass_bins[ibin+1])
#
#        hist, bins  = np.histogram( compute_fgas(x[select], y[select]), bins = fgas_bins)
#        A = 1.0
#        if norm == 'fraction':
#            A = 1.0 / np.max([1.0, np.sum(hist)])
#        plot_histogram(ax[axind], bins, hist * A, label = 'illustris', lw = line_width)
#        ax[axind].plot(bins[1:], hist, label = 'illustris', drawstyle = 'steps-post', lw = line_width)

#        x      = data['SAM']['mstar']
#        y      = data['SAM']['mcold']
#        select = (x > mass_bins[ibin]) * (x < mass_bins[ibin+1])
#
#        hist, bins = np.histogram( compute_fgas(x[select], y[select]), bins = fgas_bins)
#        A = 1.0
#        if norm == 'fraction':
#            A = 1.0 / np.max([np.sum(hist),1.0])
#        plot_histogram(ax[axind], bins, hist * A, label = 'SAM', lw = line_width)
#        ax[axind].plot(bins, hist, label = 'SAM', drawstyle = 'steps-post', lw = line_width)
#
#        x = data['MUFASA']['log_Mstar']
#        y = data['MUFASA']['log_Mcold']
#        select = (x > mass_bins[ibin]) * (x < mass_bins[ibin+1])
#        hist, bins = np.histogram(compute_fgas(x[select],y[select]), bins = fgas_bins)
#        A = 1.0
#        if norm == 'fraction':
#            A = 1.0 / np.max([np.sum(hist),1.0])
#        plot_histogram(ax[axind], bins, hist * A, label = "MUFASA", lw = line_width)
#
#        x = np.log10( data['Brooks']['Mstar'] )
#        y = np.log10( data['Brooks']['HI_Mass'] )
#        select = (x > mass_bins[ibin]) * (x < mass_bins[ibin+1])
#        hist, bins = np.histogram(compute_fgas(x[select],y[select]), bins = fgas_bins)
#        A = 1.0
#        if norm == 'fraction':
#            A = 1.0 / np.max([np.sum(hist),1.0])
#        plot_histogram(ax[axind], bins, hist * A, label = "Brooks", lw = line_width) 

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


def plot_fgas_ssfr_histograms(ssfr_bins = 10**(np.array([-13,-12,-11,-10,-9,-8])),
                              fgas_bins = np.arange(0,1.05,0.05) , norm = 'fraction'):
    _coldict = {5 : (2,3), 6 : (2,3)}
    nrow, ncol = _coldict[np.size(ssfr_bins)]

    fig, ax = plt.subplots(nrow,ncol)

    # now loop through and plot data
    _ssfr_bins = np.zeros(np.size(ssfr_bins)+1)
    _ssfr_bins[1:] = ssfr_bins

    _ssfr_bins[0] = -np.inf
    ssfr_bins = 1.0 * _ssfr_bins
    axi,axj = 0,0
    for ibin in np.arange(1, np.size(ssfr_bins)):
        axind = (axi, axj)

        x      = data['illustris']['sSFR_1Gyr']
        print x
        y      = data['illustris']['log_MHI']
        select = (x > ssfr_bins[ibin-1]) * (x<ssfr_bins[ibin])

        hist, bins  = np.histogram( compute_fgas(data['illustris']['log_Mstar'][select], 
                                                 y[select]), bins = fgas_bins)
        if norm == 'fraction':
            A = 1.0 / np.max([1.0, np.sum(hist)])
        plot_histogram(ax[axind], bins, hist * A, label = 'illustris', lw = line_width)

#        ax[axind].plot(bins[1:], hist, label = 'illustris', drawstyle = 'steps-post', lw = line_width)

#        x      = data['SAM']['mstar']
#        y      = data['SAM']['mcold']
#        select = (x > ssfr_bins[ibin]) * (x < ssfr_bins[ibin+1])
#
#        hist, bins = np.histogram( compute_fgas(x[select], y[select]), bins = fgas_bins)
#        if norm == 'fraction':
#            A = 1.0 / np.max([np.sum(hist),1.0])
#        plot_histogram(ax[axind], bins, hist * A, label = 'SAM', lw = line_width)

#        ax[axind].plot(bins, hist, label = 'SAM', drawstyle = 'steps-post', lw = line_width)

        ax[axind].set_xlabel(r'f$_{\rm gas}$')
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

    plot_fgas_histograms()
    plot_fgas_histograms(fgas_bins = np.arange(-4, 0.1, 0.1) , log_fgas = True)
    plot_fgas_ssfr_histograms()

