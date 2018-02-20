import numpy as np
import math

def randomise_data(hist,sigma):
    perbin = lambda mu: np.random.normal(mu,sigma/math.sqrt(mu))
    randomised = np.vectorize(perbin)
    return randomised(hist)

def make_variations(hist, slope):
#    print "got hist and slope", hist, slope
    if len(hist) == 1:
        hist[0] = 1
    thebins = np.arange(len(hist))
    if slope < 0:
        thebins = thebins[::-1]
    gradient = lambda b, s: b*s
    sloped = np.vectorize(gradient)
    array_up = sloped(thebins, slope)*(1/math.sqrt(hist.sum())) + hist
    array_down = hist - sloped(thebins, slope)*(1/math.sqrt(hist.sum()))
#    print "got array up",array_up
#    print "and array down",array_down
    return array_up, array_down

def build_arrays(nbins = 100, events = 1000, nsamples=1):
    allhists = {}
    nevents = events *.99/nbins
    signal_hist = (0.01*nevents/nbins)*np.arange(nbins)
    data = signal_hist
    allhists['signal_hist'] = signal_hist
    backgrounds = []
    for s in range(nsamples):
        thescale = nevents/nsamples
        print 'the scale is {}'.format(thescale)
        print 'with nevents',nevents
        print 'with bins',nbins
        print 'and ',nsamples
        background_hist = (nevents/nsamples)*np.ones(nbins)
        data = data+background_hist
        allhists['sample_{}_hist'.format(s)] = background_hist
        backgrounds.append(background_hist)
    allhists['data_hist'] = randomise_data( data, 2.)
                   
    slopes = np.linspace(2.,-2.,201)
    for nuis in range(201):
        sig_nuis_up, sig_nuis_down = make_variations(signal_hist,slopes[nuis])
        allhists['signal_up_{}'.format(slopes[nuis])] = sig_nuis_up
        allhists['signal_down_{}'.format(slopes[nuis])] = sig_nuis_down
        nsample = 0
        for background_hist in backgrounds:
            bkg_nuis_up, bkg_nuis_down = make_variations(background_hist,slopes[nuis])
            allhists['sample_{}_up_{}'.format(nsample, slopes[nuis])] = bkg_nuis_up
            allhists['sample_{}_down_{}'.format(nsample, slopes[nuis])] = bkg_nuis_down
            nsample += 1
    return allhists
