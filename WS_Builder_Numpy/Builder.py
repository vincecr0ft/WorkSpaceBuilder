import ROOT
import numpy as np
import root_numpy as rumpy
import math
ROOT.RooMsgService.instance().setGlobalKillBelow(5)

def randomise_data(hist,sigma):
    contents = rumpy.hist2array(hist)
    perbin = lambda mu: np.random.normal(mu,sigma/math.sqrt(mu))
    randomised = np.vectorize(perbin)
    return rumpy.array2hist(randomised(contents),hist)

def make_variations(hist, slope):
    np_up = hist.Clone()
    np_down = hist.Clone()
    contents = rumpy.hist2array(hist)
    if len(contents) == 1:
        contents[0] = 1
    thebins = np.arange(len(contents))
    if slope < 0:
        thebins = thebins[::-1]
    gradient = lambda b, s: b*s
    sloped = np.vectorize(gradient)
    array_up = sloped(thebins, slope)*(1/math.sqrt(contents.sum())) + contents
    array_down = contents - sloped(thebins, slope)*(1/math.sqrt(contents.sum()))
    return rumpy.array2hist(array_up,np_up), rumpy.array2hist(array_down,np_down)

def get_workspace(nchannels = 1, nsamples =1, events = 1000, nbins = 1, nnps = 0):
    nevents = events *.9/nbins 
    meas = ROOT.RooStats.HistFactory.Measurement( "meas", "meas" )
    meas.SetPOI( "SignalStrength" )
    meas.SetLumi( 1.0 )
    meas.SetLumiRelErr( 0.02 )
    meas.AddConstantParam( "Lumi" )

    for channel in range(nchannels):
        for sample in range(nsamples):
            chan = ROOT.RooStats.HistFactory.Channel( "Region{}Sample{}".format(channel,sample) )
            chan.SetStatErrorConfig(0.05, "Poisson")
            if sample < 1:
                data_hist = ROOT.TH1D("{}_observed{}".format("signal", channel),"observed", nbins, 0, nbins)
                signal_hist = ROOT.TH1D("{}_above_expected{}".format("signal", channel),"above_expected", nbins, 0, nbins)
                model_hists = [ROOT.TH1D("expected{}sample{}".format(sample, channel),"expected", nbins, 0, nbins) for sample in range(nsamples)]
                rumpy.array2hist((0.1*nevents/nbins)*np.arange(nbins),signal_hist)
                data_hist.Add(signal_hist)
                for m in model_hists:
                    rumpy.array2hist((nevents/nbins)*np.ones(nbins),m)
                    data_hist.Add(m)
                    
                data_hist = randomise_data( data_hist, 2.)
                chan.SetData( data_hist )
                
                models = []
                for s in range(nsamples):
                    model = ROOT.RooStats.HistFactory.Sample( "channel{}model{}".format(channel,s) )
                    model.SetNormalizeByTheory( False )
                    model.SetHisto( model_hists[s] )
                    model.ActivateStatError()
                    models.append(model)
                signal = ROOT.RooStats.HistFactory.Sample( "signal" )
                signal.ActivateStatError()
                signal.SetNormalizeByTheory( False )
                signal.SetHisto( signal_hist )
                signal.AddNormFactor( "SignalStrength", 1, 0, 3)
            else:
                data_hist = ROOT.TH1D("CR{}observed{}".format(sample,channel),"observed", 1, 0, 1)
                signal_hist = ROOT.TH1D("CR{}above_expected{}".format(sample,channel),"above_expected", 1, 0, 1)
                model_hists = [ROOT.TH1D("CR{}expected{}sample{}".format(sample,channel,s),"expected", 1, 0, 1) for s in range(nsamples)]

                signal_hist.SetBinContent(1,0.01*nevents)
                for m in range(len(model_hists)):
                    if m == channel:
                        model_hists[m].SetBinContent(1,0.9*nevents)
                        data_hist.Add(model_hists[m])
                    else:
                        model_hists[m].SetBinContent(1,0.1*nevents)
                        data_hist.Add(model_hists[m])

                                                                                                                            
                data_hist = randomise_data( data_hist, 2.)        
                chan.SetData( data_hist )
                
                models = []
                for sample in range(nsamples):
                    model = ROOT.RooStats.HistFactory.Sample( "model{}".format(sample) )
                    model.SetNormalizeByTheory( False )
                    model.SetHisto( model_hists[sample] )
                    model.ActivateStatError()
                    models.append(model)
                signal = ROOT.RooStats.HistFactory.Sample( "signal" )
                signal.SetNormalizeByTheory( False )
                print "this signal hist is ",signal_hist
                signal.SetHisto( signal_hist )
                signal.ActivateStatError()
                signal.AddNormFactor( "SignalStrength", 1, 0, 3)

            
            if nnps > 0:
                slopes = np.linspace(1.,-1.,nnps)
                for nuis in range(nnps):
                    print "this uncertainty is {}".format(nuis)
                    uncertainty_up   = nevents * 1.1
                    uncertainty_down = nevents * 0.9
                    signal.AddOverallSys( "signal_norm_uncertainty_{}".format(nuis),  uncertainty_down*.1, uncertainty_up*.1 )
                    for model in models:
                        model.AddOverallSys( "background_norm_uncertainty_{}".format(nuis),  uncertainty_down, uncertainty_up )
                
                    
                    sig_nuis_up, sig_nuis_down = make_variations(signal_hist,slopes[nuis]) 
                    signal_shape = ROOT.RooStats.HistFactory.HistoSys("signal_shape_{}".format(nuis))
                    signal_shape.SetHistoHigh( sig_nuis_up )
                    signal_shape.SetHistoLow( sig_nuis_down )
                    signal.AddHistoSys( signal_shape )
                    for bkg in range(nsamples):
                        background_shape = ROOT.RooStats.HistFactory.HistoSys("background_{}_shape_{}".format(bkg,nuis))
                        bkg_nuis_up, bkg_nuis_down = make_variations(model_hists[bkg],slopes[nuis])
                        background_shape.SetHistoHigh( bkg_nuis_up )
                        background_shape.SetHistoLow( bkg_nuis_down )
                        models[bkg].AddHistoSys( background_shape )

            for model in models:
                chan.AddSample( model )
            chan.AddSample( signal )
            meas.AddChannel( chan )


    hist2workspace = ROOT.RooStats.HistFactory.HistoToWorkspaceFactoryFast(meas)
    if nchannels < 1 and nsamples < 1:
        ws = hist2workspace.MakeSingleChannelModel( meas, chan )
    else:
        ws = hist2workspace.MakeCombinedModel(meas)
    iter = ws.components().fwdIterator()
    arg = iter.next()
    while arg:
        if "RooRealSum" in str(arg.IsA()):
            arg.setAttribute("BinnedLikelihood")
        arg = iter.next()
    return ws
