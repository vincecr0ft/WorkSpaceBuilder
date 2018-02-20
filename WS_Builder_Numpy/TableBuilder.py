import ROOT
import math
import numpy as np
import root_numpy as rumpy
from prettytable import PrettyTable

class Model(object):
    def __init__(self):
        self.n_channels = 0
        self.channels = {}
        
        self.n_samples = 0
        self.samples = {}
        
        self.n_nps = 0
        self.nps = {}

    def add_nuisance_pars(self, nps):
        for par in nps:
            name = 'nuis_{}'.format(self.n_nps)
            self.n_nps += 1
            self.nps[name] = par

    def add_channel(self, n_bins):
        name = 'channel_{}'.format(self.n_channels)
        self.n_channels += 1
        self.channels[name] = self.Channel(n_bins)

    class Channel:
        def __init__(self, n_bins):
            self.n_bins = n_bins
            self.samples = {}

        def add_sample(self, name, n_events, signal, shape, histosys):
            hist = self.build_array(self.n_bins, n_events, signal)
            self.samples[name] = {
                'n_events':n_events,
                'signal':signal,
                'nominal': hist
            }
            for sys_name,value in histosys:
                self.samples[name][sys_name] = (self.make_variations(hist, value))

        def build_array(self, n_bins = 10, n_events = 100, signal=False):
            if signal:
                slope = np.arange(1,n_bins+1,dtype=float)
                return np.divide(slope,np.sum(slope))*n_events
            else:
                return np.divide(np.ones(n_bins),n_bins)*n_events
            
        def make_variations(self, hist, slope):
            if len(hist) == 1:
                hist[0] = 1
            thebins = np.arange(len(hist))
            if slope < 0:
                thebins = thebins[::-1]
            slope = math.fabs(slope)
            gradient = lambda b, s: b*s
            sloped = np.vectorize(gradient)
            slope_vec = 1.+sloped(thebins, slope)
            p_sloped = np.divide(slope_vec,np.sum(slope_vec))
            n = np.sum(hist) 
            factor = 1./math.sqrt(n)
            ep_sloped = n*p_sloped*factor
            array_up = hist + ep_sloped
            array_down = hist - ep_sloped
            return array_up, array_down
        
    def add_sample(self, channels, events, signal=False, shape='linear', histosys=None):
        assert len(channels)==len(events)
        sample_name = 'sample_{}'.format(self.n_samples)
        self.n_samples += 1
        for chan,evt in zip(channels,events):
            if chan > self.n_channels:
                return 'channel {} does not exist'.format(chan)
            else:
                nps = []
                chan_name = 'channel_{}'.format(chan)
                if isinstance(histosys,list):
                    for sys in histosys:
                        nuis_name = 'nuis_{}'.format(sys)
                        if nuis_name in self.nps:
                            nps.append((nuis_name, self.nps[nuis_name]))
                elif isinstance(histosys,tuple):
                    whichsys, whichchan = histosys
                    assert len(whichsys) == len(whichchan)
                    for i, sys in enumerate(whichsys):
                        nuis_name = 'nuis_{}'.format(sys)
                        if chan in whichchan[i] and nuis_name in self.nps:
                            nps.append((nuis_name, self.nps[nuis_name]))
                        #    print 'syst {} applied to sample {} in channels {}'.format(sys,sample_name,whichchan[i])
                self.channels[chan_name].add_sample(sample_name, evt, signal, shape, nps)
                
                
    def BuildWorkspace(self, name = 'WorkspaceName'):
        meas = ROOT.RooStats.HistFactory.Measurement( name, name )
        meas.SetPOI( "SignalStrength" )
        meas.SetLumi( 1.0 )
        meas.SetLumiRelErr( 0.02 )
        meas.AddConstantParam( "Lumi" )
        for channel in self.channels:
            c = self.channels[channel]
            chan = ROOT.RooStats.HistFactory.Channel( channel )
            chan.SetStatErrorConfig(0.05, "Poisson")
            data = None
            for sample in c.samples:
                s = c.samples[sample]
                thisname = '{}_{}'.format(channel,sample)
                h1 = ROOT.TH1D(thisname, thisname, c.n_bins, 0, c.n_bins)
                sample_hist = rumpy.array2hist(s['nominal'], h1)
                
                samp = ROOT.RooStats.HistFactory.Sample( thisname )
                samp.SetNormalizeByTheory( False )
                samp.SetHisto( sample_hist )
                samp.ActivateStatError()
                if s['signal'] == True:
                    samp.AddNormFactor( "SignalStrength", 1, 0, 3)

                
                for distributions in s:
                    if 'nuis' in distributions:
                        shape = ROOT.RooStats.HistFactory.HistoSys(distributions)
                        up, down = s[distributions]
                        thisname += distributions
                        hup = ROOT.TH1D(thisname+'up', thisname+'up', c.n_bins, 0, c.n_bins)
                        hdown = ROOT.TH1D(thisname+'down', thisname+'down', c.n_bins, 0, c.n_bins)
                        shape.SetHistoHigh( rumpy.array2hist(up, hup) )
                        shape.SetHistoLow( rumpy.array2hist(down, hdown) )
                        samp.AddHistoSys( shape )
                chan.AddSample(samp)
                
                if not data is None:
                    data += np.asarray(s['nominal'])
                else:
                    data = np.asarray(s['nominal'])
            hdata = ROOT.TH1D('{}_data'.format(channel), '{}_data'.format(channel), c.n_bins, 0, c.n_bins)        
            chan.SetData( rumpy.array2hist(data, hdata))
            meas.AddChannel(chan)

            
        hist2workspace = ROOT.RooStats.HistFactory.HistoToWorkspaceFactoryFast(meas)
        if self.n_channels < 1:
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

    def Print(self):
        x = PrettyTable()
        x.field_names = ["Channel","sample","n bins","n events","signal?"]+self.nps.keys()
        for c_name in sorted(self.channels):
            chan = self.channels[c_name]
            for s_name in sorted(chan.samples):
                np_applied = ['X' if nuis in chan.samples[s_name] else '' for nuis in self.nps]
                x.add_row([
                    c_name,
                    s_name,
                    chan.n_bins,
                    chan.samples[s_name]['n_events'],
                    chan.samples[s_name]['signal']]
                    +np_applied)
        print x
