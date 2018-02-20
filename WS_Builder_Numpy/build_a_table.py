from TableBuilder import Model
import numpy as np

def gen_nps(n_nps):
    return np.random.choice(np.linspace(-2,2,201),n_nps)

m = Model()
nps = gen_nps(5)
m.add_nuisance_pars(nps)
m.add_channel(10)
m.add_channel(100)
m.add_channel(1)
m.add_sample([0,2],[10,500],signal=True,shape='linear',histosys=([1,3,4],[[0],[0,2],[0,2]]))
m.add_sample([0,1,2],[1000,500,10000],signal=False,shape='linear',histosys=[0,2,4])
m.add_sample([1,2],[750,5000],signal=False,shape='linear',histosys=([0,2,3],[[1],[1,2],[2]]))


w = m.BuildWorkspace()
m.Print()


import ROOT
mc = w.obj("ModelConfig")
pdf = mc.GetPdf() 
data = w.data("obsData")
#x = w.var("obs_x_channel_0")
mc.LoadSnapshot()

sbModel = w.obj("ModelConfig")
sbModel.SetName("S+B Model")
poi = sbModel.GetParametersOfInterest().first()
poi.setVal(1)
sbModel.SetSnapshot(ROOT.RooArgSet(poi))
bModel = sbModel.Clone()
bModel.SetName("B Model")
poi.setVal(0)
bModel.SetSnapshot(ROOT.RooArgSet(poi))
ac = ROOT.RooStats.AsymptoticCalculator(data, sbModel, bModel)
ac.SetOneSidedDiscovery(True)

asResult = ac.GetHypoTest()
asResult.Print()
