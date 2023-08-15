from ctypes.wintypes import HINSTANCE
from math import sin, cos
from ROOT import TH1F, TH2F, TMath, TFile
from datetime import datetime
from distributionbuilder import DistributionBuilder, iter

iter = 1

class SurrogateDistributionBuilder(DistributionBuilder):
    def __init__(self, histname_suffix, events, bins):
        super().__init__(histname_suffix, events, bins)
        self.base_hists = self.buildFromEvents()

    def reweightBaseHists(self) -> list[TH2F]:
        global iter
        hists: list[TH2F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH2F(histname,histname,20,-1,1,36,0,2*TMath.Pi())
            hists.append(newhist)
        if not isinstance(newhist, TH2F):
            print("Error, newHist is: ", newhist)
            
        if iter % 10 == 0:
            print("Before reweighting hists", datetime.now().strftime("%H:%M:%S"))
        for ihist in range(len(hists)):

            hist = hists[ihist]
            baseHist = self.base_hists[0][ihist]
            for ix in range(hist.GetNbinsX()):
                for iy in range(hist.GetNbinsY()):
                    bx = ix+1
                    by = iy+1
                    content = baseHist.GetBinContent(bx, by)
                    phi = hist.GetYaxis().GetBinCenter(by)
                    cosTheta = hist.GetXaxis().GetBinCenter(bx)
                    theta = TMath.ACos(cosTheta)

                    weight = self.calcWeight(theta, phi)
                    hist.SetBinContent(bx, by, content*weight)
                   # if bx == 11:
                   #     print(bx, by, cosTheta, theta, phi, weight)
                    
            
            hist.Scale(1./hist.Integral())
        if iter % 10 == 0:
            print("After reweighting hists", datetime.now().strftime("%H:%M:%S"))
        
        hmass = self.base_hists[1][0]
        hz = self.base_hists[2][0]
            
        iter = iter + 1
        return [hists, [hmass], [hz]]

    def getHists(self) -> list[TH2F]:
        return self.reweightBaseHists()