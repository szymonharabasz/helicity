from ctypes.wintypes import HINSTANCE
from math import sin, cos
from ROOT import TH1F, TMath, TFile
from datetime import datetime
from distributionbuilder_1d import DistributionBuilder_1d, iter

iter = 1

class SurrogateDistributionBuilder_1d(DistributionBuilder_1d):
    def __init__(self, histname_suffix, events, bins):
        super().__init__(histname_suffix, events, bins)
        self.base_hists = self.buildFromEvents()

    def reweightBaseHists(self) -> list[TH1F]:
        global iter
        hists: list[TH1F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH1F(histname,histname,20,-1,1)
            hists.append(newhist)
        if not isinstance(newhist, TH1F):
            print("Error, newHist is: ", newhist)
            
        for ihist in range(len(hists)):

            hist = hists[ihist]
            baseHist = self.base_hists[0][ihist]
            for ix in range(hist.GetNbinsX()):
                bx = ix+1
                content = baseHist.GetBinContent(bx)
                error = baseHist.GetBinError(bx)
                cosTheta = hist.GetXaxis().GetBinCenter(bx)
                theta = TMath.ACos(cosTheta)

                weight = self.calcWeight(theta)
                hist.SetBinContent(bx, content*weight)
                hist.SetBinError(bx, error*weight)
                    
            
            if hist.Integral() > 0:
                hist.Scale(1./hist.Integral())
        
        hmass = self.base_hists[1][0]
        hz = self.base_hists[2][0]
            
        iter = iter + 1
        return [hists, [hmass], [hz]]

    def getHists(self) -> list[TH1F]:
        return self.reweightBaseHists()