from ctypes.wintypes import HINSTANCE
from math import sin, cos
from ROOT import TH1F, TMath, TFile
from tqdm.notebook import tqdm
from datetime import datetime

iter = 0

class DistributionBuilder_1d:
    def __init__(self, histname_suffix, events, bins):
        self.histname_suffix = histname_suffix
        self.events = events
        self.bins = bins
        self.setParameters(0.0)

    def setParameters(self, lambda_theta):
        self.lambda_theta = lambda_theta

    def calcWeight(self, theta):
        return 1 + self.lambda_theta * cos(theta)**2

    def binIndex(self, event, bins):
        index = 0
        for bin in bins:
            if event.mass >= bin.m_min and event.mass < bin.m_max and event.z >= bin.z_min and event.z < bin.z_max:
                return index
            index = index + 1
        return len(bins)

    def buildFromEvents(self) -> list[TH1F]:
        global iter
        print("iter", iter)
        hists: list[TH1F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH1F(histname,histname,20,-1,1)
            hists.append(newhist)
        if not isinstance(newhist, TH1F):
            print("Error, newHist is: ", newhist)
            
        hmass = TH1F("hmass" + self.histname_suffix,"hmass" + self.histname_suffix,100,0,1000)
        hz = TH1F("hz" + self.histname_suffix,"hz" + self.histname_suffix,100,-1,1)
        ievent = 0
        nevents = len(self.events)
        print("Before processing events", datetime.now().strftime("%H:%M:%S"))
        for event in self.events:
            ievent = ievent + 1
            binIndex = self.binIndex(event, self.bins)

            if binIndex >= 0 and binIndex < len(hists):
                hist = hists[binIndex]

                weight = self.calcWeight(event.theta)
                hist.Fill(cos(event.theta), weight)
                hmass.Fill(event.mass)
                hz.Fill(event.z)
        print("After processing events", datetime.now().strftime("%H:%M:%S"))
        for hist in [*hists, hmass, hz]:
            if hist.Integral() > 0:
                hist.Scale(1./hist.Integral())

        iter = iter + 1
        return [hists, [hmass], [hz]]

    def getHists(self) -> list[TH1F]:
        return self.buildFromEvents()