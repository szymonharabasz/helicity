from ctypes.wintypes import HINSTANCE
from math import sin, cos
from ROOT import TH1F, TMath, TFile
from tqdm.notebook import tqdm
from datetime import datetime
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

    def buildFromEvents(self):
        global iter

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}]: iter ", iter)
        hists: list[TH1F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH1F(histname,histname,20,-1,1)
            hists.append(newhist)
        if not isinstance(newhist, TH1F):
            print("Error, newHist is: ", newhist)
            
        hmassLowM = TH1F("hmassLowM" + self.histname_suffix,"hmassLowM" + self.histname_suffix,100,0,1000)
        hzLowM = TH1F("hzLowM" + self.histname_suffix,"hzLowM" + self.histname_suffix,100,-1,1)
        hyLowM = TH1F("hyLowM" + self.histname_suffix,"hyLowM" + self.histname_suffix,100,-2,2)
        hmassHigM = TH1F("hmassHigM" + self.histname_suffix, "hmassHigM" + self.histname_suffix, 100, 0, 1000)
        hzHigM = TH1F("hzHigM" + self.histname_suffix, "hzHigM" + self.histname_suffix, 100, -1, 1)
        hyHigM = TH1F("hyHigM" + self.histname_suffix,"hyHigM" + self.histname_suffix,100,-2,2)
        ievent = 0
        nevents = len(self.events)

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}] Before processing events")
        print(f"Num events {nevents}")
        print(F"Num non-null events {len([evt for evt in self.events if evt is not None])}")
        for event in self.events:
            ievent = ievent + 1
            binIndex = self.binIndex(event, self.bins)

            if binIndex >= 0 and binIndex < len(hists):
                hist = hists[binIndex]

                weight = self.calcWeight(event.theta)
                hist.Fill(cos(event.theta), weight)
                if binIndex < 3 and event.mass < self.bins[2].m_max:
                    hmassLowM.Fill(event.mass)
                    hzLowM.Fill(event.z)
                    hyLowM.Fill(event.y)
                if binIndex >= 3 and event.mass >= self.bins[2].m_max:
                    hmassHigM.Fill(event.mass)
                    hzHigM.Fill(event.z)
                    hyHigM.Fill(event.y)

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}] After processing events")
        for hist in [*hists, hmassLowM, hzLowM, hyLowM]:
            if hist.Integral() > 0:
                pass
               # hist.Scale(1./hist.Integral())

        iter = iter + 1
        result = [hists, [hmassLowM, hmassHigM], [hzLowM, hzHigM], [hyLowM, hyHigM]]
        print(f"#1 result length {len(result)}")
        return result

    def getHists(self) -> list[TH1F]:
        result = self.buildFromEvents()
        print(f"#2 result length {len(result)}")
        return result