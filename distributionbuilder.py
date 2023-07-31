from ctypes.wintypes import HINSTANCE
from math import sin, cos
from ROOT import TH1F, TH2F, TMath, TFile
from tqdm.notebook import tqdm
from datetime import datetime

iter = 0

class DistributionBuilder:
    def __init__(self, histname_suffix):
        self.histname_suffix = histname_suffix
       # self.outfile = TFile("hists.root","RECREATE")

    def setParameters(self, lambda_theta, lambda_phi, lambda_theta_phi, lambda_phi_perp = 0, lambda_theta_phi_perp = 0):
        self.lambda_theta = lambda_theta
        self.lambda_phi = lambda_phi
        self.lambda_theta_phi = lambda_theta_phi
        self.lambda_phi_perp = lambda_phi_perp
        self.lambda_theta_phi_perp = lambda_theta_phi_perp

    def calcWeight(self, theta, phi):
        return 1 + self.lambda_theta * cos(theta)**2 + self.lambda_theta_phi * sin(2*theta)*cos(phi) + self.lambda_phi * sin(theta)**2 * cos(2*phi) \
            + self.lambda_theta_phi_perp * sin(2*theta) * sin(phi) + self.lambda_phi_perp * sin(theta)**2 * sin(2*phi)

    def binIndex(self, event, bins):
        index = 0
        for bin in bins:
            if event.mass >= bin.m_min and event.mass < bin.m_max and event.z >= bin.z_min and event.z < bin.z_max:
                return index
            index = index + 1
        return len(bins)

    def buildFromEvents(self, events, bins) -> list[TH2F]:
        global iter
        hists: list[TH2F] = []
        for bin in bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH2F(histname,histname,20,-1,1,36,0,2*TMath.Pi())
            hists.append(newhist)
        if not isinstance(newhist, TH2F):
            print("Error, newHist is: ", newhist)
            
        hmass = TH1F("hmass" + self.histname_suffix,"hmass" + self.histname_suffix,100,0,1000)
        hz = TH1F("hz" + self.histname_suffix,"hz" + self.histname_suffix,100,-1,1)
        ievent = 0
        nevents = len(events)
        print("Before processing events", datetime.now().strftime("%H:%M:%S"))
        for event in tqdm(events):
            if ievent % 100000 == 0:
                pass
               # print("Processing event ", ievent, " out of ", nevents)
            ievent = ievent + 1
            binIndex = self.binIndex(event, bins)

            if binIndex >= 0 and binIndex < len(hists):
                hist = hists[binIndex]

                weight = self.calcWeight(event.theta, event.phi)
               # weight = event.weight * self.calcWeight(event.theta, event.phi)
               # print(weight, event.weight, self.calcWeight(event.theta, event.phi))
               # hist.Fill(cos(event.theta), event.phi, weight)
               # print("weight, cos: ", weight, cos(event.theta), TMath.Cos(event.theta))
                hist.Fill(cos(event.theta), event.phi, weight)
                hmass.Fill(event.mass)
                hz.Fill(event.z)
        print("After processing events", datetime.now().strftime("%H:%M:%S"))
        for hist in [*hists, hmass, hz]:
            hist.Scale(1./hist.Integral())
           # hist.Write()

        iter = iter + 1
        return [hists, [hmass], [hz]]