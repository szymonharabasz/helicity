from math import sin, cos
from ROOT import TH1F, TH2F, TMath

class DistributionBuilder:
    def setParameters(self, lambda_theta, lambda_phi, lambda_theta_phi, lambda_phi_perp = 0, lambda_theta_phi_perp = 0):
        self.lambda_theta = lambda_theta
        self.lambda_phi = lambda_phi
        self.lambda_theta_phi = lambda_theta_phi
        self.lambda_phi_perp = lambda_phi_perp
        self.lambda_theta_phi_perp = lambda_theta_phi_perp

    def calcWeight(self, theta, phi):
        return 1 + self.lambda_theta + cos(theta)**2 + self.lambda_theta_phi * sin(2*theta)*cos(phi) + self.lambda_phi * sin(theta)**2 * cos(2*phi) \
            + self.lambda_theta_phi_perp * sin(2*theta) * sin(phi) + self.lambda_phi_perp * sin(theta)**2 * sin(2*phi)

    def binIndex(self, event, bins):
        index = 0
        for bin in bins:
            if event.mass >= bin.m_min and event.mass < bin.m_max and event.z >= bin.z_min and event.z < bin.z_max:
                return index
            index = index + 1
        return len(bins)

    def buildFromEvents(self, events, bins):
        hists = [
            TH2F("hist%s" % bin.suffix(),"hist%s" % bin.suffix(),20,-1,1,36,0,2*TMath.Pi()) for bin in bins
        ]
        hmass = TH1F("hmass","hmass",100,0,1000)
        hz = TH1F("hz","hz",100,-1,1)
        for event in events:
            binIndex = self.binIndex(event, bins)

            if binIndex >= 0 and binIndex < len(hists):
                hist = hists[binIndex]

                weight = self.calcWeight(event.theta, event.phi)
               # weight = event.weight * self.calcWeight(event.theta, event.phi)
               # print(weight, event.weight, self.calcWeight(event.theta, event.phi))
               # hist.Fill(cos(event.theta), event.phi, weight)
                hist.Fill(cos(event.theta), event.phi, weight)
                hmass.Fill(event.mass)
                hz.Fill(event.z)
        return (hists, hmass, hz)