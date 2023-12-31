from math import sin, cos

from ROOT import TH2F, TMath

from distributionbuilder_base import DistributionBuilderBase

iter = 0

class DistributionBuilder(DistributionBuilderBase):
    def __init__(self, histname_suffix, events, bins):
        super().__init__(histname_suffix, events, bins)
        self.setParameters(0.0, 0.0, 0.0)
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

    def add_angular(self):
        self.angular_hists: list[TH2F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH2F(histname,histname,20,-1,1,36,0,2*TMath.Pi())
            self.angular_hists.append(newhist)
        if not isinstance(newhist, TH2F):
            print("Error, newHist is: ", newhist)

    def fill_angular(self, binIndex, event):
        hist = self.angular_hists[binIndex]

        weight = self.calcWeight(event.theta, event.phi)
        hist.Fill(cos(event.theta), event.phi, weight)
