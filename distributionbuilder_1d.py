from math import cos

from ROOT import TH1F

from distributionbuilder_base import DistributionBuilderBase, iter


class DistributionBuilder_1d(DistributionBuilderBase):
    def __init__(self, histname_suffix, events, bins):
        super().__init__(histname_suffix, events, bins)
        self.setParameters(0.0)

    def setParameters(self, lambda_theta):
        self.lambda_theta = lambda_theta

    def calcWeight(self, theta):
        return 1 + self.lambda_theta * cos(theta) ** 2

    def add_angular(self):
        self.angular_hists: list[TH1F] = []
        for bin in self.bins:
            histname = "hist%s%s_iter%i" % (bin.suffix(), self.histname_suffix, iter)
            newhist = TH1F(histname, histname, 20, -1, 1)
            self.angular_hists.append(newhist)
        if not isinstance(newhist, TH1F):
            print("Error, newHist is: ", newhist)

    def fill_angular(self, binIndex, event):
        hist = self.angular_hists[binIndex]

        weight = self.calcWeight(event.theta)
        hist.Fill(cos(event.theta), weight)
