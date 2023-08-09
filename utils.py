import math
from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from surrogatedistributionbuilder import SurrogateDistributionBuilder
from ROOT import TCanvas, TFile

def calcDiff(histMC, histData, bx, by):
    contentMC = histMC.GetBinContent(bx+1,by+1)
    contentData = histData.GetBinContent(bx+1,by+1)
    errorData = histData.GetBinError(bx+1,by+1)
    if contentMC > 0 and contentData > 0:
        return pow((contentData - contentMC)/errorData, 2)
    else:
        return 0.0

def calcOneChi2(histMC, histData):
    chi2 = 0.
    ndf = 0
    for bx in range(histMC.GetXaxis().GetNbins()):
        for by in range(histMC.GetYaxis().GetNbins()):
            diff = calcDiff(histMC, histData, bx, by)
            if diff > 0:
                ndf = ndf + 1
                chi2 = chi2 + diff
    return chi2, ndf

def diffHist(histMC, histData):
    name = histMC.GetName() + "_diff"
    hdiff = histMC.Clone(name)
    for bx in range(histMC.GetXaxis().GetNbins()):
        for by in range(histMC.GetYaxis().GetNbins()):
            diff = calcDiff(histMC, histData, bx, by)
            hdiff.SetBinContent(bx, by, diff)
    return hdiff

def calcAllChi2(histsMC, histsData):
    for i in range(len(histsMC[0])):
        histMC = histsMC[0][i]
        histData = histsData[0][i]

        chi2, ndf = calcOneChi2(histMC, histData)
        print("chi2, ndf, chi2/ndf", chi2, ndf, chi2/ndf)
        yield chi2, ndf

class HistMaker:
    def __init__(self, filename, name_suffix, bins):
        self.reader = SpdmeEventsReader(filename)
        # self.builder = DistributionBuilder(name_suffix, reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder(name_suffix, self.reader.getEvents(), bins)
    def makeHists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
        self.builder.setParameters(lambda_theta, lambda_phi, lambda_theta_phi)
        return self.builder.getHists()


def makeRootPlots(histsMC, histsData):
    c = TCanvas("c","c",100,100,1200,1600)
    c.Divide(3,4)

    totentries = 0
    file = TFile("out.root","RECREATE")
    for i in range(len(histsMC[0])):
        histMC = histsMC[0][i]
        c.cd(i+1)
        histMC.Draw("COLZ")
        totentries = totentries + histMC.GetEntries()
        histMC.Write()

        histData = histsData[0][i]
        histData.Write()

    histsMC[1].Write()
    histsMC[2].Write()
    histsData[1].Write()
    histsData[2].Write()
    c.SaveAs("c.gif")
    print("Total hist entries: ", totentries)