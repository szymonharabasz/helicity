from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from ROOT import TCanvas, TFile

def calcOneChi2(histMC, histData):
    chi2 = 0.
    ndf = 0
    for bx in range(histMC.GetXaxis().GetNbins()):
        for by in range(histMC.GetYaxis().GetNbins()):
            contentMC = histMC.GetBinContent(bx+1,by+1)
            contentData = histData.GetBinContent(bx+1,by+1)
            errorData = histData.GetBinError(bx+1,by+1)
           # print("contents and error: ", contentMC, contentData, errorData)
            if contentMC > 0 and contentData > 0:
                ndf = ndf + 1
                chi2 = chi2 + pow((contentData - contentMC)/contentData, 2)

    return chi2, ndf


def calcAllChi2(histsMC, histsData):
    for i in range(len(histsMC[0])):
        histMC = histsMC[0][i]
        histData = histsData[0][i]

        chi2, ndf = calcOneChi2(histMC, histData)
        print("chi2, ndf, chi2/ndf", chi2, ndf, chi2/ndf)
        yield chi2, ndf

def makeHists(filename, name_suffix, bins, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
    reader = SpdmeEventsReader(filename)
    builder = DistributionBuilder(name_suffix)
    builder.setParameters(lambda_theta, lambda_phi, lambda_theta_phi)
    hists = builder.buildFromEvents(reader.getEvents(), bins)
    return hists


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