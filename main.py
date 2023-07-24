from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from bins import Bins
from ROOT import TCanvas, TFile
import yaml

with open("ranges.yml","r") as stream:
    try:
        borders = yaml.safe_load(stream)
        m_borders = borders["m-borders"]
        z_borders = borders["z-borders"]

        bins = []
        for i in range(len(m_borders)-1):
            for j in range(len(z_borders)-1):
                bins.append(Bins(m_borders[i], m_borders[i+1], z_borders[j], z_borders[j+1]))

#        for bin in bins:
#            print(bin)
#            print(bin.suffix())
    except yaml.YAMLError as exc:
        print(exc)

readerMC = SpdmeEventsReader()
readerMC.readEvents("medium_isotropic_eff_ag1230ag_nn_9deg.dat")
readerData = SpdmeEventsReader()
readerData.readEvents("apr12_diele_088_090_ag123ag_2500A_accepted_nn_2.dat")
builderMC = DistributionBuilder("_MC")
builderMC.setParameters(0.0, 0.0, 0.0)
histsMC = builderMC.buildFromEvents(readerMC.events, bins)
builderData = DistributionBuilder("_data")
builderData.setParameters(0.0, 0.0, 0.0)
histsData = builderData.buildFromEvents(readerData.events, bins)
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
    histMC.Write()

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
                chi2 = chi2 + pow((contentData - contentMC)/errorData, 2)

    print("chi2, ndf, chi2/ndf", chi2, ndf, chi2/ndf)

histsMC[1].Write()
histsMC[2].Write()
histsData[1].Write()
histsData[2].Write()
c.SaveAs("c.gif")
print("Total hist entries: ", totentries)