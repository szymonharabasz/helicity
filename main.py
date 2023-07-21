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

        for bin in bins:
            print(bin)
            print(bin.suffix())
    except yaml.YAMLError as exc:
        print(exc)

reader = SpdmeEventsReader()
reader.readEvents("medium_isotropic_4pi_new_less.dat")
#reader.readEvents("medium_bt_mean_embData_gen5_noRichQaCorr_richQaCI_beta2sigma_acc100mev_fakerej_notIgnoreInnerMDCmomOrg_fiducialVol_isotropic.rootat")
#reader.readEvents("medium_isotropic_eff_ag1230ag_nn_9deg.dat")
#reader.readEvents("apr12_diele_088_090_ag123ag_2500A_accepted_nn_2.dat")
event = reader.events[10]
print(event.getPx(0) + event.getPx(1))
print(event.getPy(0) + event.getPy(1))
builder = DistributionBuilder()
builder.setParameters(0.0, 0.0, 0.0)
hists = builder.buildFromEvents(reader.events, bins)
c = TCanvas("c","c",100,100,1200,1600)
c.Divide(3,4)
ic = 1
totentries = 0
file = TFile("out.root","RECREATE")
for hist in hists[0]:
    c.cd(ic)
    hist.Draw("COLZ")
    ic = ic + 1
    totentries = totentries + hist.GetEntries()
    hist.Write()
hists[1].Write()
hists[2].Write()
c.SaveAs("c.gif")
print("Total hist entries: ", totentries)