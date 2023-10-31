import math
from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from surrogatedistributionbuilder import SurrogateDistributionBuilder
from distributionbuilder_1d import DistributionBuilder_1d
from surrogatedistributionbuilder_1d import SurrogateDistributionBuilder_1d
from ROOT import TCanvas, TFile, TPaveText, gStyle, TMath, TString, TH1F

def calcDiff(histMC, histData, bx, by):
    contentMC = histMC.GetBinContent(bx+1,by+1) / histMC.Integral()
    contentData = histData.GetBinContent(bx+1,by+1) / histData.Integral()
    errorData = histData.GetBinError(bx+1,by+1) / histData.Integral()
    if contentMC > 0 and contentData > 0:
       # return pow((contentData - contentMC)/errorData, 2)
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
            hdiff.SetBinContent(bx+1, by+1, diff)
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
        #self.builder = DistributionBuilder(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder(name_suffix, self.reader.getEvents(), bins)
    def makeHists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
        self.builder.setParameters(lambda_theta, lambda_phi, lambda_theta_phi)
        return self.builder.getHists()

class HistMaker_1d:
    def __init__(self, filename, name_suffix, bins):
        self.reader = SpdmeEventsReader(filename)
        #self.builder = DistributionBuilder_1d(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder_1d(name_suffix, self.reader.getEvents(), bins)
    def makeHists(self, lambda_theta=0):
        self.builder.setParameters(lambda_theta)
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

def setOPT_text(cut_desc, x1, y1, x2, y2, color, textSize):

    description = TPaveText(x1, y1, x2, y2, "NDC")
    description.SetLineWidth(0)
    description.AddText(cut_desc)
    description.SetTextSize(textSize)
    description.SetBorderSize(0)
    description.SetTextFont(62)
    description.SetTextColor(color)
    description.SetFillColor(0)
    description.SetFillStyle(0)
    description.Draw()
    return description

def setTH1(hist, xAxisTitle, yAxisTitle, Ndevision, marker_style, marker_size, color):

    hist.GetXaxis().SetTitle(xAxisTitle)
    hist.GetYaxis().SetTitle(yAxisTitle)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetTitleSize(0.06)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetXaxis().SetNdivisions(Ndevision)
    hist.GetYaxis().SetTitleOffset(1.5)
    hist.GetXaxis().SetTitleOffset(0.9)

    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.SetMarkerStyle(marker_style)
    hist.SetMarkerSize(marker_size)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetLineWidth(2)

def setPad(canvas):
    gStyle.SetLineStyleString(22,"80 18 12 18 12 12")
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetBorderSize(0)
    canvas.SetTickx()
    canvas.SetTicky()
    canvas.SetFrameLineWidth(2)
    canvas.SetFrameBorderMode(0)
    canvas.SetFrameBorderSize(0)
    canvas.SetLeftMargin(0.17)
    canvas.SetRightMargin(0.07)
    canvas.SetTopMargin(0.074)
    canvas.SetBottomMargin(0.165)
    canvas.Range(-194.483,-10.3682,1041.38,-2.08469)

def geomAvg(h1, h2, minratio):
    name = TString(h1.GetName())
    name.ReplaceAll("NN", "Avg")
    name.ReplaceAll("PP", "Avg")
    nbins = h1.GetNbinsX()
    lowEdge = h1.GetXaxis().GetBinLowEdge(1)
    upEdge = h1.GetXaxis().GetBinUpEdge(nbins)
    result = TH1F(name.Data(), name.Data(), h1.GetNbinsX(), lowEdge, upEdge)
    for bin_x in range(1, h1.GetNbinsX()):
        for bin_y in range(1, h1.GetNbinsY()):
            for bin_z in range(1, h1.GetNbinsZ()):
                ratio1 = h1.GetBinContent(bin_x, bin_y, bin_z) /h2.GetBinContent(bin_x, bin_y, bin_z)
                ratio2 = h2.GetBinContent(bin_x, bin_y, bin_z) /h1.GetBinContent(bin_x, bin_y, bin_z)
                content = 2 * TMath.Sqrt(h1.GetBinContent(bin_x, bin_y, bin_z) * h2.GetBinContent(bin_x, bin_y, bin_z))
                x = h1.GetBinContent(bin_x, bin_y, bin_z)
                y = h2.GetBinContent(bin_x, bin_y, bin_z)
                sx = h1.GetBinError(bin_x, bin_y, bin_z)
                sy = h2.GetBinError(bin_x, bin_y, bin_z)
                dfdx2 = y / x
                dfdy2 = x / y
                error = TMath.Sqrt(dfdx2 * sx * sx + dfdy2 * sy * sy)
                if ratio1 <= minratio or ratio2 <= minratio or not TMath.Finite(ratio1) or not TMath.Finite(ratio2):
                    content = h1.GetBinContent(bin_x, bin_y, bin_z) + h2.GetBinContent(bin_x, bin_y, bin_z);
                    error = h1.GetBinError(bin_x, bin_y, bin_z) + h2.GetBinError(bin_x, bin_y, bin_z);

                result.SetBinContent(bin_x, bin_y, bin_z, content);
                result.SetBinError(bin_x, bin_y, bin_z, error);
    return result

def geomAvg1d(h1, h2, minratio):
    name = TString(h1.GetName())
    name.ReplaceAll("NN", "Avg")
    name.ReplaceAll("PP", "Avg")
    nbins = h1.GetNbinsX()
    lowEdge = h1.GetXaxis().GetBinLowEdge(1)
    upEdge = h1.GetXaxis().GetBinUpEdge(nbins)
    result = TH1F(name.Data(), name.Data(), h1.GetNbinsX(), lowEdge, upEdge)
    for bin_x in range(1, 1 + h1.GetNbinsX()):
        try:
            ratio1 = h1.GetBinContent(bin_x) /h2.GetBinContent(bin_x)
        except ZeroDivisionError:
            ratio1 = 0
        try:
            ratio2 = h2.GetBinContent(bin_x) /h1.GetBinContent(bin_x)
        except ZeroDivisionError:
            ratio2 = 0

        if ratio1 <= minratio or ratio2 <= minratio or not TMath.Finite(ratio1) or not TMath.Finite(ratio2):
            content = h1.GetBinContent(bin_x) + h2.GetBinContent(bin_x)
            error = h1.GetBinError(bin_x) + h2.GetBinError(bin_x)
        else:
            content = 2 * TMath.Sqrt(h1.GetBinContent(bin_x) * h2.GetBinContent(bin_x))
            x = h1.GetBinContent(bin_x)
            y = h2.GetBinContent(bin_x)
            sx = h1.GetBinError(bin_x)
            sy = h2.GetBinError(bin_x)
            dfdx2 = y / x
            dfdy2 = x / y
            error = TMath.Sqrt(dfdx2 * sx * sx + dfdy2 * sy * sy)

        result.SetBinContent(bin_x, content);
        result.SetBinError(bin_x, error);
    return result
