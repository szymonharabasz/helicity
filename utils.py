import math
from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from surrogatedistributionbuilder import SurrogateDistributionBuilder
from distributionbuilder_1d import DistributionBuilder_1d
from surrogatedistributionbuilder_1d import SurrogateDistributionBuilder_1d
from ROOT import TCanvas, TFile, TPaveText, gStyle, TMath, TString, TH1F


def calc_diff(hist_mc, hist_data, bx, by):
    content_mc = hist_mc.GetBinContent(bx + 1, by + 1) / hist_mc.Integral()
    content_data = hist_data.GetBinContent(bx + 1, by + 1) / hist_data.Integral()
    error_data = hist_data.GetBinError(bx + 1, by + 1) / hist_data.Integral()
    if content_mc > 0 and content_data > 0:
        # return pow((content_data - content_mc)/error_data, 2)
        return pow((content_data - content_mc) / error_data, 2)
    else:
        return 0.0


def calc_one_chi2(hist_mc, hist_data):
    chi2 = 0.
    ndf = 0
    for bx in range(hist_mc.GetXaxis().GetNbins()):
        for by in range(hist_mc.GetYaxis().GetNbins()):
            diff = calc_diff(hist_mc, hist_data, bx, by)
            if diff > 0:
                ndf = ndf + 1
                chi2 = chi2 + diff
    return chi2, ndf


def diff_hist(hist_mc, hist_data):
    name = hist_mc.GetName() + "_diff"
    hdiff = hist_mc.Clone(name)
    for bx in range(hist_mc.GetXaxis().GetNbins()):
        for by in range(hist_mc.GetYaxis().GetNbins()):
            diff = calc_diff(hist_mc, hist_data, bx, by)
            hdiff.SetBinContent(bx + 1, by + 1, diff)
    return hdiff


def calc_all_chi2(hists_mc, hists_data):
    for i in range(len(hists_mc[0])):
        hist_mc = hists_mc[0][i]
        hist_data = hists_data[0][i]

        chi2, ndf = calc_one_chi2(hist_mc, hist_data)
        print("chi2, ndf, chi2/ndf", chi2, ndf, chi2 / ndf)
        yield chi2, ndf


class HistMaker:
    def __init__(self, filename, name_suffix, bins, frame, ekin):
        self.reader = SpdmeEventsReader(filename, frame, ekin)
        # self.builder = DistributionBuilder(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder(name_suffix, self.reader.get_events(), bins)

    def make_hists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
        self.builder.setParameters(lambda_theta, lambda_phi, lambda_theta_phi)
        return self.builder.getHists()


class HistMaker1d:
    def __init__(self, filename, name_suffix, bins, frame, ekin):
        self.reader = SpdmeEventsReader(filename, frame, ekin)
        # self.builder = DistributionBuilder_1d(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder_1d(name_suffix, self.reader.get_events(), bins)

    def make_hists(self, lambda_theta=0):
        self.builder.setParameters(lambda_theta)
        return self.builder.getHists()


def make_root_plots(hists_mc, hists_data):
    c = TCanvas("c", "c", 100, 100, 1200, 1600)
    c.Divide(3, 4)

    totentries = 0
    file = TFile("out.root", "RECREATE")
    file.cd()
    for i in range(len(hists_mc[0])):
        hist_mc = hists_mc[0][i]
        c.cd(i + 1)
        hist_mc.Draw("COLZ")
        totentries = totentries + hist_mc.GetEntries()
        hist_mc.Write()

        hist_data = hists_data[0][i]
        hist_data.Write()

    hists_mc[1].Write()
    hists_mc[2].Write()
    hists_data[1].Write()
    hists_data[2].Write()
    c.SaveAs("c.gif")
    print("Total hist entries: ", totentries)


def set_opt_text(cut_desc, x1, y1, x2, y2, color, text_size):
    description = TPaveText(x1, y1, x2, y2, "NDC")
    description.SetLineWidth(0)
    description.AddText(cut_desc)
    description.SetTextSize(text_size)
    description.SetBorderSize(0)
    description.SetTextFont(62)
    description.SetTextColor(color)
    description.SetFillColor(0)
    description.SetFillStyle(0)
    description.Draw()
    return description


def set_th1(hist, x_axis_title, y_axis_title, ndevision, marker_style, marker_size, color):
    hist.GetXaxis().SetTitle(x_axis_title)
    hist.GetYaxis().SetTitle(y_axis_title)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetTitleSize(0.06)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetXaxis().SetNdivisions(ndevision)
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


def set_pad(canvas):
    gStyle.SetLineStyleString(22, "80 18 12 18 12 12")
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
    canvas.Range(-194.483, -10.3682, 1041.38, -2.08469)


def geom_avg(h1, h2, minratio):
    name = TString(h1.GetName())
    name.ReplaceAll("NN", "Avg")
    name.ReplaceAll("PP", "Avg")
    nbins = h1.GetNbinsX()
    low_edge = h1.GetXaxis().GetBinLowEdge(1)
    up_edge = h1.GetXaxis().GetBinUpEdge(nbins)
    result = TH1F(name.Data(), name.Data(), h1.GetNbinsX(), low_edge, up_edge)
    for bin_x in range(1, h1.GetNbinsX()):
        for bin_y in range(1, h1.GetNbinsY()):
            for bin_z in range(1, h1.GetNbinsZ()):
                ratio1 = h1.GetBinContent(bin_x, bin_y, bin_z) / h2.GetBinContent(bin_x, bin_y, bin_z)
                ratio2 = h2.GetBinContent(bin_x, bin_y, bin_z) / h1.GetBinContent(bin_x, bin_y, bin_z)
                content = 2 * TMath.Sqrt(h1.GetBinContent(bin_x, bin_y, bin_z) * h2.GetBinContent(bin_x, bin_y, bin_z))
                x = h1.GetBinContent(bin_x, bin_y, bin_z)
                y = h2.GetBinContent(bin_x, bin_y, bin_z)
                sx = h1.GetBinError(bin_x, bin_y, bin_z)
                sy = h2.GetBinError(bin_x, bin_y, bin_z)
                dfdx2 = y / x
                dfdy2 = x / y
                error = TMath.Sqrt(dfdx2 * sx * sx + dfdy2 * sy * sy)
                if ratio1 <= minratio or ratio2 <= minratio or not TMath.Finite(ratio1) or not TMath.Finite(ratio2):
                    content = h1.GetBinContent(bin_x, bin_y, bin_z) + h2.GetBinContent(bin_x, bin_y, bin_z)
                    error = h1.GetBinError(bin_x, bin_y, bin_z) + h2.GetBinError(bin_x, bin_y, bin_z)

                result.SetBinContent(bin_x, bin_y, bin_z, content)
                result.SetBinError(bin_x, bin_y, bin_z, error)
    return result


def geom_avg1d(h1, h2, minratio):
    name = TString(h1.GetName())
    name.ReplaceAll("NN", "Avg")
    name.ReplaceAll("PP", "Avg")
    nbins = h1.GetNbinsX()
    low_edge = h1.GetXaxis().GetBinLowEdge(1)
    up_edge = h1.GetXaxis().GetBinUpEdge(nbins)
    result = TH1F(name.Data(), name.Data(), h1.GetNbinsX(), low_edge, up_edge)
    for bin_x in range(1, 1 + h1.GetNbinsX()):
        try:
            ratio1 = h1.GetBinContent(bin_x) / h2.GetBinContent(bin_x)
        except ZeroDivisionError:
            ratio1 = 0
        try:
            ratio2 = h2.GetBinContent(bin_x) / h1.GetBinContent(bin_x)
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

        result.SetBinContent(bin_x, content)
        result.SetBinError(bin_x, error)
    return result


def ratio_err(val_a, val_b, err_a, err_b):
    err2 = math.pow(err_a / val_b, 2.) + math.pow(val_a * err_b / (val_b * val_b), 2.)
    err = math.sqrt(err2)
    return err


def symmetrize(hist):
    n = hist.GetNbinsX()
    i = n // 2 + 1
    while i <= n:
        j = n - i + 1
        ci = hist.GetBinContent(i)
        cj = hist.GetBinContent(j)
        ei = hist.GetBinError(i)
        ej = hist.GetBinError(j)

        cavg = 0.5 * (ci + cj)
        eavg = 0.5 * (ei + ej)

        hist.SetBinContent(i, cavg)
        hist.SetBinContent(j, cavg)
        hist.SetBinError(i, eavg)
        hist.SetBinError(j, eavg)

        i = i + 1
