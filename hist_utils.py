import math

from ROOT import TCanvas, TFile, TMath, TString, TH1F, TH2F

from eventsreader import SpdmeEventsReader, PickleEventsReader
from surrogatedistributionbuilder import SurrogateDistributionBuilder
from surrogatedistributionbuilder_1d import SurrogateDistributionBuilder_1d
from bins import Bins
import os.path
import torch

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
        filename_pickle = filename
        filename_pickle = filename_pickle.replace(".dat",".pkl")
        if os.path.isfile(filename_pickle):
            print("Pickle file exists")
            self.reader = PickleEventsReader(filename_pickle, frame, ekin)
        else:
            print("Picke file doesn't exist")
            self.reader = SpdmeEventsReader(filename, frame, ekin)
        # self.reader = SpdmeEventsReader(filename, frame, ekin)
        # self.builder = DistributionBuilder(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder(name_suffix, self.reader.get_events(), bins)

    def make_hists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
        self.builder.setParameters(lambda_theta, lambda_phi, lambda_theta_phi)
        return self.builder.getHists()


class HistMaker1d:
    def __init__(self, filename, name_suffix, bins, frame, ekin):
        filename_pickle = filename
        filename_pickle = filename_pickle.replace(".dat", ".pkl")
        if os.path.isfile(filename_pickle):
            print("Pickle file exists")
            self.reader = PickleEventsReader(filename_pickle, frame, ekin)
        else:
            print("Picke file doesn't exist")
            self.reader = SpdmeEventsReader(filename, frame, ekin)
        # self.builder = DistributionBuilder_1d(name_suffix, self.reader.getEvents(), bins)
        self.builder = SurrogateDistributionBuilder_1d(name_suffix, self.reader.get_events(), bins)

    def make_hists(self, lambda_theta=0):
        self.builder.setParameters(lambda_theta)
        return self.builder.getHists()


class CombinedHistMaker:
    def __init__(self, hist_maker_1, hist_maker_2, fraction):
        self.hist_maker_1 = hist_maker_1
        self.hist_maker_2 = hist_maker_2
        self.fraction = fraction

    def make_hists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):
        if isinstance(self.hist_maker_1, HistMaker) and isinstance(self.hist_maker_2, HistMaker):
            hists_1 = self.hist_maker_1.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_2 = self.hist_maker_2.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
        else:
            hists_1 = self.hist_maker_1.make_hists(lambda_theta)
            hists_2 = self.hist_maker_2.make_hists(lambda_theta)
        for h1, h2 in zip(hists_1, hists_2):
            if not isinstance(h1, list):
                h1.Add(h2, self.fraction * h1.Integral() / h2.Integral())
            else:
                for hh1, hh2 in zip(h1, h2):
                    hh1.Add(hh2, self.fraction * hh1.Integral() / hh2.Integral())
        return hists_1


class SignalHistMaker:
    def __init__(self, hist_maker_np, hist_maker_pp, hist_maker_nn, hist_maker_np_mix, hist_maker_pp_mix, hist_maker_nn_mix):
        self.hist_maker_np = hist_maker_np
        self.hist_maker_pp = hist_maker_pp
        self.hist_maker_nn = hist_maker_nn
        self.hist_maker_np_mix = hist_maker_np_mix
        self.hist_maker_pp_mix = hist_maker_pp_mix
        self.hist_maker_nn_mix = hist_maker_nn_mix

    def make_hists(self, lambda_theta=0, lambda_phi=0, lambda_theta_phi=0):

        hists_avg = []
        hists_k = []

        def subtr_bgr(hnp, hpp, hnn, hnp_mix, hpp_mix, hnn_mix):
            fnc_avg = geom_avg1d if isinstance(hpp, TH1F) else geom_avg2d

            havg = fnc_avg(hpp, hnn, 0.2)
            havg_mix = fnc_avg(hpp_mix, hnn_mix, 0.2)
            hnp_mix.Divide(havg_mix)
            hnp_mix.Scale(2.0)
            havg.Multiply(hnp_mix)

            hnp.Add(havg, -1.0)
            return havg, hnp_mix
           # for bin in range(1, hnp.GetNbinsX() + 1):
           #     hnp.SetBinContent(bin,havg.GetBinContent(bin))
           #     hnp.SetBinError(bin,havg.GetBinError(bin))

        if isinstance(self.hist_maker_np, HistMaker) and isinstance(self.hist_maker_pp, HistMaker) and isinstance(self.hist_maker_pp, HistMaker):
            hists_np = self.hist_maker_np.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_pp = self.hist_maker_pp.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_nn = self.hist_maker_nn.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_np_mix = self.hist_maker_np_mix.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_pp_mix = self.hist_maker_pp_mix.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
            hists_nn_mix = self.hist_maker_nn_mix.make_hists(lambda_theta, lambda_phi, lambda_theta_phi)
        else:
            hists_np = self.hist_maker_np.make_hists(lambda_theta)
            hists_pp = self.hist_maker_pp.make_hists(lambda_theta)
            hists_nn = self.hist_maker_nn.make_hists(lambda_theta)
            hists_np_mix = self.hist_maker_np_mix.make_hists(lambda_theta)
            hists_pp_mix = self.hist_maker_pp_mix.make_hists(lambda_theta)
            hists_nn_mix = self.hist_maker_nn_mix.make_hists(lambda_theta)
        for hnp, hpp, hnn, hnp_mix, hpp_mix, hnn_mix in zip(hists_np, hists_pp, hists_nn, hists_np_mix, hists_pp_mix, hists_nn_mix):
            if not isinstance(hnp, list):
                hist_avg, hist_k = subtr_bgr(hnp, hpp, hnn, hnp_mix, hpp_mix, hnn_mix)
                hists_avg.append(hist_avg)
                hists_k.append(hist_k)
            else:
                hists_avg.append([])
                hists_k.append([])
                for hhnp, hhpp, hhnn, hhnp_mix, hhpp_mix, hhnn_mix in zip(hnp, hpp, hnn, hnp_mix, hpp_mix, hnn_mix):
                    hist_avg, hist_k = subtr_bgr(hhnp, hhpp, hhnn, hhnp_mix, hhpp_mix, hhnn_mix)
                    hists_avg[-1].append(hist_avg)
                    hists_k[-1].append(hist_k)
        return hists_np, hists_avg, hists_k

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

def geom_avg2d(h1, h2, minratio):
    name = TString(h1.GetName())
    name.ReplaceAll("NN", "Avg")
    name.ReplaceAll("PP", "Avg")
    nbins = h1.GetNbinsX()
    low_edge_x = h1.GetXaxis().GetBinLowEdge(1)
    up_edge_x = h1.GetXaxis().GetBinUpEdge(nbins)
    low_edge_y = h1.GetYaxis().GetBinLowEdge(1)
    up_edge_y = h1.GetYaxis().GetBinUpEdge(nbins)
    result = TH2F(name.Data(), name.Data(), h1.GetNbinsX(), low_edge_x, up_edge_x, h1.GetNbinsY(), low_edge_y, up_edge_y)
    for bin_x in range(1, h1.GetNbinsX()):
        for bin_y in range(1, h1.GetNbinsY()):
            content1 = h1.GetBinContent(bin_x, bin_y)
            content2 = h2.GetBinContent(bin_x, bin_y)
            if content1 != 0 and content2 != 0:
                ratio1 = content1 / content2
                ratio2 = content2 / content1
                content = 2 * TMath.Sqrt( content1 * content2 )
                x = content1
                y = content2
                sx = h1.GetBinError(bin_x, bin_y)
                sy = h2.GetBinError(bin_x, bin_y)
                dfdx2 = y / x
                dfdy2 = x / y
                error = TMath.Sqrt(dfdx2 * sx * sx + dfdy2 * sy * sy)
                if ratio1 <= minratio or ratio2 <= minratio or not TMath.Finite(ratio1) or not TMath.Finite(ratio2):
                    content = h1.GetBinContent(bin_x, bin_y) + h2.GetBinContent(bin_x, bin_y)
                    error = h1.GetBinError(bin_x, bin_y) + h2.GetBinError(bin_x, bin_y)

                result.SetBinContent(bin_x, bin_y, content)
                result.SetBinError(bin_x, bin_y, error)
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


def hist_to_tensor(hist):
    def hist_to_tensor_1d(hist):
        n = hist.GetNbinsX()
        columns = []
        for b in range(1, n + 1):
            t = torch.tensor([
                hist.GetBinCenter(b),
                hist.GetBinContent(b),
                hist.GetBinError(b)
            ])
            columns.append(t)
        result = torch.stack(columns).transpose(0, 1)
        return result

    def hist_to_tensor_2d(hist):
        nx = hist.GetNbinsX()
        ny = hist.GetNbinsY()
        columns = []
        for bx in range(1, nx + 1):
            for by in range(1, ny + 1):
                t = torch.tensor([
                    hist.GetXaxis().GetBinCenter(bx),
                    hist.GetXaxis().GetBinCenter(by),
                    hist.GetBinContent(bx, by),
                    hist.GetBinError(bx, by)
                ])
                columns.append(t)
        result = torch.stack(columns).transpose(0, 1)
        return result

    if isinstance(hist, TH1F):
        return hist_to_tensor_1d(hist)
    elif isinstance(hist, TH2F):
        return hist_to_tensor_2d(hist)
    else:
        print("hist is none of these")
        return None

