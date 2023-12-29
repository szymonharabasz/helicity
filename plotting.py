import math

from hist_template import set_pad, set_th1
from hist_utils import diff_hist
import matplotlib.pyplot as plt
from ROOT import TCanvas, TH2F
from hist_template import set_opt_text
from hist_utils import ratio_err

axis_title = r"$cos(\theta_e^{\gamma*})$"


def oneplot(xs, ax, tensor, title):
    # pos = ax.plot(xs, tensor, ymin=0, ymax=tensor.max())
    pos = ax.plot(xs, tensor)
    # ax.set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
    ax.set_title(title)
    ax.set_xlabel(axis_title)


can_cmp_ind = 0


def plot_comparison(can, pad_nr1, pad_nr2, hist_mc, hist_data, hist_index, pull_title, bins):
    global can_cmp_ind

    curr_bin = bins[hist_index]
    title = (f"{curr_bin.m_min} < #it{{M}}_{{ee}} < {curr_bin.m_max}, {curr_bin.z_min} < cos(#theta^{{CM}}_{{"
             f"#gamma*}}) < {curr_bin.z_max}")

    plot_one_comparison(can, pad_nr1, hist_mc, hist_data, title, "cos(#theta_{e}^{#gamma*})")

    pad = can.cd(pad_nr2)
    set_pad(pad)
    pad.SetRightMargin(0.16)
    hdiff = diff_hist(hist_mc, hist_data)
    hdiff.SetTitle(pull_title)
    set_th1(hdiff, hdiff.GetXaxis().GetTitle(), "Pull value", 505, 20, 0.8, 2)
    hdiff.Draw("HIST")
    can.Update()
    can.Modified()
    can.Update()

    return hdiff


def plot_one_comparison(can, pad_nr, hist_mc, hist_data, title, xtitle):
    pad = can.cd(pad_nr)
    set_pad(pad)
    pad.SetRightMargin(0.16)
    hist_data.GetXaxis().SetTitle(xtitle)
    hist_data.SetTitle(title)
    set_th1(hist_data, hist_data.GetXaxis().GetTitle(), f"dN/d{hist_data.GetXaxis().GetTitle()} (a.u.)",
            505, 20, 0.8, 1)
    hist_data.Draw()
    hist_mc.Scale(hist_data.Integral() / hist_mc.Integral())
    hist_mc.SetLineColor(2)
    hist_mc.Draw("SAMEHIST")


def plot_comparison_2d(can, pad_nr1, pad_nr2, pad_nr3, hist_mc, hist_data, hist_index, pull_title, bins):
    global can_cmp_ind

    curr_bin = bins[hist_index]
    title = (f"{curr_bin.m_min} < #it{{M}}_{{ee}} < {curr_bin.m_max}, {curr_bin.z_min} < cos(#theta^{{CM}}_{{"
             f"#gamma*}}) < {curr_bin.z_max}")

    hist_data_px = hist_data.ProjectionX()
    hist_mc_px = hist_mc.ProjectionX()
    plot_one_comparison(can, pad_nr1, hist_mc_px, hist_data_px, title, "cos(#theta_{e}^{#gamma*})")
    hist_data_py = hist_data.ProjectionY()
    hist_mc_py = hist_mc.ProjectionY()
    plot_one_comparison(can, pad_nr2, hist_mc_py, hist_data_py, title, "#phi_{e}^{#gamma*}")

    pad = can.cd(pad_nr3)
    set_pad(pad)
    pad.SetRightMargin(0.16)
    hdiff = diff_hist(hist_mc, hist_data)
    hdiff.SetTitle(pull_title)
    set_th1(hdiff, hdiff.GetXaxis().GetTitle(), hdiff.GetYaxis().GetTitle(), 505, 20, 0.8, 2)
    hdiff.Draw("COLZ")
    can.Update()
    can.Modified()
    can.Update()

    return hdiff


def xAxisProperties(histMC, histData):
    n = 0
    meanX2 = 0
    for i, (c1, c2) in enumerate(zip(histMC, histData)):
        if c1 != 0 and c2 != 0:
            n = n + 1
            meanX2 = meanX2 + math.pow(histMC.GetXaxis().GetBinCenter(i + 1), 2)
    if n > 0:
        meanX2 = meanX2 / n
    varX2 = 0
    sigma2 = 0
    for i, (c1, c2) in enumerate(zip(histMC, histData)):
        if c1 != 0 and c2 != 0:
            center = histData.GetXaxis().GetBinCenter(i + 1)
            varX2 = varX2 + math.pow(math.pow(center, 2) - meanX2, 2)
            sigma2 = sigma2 + math.pow((c2 - c1) / c1, 2)
    if n > 2:
        sigma2 = sigma2 / (n - 2)
    return n, meanX2, varX2, sigma2

    def bin_index(x, min, max):
        return int((x - min) / (max - min) * 101)

def plot_losses(losses_all, range_used):
    fig, ax = plt.subplots(nrows=4, ncols=3)
    fig.tight_layout()
    fig.set_figheight(20)
    fig.set_figwidth(15)
    for HIST_INDEX in range_used:
        # print(losses_all[HIST_INDEX])
        ax[HIST_INDEX // 3][HIST_INDEX % 3].plot(losses_all[HIST_INDEX])
        ax[HIST_INDEX // 3][HIST_INDEX % 3].set_xlabel("Epoch")
        ax[HIST_INDEX // 3][HIST_INDEX % 3].set_ylabel("Loss ($\chi^2$)")
        ax[HIST_INDEX // 3][HIST_INDEX % 3].set_xscale("log")
        ax[HIST_INDEX // 3][HIST_INDEX % 3].set_yscale("log")
    return fig, ax

canvases = []
hdiffs = []
hmodels = []
paveTexts = []
def show_results(sign, dir_name, range_used, parameters_all, fn_get_hist_maker_mc, bins, hists_data, analyse_3d):
    with open(f'{dir_name}/results_{sign}.txt', 'w') as fout:
        for HIST_INDEX in range_used:

            # ax = plt.axes()
            # fig, ax = plt.subplots(nrows=1, ncols=1)
            lambda_theta = parameters_all[HIST_INDEX][0].item()
            if analyse_3d:
                lambda_phi = parameters_all[HIST_INDEX][1].item()
                lambda_theta_phi = parameters_all[HIST_INDEX][2].item()
            #norm = parameters_all[HIST_INDEX][1].item()
            best_hists_mc = fn_get_hist_maker_mc(sign, HIST_INDEX).make_hists(lambda_theta)
            hmodels.append(best_hists_mc[0][HIST_INDEX])

            if HIST_INDEX % 3 == 0:
                can1 = TCanvas(f"can_cmp_{HIST_INDEX}", "can", 900, 600)
                if analyse_3d:
                    can1 = TCanvas(f"can_cmp_{HIST_INDEX}", "can", 900, 900)
                    can1.Divide(3, 3)
                else:
                    can1 = TCanvas(f"can_cmp_{HIST_INDEX}", "can", 900, 600)
                    can1.Divide(3, 2)
                can1.Draw()
                canvases.append(can1)

            if analyse_3d:
                hdiff1 = plot_comparison_2d(can1, HIST_INDEX % 3 + 1, HIST_INDEX % 3 + 4, HIST_INDEX % 3 + 7, best_hists_mc[0][HIST_INDEX],
                                        hists_data[0][HIST_INDEX], HIST_INDEX, "Residuals", bins)
            else:
                hdiff1 = plot_comparison(can1, HIST_INDEX % 3 + 1, HIST_INDEX % 3 + 4, best_hists_mc[0][HIST_INDEX],
                                         hists_data[0][HIST_INDEX], HIST_INDEX, "Residuals", bins)

            hdiffs.append(hdiff1)

            n, meanX2, varX2, sigma2 = xAxisProperties(best_hists_mc[0][HIST_INDEX], hists_data[0][HIST_INDEX])
            errB0 = math.sqrt(sigma2 * (1 / n + meanX2 * meanX2 / varX2))
            errB1 = math.sqrt(sigma2 / varX2)
            ratio_error = ratio_err(lambda_theta, 1, errB1, errB0)

            can1.cd(HIST_INDEX % 3 + 1)

            caption = f"#lambda_{{#theta}} = {lambda_theta:.2f} #pm {ratio_error:.2f}"
            if HIST_INDEX < 6:
                paveText = set_opt_text(caption, 0.25, 0.26, 0.675, 0.38, 2, 0.04)
            else:
                paveText = set_opt_text(caption, 0.25, 0.76, 0.675, 0.88, 2, 0.04)
            #paveText.AddText(f"Norm = {norm:.5f}")
            if analyse_3d:
                paveText.AddText(f"#lambda_{{#phi}} = {lambda_phi:.2f}")
                paveText.AddText(f"#lambda_{{#theta#phi}} = {lambda_theta_phi:.2f}")
            paveTexts.append(paveText)

            if HIST_INDEX % 3 == 2:
                can1.SaveAs(f"{dir_name}/comparison_{HIST_INDEX}.gif")

            if analyse_3d:
                try:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}, lambda_phi = {lambda_phi:.4f}, lambda_theta_phi = {lambda_theta_phi:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}, lambda_phi = {lambda_phi:.4f}, lambda_theta_phi = {lambda_theta_phi:.4f}",
                          file=fout)
                except:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}", file=fout)
            else:
                try:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}", file=fout)
                except:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}", file=fout)
