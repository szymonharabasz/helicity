import math

from error_matrix import errors_1d_hists, errors_3d_hists
from hist_template import set_pad, set_th1
from hist_utils import diff_hist
import matplotlib.pyplot as plt
from ROOT import TCanvas, TLegend
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


def x_axis_properties(hist_mc, hist_data):
    n = 0
    mean_x2 = 0
    for i, (c1, c2) in enumerate(zip(hist_mc, hist_data)):
        if c1 != 0 and c2 != 0:
            n = n + 1
            mean_x2 = mean_x2 + math.pow(hist_mc.GetXaxis().GetBinCenter(i + 1), 2)
    if n > 0:
        mean_x2 = mean_x2 / n
    var_x2 = 0
    sigma2 = 0
    for i, (c1, c2) in enumerate(zip(hist_mc, hist_data)):
        if c1 != 0 and c2 != 0:
            center = hist_data.GetXaxis().GetBinCenter(i + 1)
            var_x2 = var_x2 + math.pow(math.pow(center, 2) - mean_x2, 2)
            sigma2 = sigma2 + math.pow((c2 - c1) / c1, 2)
    if n > 2:
        sigma2 = sigma2 / (n - 2)
    return n, mean_x2, var_x2, sigma2


def dist_properties(hist_mc, hist_data):
    n = 0
    mean_cos2_theta = 0
    mean_sin_2theta_cos_phi = 0
    mean_sin2_theta_cos_2phi = 0
    nx = hist_data.GetNbinsX()
    ny = hist_data.GetNbinsY()
    for bx in range(1, nx + 1):
        for by in range(1, ny + 1):
            c1 = hist_mc.GetBinContent(bx, by)
            c2 = hist_data.GetBinContent(bx, by)
            if c1 != 0 and c2 != 0:
                n = n + 1
                cos_theta = hist_mc.GetXaxis().GetBinCenter(bx)
                theta = math.acos(cos_theta)
                phi = hist_mc.GetYaxis().GetBinCenter(by)
                mean_cos2_theta = mean_cos2_theta + cos_theta ** 2
                mean_sin_2theta_cos_phi = mean_sin_2theta_cos_phi + math.sin(2*theta) * math.cos(phi)
                mean_sin2_theta_cos_2phi = mean_sin2_theta_cos_2phi + math.sin(theta) ** 2 * math.cos(2*phi)
    if n > 0:
        mean_cos2_theta = mean_cos2_theta / n
        mean_sin_2theta_cos_phi = mean_sin_2theta_cos_phi / n
        mean_sin2_theta_cos_2phi = mean_sin2_theta_cos_2phi / n
    var_cos2_theta = 0
    var_sin_2theta_cos_phi = 0
    var_sin2_theta_cos_2phi = 0
    sigma2 = 0
    for bx in range(1, nx + 1):
        for by in range(1, ny + 1):
            c1 = hist_mc.GetBinContent(bx, by)
            c2 = hist_data.GetBinContent(bx, by)
            if c1 != 0 and c2 != 0:
                cos_theta = hist_mc.GetXaxis().GetBinCenter(bx)
                theta = math.acos(cos_theta)
                phi = hist_mc.GetYaxis().GetBinCenter(by)
                var_cos2_theta = var_cos2_theta + (cos_theta ** 2 - mean_cos2_theta) ** 2
                var_sin_2theta_cos_phi = var_sin_2theta_cos_phi + (math.sin(2*theta) * math.cos(phi) - mean_sin_2theta_cos_phi) ** 2
                var_sin2_theta_cos_2phi = var_sin2_theta_cos_2phi + (math.sin(theta) ** 2 * math.cos(2*phi) - mean_sin2_theta_cos_2phi) ** 2
                sigma2 = sigma2 + math.pow((c2 - c1) / c1, 2)
    if n > 2:
        sigma2 = sigma2 / (n - 2)
    return n, mean_cos2_theta, var_cos2_theta, mean_sin_2theta_cos_phi, var_sin_2theta_cos_phi, mean_sin2_theta_cos_2phi, var_sin2_theta_cos_2phi, sigma2


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
            # norm = parameters_all[HIST_INDEX][1].item()
            best_hists_mc = fn_get_hist_maker_mc(sign, HIST_INDEX).make_hists(lambda_theta)
            hmodels.append(best_hists_mc[0][HIST_INDEX])

            if HIST_INDEX % 3 == 0:
                can1 = TCanvas(f"can_cmp_{HIST_INDEX}_{sign}", "can", 900, 600)
                if analyse_3d:
                    can1 = TCanvas(f"can_cmp_{HIST_INDEX}_{sign}", "can", 900, 900)
                    can1.Divide(3, 3)
                else:
                    can1 = TCanvas(f"can_cmp_{HIST_INDEX}_{sign}", "can", 900, 600)
                    can1.Divide(3, 2)
                can1.Draw()
                canvases.append(can1)

            if analyse_3d:
                hdiff1 = plot_comparison_2d(can1, HIST_INDEX % 3 + 1, HIST_INDEX % 3 + 4, HIST_INDEX % 3 + 7,
                                            best_hists_mc[0][HIST_INDEX],
                                            hists_data[0][HIST_INDEX], HIST_INDEX, "Residuals", bins)
            else:
                hdiff1 = plot_comparison(can1, HIST_INDEX % 3 + 1, HIST_INDEX % 3 + 4, best_hists_mc[0][HIST_INDEX],
                                         hists_data[0][HIST_INDEX], HIST_INDEX, "Residuals", bins)

            hdiffs.append(hdiff1)

            if analyse_3d:
                n, mean1, var1, mean2, var2, mean3, var3, sigma2 = dist_properties(best_hists_mc[0][HIST_INDEX], hists_data[0][HIST_INDEX])
                errB0 = math.sqrt(sigma2 * (1 / n + mean1 ** 2 / var1 + mean2 ** 2 / var2 + mean3 ** 2 / var3))
                errB1 = math.sqrt(sigma2 / var1)
                errB2 = math.sqrt(sigma2 / var2)
                errB3 = math.sqrt(sigma2 / var3)
                ratio_error1 = ratio_err(lambda_theta, 1, errB1, errB0)
                ratio_error2 = ratio_err(lambda_theta_phi, 1, errB2, errB0)
                ratio_error3 = ratio_err(lambda_phi, 1, errB3, errB0)
                print("Maybe covariance matrix:")
                hists_mc_null = fn_get_hist_maker_mc(sign, HIST_INDEX).make_hists(0.0)
                err_matr = errors_3d_hists(hists_data[0][HIST_INDEX], hists_mc_null[0][HIST_INDEX])
                print(err_matr)
                print("Compare:")
                errB0new = math.sqrt(err_matr[0, 0].item())
                errB1new = math.sqrt(err_matr[1, 1].item())
                errB2new = math.sqrt(err_matr[2, 2].item())
                errB3new = math.sqrt(err_matr[3, 3].item())
                ratio_new1 = ratio_err(lambda_theta, 1, errB1new, errB0new)
                ratio_new2 = ratio_err(lambda_theta_phi, 1, errB2new, errB0new)
                ratio_new3 = ratio_err(lambda_phi, 1, errB3new, errB0new)
                print(f"Old: {ratio_error1}, new: {ratio_new1}, new/old: {ratio_new1 / ratio_error1}")
                print(f"Old: {ratio_error2}, new: {ratio_new2}, new/old: {ratio_new2 / ratio_error2}")
                print(f"Old: {ratio_error3}, new: {ratio_new3}, new/old: {ratio_new3 / ratio_error3}")
                can1.cd(HIST_INDEX % 3 + 1)

                caption = f"#lambda_{{#theta}} = {lambda_theta:.2f} #pm {ratio_new1:.2f}"
                if HIST_INDEX < 6:
                    paveText = set_opt_text(caption, 0.25, 0.26, 0.675, 0.38, 2, 0.04)
                else:
                    paveText = set_opt_text(caption, 0.25, 0.76, 0.675, 0.88, 2, 0.04)
                # paveText.AddText(f"Norm = {norm:.5f}")
                paveText.AddText(f"#lambda_{{#phi}} = {lambda_phi:.2f} #pm {ratio_new3:.2f}")
                paveText.AddText(f"#lambda_{{#theta#phi}} = {lambda_theta_phi:.2f} #pm {ratio_new2:.2f}")
                paveTexts.append(paveText)

                if HIST_INDEX % 3 == 2:
                    can1.SaveAs(f"{dir_name}/comparison_{HIST_INDEX}.gif")



                try:
                    print(
                        f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_new1:.4f}, lambda_phi = {lambda_phi:.4f} +- {ratio_new3:.4f}, lambda_theta_phi = {lambda_theta_phi:.4f} +- {ratio_new2:.4f}")
                    print(
                        f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_new1:.4f}, lambda_phi = {lambda_phi:.4f} +- {ratio_new3:.4f}, lambda_theta_phi = {lambda_theta_phi:.4f} +- {ratio_new2:.4f}",
                        file=fout)
                except:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}", file=fout)
            else:
                n, meanX2, varX2, sigma2 = x_axis_properties(best_hists_mc[0][HIST_INDEX], hists_data[0][HIST_INDEX])
                errB0 = math.sqrt(sigma2 * (1 / n + meanX2 * meanX2 / varX2))
                errB1 = math.sqrt(sigma2 / varX2)
                ratio_error = ratio_err(lambda_theta, 1, errB1, errB0)

                can1.cd(HIST_INDEX % 3 + 1)

                caption = f"#lambda_{{#theta}} = {lambda_theta:.2f} #pm {ratio_error:.2f}"
                if HIST_INDEX < 6:
                    paveText = set_opt_text(caption, 0.25, 0.26, 0.675, 0.38, 2, 0.04)
                else:
                    paveText = set_opt_text(caption, 0.25, 0.76, 0.675, 0.88, 2, 0.04)
                # paveText.AddText(f"Norm = {norm:.5f}")
                paveTexts.append(paveText)

                if HIST_INDEX % 3 == 2:
                    can1.SaveAs(f"{dir_name}/comparison_{HIST_INDEX}.gif")

                print("Maybe covariance matrix:")
                hists_mc_null = fn_get_hist_maker_mc(sign, HIST_INDEX).make_hists(0.0)
                err_matr = errors_1d_hists(hists_data[0][HIST_INDEX], hists_mc_null[0][HIST_INDEX])
                print(err_matr)
                print("Compare:")
                errB0new = math.sqrt(err_matr[0, 0].item())
                errB1new = math.sqrt(err_matr[1, 1].item())
                ratio_new = ratio_err(lambda_theta, 1, errB1new, errB0new)
                print(f"Old: {errB0}, new: {errB0new}")
                print(f"Old: {errB1}, new: {errB1new}")
                print(f"Old: {ratio_error}, new: {ratio_new}, new/old: {ratio_new / ratio_error}")

                try:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f} +- {ratio_error:.4f}",
                          file=fout)
                except:
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}")
                    print(f"{HIST_INDEX}. Final result: lambda_theta = {lambda_theta:.4f}", file=fout)

            can1.SaveAs(f"{dir_name}/{can1.GetName()}.gif")
            
def show_mass_z(hists_data, histMakerMC_pi0, histMakerMC_rho, histMakerMC_mix, event_mixing, fraction, dir_name, sign):
    hmodelLowM_rho = histMakerMC_mix.make_hists(0.0) if event_mixing else histMakerMC_rho.make_hists(0.0)
    hmodelLowM_rho[2][0].SetLineColor(8)
    hmodelLowM_rho[1][0].SetLineColor(8)

    hmodelLowM = histMakerMC_pi0.make_hists(0.0)

    hmodelLowM[2][0].Scale(1.0 / hmodelLowM[2][0].Integral())
    hmodelLowM[1][0].Scale(1.0 / hmodelLowM[1][0].Integral())
    hmodelLowM_rho[2][0].Scale(1.0 / hmodelLowM_rho[2][0].Integral())
    hmodelLowM_rho[1][0].Scale(1.0 / hmodelLowM_rho[1][0].Integral())

    hmodelLowM[2][0].Add(hmodelLowM_rho[2][0], fraction)
    hmodelLowM[1][0].Add(hmodelLowM_rho[1][0], fraction)

    hmodelLowM[2][0].SetLineColor(2)
    hmodelLowM[1][0].SetLineColor(2)
    hmodelHigM = histMakerMC_rho.make_hists(0.0)
    hmodelHigM[2][1].SetLineColor(2)
    hmodelHigM[1][1].SetLineColor(2)

    can = TCanvas(f"can_mass_z_{sign}", "can", 800, 800)
    can.Divide(2, 2)
    can.Draw()

    pad = can.cd(1)
    set_pad(pad)
    dataScale = 1. / hists_data[2][0].Integral()
    hists_data[2][0].Scale(dataScale)
    hmodelLowM[2][0].Scale(1. / hmodelLowM[2][0].Integral())
    hmodelLowM_rho[2][0].Scale(1. / hmodelLowM_rho[2][0].Integral())
    hists_data[2][0].GetXaxis().SetTitle("cos(#theta^{CM}_{#gamma*})")
    hists_data[2][0].SetTitle("Masses below #pi^{0}")

    set_th1(hists_data[2][0], hists_data[2][0].GetXaxis().GetTitle(),
            f"dN/d{hists_data[2][0].GetXaxis().GetTitle()} (a.u.)", 505, 20, 0.8, 1)
    hists_data[2][0].Draw()
    hmodelLowM[2][0].Draw("SAMEHIST")
    hmodelLowM_rho[2][0].Draw("SAMEHIST")

    pad = can.cd(2)
    pad.SetLogy()
    set_pad(pad)
    hists_data[1][0].Scale(1. / hists_data[1][0].Integral())
    hmodelLowM[1][0].Scale(1. / hmodelLowM[1][0].Integral())
    hmodelLowM_rho[1][0].Scale(1. / hmodelLowM_rho[1][0].Integral())
    hists_data[1][0].SetTitle("Masses below #pi^{0}")
    hists_data[1][0].GetXaxis().SetTitle("#it{M}_{ee} (GeV/#it{c}^{2})")

    set_th1(hists_data[1][0], hists_data[1][0].GetXaxis().GetTitle(),
            f"dN/d{hists_data[1][0].GetXaxis().GetTitle()} (a.u.)", 505, 20, 0.8, 1)
    hists_data[1][0].Draw("HIST")
    hmodelLowM[1][0].Draw("SAMEHIST")
    hmodelLowM_rho[1][0].Draw("SAMEHIST")

    legend = TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(hists_data[1][0], "exp", "pl")
    legend.AddEntry(hmodelLowM[1][0], f"#pi^{{0}} + {fraction}#rho", "l")
    legend.AddEntry(hmodelLowM_rho[1][0], "#rho", "l")
    legend.Draw()

    pad = can.cd(3)
    set_pad(pad)
    hists_data[2][1].Scale(1. / hists_data[2][1].Integral())
    hmodelHigM[2][1].Scale(1. / hmodelHigM[2][1].Integral())
    hmodelLowM_rho[2][1].Scale(1. / hmodelLowM_rho[2][1].Integral())
    hists_data[2][1].GetXaxis().SetTitle("cos(#theta^{CM}_{#gamma*})")
    hists_data[2][1].SetTitle("Masses above #pi^{0}")

    set_th1(hists_data[2][1], hists_data[2][1].GetXaxis().GetTitle(),
            f"dN/d{hists_data[2][1].GetXaxis().GetTitle()} (a.u.)", 505, 20, 0.8, 1)
    hists_data[2][1].Draw()
    hmodelHigM[2][1].Draw("SAMEHIST")
    # hmodelLowM_rho[2][1].Draw("SAMEHIST")

    pad = can.cd(4)
    pad.SetLogy()
    set_pad(pad)
    hists_data[1][1].Scale(1. / hists_data[1][1].Integral())
    hmodelHigM[1][1].Scale(1. / hmodelHigM[1][1].Integral())
    hmodelLowM_rho[1][1].Scale(1. / hmodelLowM_rho[1][1].Integral())
    hists_data[1][1].GetXaxis().SetTitle("#it{M}_{ee} (GeV/#it{c}^{2})")
    hists_data[1][1].SetTitle("Masses above #pi^{0}")

    set_th1(hists_data[1][1], hists_data[1][1].GetXaxis().GetTitle(),
            f"dN/d{hists_data[1][1].GetXaxis().GetTitle()} (a.u.)", 505, 20, 0.8, 1)
    hists_data[1][1].Draw()
    hmodelHigM[1][1].Draw("SAMEHIST")
    # hmodelLowM_rho[1][1].Draw("SAMEHIST")

    can.SaveAs(f"{dir_name}/{can.GetName()}.gif")
    
    return can