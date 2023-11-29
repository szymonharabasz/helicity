import math

from hist_template import set_pad, set_th1
from hist_utils import diff_hist

axis_title = r"$cos(\theta_e^{\gamma*})$"


def oneplot(xs, ax, tensor, title):
    # pos = ax.plot(xs, tensor, ymin=0, ymax=tensor.max())
    pos = ax.plot(xs, tensor)
    # ax.set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
    ax.set_title(title)
    ax.set_xlabel(axis_title)


can_cmp_ind = 0


def plotComparison(can, pad_nr1, pad_nr2, hist_mc, histData, hist_index, pull_title, bins):
    # print(f"PLOTTING: {histMC.GetName()} and {histData.GetName()}")
    global can_cmp_ind

    curr_bin = bins[hist_index]
    title = (f"{curr_bin.m_min} < #it{{M}}_{{ee}} < {curr_bin.m_max}, {curr_bin.z_min} < cos(#theta^{{CM}}_{{"
             f"#gamma*}}) < {curr_bin.z_max}")

    pad = can.cd(pad_nr1)
    set_pad(pad)
    pad.SetRightMargin(0.16)
    histData.GetXaxis().SetTitle("cos(#theta_{e}^{#gamma*})")
    histData.SetTitle(title)
    set_th1(histData, histData.GetXaxis().GetTitle(), f"dN/d{histData.GetXaxis().GetTitle()} (a.u.)",
            505, 20, 0.8, 1)
    histData.Draw()
    hist_mc.Scale(histData.Integral() / hist_mc.Integral())
    hist_mc.SetLineColor(2)
    hist_mc.Draw("SAMEHIST")
    pad = can.cd(pad_nr2)
    set_pad(pad)
    pad.SetRightMargin(0.16)
    hdiff = diff_hist(hist_mc, histData)
    hdiff.SetTitle(pull_title)
    set_th1(hdiff, hdiff.GetXaxis().GetTitle(), "Pull value", 505, 20, 0.8, 2)
    hdiff.Draw("HIST")
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
            meanX2 = meanX2 + math.pow(histMC.GetBinCenter(i + 1), 2)
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

