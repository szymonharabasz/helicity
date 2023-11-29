from ROOT import TPaveText, gStyle


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

