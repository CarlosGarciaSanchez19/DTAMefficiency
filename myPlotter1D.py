import ROOT
import os
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

status = "Simulation Preliminary"


def plot1DEff(pathToFile, outputPlotDir, plotting, plotscaffold, chamber, processing, effpath, fil, quality):

    src = os.path.realpath(pathToFile)

    if plotting['energy'] is None:
        title = chamber + ", " + processing + ", " + plotting['sample'] + ", " + plotting['PU'] + " PU"
    else:
        title = chamber + ", " + plotting['energy'] + ", " + plotting['PU'] + " PU"

    f = ROOT.TFile.Open(src, 'read')
    h_matched = f.Get("h" + plotscaffold + "_" + chamber + "_" + processing + "_matched") # loop over chambers
    h_matched.SetDirectory(0)
    h_total = f.Get("h" + plotscaffold + "_" + chamber + "_" + processing + "_total")
    h_total.SetDirectory(0)
    f.Close()

    lowlimityaxis = plotting['lowlimityaxis']
    highlimityaxis = plotting['highlimityaxis']
    lowlimitxaxis = plotting['lowlimitxaxis']
    highlimitxaxis = plotting['highlimitxaxis']
    yaxistitle = plotting['yaxistitle']  # e.g. 'Efficiency (adim.)'
    xaxistitle = plotting['xaxistitle']  # e.g. r"$\eta$"
    yaxistitleoffset = plotting['yaxistitleoffset']  # e.g. 1.5

    eff = ROOT.TEfficiency(h_matched, h_total)

    canvas = ROOT.TCanvas('canvas', '', 500, 500)

    eff.SetTitle("")
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(0.4)
    eff.SetMarkerColor(ROOT.kBlue)
    eff.SetLineColor(ROOT.kBlue)
    eff.Draw("APZ")

    # Force ROOT to draw and process the canvas before accessing painted graph
    canvas.Update()

    # Now get the painted TGraph and style it
    g = eff.GetPaintedGraph()
    if g:
        g.SetMarkerStyle(20)
        g.SetMarkerSize(0.4)
        g.SetMarkerColor(ROOT.kBlue)
        g.SetLineColor(ROOT.kBlue)

        # Set axis titles and limits on the graph's axes
        g.GetXaxis().SetTitle(xaxistitle)
        g.GetYaxis().SetTitle(yaxistitle)
        g.GetYaxis().SetTitleOffset(yaxistitleoffset)
        g.GetXaxis().SetLimits(lowlimitxaxis, highlimitxaxis)  # For graph
        g.SetMinimum(lowlimityaxis)  # For y-axis limits in graph
        g.SetMaximum(highlimityaxis)

    # Needed to apply the above changes
    canvas.Modified()
    canvas.Update()

    firsttex = ROOT.TLatex()
    firsttex.SetTextSize(0.03)
    firsttex.SetTextFont(42)
    firsttex.DrawLatexNDC(0.04,0.91,"#scale[1.25]{       #bf{CMS}} " + "#it{" + status + "}")
    firsttex.Draw("same")

    secondtext = ROOT.TLatex()
    secondtext.SetTextSize(0.035)
    secondtext.SetTextAlign(31)
    secondtext.DrawLatexNDC(0.90, 0.91, title)
    secondtext.Draw("same")


    f = ROOT.TFile(outputPlotDir + "/" + plotscaffold + "_" + fil + "_" + chamber + "_" + processing + "_" + quality + ".root", "RECREATE")
    eff.Write("eff")
    f.Close()
    canvas.SaveAs(outputPlotDir + "/" + plotscaffold + "_" + fil + "_" + chamber + "_" + processing + "_" + quality + ".pdf")

