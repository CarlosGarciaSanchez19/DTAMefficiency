import ROOT
import array
import os
from copy import deepcopy
import CMS_lumi
import sys
from lookForPU import lookForPU
from markerColors import markerColors
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

status = "Simulation Preliminary"


def plot2DEff(pathToFile, outputPlotDir, plotting, plotscaffold, chamber, processing, effpath, fil, quality):

    src = os.path.realpath(pathToFile)

    if plotting['energy'] == '':
        title = chamber + ", " + processing + ", " + plotting['sample'] + ", " + plotting['PU'] + " PU"
    else:
        title = "PU " + plotting['PU'] + " (" + plotting['energy'] + ")"

    f = ROOT.TFile.Open(src, 'read')
    h_matched = f.Get("h" + plotscaffold + "_" + chamber + "_" + processing + "_matched") # loop over chambers
    h_matched.SetDirectory(0)
    h_total = f.Get("h" + plotscaffold + "_" + chamber + "_" + processing + "_total")
    h_total.SetDirectory(0)
    f.Close()

    # Clone h_total to create h_eff
    h_eff = h_total.Clone("h_eff")
    h_eff.Reset("ICES")
    for i in range(1, h_eff.GetNbinsX() + 1):
        for j in range(1, h_eff.GetNbinsY() + 1):
            num = h_matched.GetBinContent(i, j)
            den = h_total.GetBinContent(i, j)
            if num == 0 and den != 0:
                h_eff.SetBinContent(i, j, 0.0001)
            elif den != 0:
                h_eff.SetBinContent(i, j, num / den)
            else:
                h_eff.SetBinContent(i, j, 0)

    lowlimityaxis = plotting['lowlimityaxis']
    highlimityaxis = plotting['highlimityaxis']
    lowlimitxaxis = plotting['lowlimitxaxis']
    highlimitxaxis = plotting['highlimitxaxis']
    yaxistitle = plotting['yaxistitle']  # e.g. r"$\phi$"
    xaxistitle = plotting['xaxistitle']  # e.g. r"$\eta$"
    yaxistitleoffset = plotting['yaxistitleoffset']  # e.g. 1.5

    canvas = ROOT.TCanvas('canvas', '', 500, 500)

    canvas.SetRightMargin(0.15)

    nRGBs = 3
    nColors = 255

    # Puntos de control
    stops = array.array('d', [0.00, 0.95, 1.00])
    red   = array.array('d', [1.00, 1.00, 0.00])
    green = array.array('d', [0.00, 1.00, 1.00])
    blue  = array.array('d', [0.00, 0.00, 0.00])

    # Crea la tabla de colores personalizados
    FstColorIdx = ROOT.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, nColors)
    ROOT.gStyle.SetNumberContours(nColors)

    palette = array.array('i', [FstColorIdx + i for i in range(nColors)])
    ROOT.gStyle.SetPalette(nColors, palette)
    # ROOT.TColor.InvertPalette();

    h_eff.SetContour(75)

    h_eff.SetTitle("")
    h_eff.GetXaxis().SetTitle(xaxistitle)
    h_eff.GetYaxis().SetTitle(yaxistitle)
    h_eff.GetXaxis().SetTitleSize(0.045)
    h_eff.GetYaxis().SetTitleSize(0.045)
    h_eff.GetYaxis().SetTitleOffset(yaxistitleoffset)
    h_eff.GetYaxis().SetRangeUser(lowlimityaxis, highlimityaxis)
    h_eff.GetXaxis().SetRangeUser(lowlimitxaxis, highlimitxaxis)
    h_eff.GetZaxis().SetTitle("Efficency")
    h_eff.GetZaxis().SetTitleOffset(1.5)

    h_eff.SetMinimum(plotting['zmin'])
    h_eff.SetMaximum(plotting['zmax'])
    h_eff.Draw("COLZ")

    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    firsttex = ROOT.TLatex()
    firsttex.SetTextSize(0.035)
    firsttex.SetTextFont(42)
    firsttex.DrawLatexNDC(0.04,0.91,"#scale[1.25]{      #bf{CMS}} " + "#it{" + status + "}")
    firsttex.Draw("same")

    secondtext = ROOT.TLatex()
    secondtext.SetTextSize(0.035)
    secondtext.SetTextAlign(31)
    secondtext.DrawLatexNDC(0.90, 0.91, title)
    secondtext.Draw("same")

    thirdtext = ROOT.TLatex()
    thirdtext.SetTextSize(0.035)
    thirdtext.SetTextAlign(31)
    thirdtext.DrawLatexNDC(0.20, 0.85, chamber)
    thirdtext.Draw("same")

    canvas.SaveAs(outputPlotDir + "/" + plotting['OutputFilename'] + "_" + fil + "_" + chamber + "_" + processing + "_" + quality + ".pdf")
    f = ROOT.TFile.Open(outputPlotDir + "/" + plotting['OutputFilename'] + "_" + fil + "_" + chamber + "_" + processing + "_" + quality + ".root", 'RECREATE')
    h_eff.Write()
    f.Close()

