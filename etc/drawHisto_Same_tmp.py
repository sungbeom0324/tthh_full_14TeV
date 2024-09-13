import ROOT
import numpy as np
import os
ROOT.gStyle.SetOptStat(0)

# RDF
PRE = "FULL_1218_14TeV"
Tree = "Delphes"

# Input
indir = "./samples1/" 
tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh_di.root")
tthbb  = ROOT.RDataFrame(Tree, indir + PRE + "_tthbb_di.root")
ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb_di.root")
ttbbcc = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbcc_di.root")
ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb_di.root")
dfs = {"_tthh": tthh, "_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb, "_tthbb": tthbb}

# OutDirectory. Modify!!
outdir = "./plots/" + PRE + "/Same/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHisto(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}
        legend = ROOT.TLegend(0.6, 0.855, 0.85, 0.705)
        hist_title = hist_name.replace("_size", " multiplicity")
        hist_title = hist_title.replace("_", " ")
        hist_title = "m_{H} (min #chi^{2})" # MODIFY # 

        ymax, color = 0, 1
        for df_name, df in dfs.items():
            nbin, xmin, xmax = 25, 0, 250
            _xmax = df.Max(hist_name).GetValue()
#            if _xmax < 100: xmax = 100
#            if _xmax < 20: nbin, xmin, xmax = 9, 0, 9 #int(_xmax+2), int(_xmax+2)
#            if _xmax < 5: nbin, xmin, xmax = 10, -5, 5
#            if _xmax <= 1: nbin, xmin, xmax = 50, 0, 0.006
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            if ymax < h.GetMaximum(): ymax = h.GetMaximum()
            h.GetXaxis().SetTitle("Higgs Mass [GeV]")
            h.GetYaxis().SetTitle("Normalized Entries")
            h.GetYaxis().SetTitleOffset(1.6)           
            h.SetLineColor(color)
            color+=1
            if color in [5] :
                color += 1
            h.SetLineWidth(2)
            legend.AddEntry(h.GetValue(), df_name, "l")
            hist_dict[hist_name + "_" + df_name] = h

        first = True
        for _tmp, h in hist_dict.items():
            h.SetMaximum(ymax * 1.1)
            if first:
                h.SetFillColorAlpha(ROOT.kGray, 0.8)
                h.DrawNormalized("hist E")
                first = False
            else: h.DrawNormalized("same E")
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetTextSize(0.025)  # Adjust the text size as needed
        latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")  # Modify the text and position

        canvas.Print(outdir + PRE + "_" + hist_name + flag + ".pdf")
        canvas.Clear()


# Histogram Config #
hists_S0 = [
"close_Higgs_mass"
#    "GenHiggs_mass", 
#    "GenHiggs_bdr",
#    "GenHiggs_Ht",
#    "GenHiggs_dEta",
#    "GenHiggs_dPhi",
#    "GenHiggs_mbmb",
]

hists_S1 = hists_S0 + [
]

hists_S2 = hists_S1 + [
]

hists_S3 = hists_S2 + [
#    "bJet1_pt", "bJet2_pt", "bJet3_pt"
]


# Draw #
drawHisto(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton
#for dfname, df in dfs.items():
#    df = df.Filter("GenbJet_size > 3")
#           .Define("GenbJet1_pt", "GenbJet_pt[0]")\
#           .Define("GenbJet2_pt", "GenbJet_pt[1]")\
#           .Define("GenbJet3_pt", "GenbJet_pt[2]")\
#           .Define("GenbJet4_pt", "GenbJet_pt[3]")
#    dfs[dfname] = df
#drawHisto(hists_S1, dfs, "_S1")

# ES2 : Jet
#for dfname, df in dfs.items():
#    df = df.Filter("Jet_size >= 5")\
#           .Define("Jet1_pt", "Jet_pt[0]")\
#           .Define("Jet2_pt", "Jet_pt[1]")\
#           .Define("Jet3_pt", "Jet_pt[2]")\
#           .Define("Jet4_pt", "Jet_pt[3]")\
#           .Define("Jet5_pt", "Jet_pt[4]")
#    dfs[dfname] = df
#drawHisto(hists_S2, dfs, "_S2")

# ES3 : bJet
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHisto(hists_S3, dfs, "_S3")

