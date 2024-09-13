# conda activate py36
import ROOT
import numpy as np
import os

# Criteria
PRE = "GEN_0115_14TeV"
Tree = "Delphes"

# Input
indir = "./samples1/" 
#tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh_di.root")
ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb_di.root")
ttbbcc = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbcc_di.root")
ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb_di.root")
dfs = {"_ttbb": ttbb} #"_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb }

# Output
outdir = "./plots/"+ PRE + "/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHistoSingle(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}

        for df_name, df in dfs.items():
            nbin, xmin, xmax = 40, 0, 400
            _xmax = df.Max(hist_name).GetValue()
            if _xmax < 100: xmax = 100
            if _xmax < 20: nbin, xmin, xmax = 10, 0, 10 #int(_xmax+2), int(_xmax+2)
            if _xmax < 1: nbin, xmax = 10 ,1
            if _xmax < 0: nbin, xmin, xmax = 20, -4, 4
            hist_title = hist_name.replace("_size", " multiplicity")
            hist_title = hist_title.replace("_", " ")
            hist_title = hist_title + df_name
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            h.GetXaxis().SetTitle(hist_title)
            h.GetYaxis().SetTitle("Normalized Entries")
            h.GetYaxis().SetTitleOffset(1.5)                           
            h.SetLineWidth(2)
            hist_dict[hist_name + df_name] = h

            latex = ROOT.TLatex()
            latex.SetTextSize(0.035)  # Adjust the text size as needed
            latex.DrawLatexNDC(0.5, 0.6, "nGenbQ == 6")  # Modify the text and position

        for _tmp, h in hist_dict.items():
            h.DrawNormalized("hist")
            canvas.Print(outdir + PRE + "_" + _tmp + flag + ".pdf")
            canvas.Clear()

# Histogram Config #
hists_S0 = [
 "nGenbQ", "bJet_size"
]

hists_S1 = hists_S0 + [
]

hists_S2 = hists_S1 + [
]

hists_S3 = hists_S2 + [
]


# Draw #
drawHistoSingle(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton 
for dfname, df in dfs.items():
    df = df.Filter("Lep_size == 2 && Lep1_pt > 17 && Lep2_pt > 10")\
           .Filter("bJet_size >= 5 ")\
           .Filter("j_ht > 300")
    dfs[dfname] = df
drawHistoSingle(hists_S1, dfs, "_S1")

# ES2 : Jet
for dfname, df in dfs.items():
    df = df.Filter("nGenbQ==4")
    dfs[dfname] = df
drawHistoSingle(hists_S2, dfs, "_S2")

# ES3 : bJet
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHistoSingle(hists_S3, dfs, "_S3")
