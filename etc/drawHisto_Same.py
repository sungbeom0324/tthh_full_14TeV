import ROOT
import numpy as np
import os
ROOT.gStyle.SetOptStat(0)

# RDF # modify # 
PRE = "Gen_0220"
Tree = "Delphes"

# In #
indir = "./samples1/" 
tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh.root")
tthbb  = ROOT.RDataFrame(Tree, indir + PRE + "_tthbb.root")
ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb.root")
ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb.root")
#dfs = {"_tthh": tthh, "_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb, "_tthbb": tthbb}
dfs = {"_tthh": tthh, "_tthbb": tthbb}

# Out #
outdir = "./plots/" + PRE + "/Same/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHisto(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}
        legend = ROOT.TLegend(0.6, 0.855, 0.85, 0.705) # Upper Right #
#        legend = ROOT.TLegend(0.38, 0.855, 0.62, 0.705) # Center #
#        legend = ROOT.TLegend(0.26, 0.855, 0.50, 0.705) # Upper Left #
        hist_title = hist_name.replace("_size", " multiplicity")
        hist_title = hist_title.replace("_", " ")

        ymax, color = 0, 1
        for df_name, df in dfs.items():
            nbin, xmin, xmax = 450, 0, 450 # modify #
            _xmax = df.Max(hist_name).GetValue()
            if _xmax < 100: xmax = 100
            if _xmax < 20: nbin, xmin, xmax = 14, -7, 7 #int(_xmax+2), int(_xmax+2)
            if _xmax < 5: nbin, xmin, xmax = 4, 0, 4
            if _xmax <= 1: nbin, xmin, xmax = 50, 0, 0.006
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            if ymax < h.GetMaximum(): ymax = h.GetMaximum()
            h.GetXaxis().SetTitle("Higgs mass [GeV]") # modify #
            h.GetYaxis().SetTitle("Normalized Entries") # modify #
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
            h.SetMaximum(ymax * 1.2)
            if first:
                h.SetFillColorAlpha(ROOT.kGray, 0.8)
                h.DrawNormalized("hist")
                first = False
            else: h.DrawNormalized("same")
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetTextSize(0.025)  # Adjust the text size as needed
        latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")  # Modify the text and position
        latex.SetTextSize(0.035)  # Adjust the text size as needed
#        latex.DrawLatexNDC(0.5, 0.6, "nGenbQ == 6")  # Modify the text and position
        

        canvas.Print(outdir + PRE + "_" + hist_name + flag + "_higgs.pdf") # modify #
        canvas.Clear()


# Histogram Config #
hists_S0 = [
#    "nTop",
#    "nW", 
#    "nGenbQ",
#    "GenbQuark_pt",
#    "GenbQuark_eta",
#    "nHiggs",
#    "W_daughters",
#    "Tau_daughters",
#    "Lep_size",
#    "Higgs_daughters", 
#    "nHiggs",
#    "Q_Higgs_mass",    
    "GenHiggs_mass"
]

hists_S1 = hists_S0 + [
]

hists_S2 = hists_S1 + [
]

hists_S3 = hists_S2 + [
]


# Draw #
drawHisto(hists_S0, dfs, "_S0")

############## Draw After Event Selection  ##################

# ES1
#for dfname, df in dfs.items():
#    df = df.Filter("nGenbQ == 6")\
#           .Filter("Lep_size == 2 && Lep1_pt > 18 && Lep2_pt > 10")\
#           .Filter("j_ht > 300")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")\
#           .Define("bJet4_pt", "bJet_pt[3]")\
#           .Define("bJet5_pt", "bJet_pt[4]")
#    dfs[dfname] = df
#drawHisto(hists_S1, dfs, "_S1")

# ES2
#for dfname, df in dfs.items():
#    df = df.Filter("Jet_size >= 5")\
#           .Define("Jet1_pt", "Jet_pt[0]")\
#           .Define("Jet2_pt", "Jet_pt[1]")\
#           .Define("Jet3_pt", "Jet_pt[2]")\
#           .Define("Jet4_pt", "Jet_pt[3]")\
#           .Define("Jet5_pt", "Jet_pt[4]")
#    dfs[dfname] = df
#drawHisto(hists_S2, dfs, "_S2")

# ES3
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHisto(hists_S3, dfs, "_S3")

