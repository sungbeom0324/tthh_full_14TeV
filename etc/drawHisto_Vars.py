import ROOT
import numpy as np
import os
ROOT.gStyle.SetOptStat(0)

# Prefix and Tree
PRE = "GEN_1225_14TeV"
Tree = "Delphes"

# Input
indir = "./samples1/" 
tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh_di.root")
#ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb_di.root")
#ttbbcc = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbcc_di.root")
#ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb_di.root")
dfs = {"_tthh": tthh} #,"_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb }

# Output
outdir = "./plots/" + PRE + "/Vars/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHisto(branches, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for df_name, df in dfs.items():
        legend = ROOT.TLegend(0.68, 0.87, 0.85, 0.8)
        hist_title = "#Delta R_{bb from Higgs}" # MODIFY #

        ymax, color = 0, 1
        hist_dict = {}
        for branch in branches:
            nbin, xmin, xmax = 40, 0, 4 # MODIFY #
            h = df.Histo1D(ROOT.RDF.TH1DModel(branch, hist_title, nbin, xmin, xmax), branch)
            if ymax < h.GetMaximum(): ymax = h.GetMaximum() 
            h.GetXaxis().SetTitle("#Delta R_{bb}")               
            h.GetYaxis().SetTitle("Normalized entries")     
            h.GetYaxis().SetTitleOffset(1.8)                
            h.SetLineColor(color)                           
            color+=3                                        
            if color in [5] :
                color += 1
            h.SetLineWidth(2)                               
            legend_names = ["#Delta R_{bfh}", "#Delta R_{bb}"] # MODIFY #
            if branch == "Genbfh_dr": legend_name = legend_names[0]
            if branch == "GenbJet_dr": legend_name = legend_names[1]
            #legend_name = df_name + " " + branch.replace("_", " ")
            legend.AddEntry(h.GetValue(), legend_name, "l")
            hist_dict[branch + "_" + df_name] = h 
        first = True
        for _tmp, h in hist_dict.items():
            h.SetMaximum(ymax * 0.2)
            if first:
                h.DrawNormalized("E hist")
                first = False
            else: h.DrawNormalized("E same")
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetTextSize(0.025)  # Adjust the text size as needed
        latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")  # Modify the text and position

        canvas.Print(outdir + PRE + "_" + "dr_Comparison" + flag + ".pdf") # MODIFY # 
        canvas.Clear()

 
# Branch
hists_S0 = [
   "Genbfh_dr", "GenbJet_dr"
]
hists_S1 = hists_S0 + [
]

drawHisto(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton
'''
for dfname, df in dfs.items():
    df = df.Filter("Genb1JetFromHiggs1_pt>=0")\
           .Filter("Genb2JetFromHiggs1_pt>=0")\
           .Filter("Genb1JetFromHiggs2_pt>=0")\
           .Filter("Genb2JetFromHiggs2_pt>=0")\
           .Define("Asymmetry_bfh1", "Genb1JetFromHiggs1_pt-Genb2JetFromHiggs1_pt")\
           .Define("Asymmetry_bfh2", "Genb1JetFromHiggs2_pt-Genb2JetFromHiggs2_pt")
dfs[dfname] = df
drawHisto(hists_S1, dfs, "_S1")
'''
