# conda activate py36
import ROOT
import numpy as np
import os

# Criteria
PRE = "CF_1223_14TeV"
Tree = "Delphes"
fxmin, fxmax = 50, 200
fit_function = ROOT.TF1("fit_function", "gaus", fxmin, fxmax)
#fit_function = ROOT.TF1("fit_function", "[0] / (2 * TMath::Pi()) * ([1] / ( (x - [2])**2 + ([1]**2) / 4 ))", fxmin, fxmax)
#fit_function.SetParameters(1, 1, 125)

# Input
indir = "./samples2/" 
tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh_di.root")
#ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb_di.root")
#ttbbcc = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbcc_di.root")
#ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb_di.root")
dfs = {"_tthh": tthh} #"_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb }

# Output
outdir = "./plots/"+ PRE + "/Fitting/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHistoSingle(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}

        for df_name, df in dfs.items():
            nbin, xmin, xmax = 80, 0, 400
            _xmax = df.Max(hist_name).GetValue()
            if _xmax < 100: xmax = 100
            if _xmax < 20: nbin, xmax = 60, 6 #int(_xmax+2), int(_xmax+2)
            if _xmax < 0: nbin, xmin, xmax = 20, -4, 4
            hist_title = hist_name.replace("_size", " multiplicity")
            hist_title = hist_title.replace("_", " ")
            hist_title = hist_title + df_name
            hist_title = "Gen Higgs mass"
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            h.GetXaxis().SetTitle('Higgs Mass m_{bb} [GeV]')
            h.GetYaxis().SetTitle("Entries")
            h.GetYaxis().SetTitleOffset(1.5)                           
            h.SetLineWidth(2)
            hist_dict[hist_name + df_name] = h

        for _tmp, h in hist_dict.items():
            fit_result = h.Fit(fit_function, "S")
            fit_status = fit_result.Status()
            fit_parameters = fit_function.GetParameters()
            chi_square = fit_result.Chi2()
            ndf = fit_result.Ndf()

            h.Draw("hist")
#            line = ROOT.TLine(0.4, 0, 0.4, 375)
#            line.SetLineColor(ROOT.kRed)
#            line.SetLineWidth(2)
#            line.Draw("same")

            fit_function.Draw("same")
            stats = ROOT.TPaveStats(0.6, 0.7, 0.9, 0.5)

            latex = ROOT.TLatex()
            latex.SetTextSize(0.025)  # Adjust the text size as needed
            latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")  # Modify the text and position

#            pave = ROOT.TPaveText(0.65, 0.45, 0.85, 0.6, "NDC")
#            pave.SetFillColor(0)
#            pave.SetTextSize(0.04)
#            pave.AddText("Fit Results:")
#            pave.AddText("Status: {}".format(fit_status))
#            pave.AddText("Mean: {:.2f}".format(fit_function.GetParameter(1)))
#            pave.AddText("Sigma: {:.2f}".format(fit_function.GetParameter(2)))
#            pave.AddText("Chi^2/NDF: {:.2f} / {}".format(chi_square, ndf))
#            pave.Draw("same")
            canvas.Print(outdir + PRE + "_" + _tmp + flag + ".pdf")
            canvas.Clear()

# Histogram Config #
hists_S0 = [
     "GenHiggs_mass"

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
#for dfname, df in dfs.items():
#    df = df.Filter("GenHiggs1_mass>0 & GenHiggs2_mass>0")
#    dfs[dfname] = df
#drawHistoSingle(hists_S1, dfs, "_S1")

# ES2 : Jet
#for dfname, df in dfs.items():
#    df = df.Filter("Jet_size >= 5")\
#           .Define("Jet1_pt", "Jet_pt[0]")\
#           .Define("Jet2_pt", "Jet_pt[1]")\
#           .Define("Jet3_pt", "Jet_pt[2]")\
#           .Define("Jet4_pt", "Jet_pt[3]")\
#           .Define("Jet5_pt", "Jet_pt[4]")
#    dfs[dfname] = df
#drawHistoSingle(hists_S2, dfs, "_S2")

# ES3 : bJet
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHistoSingle(hists_S3, dfs, "_S3")
