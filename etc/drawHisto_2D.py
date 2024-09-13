import ROOT
import numpy as np

df = ROOT.RDataFrame("Delphes", "./samples1/Gen_0220_ttbb.root")
canvas = ROOT.TCanvas("c1", "c1")

#hist = df.Histo2D(ROOT.RDF.TH2DModel("bJet pt vs eta", "bJet pt vs eta", 40, 0, 400, 30, -3, 3), "bJet_pt", "bJet_eta")
#hist.GetXaxis().SetTitle("p_{T} [GeV]")
#hist.GetYaxis().SetTitle("MET")

hist = df.Histo2D(ROOT.RDF.TH2DModel("Jet Ht vs MET", "bJet pt vs eta", 40, 0, 400, 100, 0, 100), "j_ht", "MET_E")
hist.GetXaxis().SetTitle("p_{T} [GeV]")
hist.GetYaxis().SetTitle("MET")
hist.Draw("colz")
canvas.Print("hist.pdf")
