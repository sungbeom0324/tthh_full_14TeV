import matplotlib.pyplot as plt
import numpy as np
import ROOT

df = ROOT.RDataFrame("Delphes", "tag_1_delphes_events.root")

# Define new variables
df = df.Define("bJet_size", "Sum(Jet.BTag)")

# Event Selections
df = df.Filter("bJet_size>=2")
df_lep1 = df.Filter("Electron_size + Muon_size == 1")
df_lep2 = df.Filter("Electron_size + Muon_size >= 2")

h_njet = df.Histo1D(ROOT.RDF.TH1DModel("Jet_size", "Jet_size", 15, 0, 15), "Jet_size")
h_nbjet = df.Histo1D(ROOT.RDF.TH1DModel("bJet_size", "bJet_size", 10, 0, 10), "bJet_size")

h_lep1_nele = df_lep1.Histo1D(ROOT.RDF.TH1DModel("Electron_size", "lep1_Electron_size", 5, 0, 5), "Electron_size")
h_lep1_nmu = df_lep1.Histo1D(ROOT.RDF.TH1DModel("Muon_size", "lep1_Muon_size", 5, 0, 5), "Muon_size")
h_lep2_nele = df_lep2.Histo1D(ROOT.RDF.TH1DModel("Electron_size", "lep2_Electron_size", 5, 0, 5), "Electron_size")
h_lep2_nmu = df_lep2.Histo1D(ROOT.RDF.TH1DModel("Muon_size", "lep2_Muon_size", 5, 0, 5), "Muon_size")


canvas = ROOT.TCanvas("c", "canvas", 800, 600)
canvas.Print("output.pdf[")
canvas.Clear()

h_njet.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()
h_nbjet.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()

h_lep1_nele.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()
h_lep2_nele.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()
h_lep1_nmu.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()
h_lep2_nmu.Draw("hist")
canvas.Print("output.pdf")
canvas.Clear()


canvas.Print("output.pdf]")
