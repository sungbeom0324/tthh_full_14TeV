# Draw histograms with modules. conda activate py36 
from utils.drawHistoModules import *

# in a file.
infile = "/home/stiger97/github/tthh_full_14TeV/skimmed/B_0907_tthh.root"
tree = "Delphes" # modify #

# in multiple files.
indir = "./skimmed/"
PRE = "os2l"

# outdir = "./plots/" + PRE + "/XX/"

# examples #
#drawHisto2D(infile, tree, "nbJet vs MET QCD", "nbJet", "MET [GeV]", "bJet_size", 10,0,10, "MET_E", 100,0,200, "CF_0308")
#drawHistoSame(indir, tree, "GenbQ_pT", "pT [GeV]", "Normalized Entries", "GenbQuark_pt", 100, 0, 300, "Gen_0220", lepo=1)
#----------# 

# Single #
#drawHistoSingle(infile, tree, "Category", "Category (bJet pair from Higgs)", "Entries", "bCat_higgs5_2Mat", 11, 0, 11, PRE)
#drawHistoSingle(infile, tree, "Categories", "Category (bJet pair from Higgs)", "Entries", "bCat_higgs5_2Mat_multi", 11, 0, 11, PRE)
#drawHistoSingle(infile, tree, "bJets from Higgs", "Number of bJet from Higgs", "Entries", "nMatched_bFromHiggs", 5, 0, 5, PRE)

# GEN #
# mu
#drawHistoSame(indir, tree, "Gen Muon multiplicity", "N_{\mu}", "Normalized Entries", "nGenMuon", 10, 0, 10, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen Muon p_{T}", "p_{T, \mu}", "Normalized Entries", "GenMuon_pt", 15, 0, 150, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen Muon \eta", "\eta_{\mu}", "Normalized Entries", "GenMuon_eta", 16, -4, 4, PRE, lepo=1)
# electron
#drawHistoSame(indir, tree, "Gen Electron multiplicity", "N_{e}", "Normalized Entries", "nGenElectron", 5, 0, 5, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen Electron p_{T}", "p_{T, e}", "Normalized Entries", "GenElectron_pt", 10, 0, 100, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen Electron \eta", "\eta_{e}", "Normalized Entries", "GenElectron_eta", 16, -4, 4, PRE, lepo=1)
# bJet
#drawHistoSame(indir, tree, "Gen bJet multiplicity", "N_{bJet}", "Normalized Entries", "GenbJet_size", 9, 0, 9, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen bJet H_{T}", "H_{T, bJet}", "Normalized Entries", "GenbJet_ht", 20, 0, 1000, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen bJet \eta", "\eta_{b}", "Normalized Entries", "GenbQuark_eta", 16, -4, 4, PRE, lepo=1)
#drawHistoSingle(infile, tree, "Matchability", "Category of b-tagged jet pair from (leading) Higgs", "Normalized Entries", "bCat_higgs5_2Mat", 11, 0, 11, PRE, lepo=3)
# Higgs
#drawHistoSingle(indir, tree, "Gen Higgs p_{T}", "p_{T, higgs}", "Normalized Entries", "GenHiggs_pt", 20, 0, 200, PRE, lepo=1)
#drawHistoSingle(indir, tree, "Gen Higgs \eta", "\eta_{higgs}", "Normalized Entries", "GenHiggs_eta", 16, -4, 4, PRE, lepo=1)
#drawHistoSingle(indir, tree, "Gen Higgs dR", "\DeltaR_{hh}", "Normalized Entries", "Higgs_Seperation", 12, 0, 6, PRE, lepo=3)
#drawHisto2D(infile, tree, "pt vs eta Higgs", "p_{T, h}", "\eta [GeV]", "GenHiggs_pt", 20,0,200, "GenHiggs_eta", 12,-3,3, "GEN_INC")
#drawHisto2D(infile, tree, "eta vs phi Higgs", "\eta [GeV]", "\phi", "GenHiggs_eta", 12,-3,3, "GenHiggs_phi", 12,-3,3, "GEN_INC")
# b from Higgs
#drawHistoSingle_2Branches(infile, tree, "dR compare bfh vs bb", "\DeltaR_{bb}", "Normalized Entries", "GenbJet_dr", "Genbfh_dr", 12, 0, 6, PRE, lepo=3)
#drawHistoSame(indir, tree, "Gen bJet \DeltaR", "\DeltaR_{bb}", "Normalized Entries", "GenbJet_dr", 24, 0, 6, PRE, lepo=1)
#drawHistoSingle(infile, tree, "Gen bJet form Higgs \DeltaR", "\DeltaR_{bfh}", "Normalized Entries", "Genbfh_dr", 24, 0, 6, PRE, lepo=3)
#drawHistoSingle(infile, tree, "Higgs daughters", "PID", "Normalized Entries", "Higgs_daughters", 60,-30,30, PRE, lepo=3)


# Subplot Multiplicity.
#drawHistoSame_Sub(indir, tree, "Gen bQ multiplicity sub", "N_{b}", "Normalized Entries", "nGenbQ", 10, 0, 10, PRE, lepo=1)


#drawHistoSame(indir, tree, "Gen bQuark multiplicity", "Number of Gen bQuarks", "Normalized Entries", "nGenbQ", 10, 0, 10, PRE, lepo=1)
#drawHistoSame(indir, tree, "bJet p^{T}", "bJet p^{T} [GeV]", "Normalized Entries", "bJet_pt", 50, 0, 500, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen b-jet multiplicity", "Number of Gen b-jets", "Normalized Entries", "GenbJet_size", 10, 0, 10, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen b-Jet H_{T}", "Gen b-Jet H_{T} [GeV]", "Normalized Entries", "GenbJet_ht", 15, 0, 1500, PRE, lepo=1)
#drawHistoSame(indir, tree, "Max Gen b-Jet p_{T}", "p_{T} [GeV]", "Normalized Entries", "GenbJet_pt_max", 10, 0, 500, PRE, lepo=1)
#drawHistoSame(indir, tree, "\DeltaR between b-tagged Jets", "Gen \DeltaR_{bb}", "Normalized Entries", "GenbJet_dr", 14, 0, 7, PRE, lepo=1)
#drawHistoSame(indir, tree, "\DeltaR", "Gen bJet from identical Higgs \DeltaR_{bb}", "Normalized Entries", "Genbfh_dr", 40, 0.1, 4.1, PRE, lepo=1)
#drawHistoSame(indir, tree, "Min \DeltaR between b-tagged Jets", "Gen \DeltaR^{min}_{bb}", "Normalized Entries", "GenbJet_dr_min", 10, 0, 5, PRE, lepo=1)
#drawHistoSame(indir, tree, "2nd Min \DeltaR between b-tagged Jets", "Gen \DeltaR^{2nd min}_{bb}", "Normalized Entries", "GenbJet_dr_2nd_min", 10, 0, 5, PRE, lepo=1)
#drawHistoSame(indir, tree, "3rd Min \DeltaR between b-tagged Jets", "Gen \DeltaR^{3rd min}_{bb}", "Normalized Entries", "GenbJet_dr_3rd_min", 10, 0, 5, PRE, lepo=1)
#drawHistoSame(indir, tree, "Gen Higgs mass (reco from Gen b-Jet pair)", "Gen m_{bb from Higgs} [GeV]", "Normalized Entries", "GenHiggs_mass", 30, 0, 300, PRE, lepo=1)


# CF, Cutflow Distributions #
#drawHistoSame_Sub(indir, tree, "Lepton multiplicity", "N_{lep}", "Normalized Entries", "Lep_size", 5, 0, 5, PRE, lepo=1)
#drawHistoSame_Sub(indir, tree, "b-jet multiplicity", "Number of b-tagged jets", "Normalized Entries", "bJet_size", 10, 0, 10, PRE, lepo=1)
#drawHistoSame_Sub(indir, tree, "Scalar Sum of Jet p_{T}", "H_{T} [GeV]", "Normalized Entries", "j_ht", 20, 0, 2000, PRE, lepo=1)
#drawHistoSame(indir, tree, "Lepton multiplicity", "N_{lep}", "Normalized Entries", "Lep_size", 5, 0, 5, PRE, lepo=1)
#drawHistoSame(indir, tree, "SS or OS", "l1 x l2 charge", "Normalized Entries", "SS_OS_DL", 5, -2, 3, PRE, lepo=1)
#drawHistoSame(indir, tree, "b-jet multiplicity", "Number of b-tagged jets", "Normalized Entries", "bJet_size", 10, 0, 10, PRE, lepo=1)
#drawHistoSame(indir, tree, "Scalar Sum of Jet p_{T}", "H_{T} [GeV]", "Normalized Entries", "j_ht", 16, 0, 1600, PRE, lepo=1)


# After Selecton, Kinematic Distribution after Event Selection #
# XS1
#drawHistoSame_Sub(indir, tree, "Leading lepton p_{T}", "p_{T} of Leading lepton [GeV]", "Normalized Entries", "Lep1_pt", 20, 0, 400, PRE, "OS1", lepo=1)
#drawHistoSame_Sub(indir, tree, "SubLeading lepton p_{T}", "p_{T} of Sub-leading lepton [GeV]", "Normalized Entries", "Lep2_pt", 10, 0, 200,  PRE, "OS1", lepo=1)
#drawHistoSame_Sub(indir, tree, "Leading lepton eta", "\eta of Leading lepton", "Normalized Entries", "Lep1_eta", 16, -4, 4, PRE, "OS1", lepo=1)
#drawHistoSame_Sub(indir, tree, "SubLeading lepton eta", "\eta of Sub-leading lepton", "Normalized Entries", "Lep2_eta", 16, -4, 4, PRE, "OS1", lepo=1)
#drawHistoSame_Sub(indir, tree, "Min Chi Higgs mass", "Minimum \chi^2 m_{h1}", "Normalized Entries", "chi_Higgs_m", 25, 0, 250, PRE, "OS3", lepo=1)
drawHistoSame(indir, tree, "Min Chi Higgs1 mass", "Minimum \chi^2 m_{h1}", "Normalized Entries", "chi_Higgs1_m", 25, 0, 250, PRE, lepo=1)
drawHistoSame(indir, tree, "Min Chi Higgs2 mass", "Minimum \chi^2 m_{h2}", "Normalized Entries", "chi_Higgs2_m", 25, 0, 250, PRE, lepo=1)
drawHistoSame(indir, tree, "Final Chi", "\chi_{HH}", "Normalized Entries", "Final_Chi2", 10, 0, 100, PRE, lepo=1)

# XS2
drawHistoSame_Sub(indir, tree, "Leading bJet p_{T}", "p_{T} of Leading b-tagged jet [GeV]", "Normalized Entries", "bJet1_pt", 30, 0, 600, PRE, "SS3", lepo=1)
drawHistoSame_Sub(indir, tree, "Leading bJet eta", "\eta of Leading b-tagged jet", "Normalized Entries", "bJet1_eta", 16, -4, 4, PRE, "SS3", lepo=1)

# XS3
#drawHistoSame_Sub(indir, tree, "Max bJet dR", "Maximum \DeltaR of b-tagged jets", "Normalized Entries", "bb_max_dr", 12, 0, 6, PRE, "OS3", lepo=1)
#drawHistoSame_Sub(indir, tree, "Min bJet dR", "Minimum \DeltaR of b-tagged jets", "Normalized Entries", "bb_min_dr", 12, 0, 3, PRE, "OS3", lepo=1)



# DNN_kinematic.root, Newly Defined ones #
#drawHistoSame_Single_Sub(infile, tree, "di-Higgs mass (DNN)", "m_{HH} [GeV]", "Normalized Entries", "higgs_mass_sum", 30, 0, 1500, "FULL_0522", lepo=1)
#drawHistoSame_Single_Sub(infile, tree, "Higgs mass (DNN)", "m_{H1} [GeV]", "Normalized Entries", "higgs_mass", 25, 0, 250, "DI_OSDL", lepo=1)
#drawHistoSame_Single_Sub(infile, tree, "2^{nd} Higgs mass (DNN)", "m_{H2} [GeV]", "Normalized Entries", "higgs_mass_sub", 25, 0, 500, "DI_OSDL", lepo=1)
#drawHistoSame_Single_Sub(infile, tree, "\DeltaR of b-Jets from Higgs1", "\DeltaR^{bfh}", "Normalized Entries", "bfh_dr", 35, 0, 7, "DI_OSDL", lepo=1)
#drawHistoSame_Single_Sub(infile, tree, "H_{T} of Jets from Higgs1", "H^{from Higgs}_{T} [GeV]", "Normalized Entries", "bfh_Ht", 28, 0, 700, "DI_OSDL", lepo=1)

# DNN2_result.root 
#drawHistoStack_Single(infile, tree, "DNN Discriminator", "DNN", "Normalized Entries", "DNN", 20, 0, 1, "DI_OSDL", lepo=1)
#drawHistoStack_Single(infile, tree, "tthh", "DNN G1 prediction score", "Number of events", "G1", 20, 0, 1, "DI_OSDL", lepo=1)
#drawHistoStack_Single(infile, tree, "tthbb", "DNN G2 prediction score", "Number of events", "G2", 20, 0, 1, "DI_OSDL", lepo=1)
#drawHistoStack_Single(infile, tree, "ttbb", "DNN G3 prediction score", "Number of events", "G3", 20, 0, 1, "DI_OSDL", lepo=1)
#drawHistoStack_Single(infile, tree, "ttbbbb", "DNN G4 prediction score", "Number of events", "G4", 20, 0, 1, "DI_OSDL", lepo=1)
