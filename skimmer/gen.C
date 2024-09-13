#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "../classes/DelphesClasses.h"
#include "../external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TMath.h>

#include "../utils/utility2.h"

// modify //
void gen(std::string channel, std::string outdir="./skimmed/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/stiger97/HLLHC_tthh/DATA/inc/"+channel+"/Events/tag_1_*.root";
    std::cout << infile << std::endl;
    std::cout << outdir << std::endl;
    auto treename = "Delphes";
    auto _df = ROOT::RDataFrame(treename, infile);
    
    // Constants, Basic Branches //
    auto df0 = _df.Define("T_pid", "int(6)")
                  .Define("H_pid", "int(25)")
                  .Define("g_pid", "int(21)")
                  .Define("W_pid", "int(24)")
                  .Define("b_pid", "int(5)")
                  .Define("c_pid", "int(4)")
                  .Define("e_pid", "int(11)")
                  .Define("mu_pid", "int(13)")
                  .Define("int0", "int(0)").Define("int1", "int(1)").Define("int2", "int(2)")
                  .Define("float0", "float(0)")
                  .Define("drmax1", "float(0.15)").Define("drmax2", "float(0.4)")

                  .Define("ParticlePID", {"Particle.PID"})
                  .Define("ParticlePT", {"Particle.PT"})
                  .Define("D1", {"Particle.D1"})
                  .Define("D2", {"Particle.D2"})
                  .Define("M1", {"Particle.M1"})
                  .Define("M2", {"Particle.M2"})
                  .Define("GenMissingET_met", "GenMissingET.MET")
                  .Define("GenMissingET_eta", "GenMissingET.Eta")
                  .Define("GenMissingET_phi", "GenMissingET.Phi");
                  
    // Gen and Matching //
    auto df1 = df0.Define("isLast", ::isLast, {"Particle.PID", "Particle.D1", "Particle.D2"})
                  .Define("Top", "abs(Particle.PID) == 6 && isLast").Define("nTop", "Sum(Top)")
                  .Define("Higgs", "abs(Particle.PID) == 25 && isLast").Define("nHiggs", "Sum(Higgs)")
                  .Define("fin_Higgs_idx", ::findidx, {"ParticlePID", "Higgs"})
                  .Define("Higgs_daughters", ::daughterPID, {"ParticlePID", "D1", "D2", "fin_Higgs_idx"})
                  .Define("W", "abs(Particle.PID) == 24 && isLast").Define("nW", "Sum(W)")
                  .Define("fin_W_idx", ::findidx, {"ParticlePID", "W"})
                  .Define("W_daughters", ::daughterPID, {"ParticlePID", "D1", "D2", "fin_W_idx"})
                  .Define("Tau", "abs(Particle.PID) == 15 && isLast")
                  .Define("fin_Tau_idx", ::findidx, {"Particle.PID", "Tau"})
                  .Define("Tau_daughters", ::daughterPID, {"ParticlePID", "D1", "D2", "fin_Tau_idx"})
                  .Define("Lep", "abs(Particle.PID) == 11 || abs(Particle.PID) == 13 && isLast").Define("nLep", "Sum(Lep)")
                  .Define("fin_Lep_idx", ::findidx, {"ParticlePID", "Lep"})
                  .Define("Lep_mothers", ::motherPID, {"ParticlePID", "M1", "M2", "fin_Lep_idx"})
                  .Define("GenMuon", "abs(Particle.PID) == 13 && isLast").Define("nGenMuon", "Sum(GenMuon)")
                  .Define("GenElectron", "abs(Particle.PID) == 11 && isLast").Define("nGenElectron", "Sum(GenElectron)")
                  .Define("GenbQuark", "abs(Particle.PID) == 5 && isLast").Define("nGenbQ", "Sum(GenbQuark)")
                  .Define("GencQuark", "abs(Particle.PID) == 4 && isLast").Define("nGencQ", "Sum(GencQuark)")
                  .Define("GenbQuark_pt", "Particle.PT[GenbQuark]")
                  .Define("GenbQuark_eta", "Particle.Eta[GenbQuark]")
                  .Define("GenbQuark_phi", "Particle.Phi[GenbQuark]")
                  .Define("GenbQuark_mass", "Particle.Mass[GenbQuark]")
                  .Define("GenMuon_pt", "Particle.PT[GenMuon]")
                  .Define("GenMuon_eta", "Particle.Eta[GenMuon]")
                  .Define("GenElectron_pt", "Particle.PT[GenElectron]")
                  .Define("GenElectron_eta", "Particle.Eta[GenElectron]")
                  .Define("GenHiggs_pt", "Particle.PT[Higgs]")
                  .Define("GenHiggs_eta", "Particle.Eta[Higgs]")
                  .Define("GenHiggs_phi", "Particle.Phi[Higgs]")

                  // Gen bJet
                  .Define("GenJetBTag", ::dRMatching, {"GenbQuark", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenbJet_size", "Sum(GenJetBTag)")
                  .Define("GenbJet_pt", "GenJet.PT[GenJetBTag]")
                  .Define("GenbJet_pt_sorted", ::Sort, {"GenbJet_pt"})
                  .Define("GenbJet_pt_max", "GenbJet_pt_sorted[0]")
                  .Define("GenbJet_eta", "GenJet.Eta[GenJetBTag]")
                  .Define("GenbJet_phi", "GenJet.Phi[GenJetBTag]")
                  .Define("GenbJet_mass", "GenJet.Mass[GenJetBTag]")
                  .Define("GenbJet_dr", ::all_dR, {"GenbJet_pt", "GenbJet_eta", "GenbJet_phi", "GenbJet_mass"})
//                  .Define("GenbJet_dr_min", ::Min_element, {"GenbJet_dr"})
                  .Define("GenbJet_dr_sorted", ::Sort_up, {"GenbJet_dr"})
                  .Define("GenbJet_dr_min", "GenbJet_dr_sorted[0]")
                  .Define("GenbJet_dr_2nd_min", "GenbJet_dr_sorted[1]")
                  .Define("GenbJet_dr_3rd_min", "GenbJet_dr_sorted[2]")
                  .Define("GenbJet_ht", ::Ht, {"GenbJet_pt"})

                  // Find Last Particles
                  .Define("FinalGenPart_idx", ::FinalParticle_idx, {"Particle.PID", "Particle.PT", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Top", "Higgs"})
                  .Define("Top1_idx", "FinalGenPart_idx[0]")
                  .Define("GenbFromTop1_idx", "FinalGenPart_idx[1]")
                  .Define("Top2_idx", "FinalGenPart_idx[2]")
                  .Define("GenbFromTop2_idx", "FinalGenPart_idx[3]")
                  .Define("Higgs1_idx", "FinalGenPart_idx[4]")
                  .Define("Genb1FromHiggs1_idx", "FinalGenPart_idx[5]")
                  .Define("Genb2FromHiggs1_idx", "FinalGenPart_idx[6]")
                  .Define("Higgs2_idx", "FinalGenPart_idx[7]")
                  .Define("Genb1FromHiggs2_idx", "FinalGenPart_idx[8]")
                  .Define("Genb2FromHiggs2_idx", "FinalGenPart_idx[9]")
                  .Define("bQuarkFromTop1_pt", "Particle.PT[GenbFromTop1_idx]")
                  .Define("bQuarkFromTop2_pt", "Particle.PT[GenbFromTop2_idx]")
                  .Define("b1QuarkFromHiggs1_pt", "Particle.PT[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_pt", "Particle.PT[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_pt", "Particle.PT[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_pt", "Particle.PT[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_eta", "Particle.Eta[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_eta", "Particle.Eta[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_eta", "Particle.Eta[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_eta", "Particle.Eta[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_phi", "Particle.Phi[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_phi", "Particle.Phi[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_phi", "Particle.Phi[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_phi", "Particle.Phi[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_mass", "Particle.Mass[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_mass", "Particle.Mass[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_mass", "Particle.Mass[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_mass", "Particle.Mass[Genb2FromHiggs2_idx]")
                  .Define("Q_Higgs1_var", ::GenHiggsReco, {"b1QuarkFromHiggs1_pt", "b1QuarkFromHiggs1_eta", "b1QuarkFromHiggs1_phi", "b1QuarkFromHiggs1_mass", "b2QuarkFromHiggs1_pt", "b2QuarkFromHiggs1_eta", "b2QuarkFromHiggs1_phi", "b2QuarkFromHiggs1_mass"})
                  .Define("Q_Higgs2_var", ::GenHiggsReco, {"b1QuarkFromHiggs2_pt", "b1QuarkFromHiggs2_eta", "b1QuarkFromHiggs2_phi", "b1QuarkFromHiggs2_mass", "b2QuarkFromHiggs2_pt", "b2QuarkFromHiggs2_eta", "b2QuarkFromHiggs2_phi", "b2QuarkFromHiggs2_mass"})
                  .Define("Q_Higgs1_mass", "Q_Higgs1_var[3]")
                  .Define("Q_Higgs2_mass", "Q_Higgs2_var[3]")
                  .Define("Q_Higgs_mass", ::ConcatFloat, {"Q_Higgs1_mass", "Q_Higgs2_mass"})

                  // Higgs
                  .Define("h1_pt", "Particle.PT[Higgs1_idx]")
                  .Define("h2_pt", "Particle.PT[Higgs2_idx]")
                  .Define("h1_eta", "Particle.Eta[Higgs1_idx]")
                  .Define("h2_eta", "Particle.Eta[Higgs2_idx]")
                  .Define("h1_phi", "Particle.Phi[Higgs1_idx]")
                  .Define("h2_phi", "Particle.Phi[Higgs2_idx]")
                  .Define("h1_m", "Particle.Mass[Higgs1_idx]")
                  .Define("h2_m", "Particle.Mass[Higgs2_idx]")
                  .Define("Higgs_Seperation", ::dR2, {"h1_pt", "h1_eta", "h1_phi", "h1_m", "h2_pt", "h2_eta", "h2_phi", "h2_m"})

                  // GenJet
                  .Define("GenJet_pt", "GenJet.PT")
                  .Define("GenJet_eta", "GenJet.Eta")
                  .Define("GenJet_Avg", ::Avg, {"GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenJet_dr_avg", "GenJet_Avg[0]")
                  .Define("GenJet_dEta_avg", "GenJet_Avg[1]")
                  .Define("GenJet_dPhi_avg", "GenJet_Avg[2]")

                  // Gen bJet Matching
                  .Define("GenbJetFromTop1_idx", ::dRMatching_idx, {"GenbFromTop1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenbJetFromTop2_idx", ::dRMatching_idx, {"GenbFromTop2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb1JetFromHiggs1_idx", ::dRMatching_idx, {"Genb1FromHiggs1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb2JetFromHiggs1_idx", ::dRMatching_idx, {"Genb2FromHiggs1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb1JetFromHiggs2_idx", ::dRMatching_idx, {"Genb1FromHiggs2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb2JetFromHiggs2_idx", ::dRMatching_idx, {"Genb2FromHiggs2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("nGenOverlap", ::nOverlap, {"GenbJetFromTop1_idx", "GenbJetFromTop2_idx", "Genb1JetFromHiggs1_idx", "Genb2JetFromHiggs1_idx", "Genb1JetFromHiggs2_idx", "Genb2JetFromHiggs2_idx"})
                  .Define("GenOverlap_bt1", "nGenOverlap[0]")
                  .Define("GenOverlap_bt2", "nGenOverlap[1]")
                  .Define("GenOverlap_b1h1", "nGenOverlap[2]")
                  .Define("GenOverlap_b2h1", "nGenOverlap[3]")
                  .Define("GenOverlap_b1h2", "nGenOverlap[4]")
                  .Define("GenOverlap_b2h2", "nGenOverlap[5]")
                  .Define("GenMuddiness", "nGenOverlap[6]");

    // 4 Vector of Gen // ::idx_var gives -999 for idx = -1.  
    auto df2 = df1.Define("Genb1JetFromHiggs1_pt", ::idx_var, {"GenJet.PT", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_pt", ::idx_var, {"GenJet.PT", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_pt", ::idx_var, {"GenJet.PT", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_pt", ::idx_var, {"GenJet.PT", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_pt", ::idx_var, {"GenJet.PT", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_pt", ::idx_var, {"GenJet.PT", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_eta", ::idx_var, {"GenJet.Eta", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_eta", ::idx_var, {"GenJet.Eta", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_eta", ::idx_var, {"GenJet.Eta", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_eta", ::idx_var, {"GenJet.Eta", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_eta", ::idx_var, {"GenJet.Eta", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_eta", ::idx_var, {"GenJet.Eta", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_phi", ::idx_var, {"GenJet.Phi", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_phi", ::idx_var, {"GenJet.Phi", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_phi", ::idx_var, {"GenJet.Phi", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_phi", ::idx_var, {"GenJet.Phi", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_phi", ::idx_var, {"GenJet.Phi", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_phi", ::idx_var, {"GenJet.Phi", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_mass", ::idx_var, {"GenJet.Mass", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_mass", ::idx_var, {"GenJet.Mass", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_mass", ::idx_var, {"GenJet.Mass", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_mass", ::idx_var, {"GenJet.Mass", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_mass", ::idx_var, {"GenJet.Mass", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_mass", ::idx_var, {"GenJet.Mass", "GenbJetFromTop2_idx"})
                  .Define("GenbJet_pt_scheme", ::pt_scheme, {"Genb1JetFromHiggs1_pt", "Genb2JetFromHiggs1_pt", "Genb1JetFromHiggs2_pt", "Genb2JetFromHiggs2_pt", "GenbJetFromTop1_pt", "GenbJetFromTop2_pt"})

                  // Higgs Reco From GenJets 
                  .Define("GenHiggs1_var", ::GenHiggsReco, {"Genb1JetFromHiggs1_pt", "Genb1JetFromHiggs1_eta", "Genb1JetFromHiggs1_phi", "Genb1JetFromHiggs1_mass", "Genb2JetFromHiggs1_pt", "Genb2JetFromHiggs1_eta", "Genb2JetFromHiggs1_phi", "Genb2JetFromHiggs1_mass"})
                  .Define("GenHiggs2_var", ::GenHiggsReco, {"Genb1JetFromHiggs2_pt", "Genb1JetFromHiggs2_eta", "Genb1JetFromHiggs2_phi", "Genb1JetFromHiggs2_mass", "Genb2JetFromHiggs2_pt", "Genb2JetFromHiggs2_eta", "Genb2JetFromHiggs2_phi", "Genb2JetFromHiggs2_mass"})
                  .Define("GenHiggs1_pt", "GenHiggs1_var[0]")
                  .Define("GenHiggs1_eta", "GenHiggs1_var[1]")
                  .Define("GenHiggs1_phi", "GenHiggs1_var[2]")
                  .Define("GenHiggs1_mass", "GenHiggs1_var[3]")
                  .Define("Genbfh1_dr", "GenHiggs1_var[4]")
                  .Define("Genbfh1_Ht", "GenHiggs1_var[5]")
                  .Define("Genbfh1_dEta", "GenHiggs1_var[6]")
                  .Define("Genbfh1_dPhi", "GenHiggs1_var[7]")
                  .Define("Genbfh1_mbmb", "GenHiggs1_var[8]")
                  .Define("GenHiggs2_pt", "GenHiggs2_var[0]")
                  .Define("GenHiggs2_eta", "GenHiggs2_var[1]")
                  .Define("GenHiggs2_phi", "GenHiggs2_var[2]")
                  .Define("GenHiggs2_mass", "GenHiggs2_var[3]")
                  .Define("Genbfh2_dr", "GenHiggs2_var[4]")
                  .Define("Genbfh2_Ht", "GenHiggs2_var[5]")
                  .Define("Genbfh2_dEta", "GenHiggs2_var[6]")
                  .Define("Genbfh2_dPhi", "GenHiggs2_var[7]")
                  .Define("Genbfh2_mbmb", "GenHiggs2_var[8]")
                  .Define("GenHiggs_mass", ::ConcatFloat, {"GenHiggs1_mass", "GenHiggs2_mass"})
                  .Define("Genbfh_dr", ::ConcatFloat, {"Genbfh1_dr", "Genbfh2_dr"})
                  .Define("Genbfh_Ht", ::ConcatFloat, {"Genbfh1_Ht", "Genbfh2_Ht"})
                  .Define("Genbfh_dEta", ::ConcatFloat, {"Genbfh1_dEta", "Genbfh2_dEta"})
                  .Define("Genbfh_dPhi", ::ConcatFloat, {"Genbfh1_dPhi", "Genbfh2_dPhi"})
                  .Define("Genbfh_mbmb", ::ConcatFloat, {"Genbfh1_mbmb", "Genbfh2_mbmb"});
                 
    auto df3 = df2.Define("goodJet", "Jet.PT>=30 && abs(Jet.Eta)<3.0")
                  .Define("goodElectron", "Electron.PT>=23 && abs(Electron.Eta)<3.0")
                  .Define("goodMuon", "Muon.PT>=17 && abs(Muon.Eta)<2.8")

                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_mass", "Jet.Mass[goodJet]")
                  .Define("Jet_E", ::GetE, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
                  .Define("Jet_btag", "Jet.BTag[goodJet] == 1 || Jet.BTag[goodJet] == 3")
                  .Define("Jet_btag_delphes", "Jet.BTag[goodJet]")
                  .Define("Jet_ctag", "Jet.BTag[goodJet] != 3 && Jet.BTag[goodJet] == 2")

                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_mass", "Jet_mass[Jet_btag]")
                  .Define("bJet_E", ::GetE, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
                  .Define("bJet_size", "Sum(Jet_btag)")
                  .Define("cJet_size", "Sum(Jet_ctag)")
                  .Define("bJet_size_delphes", "Sum(Jet_btag_delphes)")

                  .Define("Muon_pt", "Muon.PT[goodMuon]")
                  .Define("Muon_eta", "Muon.Eta[goodMuon]")
                  .Define("Muon_phi", "Muon.Phi[goodMuon]")
                  .Define("Muon_t", "Muon.T[goodMuon]")
                  .Define("nMuon", "Sum(goodMuon)")
                  .Define("Muon_charge", "Muon.Charge[goodMuon]")
                  .Define("Electron_pt", "Electron.PT[goodElectron]")
                  .Define("Electron_eta", "Electron.Eta[goodElectron]")
                  .Define("Electron_phi", "Electron.Phi[goodElectron]")
                  .Define("Electron_t", "Electron.T[goodElectron]")
                  .Define("nElectron", "Sum(goodElectron)")
                  .Define("Electron_charge", "Electron.Charge[goodElectron]")
                  .Define("Lep_size", "nMuon + nElectron")
                  .Define("Lep_4vec", ::TwoLeptons, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "Muon_charge", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "Electron_charge"})
                  .Define("Lep1_pt", "Lep_4vec[0]")
                  .Define("Lep1_eta", "Lep_4vec[1]")
                  .Define("Lep1_phi", "Lep_4vec[2]")
                  .Define("Lep1_t", "Lep_4vec[3]")
                  .Define("Lep1_ch", "Lep_4vec[4]")
                  .Define("Lep2_pt", "Lep_4vec[5]")
                  .Define("Lep2_eta", "Lep_4vec[6]")
                  .Define("Lep2_phi", "Lep_4vec[7]")
                  .Define("Lep2_t", "Lep_4vec[8]")
                  .Define("Lep2_ch", "Lep_4vec[9]")
                  .Define("SS_OS_DL", "Lep1_ch*Lep2_ch")

                  .Define("MET_E", "MissingET.MET")
                  .Define("MET_Eta", "MissingET.Eta")
                  .Define("MET_Phi", "MissingET.Phi")

                  .Define("j_ht", ::Ht, {"Jet_pt"});
//    auto df3 = df2.Filter("nGenbQ == 4");

    std::initializer_list<std::string> variables = {

// Gen // 
"nGenMuon", "GenMuon_pt", "GenMuon_eta", 
"nGenElectron", "GenElectron_pt", "GenElectron_eta",
"GenJet_pt", "GenJet_eta", 

"Higgs_Seperation", 

"nTop", "nHiggs", "nW","W", 
"fin_W_idx", "W_daughters", "Tau_daughters", 
"Higgs_daughters", "GenHiggs_pt", "GenHiggs_eta", "GenHiggs_phi",
"nLep", "Lep_mothers",
"GenbJet_ht", "GenbJet_pt_max",
"GenbJet_dr", "GenbJet_dr_min", "GenbJet_dr_2nd_min", "GenbJet_dr_3rd_min", "Genbfh_dr", "GenbJet_size", 
"nGenbQ", "nGencQ", "GenbQuark_pt", "GenbQuark_eta",

"Q_Higgs_mass", "GenHiggs_mass",


"Jet_size", "bJet_size", "Lep_size", "Lep1_pt", "Lep2_pt", "j_ht", "bJet_pt", "Jet_pt", "bJet_eta",
"Lep1_ch", "Lep2_ch", "SS_OS_DL", 
"GenMissingET_met", "MET_E",
"Jet_btag", "Jet_btag_delphes", "Jet_ctag", "bJet_size_delphes", "cJet_size"

    };

    // modify // 
    df3.Snapshot(treename, outdir+ "gen_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
