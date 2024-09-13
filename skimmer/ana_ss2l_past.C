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

#include "../utils/utility.h"

// MODIFY!!
void ana_ss2l(std::string channel, std::string outdir="../skimmed/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/stiger97/HLLHC_tthh/DATA/semi/"+channel+"/Events/tag_1_*.root"; // modify //
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
                  .Define("GenJetBTag", {"GenJet.BTag"}).Define("nGenbJet", "Sum(GenJetBTag)")
                  .Define("JetBTag", {"Jet.BTag"}).Define("nbJet", "Sum(JetBTag)")
                  .Define("GenMissingET_met", "GenMissingET.MET")
                  .Define("GenMissingET_eta", "GenMissingET.Eta")
                  .Define("GenMissingET_phi", "GenMissingET.Phi");
                  
    // Gen and Matching //
    auto df1 = df0.Define("isLast", ::isLast, {"Particle.PID", "Particle.D1", "Particle.D2"})
                  .Define("Top", "abs(Particle.PID) == 6 && isLast").Define("nTop", "Sum(Top)")
                  .Define("Higgs", "abs(Particle.PID) == 25 && isLast").Define("nHiggs", "Sum(Higgs)")
                  .Define("W", "abs(Particle.PID) == 24 && isLast").Define("nW", "Sum(W)")
                  .Define("GenbQuark", "abs(Particle.PID) == 5 && isLast").Define("nGenbQ", "Sum(GenbQuark)")
                  .Define("GencQuark", "abs(Particle.PID) == 4 && isLast")
//                  .Filter("nGenbQ == 4")  // CAUTION. THIS OPTION IS ONLY TTBB

                  // Find Last Particles
                  .Define("FinGenPtc_idx", ::FinalParticle_idx, {"Particle.PID", "Particle.PT", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Top", "Higgs"})
                  .Define("Top1_idx", "FinGenPtc_idx[0]")
                  .Define("Gbft1_idx", "FinGenPtc_idx[1]")
                  .Define("Top2_idx", "FinGenPtc_idx[2]")
                  .Define("Gbft2_idx", "FinGenPtc_idx[3]")
                  .Define("Higgs1_idx", "FinGenPtc_idx[4]")
                  .Define("Gb1fh1_idx", "FinGenPtc_idx[5]")
                  .Define("Gb2fh1_idx", "FinGenPtc_idx[6]")
                  .Define("Higgs2_idx", "FinGenPtc_idx[7]")
                  .Define("Gb1fh2_idx", "FinGenPtc_idx[8]")
                  .Define("Gb2fh2_idx", "FinGenPtc_idx[9]")
                  //pt
                  .Define("h1_pt", "Particle.PT[Higgs1_idx]")
                  .Define("h2_pt", "Particle.PT[Higgs2_idx]")
                  .Define("bqft1_pt", "Particle.PT[Gbft1_idx]")
                  .Define("bqft2_pt", "Particle.PT[Gbft2_idx]")
                  .Define("b1qfh1_pt", "Particle.PT[Gb1fh1_idx]")
                  .Define("b2qfh1_pt", "Particle.PT[Gb2fh1_idx]")
                  .Define("b1qfh2_pt", "Particle.PT[Gb1fh2_idx]")
                  .Define("b2qfh2_pt", "Particle.PT[Gb2fh2_idx]")
                  //eta
                  .Define("h1_eta", "Particle.Eta[Higgs1_idx]")
                  .Define("h2_eta", "Particle.Eta[Higgs2_idx]")
                  .Define("b1qfh1_eta", "Particle.Eta[Gb1fh1_idx]")
                  .Define("b2qfh1_eta", "Particle.Eta[Gb2fh1_idx]")
                  .Define("b1qfh2_eta", "Particle.Eta[Gb1fh2_idx]")
                  .Define("b2qfh2_eta", "Particle.Eta[Gb2fh2_idx]")
                  //phi
                  .Define("h1_phi", "Particle.Phi[Higgs1_idx]")
                  .Define("h2_phi", "Particle.Phi[Higgs2_idx]")
                  .Define("b1qfh1_phi", "Particle.Phi[Gb1fh1_idx]")
                  .Define("b2qfh1_phi", "Particle.Phi[Gb2fh1_idx]")
                  .Define("b1qfh2_phi", "Particle.Phi[Gb1fh2_idx]")
                  .Define("b2qfh2_phi", "Particle.Phi[Gb2fh2_idx]")
                  //m
                  .Define("h1_m", "Particle.Mass[Higgs1_idx]")
                  .Define("h2_m", "Particle.Mass[Higgs2_idx]")
                  .Define("b1qfh1_m", "Particle.Mass[Gb1fh1_idx]")
                  .Define("b2qfh1_m", "Particle.Mass[Gb2fh1_idx]")
                  .Define("b1qfh2_m", "Particle.Mass[Gb1fh2_idx]")
                  .Define("b2qfh2_m", "Particle.Mass[Gb2fh2_idx]")

                  .Define("Q_Higgs1_var", ::GenHiggsReco, {"b1qfh1_pt", "b1qfh1_eta", "b1qfh1_phi", "b1qfh1_m", "b2qfh1_pt", "b2qfh1_eta", "b2qfh1_phi", "b2qfh1_m"})
                  .Define("Q_Higgs2_var", ::GenHiggsReco, {"b1qfh2_pt", "b1qfh2_eta", "b1qfh2_phi", "b1qfh2_m", "b2qfh2_pt", "b2qfh2_eta", "b2qfh2_phi", "b2qfh2_m"})
                  .Define("Q_Higgs1_m", "Q_Higgs1_var[3]")
                  .Define("Q_Higgs2_m", "Q_Higgs2_var[3]")
                  .Define("Q_Higgs_m", ::ConcatFloat, {"Q_Higgs1_m", "Q_Higgs2_m"})

                  //dR
                  .Define("Higgs_Seperation", ::dR2, {"h1_pt", "h1_eta", "h1_phi", "h1_m", "h2_pt", "h2_eta", "h2_phi", "h2_m"})

                  // GenJet
                  .Define("GenJet_dR", ::all_dR, {"GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenJet_Avg", ::Avg, {"GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenJet_dr_avg", "GenJet_Avg[0]")
                  .Define("GenJet_dEta_avg", "GenJet_Avg[1]")
                  .Define("GenJet_dPhi_avg", "GenJet_Avg[2]")

                  // Gen bJet Matching // Gbjft1 = Genb1JetFromTop1
                  .Define("Gbjft1_idx", ::dRMatching_idx, {"Gbft1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Gbjft2_idx", ::dRMatching_idx, {"Gbft2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Gbj1fh1_idx", ::dRMatching_idx, {"Gb1fh1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Gbj2fh1_idx", ::dRMatching_idx, {"Gb2fh1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Gbj1fh2_idx", ::dRMatching_idx, {"Gb1fh2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Gbj2fh2_idx", ::dRMatching_idx, {"Gb2fh2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("nGenOverlap", ::nOverlap, {"Gbjft1_idx", "Gbjft2_idx", "Gbj1fh1_idx", "Gbj2fh1_idx", "Gbj1fh2_idx", "Gbj2fh2_idx"})
                  .Define("GenOverlap_bt1", "nGenOverlap[0]")
                  .Define("GenOverlap_bt2", "nGenOverlap[1]")
                  .Define("GenOverlap_b1h1", "nGenOverlap[2]")
                  .Define("GenOverlap_b2h1", "nGenOverlap[3]")
                  .Define("GenOverlap_b1h2", "nGenOverlap[4]")
                  .Define("GenOverlap_b2h2", "nGenOverlap[5]")
                  .Define("GenMuddiness", "nGenOverlap[6]")

                  // Reco bJet Matching
                  .Define("bjft1_idx", ::dRMatching_idx, {"Gbjft1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bjft2_idx", ::dRMatching_idx, {"Gbjft2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bj1fh1_idx", ::dRMatching_idx, {"Gbj1fh1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bj2fh1_idx", ::dRMatching_idx, {"Gbj2fh1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bj1fh2_idx", ::dRMatching_idx, {"Gbj1fh2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bj2fh2_idx", ::dRMatching_idx, {"Gbj2fh2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})

                  .Define("nOverlap", ::nOverlap, {"bjft1_idx", "bjft2_idx", "bj1fh1_idx", "bj2fh1_idx", "bj1fh2_idx", "bj2fh2_idx"})
                  .Define("Overlap_bt1", "nOverlap[0]")
                  .Define("Overlap_bt2", "nOverlap[1]")
                  .Define("Overlap_b1h1", "nOverlap[2]")
                  .Define("Overlap_b2h1", "nOverlap[3]")
                  .Define("Overlap_b1h2", "nOverlap[4]")
                  .Define("Overlap_b2h2", "nOverlap[5]")
                  .Define("Muddiness", "nOverlap[6]")
                  
                  .Define("nMatchedbJet", ::NumberOf, {"bj1fh1_idx", "bj2fh1_idx", "bj1fh2_idx", "bj2fh2_idx", "bjft1_idx", "bjft2_idx"})
                  .Define("nMatchedbJet_FromHiggs1", "nMatchedbJet[0]")
                  .Define("nMatchedbJet_FromHiggs2", "nMatchedbJet[1]")
                  .Define("nMatchedbJet_FromTop1", "nMatchedbJet[2]")
                  .Define("nMatchedbJet_FromTop2", "nMatchedbJet[3]")
                  .Define("nMatchedbJet_all", "nMatchedbJet[4]");

    // 4 Vector of Gen // ::idx_var gives -999 for idx = -1.  
    auto df2 = df1.Define("Gbj1fh1_pt", ::idx_var, {"GenJet.PT", "Gbj1fh1_idx"})
                  .Define("Gbj2fh1_pt", ::idx_var, {"GenJet.PT", "Gbj2fh1_idx"})
                  .Define("Gbj1fh2_pt", ::idx_var, {"GenJet.PT", "Gbj1fh2_idx"})
                  .Define("Gbj2fh2_pt", ::idx_var, {"GenJet.PT", "Gbj2fh2_idx"})
                  .Define("Gbjft1_pt", ::idx_var, {"GenJet.PT", "Gbjft1_idx"})
                  .Define("Gbjft2_pt", ::idx_var, {"GenJet.PT", "Gbjft2_idx"})

                  .Define("Gbj1fh1_eta", ::idx_var, {"GenJet.Eta", "Gbj1fh1_idx"})
                  .Define("Gbj2fh1_eta", ::idx_var, {"GenJet.Eta", "Gbj2fh1_idx"})
                  .Define("Gbj1fh2_eta", ::idx_var, {"GenJet.Eta", "Gbj1fh2_idx"})
                  .Define("Gbj2fh2_eta", ::idx_var, {"GenJet.Eta", "Gbj2fh2_idx"})
                  .Define("Gbjft1_eta", ::idx_var, {"GenJet.Eta", "Gbjft1_idx"})
                  .Define("Gbjft2_eta", ::idx_var, {"GenJet.Eta", "Gbjft2_idx"})

                  .Define("Gbj1fh1_phi", ::idx_var, {"GenJet.Phi", "Gbj1fh1_idx"})
                  .Define("Gbj2fh1_phi", ::idx_var, {"GenJet.Phi", "Gbj2fh1_idx"})
                  .Define("Gbj1fh2_phi", ::idx_var, {"GenJet.Phi", "Gbj1fh2_idx"})
                  .Define("Gbj2fh2_phi", ::idx_var, {"GenJet.Phi", "Gbj2fh2_idx"})
                  .Define("Gbjft1_phi", ::idx_var, {"GenJet.Phi", "Gbjft1_idx"})
                  .Define("Gbjft2_phi", ::idx_var, {"GenJet.Phi", "Gbjft2_idx"})

                  .Define("Gbj1fh1_m", ::idx_var, {"GenJet.Mass", "Gbj1fh1_idx"})
                  .Define("Gbj2fh1_m", ::idx_var, {"GenJet.Mass", "Gbj2fh1_idx"})
                  .Define("Gbj1fh2_m", ::idx_var, {"GenJet.Mass", "Gbj1fh2_idx"})
                  .Define("Gbj2fh2_m", ::idx_var, {"GenJet.Mass", "Gbj2fh2_idx"})
                  .Define("Gbjft1_m", ::idx_var, {"GenJet.Mass", "Gbjft1_idx"})
                  .Define("Gbjft2_m", ::idx_var, {"GenJet.Mass", "Gbjft2_idx"})
                  .Define("GenbJet_pt_scheme", ::pt_scheme, {"Gbj1fh1_pt", "Gbj2fh1_pt", "Gbj1fh2_pt", "Gbj2fh2_pt", "Gbjft1_pt", "Gbjft2_pt"})

                  // Higgs Reco From GenJets 
                  .Define("GenHiggs1_var", ::GenHiggsReco, {"Gbj1fh1_pt", "Gbj1fh1_eta", "Gbj1fh1_phi", "Gbj1fh1_m", "Gbj2fh1_pt", "Gbj2fh1_eta", "Gbj2fh1_phi", "Gbj2fh1_m"})
                  .Define("GenHiggs2_var", ::GenHiggsReco, {"Gbj1fh2_pt", "Gbj1fh2_eta", "Gbj1fh2_phi", "Gbj1fh2_m", "Gbj2fh2_pt", "Gbj2fh2_eta", "Gbj2fh2_phi", "Gbj2fh2_m"})
                  .Define("GenHiggs1_pt", "GenHiggs1_var[0]")
                  .Define("GenHiggs1_eta", "GenHiggs1_var[1]")
                  .Define("GenHiggs1_phi", "GenHiggs1_var[2]")
                  .Define("GenHiggs1_m", "GenHiggs1_var[3]")
                  .Define("GenHiggs1_bdr", "GenHiggs1_var[4]")
                  .Define("GenHiggs1_Ht", "GenHiggs1_var[5]")
                  .Define("GenHiggs1_dEta", "GenHiggs1_var[6]")
                  .Define("GenHiggs1_dPhi", "GenHiggs1_var[7]")
                  .Define("GenHiggs1_mbmb", "GenHiggs1_var[8]")
                  .Define("GenHiggs2_pt", "GenHiggs2_var[0]")
                  .Define("GenHiggs2_eta", "GenHiggs2_var[1]")
                  .Define("GenHiggs2_phi", "GenHiggs2_var[2]")
                  .Define("GenHiggs2_m", "GenHiggs2_var[3]")
                  .Define("GenHiggs2_bdr", "GenHiggs2_var[4]")
                  .Define("GenHiggs2_Ht", "GenHiggs2_var[5]")
                  .Define("GenHiggs2_dEta", "GenHiggs2_var[6]")
                  .Define("GenHiggs2_dPhi", "GenHiggs2_var[7]")
                  .Define("GenHiggs2_mbmb", "GenHiggs2_var[8]")
                  .Define("GenHiggs_m", ::ConcatFloat, {"GenHiggs1_m", "GenHiggs2_m"})
                  .Define("GenHiggs_bdr", ::ConcatFloat, {"GenHiggs1_bdr", "GenHiggs2_bdr"})
                  .Define("GenHiggs_Ht", ::ConcatFloat, {"GenHiggs1_Ht", "GenHiggs2_Ht"})
                  .Define("GenHiggs_dEta", ::ConcatFloat, {"GenHiggs1_dEta", "GenHiggs2_dEta"})
                  .Define("GenHiggs_dPhi", ::ConcatFloat, {"GenHiggs1_dPhi", "GenHiggs2_dPhi"})
                  .Define("GenHiggs_mbmb", ::ConcatFloat, {"GenHiggs1_mbmb", "GenHiggs2_mbmb"});

    // Reco Matching //
    auto df3 = df2.Define("b1JetFromHiggs1_pt", ::idx_var, {"Jet.PT", "bj1fh1_idx"})
                  .Define("b2JetFromHiggs1_pt", ::idx_var, {"Jet.PT", "bj2fh1_idx"})
                  .Define("b1JetFromHiggs2_pt", ::idx_var, {"Jet.PT", "bj1fh2_idx"})
                  .Define("b2JetFromHiggs2_pt", ::idx_var, {"Jet.PT", "bj2fh2_idx"})
                  .Define("bJetFromTop1_pt", ::idx_var, {"Jet.PT", "bjft1_idx"})
                  .Define("bJetFromTop2_pt", ::idx_var, {"Jet.PT", "bjft2_idx"})

                  .Define("Muon_size_del", "Muon_size")
                  .Define("Electron_size_del", "Electron_size")
                  .Define("Jet_size_del", "Jet_size");

    // Reco // Object Selection //
    auto df4 = df3.Define("goodJet", "Jet.PT>=30 && abs(Jet.Eta)<3.0")
                  .Define("goodElectron", "Electron.PT>=23 && abs(Electron.Eta)<3.0")
                  .Define("goodMuon", "Muon.PT>=17 && abs(Muon.Eta)<2.8")

                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_m", "Jet.Mass[goodJet]")
                  .Define("Jet_E", ::GetE, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_m"})
                  .Define("Jet_btag", "Jet.BTag[goodJet]==1 || Jet.BTag[goodJet]==3")
                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_m", "Jet_m[Jet_btag]")
                  .Define("bJet_E", ::GetE, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_m"})
                  .Define("bJet_size", "Sum(Jet_btag)")

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

                  .Define("j_ht", ::Ht, {"Jet_pt"})
                  // You Must Use [bi_Tag] Before Filter // 
                  // Event Selection // 
                  .Filter("Lep_size == 2")
                  .Filter("SS_OS_DL == 1")
                  .Filter("bJet_size >= 3")
                  .Filter("j_ht > 300")

                  .Define("bJet_pt_scheme", ::pt_scheme, {"b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", "bJetFromTop1_pt", "bJetFromTop2_pt"})
                  .Define("isMatchable", ::Matchable, {"bJet_pt", "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", "bJetFromTop1_pt", "bJetFromTop2_pt"})

                  .Define("Jet1_pt", "Jet_pt[0]").Define("Jet1_eta", "Jet_eta[0]").Define("Jet1_phi", "Jet_phi[0]")
                  .Define("Jet1_m", "bJet_m[0]")
                  .Define("Jet2_pt", "Jet_pt[1]").Define("Jet2_eta", "Jet_eta[1]").Define("Jet2_phi", "Jet_phi[1]")
                  .Define("Jet2_m", "Jet_m[1]")
                  .Define("Jet3_pt", "Jet_pt[2]").Define("Jet3_eta", "Jet_eta[2]").Define("Jet3_phi", "Jet_phi[2]")
                  .Define("Jet3_m", "Jet_m[2]")
                  .Define("Jet4_pt", "Jet_pt[3]").Define("Jet4_eta", "Jet_eta[3]").Define("Jet4_phi", "Jet_phi[3]")
                  .Define("Jet4_m", "Jet_m[3]")
                  .Define("Jet5_pt", "Jet_pt[4]").Define("Jet5_eta", "Jet_eta[4]").Define("Jet5_phi", "Jet_phi[4]")
                  .Define("Jet5_m", "Jet_m[4]")
                  .Define("bJet1_pt", "bJet_pt[0]").Define("bJet1_eta", "bJet_eta[0]").Define("bJet1_phi", "bJet_phi[0]")
                  .Define("bJet1_m", "bJet_m[0]")
                  .Define("bJet2_pt", "bJet_pt[1]").Define("bJet2_eta", "bJet_eta[1]").Define("bJet2_phi", "bJet_phi[1]")
                  .Define("bJet2_m", "bJet_m[1]")
                  .Define("bJet3_pt", "bJet_pt[2]").Define("bJet3_eta", "bJet_eta[2]").Define("bJet3_phi", "bJet_phi[2]")
                  .Define("bJet3_m", "bJet_m[2]")
                  .Define("bJet4_pt", "bJet_pt[3]").Define("bJet4_eta", "bJet_eta[3]").Define("bJet4_phi", "bJet_phi[3]")
                  .Define("bJet4_m", "bJet_m[3]")
                  .Define("bJet5_pt", "bJet_pt[4]").Define("bJet5_eta", "bJet_eta[4]").Define("bJet5_phi", "bJet_phi[4]")
                  .Define("bJet5_m", "bJet_m[4]")

                  .Define("j1j2_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_m", "Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_m"})
                  .Define("j1j3_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_m", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_m"})
                  .Define("j1j4_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_m", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_m"})
                  .Define("j1j5_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_m", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_m"})
                  .Define("j2j3_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_m", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_m"})
                  .Define("j2j4_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_m", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_m"})
                  .Define("j2j5_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_m", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_m"})
                  .Define("j3j4_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_m", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_m"})
                  .Define("j3j5_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_m", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_m"})
                  .Define("j4j5_dr", ::dR2, {"Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_m", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_m"})
                  .Define("jj_dr", ::ConcatFloat_withoutSort_10, {"j1j2_dr", "j1j3_dr", "j1j4_dr", "j1j5_dr", "j2j3_dr", "j2j4_dr", "j2j5_dr", "j3j4_dr", "j3j5_dr", "j4j5_dr"}) //Be Careful

                  .Define("j_Vars", ::Vars, {"jj_dr", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_m", "Jet_E"})
                  .Define("jj_avg_dr", "j_Vars[0]")
                  .Define("jj_max_dr", "j_Vars[1]")
                  .Define("jj_min_dr", "j_Vars[2]")
                  .Define("jj_dEta_WhenMaxdR", "j_Vars[3]")
                  .Define("j_cent", "j_Vars[5]")
                  .Define("jj_max_deta", "j_Vars[6]")
                  .Define("jj_max_m", "j_Vars[7]")
                  .Define("jj_twist", "j_Vars[8]")

                  .Define("b1b2_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m", "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m"})
                  .Define("b1b3_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m"})
                  .Define("b1b4_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m"})
                  .Define("b1b5_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m"})
                  .Define("b2b3_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m"})
                  .Define("b2b4_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m"})
                  .Define("b2b5_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m"})
                  .Define("b3b4_dr", ::dR2, {"bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m"})
                  .Define("b3b5_dr", ::dR2, {"bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m"})
                  .Define("b4b5_dr", ::dR2, {"bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m"})
                  .Define("bb_dr", ::ConcatFloat_withoutSort_10, {"b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr", "b2b3_dr", "b2b4_dr", "b2b5_dr", "b3b4_dr", "b3b5_dr", "b4b5_dr"})

                  
                  .Define("b_Vars", ::Vars, {"bb_dr", "bJet_pt", "bJet_eta", "bJet_phi", "bJet_m", "bJet_E"})
                  .Define("bb_avg_dr", "b_Vars[0]")
                  .Define("bb_max_dr", "b_Vars[1]")
                  .Define("bb_min_dr", "b_Vars[2]")
                  .Define("bb_dEta_WhenMaxdR", "b_Vars[3]")
                  .Define("b_ht", "b_Vars[4]")
                  .Define("b_cent", "b_Vars[5]")
                  .Define("bb_min_deta", "b_Vars[6]")
                  .Define("bb_max_m", "b_Vars[7]")
                  .Define("bb_twist", "b_Vars[8]")


                   // bJet Back Tracing //  
                  .Define("bJetFrom", ::bJetFrom, {"bJet_pt", "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", "bJetFromTop1_pt", "bJetFromTop2_pt"})
                  .Define("b1", "bJetFrom[0]")
                  .Define("b2", "bJetFrom[1]")
                  .Define("b3", "bJetFrom[2]")
                  .Define("b4", "bJetFrom[3]")
                  .Define("b5", "bJetFrom[4]")
                  .Define("nMatched_bFromTop", "bJetFrom[5]")
                  .Define("nMatched_bFromHiggs", "bJetFrom[6]")
                  .Define("nMatched_bJet", "bJetFrom[7]")

                  // Answer Categorizations for DNN // 
                  .Define("bCat_higgs5_2Mat", ::bCat_higgs5_2Mat, {"b1", "b2", "b3", "b4", "b5"})
                  .Define("bCat_higgs5_2Mat_multi", ::bCat_higgs5_2Mat_multi, {"b1", "b2", "b3", "b4", "b5"})
                  .Define("bCat_higgs5_2Mat_1", "bCat_higgs5_2Mat_multi[0]")
                  .Define("bCat_higgs5_2Mat_2", "bCat_higgs5_2Mat_multi[1]")
                  .Define("bCat_top_1", ::bCat_top_1, {"b1", "b2", "b3", "b4", "b5"})

                  .Define("Event_shapes", ::Event_shapes, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_m"})
                  .Define("aplanarity", "Event_shapes[0]")
                  .Define("sphericity", "Event_shapes[1]");
                  

    // Reconstruction of Higgs //
    auto df5 = df4.Define("RecoHiggs", ::RecoHiggs, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_m", "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt"})
                  .Define("chi2_Higgs_pt", "RecoHiggs[0]")
                  .Define("chi2_Higgs_eta", "RecoHiggs[1]")
                  .Define("chi2_Higgs_phi", "RecoHiggs[2]")
                  .Define("chi2_Higgs_m", "RecoHiggs[3]")
                  .Define("Matched_idx1", "RecoHiggs[4]")
                  .Define("Matched_idx2", "RecoHiggs[5]")
                  .Define("Correct_Chi", "RecoHiggs[6]")
                  .Define("Chi_min", "RecoHiggs[7]");



//    auto df6 = df5.Filter("isMatchable > 0"); 

    std::initializer_list<std::string> variables = {

    "nGenbQ",

    "Q_Higgs_m", "Q_Higgs1_m", "Q_Higgs2_m",
    "Gbj1fh1_pt", "Gbj2fh1_pt", "Gbj1fh2_pt", "Gbj2fh2_pt", 
    "Gbjft1_pt", "Gbjft2_pt",
    "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", 
    "bJetFromTop1_pt", "bJetFromTop2_pt",
    "GenOverlap_bt1", "GenOverlap_bt2", "GenOverlap_b1h1", "GenOverlap_b2h1", "GenOverlap_b1h2", "GenOverlap_b2h2",
    "GenMuddiness",
    "Overlap_bt1", "Overlap_bt2", "Overlap_b1h1", "Overlap_b2h1", "Overlap_b1h2", "Overlap_b2h2",
    "Muddiness",

    "nMatchedbJet",
                     
    "GenbJet_pt_scheme", "bJetFrom", 
    "b1", "b2", "b3", "b4", "b5", "nMatched_bFromTop", "nMatched_bFromHiggs", "nMatched_bJet",

    "bCat_higgs5_2Mat", "bCat_top_1",
    "bCat_higgs5_2Mat_1", "bCat_higgs5_2Mat_2",

    "GenJet_dR",

     "Higgs_Seperation", "h1_pt", "h1_eta", "h1_phi", "h1_m", "h2_pt", "h2_eta", "h2_phi", "h2_m",
    //----------------Reco---------------------//

     "Jet_pt", "Jet_eta", "Jet_phi", "Jet_m", "Jet_size",
     "bJet_pt", "bJet_eta", "bJet_phi", "bJet_m", "bJet_size",
     "Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_m",
     "Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_m",
     "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_m",
     "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_m",
     "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_m",
     "j1j2_dr", "j1j3_dr", "j1j4_dr", "j1j5_dr", "j2j3_dr", "j2j4_dr", "j2j5_dr", 
     "j3j4_dr", "j3j5_dr", "j4j5_dr",

     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m",
     "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m",
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr", "b2b3_dr", "b2b4_dr", "b2b5_dr", "b3b4_dr", "b3b5_dr", "b4b5_dr",
     
     "Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "nMuon", "Lep_size", "Muon_charge",
     "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "nElectron", "Electron_charge",
     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t", "Lep1_ch",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t", "Lep2_ch",
     "SS_OS_DL",

     "MET_E", "MET_Eta", "MET_Phi",


    //----------------New-----------------------//

     "bb_dr", "b_Vars", "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_min_deta", "bb_max_m", "bb_twist",
    "jj_dr", "j_Vars", "jj_avg_dr", "jj_max_dr", "jj_min_dr", "j_ht", "jj_dEta_WhenMaxdR", "j_cent", "jj_max_deta", "jj_max_m", "jj_twist", 
    "chi2_Higgs_pt", "chi2_Higgs_eta", "chi2_Higgs_phi", "chi2_Higgs_m",

     "isMatchable", "Matched_idx1", "Matched_idx2", "Correct_Chi", "Chi_min", 
     "Event_shapes", "aplanarity", "sphericity",

     "Muon_size_del", "Electron_size_del", "Jet_size_del",

     "bCat_higgs5_2Mat_multi"
    };

    // MODIFY!!
    df5.Snapshot(treename, outdir+ "ss2l_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
