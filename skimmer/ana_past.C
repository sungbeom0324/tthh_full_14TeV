
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
void ana(std::string channel, std::string outdir="../skimmed/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/stiger97/HLLHC_tthh/DATA/inc/"+channel+"/Events/tag_*.root"; // modify //
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
                  .Define("isAdd", ::isAdd, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2"})
                  .Define("Top", "abs(Particle.PID) == 6 && isLast").Define("nTop", "Sum(Top)")
                  .Define("Higgs", "abs(Particle.PID) == 25 && isLast").Define("nHiggs", "Sum(Higgs)")
                  .Define("W", "abs(Particle.PID) == 24 && isLast").Define("nW", "Sum(W)")
                  .Define("GenbQuark", "abs(Particle.PID) == 5 && isLast").Define("nGenbQ", "Sum(GenbQuark)")
                  .Define("GenAddbQuark", "abs(Particle.PID) == 5 && isLast && isAdd")
                  .Define("nGenAddbQ", "Sum(GenAddbQuark)")
                  .Define("GencQuark", "abs(Particle.PID) == 4 && isLast")

//                  .Filter("nGenbQ == 4")  // CAUTION. THIS OPTION IS ONLY TTBB
                  .Define("FinGenPtc_idx", ::FinalParticle_idx, {"Particle.PID", "Particle.PT", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Top", "Higgs"})
                  .Define("Higgs1_idx", "FinGenPtc_idx[4]")
                  .Define("Higgs2_idx", "FinGenPtc_idx[7]")
                  //pt
                  .Define("h1_pt", "Particle.PT[Higgs1_idx]")
                  .Define("h2_pt", "Particle.PT[Higgs2_idx]")
                  .Define("h_pt", ::ConcatFloat_WithoutSort, {"h1_pt", "h2_pt"})
                  //eta
                  .Define("h1_eta", "Particle.Eta[Higgs1_idx]")
                  .Define("h2_eta", "Particle.Eta[Higgs2_idx]")
                  .Define("h_eta", ::ConcatFloat_WithoutSort, {"h1_eta", "h2_eta"})
                  //phi
                  .Define("h1_phi", "Particle.Phi[Higgs1_idx]")
                  .Define("h2_phi", "Particle.Phi[Higgs2_idx]")
                  .Define("h_phi", ::ConcatFloat_WithoutSort, {"h1_phi", "h2_phi"})
                  //m
                  .Define("h1_m", "Particle.Mass[Higgs1_idx]")
                  .Define("h2_m", "Particle.Mass[Higgs2_idx]")
                  .Define("h_m", ::ConcatFloat_WithoutSort, {"h1_m", "h2_m"})
                  //dR
                  .Define("Higgs_Separation", ::dR2, {"h1_pt", "h1_eta", "h1_phi", "h1_m", "h2_pt", "h2_eta", "h2_phi", "h2_m"});
    // Reco //
    auto df2 = df1.Define("goodJet", "Jet.PT>=30 && abs(Jet.Eta)<3.0")
                  .Define("goodElectron", "Electron.PT>=23 && abs(Electron.Eta)<3.0")
                  .Define("goodMuon", "Muon.PT>=17 && abs(Muon.Eta)<2.8")

                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_mass", "Jet.Mass[goodJet]")
                  .Define("Jet_E", ::GetE, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
                  .Define("Jet_btag", "Jet.BTag[goodJet]==1 || Jet.BTag[goodJet]==3")
                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_mass", "Jet_mass[Jet_btag]")
                  .Define("bJet_E", ::GetE, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
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

                  .Define("Muon_IsolationVar", {"Muon.IsolationVar"})
                  .Define("Electron_IsolationVar", {"Electron.IsolationVar"})

                  .Define("j_ht", ::Ht, {"Jet_pt"})
                  // You Must Use [bi_Tag] Before Filter // 
                  

                  .Define("Jet1_pt", "Jet_pt[0]").Define("Jet1_eta", "Jet_eta[0]").Define("Jet1_phi", "Jet_phi[0]")
                  .Define("Jet1_mass", "bJet_mass[0]")
                  .Define("Jet2_pt", "Jet_pt[1]").Define("Jet2_eta", "Jet_eta[1]").Define("Jet2_phi", "Jet_phi[1]")
                  .Define("Jet2_mass", "Jet_mass[1]")
                  .Define("Jet3_pt", "Jet_pt[2]").Define("Jet3_eta", "Jet_eta[2]").Define("Jet3_phi", "Jet_phi[2]")
                  .Define("Jet3_mass", "Jet_mass[2]")
                  .Define("Jet4_pt", "Jet_pt[3]").Define("Jet4_eta", "Jet_eta[3]").Define("Jet4_phi", "Jet_phi[3]")
                  .Define("Jet4_mass", "Jet_mass[3]")
                  .Define("Jet5_pt", "Jet_pt[4]").Define("Jet5_eta", "Jet_eta[4]").Define("Jet5_phi", "Jet_phi[4]")
                  .Define("Jet5_mass", "Jet_mass[4]")
                  .Define("bJet1_pt", "bJet_pt[0]").Define("bJet1_eta", "bJet_eta[0]").Define("bJet1_phi", "bJet_phi[0]")
                  .Define("bJet1_mass", "bJet_mass[0]")
                  .Define("bJet2_pt", "bJet_pt[1]").Define("bJet2_eta", "bJet_eta[1]").Define("bJet2_phi", "bJet_phi[1]")
                  .Define("bJet2_mass", "bJet_mass[1]")
                  .Define("bJet3_pt", "bJet_pt[2]").Define("bJet3_eta", "bJet_eta[2]").Define("bJet3_phi", "bJet_phi[2]")
                  .Define("bJet3_mass", "bJet_mass[2]")
                  .Define("bJet4_pt", "bJet_pt[3]").Define("bJet4_eta", "bJet_eta[3]").Define("bJet4_phi", "bJet_phi[3]")
                  .Define("bJet4_mass", "bJet_mass[3]")
                  .Define("bJet5_pt", "bJet_pt[4]").Define("bJet5_eta", "bJet_eta[4]").Define("bJet5_phi", "bJet_phi[4]")
                  .Define("bJet5_mass", "bJet_mass[4]")

                  .Define("j1j2_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass"})
                  .Define("j1j3_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass"})
                  .Define("j1j4_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j1j5_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j2j3_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass"})
                  .Define("j2j4_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j2j5_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j3j4_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j3j5_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j4j5_dr", ::dR2, {"Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("jj_dr", ::ConcatFloat_withoutSort_10, {"j1j2_dr", "j1j3_dr", "j1j4_dr", "j1j5_dr", "j2j3_dr", "j2j4_dr", "j2j5_dr", "j3j4_dr", "j3j5_dr", "j4j5_dr"}) //Be Careful

                  .Define("j_Vars", ::Vars, {"jj_dr", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_E"})
                  .Define("jj_avg_dr", "j_Vars[0]")
                  .Define("jj_max_dr", "j_Vars[1]")
                  .Define("jj_min_dr", "j_Vars[2]")
                  .Define("jj_dEta_WhenMaxdR", "j_Vars[3]")
                  .Define("j_cent", "j_Vars[5]")
                  .Define("jj_max_deta", "j_Vars[6]")
                  .Define("jj_max_mass", "j_Vars[7]")
                  .Define("jj_twist", "j_Vars[8]")

                  .Define("b1b2_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass"})
                  .Define("b1b3_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass"})
                  .Define("b1b4_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("b1b5_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass"})
                  .Define("b2b3_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass"})
                  .Define("b2b4_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("b2b5_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass"})
                  .Define("b3b4_dr", ::dR2, {"bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("b3b5_dr", ::dR2, {"bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass"})
                  .Define("b4b5_dr", ::dR2, {"bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass", "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass"})
                  .Define("bb_dr", ::ConcatFloat_withoutSort_10, {"b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr", "b2b3_dr", "b2b4_dr", "b2b5_dr", "b3b4_dr", "b3b5_dr", "b4b5_dr"})

                  
                  .Define("b_Vars", ::Vars, {"bb_dr", "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_E"})
                  .Define("bb_avg_dr", "b_Vars[0]")
                  .Define("bb_max_dr", "b_Vars[1]")
                  .Define("bb_min_dr", "b_Vars[2]")
                  .Define("bb_dEta_WhenMaxdR", "b_Vars[3]")
                  .Define("b_ht", "b_Vars[4]")
                  .Define("b_cent", "b_Vars[5]")
                  .Define("bb_max_deta", "b_Vars[6]")
                  .Define("bb_max_mass", "b_Vars[7]")
                  .Define("bb_twist", "b_Vars[8]");


    std::initializer_list<std::string> variables = {

    //-----------Gen--------------------//

     "isAdd", "nGenbQ", "nGenAddbQ",
     "Higgs_Separation", 
     "h1_pt", "h1_eta", "h1_phi", "h1_m", "h2_pt", "h2_eta", "h2_phi", "h2_m",
     "h_pt", "h_eta", "h_phi", "h_m", 

    //-----------Cut Flow--------------------//
    "Electron_pt", "Muon_pt",
    //
    "Lep_size", "Lep1_pt", "Lep2_pt", "Lep1_eta", "Lep2_eta", "Lep1_ch", "Lep2_ch", "SS_OS_DL", 
    "Muon_IsolationVar", "Electron_IsolationVar", "nMuon", "nElectron", "Muon_charge", "Electron_charge",
    "Jet_size", "Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", 
    "j_ht", "MET_E",
    "bJet_size", "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", 
    "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass",
    "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass",
    "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass",
    "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass",
    "bb_max_dr", "bb_min_dr", "bb_dEta_WhenMaxdR", "b_ht", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist", 

    };

    // MODIFY!!
    df2.Snapshot(treename, outdir+ "inc_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
