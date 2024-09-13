
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TMath.h>

#include "utility2.h"

// MODIFY!!
void gen_tthh_df(std::string channel, std::string outdir="./samples1/"){
    gSystem->Load("libDelphes");

//    auto infile = "/data1/users/itseyes/tthh/13.6_di_test_MET/"+channel+"/*.root";
    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/*.root";
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
                  .Define("JetBTag", {"Jet.BTag"})
                  .Define("GenMissingET_met", "GenMissingET.MET")
                  .Define("GenMissingET_eta", "GenMissingET.Eta")
                  .Define("GenMissingET_phi", "GenMissingET.Phi");
                  

    // Gen_bi //   // Quark, FromQuark, AddQuark
    auto df1 = df0.Define("isLast", ::isLast, {"Particle.PID", "Particle.D1", "Particle.D2"})
                  .Define("Top", "abs(Particle.PID) == 6 && isLast").Define("nTop", "Sum(Top)")
                  .Define("Higgs", "abs(Particle.PID) == 25 && isLast").Define("nHiggs", "Sum(Higgs)")
                  .Define("W", "abs(Particle.PID) == 24 && isLast").Define("nW", "Sum(W)")
                  .Define("GenbQuark", "abs(Particle.PID) == 5 && isLast")
                  .Define("GencQuark", "abs(Particle.PID) == 4 && isLast")

                  .Define("FinalGenPart_idx", ::FinalGenPart_idx, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Top", "Higgs"})
                  .Define("GenbFromTop1_idx", "FinalGenPart_idx[1]")
                  .Define("GenlepFromTop1_idx", "FinalGenPart_idx[2]")
                  .Define("GenbFromTop2_idx", "FinalGenPart_idx[4]")
                  .Define("GenlepFromTop2_idx", "FinalGenPart_idx[5]")

                  .Define("GenbJetFromTop1_idx", ::dRMatching_idx, {"GenbFromTop1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenbJetFromTop2_idx", ::dRMatching_idx, {"GenbFromTop2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})

                  .Define("bJetFromTop1_idx", ::dRMatching_idx, {"GenbFromTop1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bJetFromTop2_idx", ::dRMatching_idx, {"GenbFromTop2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})

                  .Define("muonFromTop1_idx", ::dRMatching_idx, {"GenlepFromTop1_idx", "drmax1", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "Muon.PT", "Muon.Eta", "Muon.Phi", "Muon.PT"})
                  .Define("elecFromTop1_idx", ::dRMatching_idx, {"GenlepFromTop1_idx", "drmax1", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "Electron.PT", "Electron.Eta", "Electron.Phi", "Electron.PT"})
                  .Define("muonFromTop2_idx", ::dRMatching_idx, {"GenlepFromTop1_idx", "drmax1", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "Muon.PT", "Muon.Eta", "Muon.Phi", "Muon.PT"})
                  .Define("elecFromTop2_idx", ::dRMatching_idx, {"GenlepFromTop1_idx", "drmax1", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "Electron.PT", "Electron.Eta", "Electron.Phi", "Electron.PT"});

    // 4 Vector of Particle //
    auto df2 = df1.Filter("GenbFromTop1_idx>=0 && GenlepFromTop1_idx>=0 && GenbFromTop2_idx>=0 && GenlepFromTop2_idx>=0")
                  .Define("GenbFromTop1_pt",  "Particle.PT[GenbFromTop1_idx]")
                  .Define("GenbFromTop1_eta",  "Particle.Eta[GenbFromTop1_idx]")
                  .Define("GenbFromTop1_phi",  "Particle.Phi[GenbFromTop1_idx]")
                  .Define("GenbFromTop1_mass",  "Particle.Mass[GenbFromTop1_idx]")
                  .Define("GenlepFromTop1_pt", "Particle.PT[GenlepFromTop1_idx]")
                  .Define("GenlepFromTop1_eta", "Particle.Eta[GenlepFromTop1_idx]")
                  .Define("GenlepFromTop1_phi", "Particle.Phi[GenlepFromTop1_idx]")
                  .Define("GenlepFromTop1_mass", "Particle.Mass[GenlepFromTop1_idx]")
                  .Define("GenbFromTop2_pt",  "Particle.PT[GenbFromTop2_idx]")
                  .Define("GenbFromTop2_eta",  "Particle.Eta[GenbFromTop2_idx]")
                  .Define("GenbFromTop2_phi",  "Particle.Phi[GenbFromTop2_idx]")
                  .Define("GenbFromTop2_mass",  "Particle.Mass[GenbFromTop2_idx]")
                  .Define("GenlepFromTop2_pt", "Particle.PT[GenlepFromTop2_idx]")
                  .Define("GenlepFromTop2_eta", "Particle.Eta[GenlepFromTop2_idx]")
                  .Define("GenlepFromTop2_phi", "Particle.Phi[GenlepFromTop2_idx]")
                  .Define("GenlepFromTop2_mass", "Particle.Mass[GenlepFromTop2_idx]")

                  .Define("GenlepFromTop_dr", ::dR2, {"GenlepFromTop1_pt", "GenlepFromTop1_eta", "GenlepFromTop1_phi", "GenlepFromTop1_mass", "GenlepFromTop2_pt", "GenlepFromTop2_eta", "GenlepFromTop2_phi", "GenlepFromTop2_mass"});
                 

    // 4 Vector of GenJet //
    auto df3 = df2.Filter("GenbJetFromTop1_idx>=0 && GenbJetFromTop2_idx>=0")
                  .Define("GenbJetFromTop1_pt",  "GenJet.PT[GenbJetFromTop1_idx]")
                  .Define("GenbJetFromTop1_eta",  "GenJet.Eta[GenbJetFromTop1_idx]")
                  .Define("GenbJetFromTop1_phi",  "GenJet.Phi[GenbJetFromTop1_idx]")
                  .Define("GenbJetFromTop1_mass",  "GenJet.Mass[GenbJetFromTop1_idx]")
                  .Define("GenbJetFromTop2_pt",  "GenJet.PT[GenbJetFromTop2_idx]")
                  .Define("GenbJetFromTop2_eta",  "GenJet.Eta[GenbJetFromTop2_idx]")
                  .Define("GenbJetFromTop2_phi",  "GenJet.Phi[GenbJetFromTop2_idx]")
                  .Define("GenbJetFromTop2_mass",  "GenJet.Mass[GenbJetFromTop2_idx]");

    auto df4 = df3.Define("GenTop1Mass", ::GenTopMass, {"GenbJetFromTop1_pt", "GenbJetFromTop1_eta", "GenbJetFromTop1_phi", "GenbJetFromTop1_mass", "GenlepFromTop1_pt", "GenlepFromTop1_eta", "GenlepFromTop1_phi", "GenlepFromTop1_mass"})
                  .Define("GenTop2Mass", ::GenTopMass, {"GenbJetFromTop2_pt", "GenbJetFromTop2_eta", "GenbJetFromTop2_phi", "GenbJetFromTop2_mass", "GenlepFromTop2_pt", "GenlepFromTop2_eta", "GenlepFromTop2_phi", "GenlepFromTop2_mass"})
                  .Define("GenTopMass", ::ConcatFloat, {"GenTop1Mass", "GenTop2Mass"});

    std::initializer_list<std::string> variables = {
                      // Basic Branches // 
    "ParticlePID", "ParticlePT", "D1", "D2",
                      // Gen Quark_bi //

                      // Gen Jet_bi //

                      // nGen Quark//

                      // nGen Jet //

                      // Gen Quark 4-vector //
    "GenbFromTop1_pt",
                      // Gen Jet 4-vector //
    "GenbJetFromTop1_pt",                    
                      // Gen dR //
                      
                      // Recoed Gen Mass //
    "GenTop1Mass", "GenTop2Mass", "GenTopMass",
                      //-------------------------------------- Reco----------------------------------------//
                      

                      //-------------------------------------- New ----------------------------------------//
    "GenlepFromTop_dr"                       
    };

    // MODIFY!!
    df4.Snapshot(treename, outdir+ "Past_di_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
