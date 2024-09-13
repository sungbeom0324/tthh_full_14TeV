
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
void reco_tthh_df(std::string channel, std::string outdir="./samples1/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/tag_1*.root";
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
                  .Define("JetBTag", {"Jet.BTag"});

    // Reco //
    auto df1 = df0.Define("goodJet", "Jet.PT>=30 && abs(Jet.Eta)<2.4")
                  .Define("goodElectron", "Electron.PT>=20 && abs(Electron.Eta)<2.4")
                  .Define("goodMuon", "Muon.PT>=20 && abs(Muon.Eta)<2.4")
                    
                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_mass", "Jet.Mass[goodJet]")
                  .Define("Jet_btag", "Jet.BTag[goodJet]")
                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_mass", "Jet_mass[Jet_btag]")
                  .Define("bJet_size", "Sum(Jet_btag)")

                  .Define("Muon_pt", "Muon.PT[goodMuon]")
                  .Define("Muon_eta", "Muon.Eta[goodMuon]")
                  .Define("Muon_phi", "Muon.Phi[goodMuon]")
                  .Define("nMuon", "Sum(goodMuon)")
                  
                  .Define("Electron_pt", "Electron.PT[goodElectron]")
                  .Define("Electron_eta", "Electron.Eta[goodElectron]")
                  .Define("Electron_phi", "Electron.Phi[goodElectron]")
                  .Define("nElectron", "Sum(goodElectron)")

                  .Define("Lep_size", "nMuon + nElectron")

                  .Filter("Lep_size == 2")
                  .Define("Lep_4vec", ::TwoLeptons, {"Muon_pt", "Muon_eta", "Muon_phi", "Electron_pt", "Electron_eta", "Electron_phi"})
                  .Define("Lep1_pt", "Lep_4vec[0]")
                  .Define("Lep1_eta", "Lep_4vec[1]")
                  .Define("Lep1_phi", "Lep_4vec[2]")
                  .Define("Lep2_pt", "Lep_4vec[3]")
                  .Define("Lep2_eta", "Lep_4vec[4]")
                  .Define("Lep2_phi", "Lep_4vec[5]");


    // Reconstruction of Higgs //
    auto df2 = df1.Define("RecoHiggs", ::RecoHiggs, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
                  .Define("Higgs_pt", "RecoHiggs[0]")
                  .Define("Higgs_eta", "RecoHiggs[1]")
                  .Define("Higgs_phi", "RecoHiggs[2]")
                  .Define("Higgs_mass", "RecoHiggs[3]");

    std::initializer_list<std::string> variables = {
                      
                      //-------------------------------------- Reco----------------------------------------//
     "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_size",
     "Muon_pt", "Muon_eta", "Muon_phi", "nMuon",
     "Electron_pt", "Electron_eta", "Electron_phi", "nElectron",
     "Lep1_pt", "Lep1_eta", "Lep1_phi",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", 
                      //-------------------------------------- New ----------------------------------------//
                        
    };

    // MODIFY!!
    df2.Snapshot(treename, outdir+ "Strange_J30L20_0.2dpT_NodR_Lep4vec_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
