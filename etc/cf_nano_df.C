
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
void cf_nano_df(std::string channel, std::string outdir="./samples1/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/common/NanoAOD/forTTHHto4b/RunIISummer20UL18NanoAODv9/TTHHTo4b_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/"+channel+"/*.root";
    std::cout << infile << std::endl;
    std::cout << outdir << std::endl;
    auto treename = "Events";

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
                  .Define("int3", "int(3)").Define("int4", "int(4)")
                  .Define("float0", "float(0)")
                  .Define("drmax1", "float(0.15)").Define("drmax2", "float(0.4)");
                  
    // Gen and Matching //
    auto df1 = df0.Define("nWlep", ::isDilep, {"GenPart_pdgId", "GenPart_genPartIdxMother"});
    
    // Reco //
    auto df2 = df1.Define("goodElectron", "Electron_pt>=10 && abs(Electron_eta)<2.8")
                  .Define("goodMuon", "Muon_pt>=10 && abs(Muon_eta)<2.8")

                  .Define("Mu_pt", "Muon_pt[goodMuon]")
                  .Define("Mu_eta", "Muon_eta[goodMuon]")
                  .Define("Mu_phi", "Muon_phi[goodMuon]")
                  .Define("Mu_mass", "Muon_mass[goodMuon]")
                  .Define("nMu", "Sum(goodMuon)")
                  .Define("Elec_pt", "Electron_pt[goodElectron]")
                  .Define("Elec_eta", "Electron_eta[goodElectron]")
                  .Define("Elec_phi", "Electron_phi[goodElectron]")
                  .Define("Elec_mass", "Electron_mass[goodElectron]")
                  .Define("nElec", "Sum(goodElectron)")
                  .Define("nLepton", "nMu + nElec");

    auto df3 = df2.Define("Mu_looseId", "Muon_looseId==1")
                  .Define("Mu_tightId", "Muon_tightId==1")
                  .Define("Mu_multiIsoId_L", "Muon_multiIsoId==1") // 1=MultiIsoLoose, 2=MultiIsoMedium
                  .Define("Mu_multiIsoId_M", "Muon_multiIsoId==2") // 1=MultiIsoLoose, 2=MultiIsoMedium
                  .Define("Mu_miniIsoId_L", "Muon_miniIsoId==1") // 1=MiniIsoLoose, 2=MiniIsoMedium, 3=Tight, 4=VeryTight
                  .Define("Mu_miniIsoId_M", "Muon_miniIsoId==2") // 1=MiniIsoLoose, 2=MiniIsoMedium, 3=Tight, 4=VeryTight
                  .Define("Mu_miniIsoId_T", "Muon_miniIsoId==3") // 1=MiniIsoLoose, 2=MiniIsoMedium, 3=Tight, 4=VeryTight
                  .Define("Mu_miniIsoId_VT", "Muon_miniIsoId==4") // 1=MiniIsoLoose, 2=MiniIsoMedium, 3=Tight, 4=VeryTight
                  .Define("Elec_MVAisoLoose", "Electron_mvaFall17V2Iso_WPL==1")
                  .Define("Elec_MVAiso90", "Electron_mvaFall17V2Iso_WP90==1")
                  .Define("Elec_MVAiso80", "Electron_mvaFall17V2Iso_WP80==1")

                  .Define("nMu_looseId", ::Multiplicity, {"Mu_looseId", "int1"})
                  .Define("nMu_tightId", ::Multiplicity, {"Mu_tightId", "int1"})
                  .Define("nMu_multiIsoId_L", ::Multiplicity, {"Mu_multiIsoId_L", "int1"})
                  .Define("nMu_multiIsoId_M", ::Multiplicity, {"Mu_multiIsoId_M", "int2"})
                  .Define("nMu_miniIsoId_L", ::Multiplicity, {"Mu_miniIsoId_L", "int1"})
                  .Define("nMu_miniIsoId_M", ::Multiplicity, {"Mu_miniIsoId_M", "int2"})
                  .Define("nMu_miniIsoId_T", ::Multiplicity, {"Mu_miniIsoId_T", "int3"})
                  .Define("nMu_miniIsoId_VT", ::Multiplicity, {"Mu_miniIsoId_VT", "int4"})
                  .Define("nElec_MVAisoLoose", ::Multiplicity, {"Elec_MVAisoLoose", "int1"})
                  .Define("nElec_MVAiso90", ::Multiplicity, {"Elec_MVAiso90", "int1"})
                  .Define("nElec_MVAiso80", ::Multiplicity, {"Elec_MVAiso80", "int1"});
                 // You Must Use [bi_Tag] Before Filter // 
                  


//    auto df6 = df5.Filter("isMatchable > 0"); 

    std::initializer_list<std::string> variables = {

    "nWlep",
    //-----------Cut Flow--------------------//
    "nMu", "nElec", "nLepton",
    "Mu_pt", "Mu_eta", "Mu_phi", "Mu_mass", 
    "Elec_pt", "Elec_eta", "Elec_phi", "Elec_mass",

    //-------------ID------------------------//
    "Mu_looseId", "Mu_tightId", "Mu_multiIsoId_L", "Mu_multiIsoId_M", "Mu_miniIsoId_L", "Mu_miniIsoId_M", "Mu_miniIsoId_T", "Mu_miniIsoId_VT",
    "Elec_MVAisoLoose", "Elec_MVAiso90", "Elec_MVAiso80",

    "nMu_looseId", "nMu_tightId", "nMu_multiIsoId_L", "nMu_multiIsoId_M", "nMu_miniIsoId_L", "nMu_miniIsoId_M", "nMu_miniIsoId_T", "nMu_miniIsoId_VT",
    "nElec_MVAisoLoose", "nElec_MVAiso90", "nElec_MVAiso80"
        
    };

    // MODIFY!!
    df3.Snapshot(treename, outdir+ "CF_0326_" + "tthh_CMS" + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
