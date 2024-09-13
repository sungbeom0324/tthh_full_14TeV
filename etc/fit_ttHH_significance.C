#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestPlot.h"

void fit_ttHH_significance() {
    // Load the ROOT file
    TFile* file = TFile::Open("path_to_your_file.root");
    
    // Get the TTree
    TTree* tree = (TTree*)file->Get("Delphes");

    // Define histograms
    TH1F* h_ttHH = new TH1F("h_ttHH", "ttHH", 100, 0, 1); // 신호 히스토그램 (ttHH)
    TH1F* h_ttbb = new TH1F("h_ttbb", "ttbb", 100, 0, 1); // 백그라운드 히스토그램 (ttbb)
    TH1F* h_ttbbbb = new TH1F("h_ttbbbb", "ttbbbb", 100, 0, 1); // 백그라운드 히스토그램 (ttbbbb)
    TH1F* h_ttHbb = new TH1F("h_ttHbb", "ttHbb", 100, 0, 1); // 백그라운드 히스토그램 (ttHbb)

    // Fill histograms from branches
    tree->Draw("DNN_ttHH >> h_ttHH", "DNN_ttHH");
    tree->Draw("DNN_ttbb >> h_ttbb", "DNN_ttbb");
    tree->Draw("DNN_ttbbbb >> h_ttbbbb", "DNN_ttbbbb");
    tree->Draw("DNN_ttHbb >> h_ttHbb", "DNN_ttHbb");

    // Define observable
    RooRealVar x("x", "DNN Score", 0, 1);

    // Convert histograms to RooDataHist
    RooDataHist dh_ttHH("dh_ttHH", "ttHH", x, Import(*h_ttHH));
    RooDataHist dh_ttbb("dh_ttbb", "ttbb", x, Import(*h_ttbb));
    RooDataHist dh_ttbbbb("dh_ttbbbb", "ttbbbb", x, Import(*h_ttbbbb));
    RooDataHist dh_ttHbb("dh_ttHbb", "ttHbb", x, Import(*h_ttHbb));

    // Convert RooDataHist to RooHistPdf
    RooHistPdf pdf_ttHH("pdf_ttHH", "ttHH PDF", x, dh_ttHH);
    RooHistPdf pdf_ttbb("pdf_ttbb", "ttbb PDF", x, dh_ttbb);
    RooHistPdf pdf_ttbbbb("pdf_ttbbbb", "ttbbbb PDF", x, dh_ttbbbb);
    RooHistPdf pdf_ttHbb("pdf_ttHbb", "ttHbb PDF", x, dh_ttHbb);

    // Define coefficients for signal and background
    RooRealVar mu("mu", "signal strength", 1, 0, 10);
    RooRealVar n_ttbb("n_ttbb", "ttbb events", h_ttbb->Integral());
    RooRealVar n_ttbbbb("n_ttbbbb", "ttbbbb events", h_ttbbbb->Integral());
    RooRealVar n_ttHbb("n_ttHbb", "ttHbb events", h_ttHbb->Integral());

    // Model: total PDF (signal + background)
    RooAddPdf model("model", "signal + background", RooArgList(pdf_ttHH, pdf_ttbb, pdf_ttbbbb, pdf_ttHbb), RooArgList(mu, n_ttbb, n_ttbbbb, n_ttHbb));

    // Generate a toy dataset based on the total model
    RooDataHist* dh_data = model.generateBinned(x, h_ttHH->GetEntries() + h_ttbb->GetEntries() + h_ttbbbb->GetEntries() + h_ttHbb->GetEntries(), RooFit::Extended());

    // Fit model to data
    RooFitResult* fitResult = model.fitTo(*dh_data, Save());

    // Profile likelihood calculator
    RooStats::ProfileLikelihoodCalculator plc(*dh_data, model, RooArgSet(mu));
    RooStats::LikelihoodInterval* interval = plc.GetInterval();

    // Print results
    std::cout << "Fit results:" << std::endl;
    fitResult->Print("v");

    std::cout << "Profile likelihood results:" << std::endl;
    std::cout << "  68% CI for mu: [" << interval->LowerLimit(mu) << ", " << interval->UpperLimit(mu) << "]" << std::endl;

    // HypoTestInverter for significance calculation
    RooStats::AsymptoticCalculator ac(*dh_data, model, model);
    RooStats::HypoTestInverter calc("calc", "", ac);
    calc.SetConfidenceLevel(0.95);
    calc.UseCLs(true);

    RooStats::HypoTestResult* hypoTestResult = calc.GetInterval();
    double significance = hypoTestResult->Significance();

    std::cout << "Significance: " << significance << " sigma" << std::endl;

    // Plot results
    RooPlot* frame = x.frame();
    dh_data->plotOn(frame);
    model.plotOn(frame);

    TCanvas* canvas = new TCanvas("canvas", "Fit results", 800, 600);
    frame->Draw();
    canvas->SaveAs("fit_results.png");

    // Clean up
    delete file;
    delete h_ttHH;
    delete h_ttbb;
    delete h_ttbbbb;
    delete h_ttHbb;
    delete dh_data;
}


