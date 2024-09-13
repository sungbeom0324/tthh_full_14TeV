import ROOT
from ROOT import RooFit, RooRealVar, RooDataHist, RooHistPdf, RooAddPdf, RooFitResult, RooPlot, RooStats

def fit_ttHH_significance():
    # Load the ROOT file
    infile = ROOT.TFile.Open("path_to_your_file.root")
    
    # Get the TTree
    tree = infile.Get("Delphes")

    # Define histograms
    h_ttHH = ROOT.TH1F("h_ttHH", "ttHH", 100, 0, 1) # 신호 히스토그램 (ttHH)
    h_ttbb = ROOT.TH1F("h_ttbb", "ttbb", 100, 0, 1) # 백그라운드 히스토그램 (ttbb)
    h_ttbbbb = ROOT.TH1F("h_ttbbbb", "ttbbbb", 100, 0, 1) # 백그라운드 히스토그램 (ttbbbb)
    h_ttHbb = ROOT.TH1F("h_ttHbb", "ttHbb", 100, 0, 1) # 백그라운드 히스토그램 (ttHbb)

    # Fill histograms from branches
    tree.Draw("DNN_ttHH >> h_ttHH", "DNN_ttHH")
    tree.Draw("DNN_ttbb >> h_ttbb", "DNN_ttbb")
    tree.Draw("DNN_ttbbbb >> h_ttbbbb", "DNN_ttbbbb")
    tree.Draw("DNN_ttHbb >> h_ttHbb", "DNN_ttHbb")

    # Define observable
    x = RooRealVar("x", "DNN Score", 0, 1)

    # Convert histograms to RooDataHist
    dh_ttHH = RooDataHist("dh_ttHH", "ttHH", RooArgList(x), h_ttHH)
    dh_ttbb = RooDataHist("dh_ttbb", "ttbb", RooArgList(x), h_ttbb)
    dh_ttbbbb = RooDataHist("dh_ttbbbb", "ttbbbb", RooArgList(x), h_ttbbbb)
    dh_ttHbb = RooDataHist("dh_ttHbb", "ttHbb", RooArgList(x), h_ttHbb)

    # Convert RooDataHist to RooHistPdf
    pdf_ttHH = RooHistPdf("pdf_ttHH", "ttHH PDF", RooArgSet(x), dh_ttHH)
    pdf_ttbb = RooHistPdf("pdf_ttbb", "ttbb PDF", RooArgSet(x), dh_ttbb)
    pdf_ttbbbb = RooHistPdf("pdf_ttbbbb", "ttbbbb PDF", RooArgSet(x), dh_ttbbbb)
    pdf_ttHbb = RooHistPdf("pdf_ttHbb", "ttHbb PDF", RooArgSet(x), dh_ttHbb)

    # Define coefficients for signal and background
    mu = RooRealVar("mu", "signal strength", 1, 0, 10)
    n_ttbb = RooRealVar("n_ttbb", "ttbb events", h_ttbb.Integral())
    n_ttbbbb = RooRealVar("n_ttbbbb", "ttbbbb events", h_ttbbbb.Integral())
    n_ttHbb = RooRealVar("n_ttHbb", "ttHbb events", h_ttHbb.Integral())

    # Model: total PDF (signal + background)
    model = RooAddPdf("model", "signal + background", RooArgList(pdf_ttHH, pdf_ttbb, pdf_ttbbbb, pdf_ttHbb), RooArgList(mu, n_ttbb, n_ttbbbb, n_ttHbb))

    # Generate a toy dataset based on the total model
    dh_data = model.generateBinned(RooArgSet(x), int(h_ttHH.GetEntries() + h_ttbb.GetEntries() + h_ttbbbb.GetEntries() + h_ttHbb.GetEntries()), RooFit.Extended())

    # Fit model to data
    fitResult = model.fitTo(dh_data, RooFit.Save())

    # Profile likelihood calculator
    plc = RooStats.ProfileLikelihoodCalculator(dh_data, model, RooArgSet(mu))
    interval = plc.GetInterval()

    # Print results
    print("Fit results:")
    fitResult.Print("v")

    print("Profile likelihood results:")
    print(f"  68% CI for mu: [{interval.LowerLimit(mu)}, {interval.UpperLimit(mu)}]")

    # HypoTestInverter for significance calculation
    ac = RooStats.AsymptoticCalculator(dh_data, model, model)
    calc = RooStats.HypoTestInverter(ac)
    calc.SetConfidenceLevel(0.95)
    calc.UseCLs(True)

    hypoTestResult = calc.GetInterval()
    significance = hypoTestResult.Significance()

    print(f"Significance: {significance} sigma")

    # Plot results
    frame = x.frame()
    dh_data.plotOn(frame)
    model.plotOn(frame)

    canvas = ROOT.TCanvas("canvas", "Fit results", 800, 600)
    frame.Draw()
    canvas.SaveAs("fit_results.png")

# Run the function
fit_ttHH_significance()

