import ROOT
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgSet, RooDataHist, RooHistPdf, RooAddPdf, RooArgList

# 1. ROOT 파일 읽기
file = ROOT.TFile("/home/stiger97/github/tthh_14TeV/DNN_result/FULL_0522/LetsFind_tthh/bCat_higgs5_2Mat/dnn2_result.root")

# TTree 가져오기
tree = file.Get("Delphes")  # tree_name을 실제 트리 이름으로 대체

# DNN 변수와 category 변수 정의
dnn_output = RooRealVar("DNN", "DNN output", 0.0, 1.0)
category = RooRealVar("category", "event category", 0, 3)

# RooDataSet 생성
data = RooDataSet("data", "data", ROOT.RooArgSet(dnn_output, category), RooFit.Import(tree))


# 2. 신호와 데이터셋 분리
# 신호와 배경을 분리하는 카테고리 필터링
signal_data = data.reduce("category==0")
tthbb_data = data.reduce("category==1")
ttbb_data = data.reduce("category==2")
ttbbbb_data = data.reduce("category==3")

# 배경 데이터셋 병합
bkg_data = tthbb_data.Clone()
bkg_data.append(ttbb_data)
bkg_data.append(ttbbbb_data)

# 3. 히스토그램 기반 PDF 생성 ("명확한 Gaussian 과 같은 분포가 보이지 않을 때 히스토그램 기반으로 정의할 수 있음.")
# 히스토그램 정의
signal_hist = ROOT.TH1F("signal_hist", "Signal Histogram", 100, 0, 1)
bkg_hist = ROOT.TH1F("bkg_hist", "Background Histogram", 100, 0, 1)

# 히스토그램에 데이터 채우기
for i in range(signal_data.numEntries()):
    dnn_output.setVal(signal_data.get(i).getRealValue("DNN"))
    signal_hist.Fill(dnn_output.getVal())

for i in range(bkg_data.numEntries()):
    dnn_output.setVal(bkg_data.get(i).getRealValue("DNN"))
    bkg_hist.Fill(dnn_output.getVal())

# 히스토그램을 RooDataHist로 변환
signal_roodatahist = RooDataHist("signal_roodatahist", "Signal RooDataHist", RooArgSet(dnn_output), signal_hist)
bkg_roodatahist = RooDataHist("bkg_roodatahist", "Background RooDataHist", RooArgSet(dnn_output), bkg_hist)

# 히스토그램 기반 PDF 생성
signal_pdf = RooHistPdf("signal_pdf", "Signal PDF", RooArgSet(dnn_output), signal_roodatahist)
bkg_pdf = RooHistPdf("bkg_pdf", "Background PDF", RooArgSet(dnn_output), bkg_roodatahist)

# 4. 모델 정의 및 피팅
nsig = RooRealVar("nsig", "number of signal events", 100, 0, 10000)
nbkg = RooRealVar("nbkg", "number of background events", 1000, 0, 100000)
model = RooAddPdf("model", "Signal + Background", RooArgList(signal_pdf, bkg_pdf), RooArgList(nsig, nbkg))

# 데이터 피팅
fit_result = model.fitTo(data, RooFit.Save())

# 5. 신호 한계 계산.
import ROOT.RooStats as RooStats

plc = RooStats.ProfileLikelihoodCalculator(data, model, ROOT.RooArgSet(nsig))
plc.SetConfidenceLevel(0.95)
interval = plc.GetInterval()

lowerLimit = interval.LowerLimit(nsig)
upperLimit = interval.UpperLimit(nsig)

print(f"95% 신뢰 수준에서 신호의 한계: [{lowerLimit}, {upperLimit}]")

