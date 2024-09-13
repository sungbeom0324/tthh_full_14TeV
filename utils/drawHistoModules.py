import ROOT
import numpy
import os

def drawHistoSingle(infile, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Same/"
    try:
        os.makedirs(outdir)
    except:
        pass

    # tthh 샘플만 사용
    tthh = ROOT.RDataFrame(tree, infile)

    canvas = ROOT.TCanvas("c", "c", 400, 400)
    canvas.SetLeftMargin(0.15)  # Y축 레이블을 위한 왼쪽 여백 조정

    # 범례 설정
    legend_position = {1: (0.72, 0.84, 0.87, 0.46), 2: (0.38, 0.855, 0.57, 0.705), 3: (0.72, 0.84, 0.87, 0.76)}
    legend = ROOT.TLegend(*legend_position.get(lepo, legend_position[3]))
    legend.SetTextSize(0.028)
    legend.SetEntrySeparation(0.20)

    # 히스토그램 생성
    h = tthh.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch)
    ymax = h.GetMaximum() 
    h.SetMaximum(1.2 * ymax)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(4)
#    h.SetLineStyle(2)  # 대시(dashed) 라인 스타일 설정
    h.SetStats(0)
    legend.AddEntry(h.GetValue(), "t\\bar{t}HH", "f")

    h.DrawNormalized("hist")

    legend.SetBorderSize(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.19, 0.84, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.67, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")

    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch + ".pdf")
    canvas.Clear()

def drawHistoSingle_2Branches(infile, tree, title, xtitle, ytitle, branch1, branch2, nbin, xmin, xmax, PRE, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Same/"
    try:
        os.makedirs(outdir)
    except:
        pass

    # tthh 샘플만 사용
    tthh = ROOT.RDataFrame(tree, infile).Filter("Genbfh_dr[0]>0")

    canvas = ROOT.TCanvas("c", "c", 400, 400)
    canvas.SetLeftMargin(0.15)  # Y축 레이블을 위한 왼쪽 여백 조정

    # 범례 설정
    legend_position = {1: (0.72, 0.84, 0.87, 0.46), 2: (0.38, 0.855, 0.57, 0.705), 3: (0.66, 0.85, 0.81, 0.77)}
    legend = ROOT.TLegend(*legend_position.get(lepo, legend_position[3]))
    legend.SetTextSize(0.028)
    legend.SetEntrySeparation(0.20)

    # 첫 번째 히스토그램 생성
    h1 = tthh.Histo1D(ROOT.RDF.TH1DModel(title + " - " + branch1, title, nbin, xmin, xmax), branch1)
    ymax1 = h1.GetMaximum()
    h1.SetMaximum(1.3 * ymax1)
    h1.GetXaxis().SetTitle(xtitle)
    h1.GetYaxis().SetTitle(ytitle)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetLineWidth(4)
    h1.SetStats(0)
    legend.AddEntry(h1.GetValue(), "b pair (all)", "l")

    # 두 번째 히스토그램 생성
    h2 = tthh.Histo1D(ROOT.RDF.TH1DModel(title + " - " + branch2, title, nbin, xmin, xmax), branch2)
#    ymax2 = h2.GetMaximum()
#    h2.SetMaximum(1.2 * ymax2)
    h2.SetLineColor(ROOT.kBlue)
    h2.SetLineWidth(4)
    h2.SetStats(0)
    legend.AddEntry(h2.GetValue(), "b pair from Higgs", "l")

    # 첫 번째 히스토그램을 그린 후 두 번째 히스토그램을 같은 캔버스에 그리기
    h1.DrawNormalized("hist")
    h2.DrawNormalized("hist same")

    legend.SetBorderSize(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.19, 0.84, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.67, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")

    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch1 + "_vs_" + branch2 + ".pdf")
    canvas.Clear() 

def drawHistoSingleSum(infile, tree, title, xtitle, ytitle, branch1, branch2, nbin, xmin, xmax, PRE):
    outdir = "./plots/" + PRE + "/"
    try : os.makedirs(outdir)
    except : pass
    df = ROOT.RDataFrame(tree, infile)
    df = df.Define("sumbranch", "{} + {}".format(branch1, branch2))     
    canvas = ROOT.TCanvas("c1", "c1")
    h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), "sumbranch")
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.Draw("colz")
    title = title.replace(" ", "_")
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")
    canvas.Print(outdir+ title + ".pdf")
    canvas.Clear()

def drawHisto2D(infile, tree, title, xtitle, ytitle, xbranch, nxbin, xmin, xmax, ybranch, nybin, ymin, ymax, PRE):
    outdir = "./plots/" + PRE + "/2D/"
    try:
        os.makedirs(outdir)
    except:
        pass
    
    df = ROOT.RDataFrame(tree, infile)
    
    # Create canvas with adjusted margin
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    canvas.SetLeftMargin(0.15)  # Adjust left margin for Y-axis label

    # Create 2D histogram
    h = df.Histo2D(ROOT.RDF.TH2DModel("", "", nxbin, xmin, xmax, nybin, ymin, ymax), xbranch, ybranch)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.SetStats(0)
    
    # Draw 2D histogram
    h.Draw("colz")

    # Add text labels to the canvas
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.19, 0.91, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.67, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")
    
    # Save the canvas to a file
    title = title.replace(" ", "_")
    canvas.Print(outdir + title + ".pdf")
    canvas.Clear()

def drawHistoSame(indir, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Same/"
    try : os.makedirs(outdir)
    except : pass
    
    tthh   = ROOT.RDataFrame(tree, indir + PRE + "_tthh.root")
    tth    = ROOT.RDataFrame(tree, indir + PRE + "_tth.root")
    ttbbh  = ROOT.RDataFrame(tree, indir + PRE + "_ttbbh.root")
    ttzh   = ROOT.RDataFrame(tree, indir + PRE + "_ttzh.root")
    ttvv   = ROOT.RDataFrame(tree, indir + PRE + "_ttvv.root")
    ttbbv  = ROOT.RDataFrame(tree, indir + PRE + "_ttbbv.root")
    ttbb   = ROOT.RDataFrame(tree, indir + PRE + "_ttbb.root")
    ttbbbb = ROOT.RDataFrame(tree, indir + PRE + "_ttbbbb.root")
    tttt   = ROOT.RDataFrame(tree, indir + PRE + "_tttt.root").Range(10000)
    ROOT.gStyle.SetPadTickX(1)  # X축의 위쪽에 tick 추가
    ROOT.gStyle.SetPadTickY(1)  # Y축의 오른쪽에 tick 추가

    dfs = {"ttHH":tthh, "ttH":tth, "ttbbH":ttbbh, "ttZH":ttzh, "ttVV":ttvv, "ttbbV":ttbbv, "ttbb":ttbb, "ttbbbb":ttbbbb, "tttt":tttt}
#    dfs = {}
#    for key, _df in _dfs.items():
#        dfs[key] = _df.Filter("Lep_size==2").Filter("SS_OS_DL==-1").Filter("bJet_size>=5").Filter("j_ht>300") # modify! #

#    dfs = {"ttHH": tthh, "ttH": tth, "ttbbH": ttbbh, "ttZH": ttzh, 
#           "ttVV": ttvv, "ttbbV": ttbbv, "ttbb": ttbb, "ttbbbb": ttbbbb, "tttt": tttt}

    legs = {
        "ttHH" : "t\\bar{t}HH",
        "ttH"  : "t\\bar{t}H", "ttbbH" : "t\\bar{t}b\\bar{b}H", "ttZH" : "t\\bar{t}ZH",
        "ttVV" : "t\\bar{t}VV", "ttbbV" : "t\\bar{t}b\\bar{b}V", 
        "ttbb" : "t\\bar{t}b\\bar{b}", "ttbbbb" : "t\\bar{t}b\\bar{b}b\\bar{b}", "tttt" : "t\\bar{t}t\\bar{t}"
    }

    colors = {
        "ttHH": ROOT.kBlack,
        "ttH" : ROOT.kGray+1, "ttbbH": ROOT.kGray + 3, "ttZH": ROOT.kGray + 5,
        "ttVV": ROOT.kGreen, "ttbbV": ROOT.kGreen + 2,
        "ttbb": ROOT.kBlue, "ttbbbb": ROOT.kCyan, 
        "tttt": ROOT.kRed
    }
    
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    canvas.SetLeftMargin(0.15)  # Adjust the left margin to make space for the y-axis label
    # lepo : xy xy
    legend_position = {1: (0.72, 0.84, 0.87, 0.46), 2: (0.38, 0.855, 0.57, 0.705), 3: (0.26, 0.855, 0.45, 0.705)}
    legend = ROOT.TLegend(*legend_position.get(lepo, legend_position[3]))
    legend.SetTextSize(0.028)  # Adjust this value to make the legend larger or smaller
    legend.SetEntrySeparation(0.20)  # Adjust the separation between legend entries
    ymax = 0
   
    hist_dict = {}
    for df_name, df in dfs.items():
        h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch)
        if h.GetEntries() == 0:
            continue
        if ymax < h.GetMaximum(): ymax = h.GetMaximum() # 실제 수가 큰게 들어가면 그게 ymax 되고 그 이후 norm.
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.SetLineColor(colors[df_name])
        h.SetLineWidth(2)
        h.SetStats(0)
        if df_name == "ttHH":
            h.SetLineStyle(2)  # Dashed line style
            h.SetLineWidth(4)  # Thick line
            legend.AddEntry(h.GetValue(), " " + legs[df_name], "f")
        else:
            legend.AddEntry(h.GetValue(), " " + legs[df_name], "f")
        hist_dict[branch + "_" + df_name] = h

    first = True
    for _tmp, h in hist_dict.items():
        if first:
            h.SetMaximum(1.2 * ymax)  # modify #
            h.DrawNormalized("hist")
            first = False
        else:
#            h.SetMaximum(1.2 * ymax)  # modify #
            h.DrawNormalized("same")

    legend.SetBorderSize(0)
    legend.Draw()
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
#    latex.DrawLatexNDC(0.17, 0.91, "Delphes Simulation")
    latex.DrawLatexNDC(0.19, 0.84, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.67, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")
    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch + ".pdf")  # modify #
    canvas.Clear()




def drawHistoSame_Single(infile, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, lepo):
    outdir = "./plots/" + PRE + "/Same/"
    try : os.makedirs(outdir)
    except : pass
    _df   = ROOT.RDataFrame(tree, infile)
    df = _df.Filter("higgs_mass>0")
    tthh = df.Filter("category==0")
    tthbb = df.Filter("category==1")
    ttbb = df.Filter("category==2")
    ttbbbb = df.Filter("category==3")
    dfs = {"tthh": tthh, "ttbbbb": ttbbbb, "ttbb": ttbb, "tthbb": tthbb}
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    ymax, color = 0, 1
    if lepo == 1 : legend = ROOT.TLegend(0.70, 0.83, 0.85, 0.68) # Upper Right #
    elif lepo == 2 : legend = ROOT.TLegend(0.38, 0.855, 0.57, 0.705) # Center #  
    else : legend = ROOT.TLegend(0.26, 0.855, 0.45, 0.705) # Upper Left #
    hist_dict = {}
    for df_name, df in dfs.items():
        h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch)
        if ymax < h.GetMaximum(): ymax = h.GetMaximum();
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetStats(0)
        if df_name == "tthh":
            h.SetFillColor(ROOT.kBlack) # 검은색으로 설정
            h.SetFillStyle(1001) # 완전히 채워진 스타일
            legend.AddEntry(h.GetValue(), " "+df_name, "f") # 면적 색칠
        else:
            if color in [5] : color += 1;
            legend.AddEntry(h.GetValue(), " "+df_name, "l")
        color += 1
        hist_dict[branch + "_" + df_name] = h    
    first = True
    for _tmp, h in hist_dict.items():
        h.SetMaximum(ymax * 1.1)
        if first:
            h.SetFillColorAlpha(ROOT.kGray, 0.8)
            h.DrawNormalized("hist")
            first = False
        else: h.DrawNormalized("same")

    legend.SetBorderSize(0)        
    legend.Draw()
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.68, 0.91, "HL-LHC (\\sqrt{s} =14TeV)")
    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch + ".pdf") # modify #
    canvas.Clear()


def drawHistoSame_Sub(indir, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, S, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Same/"
    try:
        os.makedirs(outdir)
    except:
        pass

    # 샘플 데이터프레임 로드
    tthh = ROOT.RDataFrame(tree, indir + PRE + "_tthh.root")
    tth = ROOT.RDataFrame(tree, indir + PRE + "_tth.root")
    ttbbh = ROOT.RDataFrame(tree, indir + PRE + "_ttbbh.root")
    ttzh = ROOT.RDataFrame(tree, indir + PRE + "_ttzh.root")
    ttvv = ROOT.RDataFrame(tree, indir + PRE + "_ttvv.root")
    ttbbv = ROOT.RDataFrame(tree, indir + PRE + "_ttbbv.root")
    ttbb = ROOT.RDataFrame(tree, indir + PRE + "_ttbb.root")
    ttbbbb = ROOT.RDataFrame(tree, indir + PRE + "_ttbbbb.root")
    tttt = ROOT.RDataFrame(tree, indir + PRE + "_tttt.root")#.Range(10000)

    # 데이터프레임 딕셔너리
    _dfs = {"ttHH":tthh, "ttH":tth, "ttbbH":ttbbh, "ttZH":ttzh, "ttVV":ttvv, "ttbbV":ttbbv, "ttbb":ttbb, "ttbbbb":ttbbbb, "tttt":tttt}
    dfs = {}
    for key, _df in _dfs.items():
        dfs[key] = _df.Filter("Lep_size==2").Filter("SS_OS_DL==-1").Filter("bJet_size>=5").Filter("j_ht>300") # modify! #
    dfs = {"ttHH": tthh, "ttH": tth, "ttbbH": ttbbh, "ttZH": ttzh, 
           "ttVV": ttvv, "ttbbV": ttbbv, "ttbb": ttbb, "ttbbbb": ttbbbb, "tttt": tttt}

    # 가중치 설정
    weights = {"ttHH": 1.0, "ttH": 1.0, "ttbbH": 1.0, "ttZH": 1.0, "ttVV": 1.0, "ttbbV": 1.0, "ttbb": 1.0, "ttbbbb": 1.0, "tttt": 1.0}
#    weights = {"ttHH": 1.0, "ttH": 612, "ttbbH": 2.919, "ttZH": 1.543, "ttVV": 14.62, "ttbbV": 5.011, "ttbb": 2395, "ttbbbb": 8.918, "tttt": 11.81}

    colors = {
        "ttHH": ROOT.kBlack,
        "ttH": ROOT.kGray+1, "ttbbH": ROOT.kGray + 3, "ttZH": ROOT.kGray + 5,
        "ttVV": ROOT.kGreen, "ttbbV": ROOT.kGreen + 2,
        "ttbb": ROOT.kBlue, "ttbbbb": ROOT.kCyan,
        "tttt": ROOT.kRed
    }

    ROOT.gStyle.SetPadTickX(1)  # X축 위쪽에 tick 추가
    ROOT.gStyle.SetPadTickY(1)  # Y축 오른쪽에 tick 추가

    canvas = ROOT.TCanvas("c", "c", 800, 900)
    canvas.SetLeftMargin(0.15)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.025)  # 하단 마진 제거
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.07, 1, 0.3)
    pad2.SetTopMargin(0.025)  # 상단 마진 제거
    pad2.SetBottomMargin(0.25)  # 비율 플롯의 라벨을 위한 공간 확보
    pad2.Draw()

    pad1.cd()  # 메인 플롯 그리기
    legend_position = {1: (0.72, 0.84, 0.87, 0.46), 2: (0.38, 0.855, 0.57, 0.705), 3: (0.26, 0.855, 0.45, 0.705)}
    legend = ROOT.TLegend(*legend_position.get(lepo, legend_position[3]))
    legend.SetTextSize(0.028)
    legend.SetEntrySeparation(0.20)

    hist_dict = {}
    ymax = 0
    sum_hist = None
    for df_name, df in dfs.items():
        h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch)
        if ymax < h.GetMaximum():
            ymax = h.GetMaximum()

        h.GetXaxis().SetTitle(xtitle)
        h.GetXaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle(ytitle)
        h.SetLineColor(colors[df_name])
        h.SetLineWidth(2)
        h.SetStats(0)

        if df_name == "ttHH":
            h.SetLineStyle(2)  # Dashed line style
            h.SetLineWidth(4)  # Thick line
            legend.AddEntry(h.GetValue(), " " + "t\\bar{t}HH", "f")
        else:
            legend.AddEntry(h.GetValue(), " " + df_name, "f")
            if sum_hist is None:
                sum_hist = h.Clone("weight_sum")
            else:
                sum_hist.Add(h.GetValue(), weights[df_name])  # 개별 가중치 적용
        
        hist_dict[branch + "_" + df_name] = h

    first = True
    for _tmp, h in hist_dict.items():
        if first:
            h.SetMaximum(1.4 * ymax)
            h.DrawNormalized("hist")
            first = False
        else:
            h.DrawNormalized("same")

    ratio_hist = hist_dict[branch + "_" + "ttHH"].Clone("ratio")
    sum_hist.Scale(1.0 / sum_hist.Integral())
    ratio_hist.Scale(1.0 / ratio_hist.Integral())
    ratio_hist.Divide(sum_hist)
    pad2.cd()  # 비율 플롯 그리기
    ratio_hist.GetXaxis().SetTitle(xtitle)
    ratio_hist.GetXaxis().SetTitleSize(0.11)
    ratio_hist.GetXaxis().SetLabelSize(0.1)
    ratio_hist.SetTitle("")
    ratio_hist.GetYaxis().SetTitle("Sig/Bkg")
    ratio_hist.GetYaxis().SetTitleSize(0.1)
    ratio_hist.GetYaxis().SetTitleOffset(0.25)
    ratio_hist.GetYaxis().SetLabelSize(0.05)
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMinimum(0)  # 비율 플롯의 최소값 설정
    ratio_hist.SetMaximum(2)  # 비율 플롯의 최대값 설정
    ratio_hist.Draw("ep")
    line = ROOT.TLine(ratio_hist.GetXaxis().GetXmin(), 1, ratio_hist.GetXaxis().GetXmax(), 1)  # 수평선 생성
    line.SetLineColor(ROOT.kRed)  # 선의 색상 설정 (예: 빨간색)
    line.SetLineWidth(2)  # 선의 두께 설정
    line.SetLineStyle(2)  # 선 스타일 설정 (예: 점선)
    line.Draw()  # 수평선 그리기

    pad1.cd()
    legend.SetBorderSize(0)
    legend.Draw()
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.15, 0.84, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.70, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")
    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + S + "_" + branch + ".pdf")
    canvas.Clear()

def drawHistoSame_Single_Sub(infile, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Same/"
    try:
        os.makedirs(outdir)
    except:
        pass

    df = ROOT.RDataFrame(tree, infile)

    # G1, G2, G3, G4로 필터링
    G1 = df.Filter("category == 0")
    G2 = df.Filter("category == 1")
    G3 = df.Filter("category == 2")
    G4 = df.Filter("category == 3")

    # 데이터프레임 딕셔너리
    dfs = {"G1": G1, "G2": G2, "G3": G3, "G4": G4}

    # G1도 포함하여 가중치 명시
    weights = {"G1": 1.0, "G2": 1.0, "G3": 1.0, "G4": 1.0}

    # 색상 설정
    colors = {
        "G1": ROOT.kBlack,
        "G2": ROOT.TColor.GetColorTransparent(ROOT.kBlue, 0.6),  # 50% 투명도
        "G3": ROOT.TColor.GetColorTransparent(ROOT.kRed, 0.6),    # 50% 투명도
        "G4": ROOT.TColor.GetColorTransparent(ROOT.kGreen + 2, 0.6)  # 50% 투명도
    }

    # 레전드 라벨 설정
    legs = {
        "G1": "Label for G1",
        "G2": "Label for G2",
        "G3": "Label for G3",
        "G4": "Label for G4"
    }

    ROOT.gStyle.SetPadTickX(1)  # X축 위쪽에 tick 추가
    ROOT.gStyle.SetPadTickY(1)  # Y축 오른쪽에 tick 추가

    canvas = ROOT.TCanvas("c", "c", 800, 900)
    canvas.SetLeftMargin(0.15)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.025)
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.07, 1, 0.3)
    pad2.SetTopMargin(0.025)
    pad2.SetBottomMargin(0.25)
    pad2.Draw()
    pad1.cd()

    ymax, color = 0, 1

    # 범례 위치 설정
    legend_position = {1: (0.72, 0.84, 0.87, 0.46), 2: (0.38, 0.855, 0.57, 0.705), 3: (0.26, 0.855, 0.45, 0.705)}
    legend = ROOT.TLegend(*legend_position.get(lepo, legend_position[3]))
    legend.SetTextSize(0.028)
    legend.SetEntrySeparation(0.20)

    hist_dict = {}
    sum_hist = None

    # 히스토그램 생성 및 꾸미기
    for df_name, df in dfs.items():
        h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch)
        if ymax < h.GetMaximum():
            ymax = h.GetMaximum()

        h.GetXaxis().SetTitle(xtitle)
        h.GetXaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle(ytitle)
        h.SetLineColor(colors[df_name])
        h.SetLineWidth(2)
        h.SetStats(0)

        # 범례에 라벨 추가
        legend.AddEntry(h.GetValue(), legs[df_name], "f")

        if df_name == "G1":
            h.SetLineStyle(2)  # Dashed line style
            h.SetLineWidth(4)  # Thick line
        else:
            if sum_hist is None:
                sum_hist = h.Clone("weight_sum")
            else:
                sum_hist.Add(h.GetValue(), weights[df_name])

        hist_dict[branch + "_" + df_name] = h

    first = True
    for _tmp, h in hist_dict.items():
        h.SetMaximum(1.4 * ymax)
        if first:
            h.DrawNormalized("hist")
            first = False
        else:
            h.DrawNormalized("same")

    # 비율 히스토그램 생성 및 그리기
    ratio_hist = hist_dict[branch + "_" + "G1"].Clone("ratio")
    sum_hist.Scale(1.0 / sum_hist.Integral())
    ratio_hist.Scale(1.0 / ratio_hist.Integral())
    ratio_hist.Divide(sum_hist)

    pad2.cd()
    ratio_hist.GetXaxis().SetTitle(xtitle)
    ratio_hist.GetXaxis().SetTitleSize(0.11)
    ratio_hist.GetXaxis().SetLabelSize(0.1)
    ratio_hist.SetTitle("")
    ratio_hist.GetYaxis().SetTitle("Sig/Bkg")
    ratio_hist.GetYaxis().SetTitleSize(0.1)
    ratio_hist.GetYaxis().SetTitleOffset(0.25)
    ratio_hist.GetYaxis().SetLabelSize(0.05)
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMinimum(0)
    ratio_hist.SetMaximum(2)
    ratio_hist.Draw("ep")

    line = ROOT.TLine(ratio_hist.GetXaxis().GetXmin(), 1, ratio_hist.GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()

    pad1.cd()
    legend.SetBorderSize(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.15, 0.84, "MadGraph5 Simulation")
    latex.DrawLatexNDC(0.70, 0.91, "HL-LHC, \\sqrt{s} = 14 TeV")

    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch + ".pdf")
    canvas.Clear()


def drawHistoStack_Single(infile, tree, title, xtitle, ytitle, branch, nbin, xmin, xmax, PRE, lepo):
    title = ""
    outdir = "./plots/" + PRE + "/Stack/"
    try:
        os.makedirs(outdir)
    except:
        pass

    # 데이터를 읽어들임
    df = ROOT.RDataFrame(tree, infile)

    # Scale Factor s.t. ttbb = 5368
#    weights = {"G1": 19, "G2": 2093, "G3": 919, "G4": 769} #SSDL
    weights = {"G1": 13, "G2": 1000, "G3": 4500, "G4": 413} #OSDL


    G1 = df.Filter("category==3")
    h_G1 = G1.Histo1D(("G1_hist", "G1_hist", nbin, xmin, xmax), branch)
    total_G1 = h_G1.Integral()
    SF = 1 / total_G1
    scaled_weights = {key: value * SF for key, value in weights.items()}

    # 각 카테고리에 새로운 가중치 적용
    G1 = df.Filter("category==0").Define("weight", f"{scaled_weights['G1']*300}")
    G2 = df.Filter("category==1").Define("weight", f"{scaled_weights['G2']}")
    G3 = df.Filter("category==2").Define("weight", f"{scaled_weights['G3']}")
    G4 = df.Filter("category==3").Define("weight", f"{scaled_weights['G4']}")

    dfs = {"G1": G1, "G2": G2, "G3": G3, "G4": G4}
    colors = {
        "G1": ROOT.kBlack,
        "G2": ROOT.TColor.GetColorTransparent(ROOT.kBlue, 0.6),  # 50% 투명도
        "G3": ROOT.TColor.GetColorTransparent(ROOT.kRed, 0.6),    # 50% 투명도
        "G4": ROOT.TColor.GetColorTransparent(ROOT.kGreen + 2, 0.6)  # 50% 투명도
    }

    canvas = ROOT.TCanvas("c", "c", 800, 800)
    canvas.SetLeftMargin(0.15)  # Adjust the left margin to make space for the y-axis label
#    canvas.SetLogy()  # Use logarithmic scale for the y-axis
    ymax, color = 0, 1
    if lepo == 1:
        legend = ROOT.TLegend(0.65, 0.83, 0.85, 0.68)  # Upper Right #
    elif lepo == 2:
        legend = ROOT.TLegend(0.38, 0.855, 0.57, 0.705)  # Center #
    else:
        legend = ROOT.TLegend(0.26, 0.855, 0.45, 0.705)  # Upper Left #
    legs = {"G1" : "G1", "G2" : "G2", "G3" : "G3", "G4" : "G4"}
    legend.SetTextSize(0.035)  # Adjust this value to make the legend larger or smaller
    legend.SetEntrySeparation(0.02)  # Adjust the separation between legend entries

    hist_dict = {}
    stack = ROOT.THStack("hs", title)
    for df_name, df in dfs.items():
        h = df.Histo1D(ROOT.RDF.TH1DModel(title, title, nbin, xmin, xmax), branch, "weight")
        if ymax < h.GetMaximum():
            ymax = h.GetMaximum()
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.SetLineColor(colors[df_name])

#        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetStats(0)
        if df_name != "G1":
            if color in [5]:
                color += 1
            h.SetFillColor(colors[df_name])
            legend.AddEntry(h.GetValue(), " " + legs[df_name] , "f")  # Filled area
        else:
            h.SetLineStyle(2)  # 점선 스타일 설정
            h.SetLineColor(ROOT.kBlack)  # 검은색으로 설정
            h.SetLineWidth(3)
            legend.AddEntry(h.GetValue(), " " + legs[df_name] + " (x300)", "l")  # Line
        color += 1
        hist_dict[branch + "_" + df_name] = h
        if df_name != "G1":
            stack.Add(h.GetValue())

    stack.Draw("hist")
    stack.GetXaxis().SetTitle(xtitle)
    stack.GetYaxis().SetTitle(ytitle)
    stack.SetMaximum(ymax * 1.2)

    # "tthh" 히스토그램을 점선으로 추가
    hist_dict[branch + "_G1"].Draw("hist same")

    legend.SetBorderSize(0)
    legend.Draw()
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    latex.DrawLatexNDC(0.15, 0.91, "Delphes")
    latex.DrawLatexNDC(0.55, 0.91, "HL-LHC (#sqrt{s} =14TeV, L = 3000fb^{-1})")

    title = title.replace(" ", "_")
    canvas.Print(outdir + PRE + "_" + branch + ".pdf")
    canvas.Clear()
    print(outdir + PRE + "_" + branch + ".pdf")

