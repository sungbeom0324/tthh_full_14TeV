# at your wd, python cutflow/cutflow.py
import ROOT
import numpy as np

# Criteria
PRE = "inc"
Tree = "Delphes"
print(PRE)

# Input
tthh = "skimmed/" + PRE + "_tthh.root"
tth = "skimmed/" + PRE + "_tth.root"
ttbbh = "skimmed/" + PRE + "_ttbbh.root"
ttzh = "skimmed/" + PRE + "_ttzh.root"
ttvv = "skimmed/" + PRE + "_ttvv.root"
ttbbv = "skimmed/" + PRE + "_ttbbv.root"
ttbb = "skimmed/" + PRE + "_ttbb.root"
ttbbbb = "skimmed/" + PRE + "_ttbbbb.root"
tttt = "skimmed/" + PRE + "_tttt.root"

TreeName = "Delphes"

# Luminosity [fb^-1]
L = 3000
BR = 1.0
#BR = 0.457 # Hadronic
#BR = 0.438 # Semileptonic
#BR = 0.105 # Dileptonic

# [HLLHC : Inclusive, fb]
x_tthh  = 0.948 * L * BR
x_tth   = 612 * L * BR
x_ttbbh = 2.919 * L * BR
x_ttzh  = 1.543 * L * BR
x_ttvv  = 14.62 * L * BR
x_ttbbv = 5.011 * L * BR
x_ttbb =  2661 * 0.9027 * L * BR 
x_ttbbbb =  8.918 * L * BR
x_tttt = 11.81 * L

#
crossx = {"x_tthh" : x_tthh, "x_tth" : x_tth, "x_ttbbh" : x_ttbbh, "x_ttzh" : x_ttzh, "x_ttvv" : x_ttvv, "x_ttbbv" : x_ttbbv, "x_ttbb" : x_ttbb, "x_ttbbbb" : x_ttbbbb, "x_tttt" : x_tttt}
for key, val in crossx.items():
    print(key+" : "+str(round(val/3000.,2)))

# RDF
tthh = ROOT.RDataFrame(Tree, tthh)
tth = ROOT.RDataFrame(Tree, tth)
ttbbh = ROOT.RDataFrame(Tree, ttbbh)
ttzh = ROOT.RDataFrame(Tree, ttzh)
ttvv = ROOT.RDataFrame(Tree, ttvv)
ttbbv = ROOT.RDataFrame(Tree, ttbbv)
ttbb = ROOT.RDataFrame(Tree, ttbb)
ttbbbb = ROOT.RDataFrame(Tree, ttbbbb)
tttt = ROOT.RDataFrame(Tree, tttt)

print("Calculating Acceptance and Cutflow")

# MODIFY!! E.S.
def Acceptance(df, df_name):
    Accept = []
    S0 = float(df.Count().GetValue())
    df = df.Filter("Lep_size == 2 && SS_OS_DL == 1"); S1 = float(df.Count().GetValue())
    df = df.Filter("bJet_size >= 3"); S2 = float(df.Count().GetValue())
    df = df.Filter("j_ht>300"); S3 = float(df.Count().GetValue())
    DNN = float(input(f"{df_name} DNN Efficiency SSDL = "))
    S4 = S3 * DNN
    Accept.extend([S0,S1,S2,S3,S4])
    print(Accept)
    return Accept

print("________ACCEPTANCE________")
tthh = Acceptance(tthh, "tthh")
tth = Acceptance(tth, "tth")
ttbbh = Acceptance(ttbbh, "ttbbh")
ttzh = Acceptance(ttzh, "ttzh")
ttvv = Acceptance(ttvv, "ttvv")
ttbbv = Acceptance(ttbbv, "ttbbv")
ttbb = Acceptance(ttbb, "ttbb")
ttbbbb = Acceptance(ttbbbb, "ttbbbb")
tttt = Acceptance(tttt, "tttt")

Acc = {     # Normalize Acceptance to the cross section * L.
    "tthh" : [tthh, x_tthh/tthh[0]], 
    "tth" : [tth, x_tth/tth[0]],
    "ttbbh" : [ttbbh, x_ttbbh/ttbbh[0]],
    "ttzh" : [ttzh, x_ttzh/ttzh[0]],
    "ttvv" : [ttvv, x_ttvv/ttvv[0]],
    "ttbbv" : [ttbbv, x_ttbbv/ttbbv[0]],
    "ttbb" : [ttbb, x_ttbb/ttbb[0]],
    "ttbbbb" : [ttbbbb, x_ttbbbb/ttbbbb[0]],
    "tttt" : [tttt, x_tttt/tttt[0]],
}

# Cutflow
def Cutflow(Acc):
    for key, value in Acc.items():
        value[0] = [element * value[1] for element in value[0]] # value[0] = [S0,S1,S2,S3], value[1] = Weight.
        rounded = [round(num, 2) for num in value[0]]
        print(key, rounded)
    return Acc

print("__________CUTFLOW__________")        
CF = Cutflow(Acc)

print(" ")
print("________SIGNIFICANCE________")

# Significance # Modify! # 
for i in range(0,5): # 5 Cuts + 1 No cut
    print("Significance of ES :", i)
    Significance = CF["tthh"][0][i]/np.sqrt(CF["tth"][0][i] + CF["ttbbh"][0][i] + CF["ttzh"][0][i] + CF["ttvv"][0][i] + CF["ttbbv"][0][i] +CF["ttbb"][0][i] + CF["ttbbbb"][0][i] + CF["tttt"][0][i])
    print("Significance: {:.2f}".format(Significance))
    print("----Done----")

