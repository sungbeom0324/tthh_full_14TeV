#ng histograms and save in a pdf.
# e.g. tthh
import ROOT
import numpy as np

# Criteria
PRE = "CF_0721"
Tree = "Delphes"
print(PRE)

# Input
tthh_inc = "samples1/" + PRE + "_tthh_inc.root"
tthh_had = "samples1/" + PRE + "_tthh_had.root"
tthh_semi = "samples1/" + PRE + "_tthh_semi.root"
tthh_di = "samples1/" + PRE + "_tthh_di.root"
TreeName = "Delphes"

# Luminosity [fb^-1]
L = 3000

# [HLLHC : Inclusive, fb]
x_tthh_inc  = 0.948 * L
x_tthh_had  = 0.948 * L * 0.46
x_tthh_semi  = 0.948 * L * 0.44
x_tthh_di  = 0.948 * L * 0.10


# RDF
tthh_inc = ROOT.RDataFrame(Tree, tthh_inc)
tthh_had = ROOT.RDataFrame(Tree, tthh_had)
tthh_semi = ROOT.RDataFrame(Tree, tthh_semi)
tthh_di = ROOT.RDataFrame(Tree, tthh_di)

print("Calculating Acceptance and Cutflow")

# MODIFY!! E.S.
def Acceptance(df, df_name):
    Accept = []
    S0 = float(df.Count().GetValue())
    df = df.Filter("Lep_size == 2"); S1 = float(df.Count().GetValue())
    df = df.Filter("SS_OS_DL == -1"); S2 = float(df.Count().GetValue())
    df = df.Filter("bJet_size >= 5"); S3 = float(df.Count().GetValue())
    DNN = float(input(f"{df_name} DNN Efficiency OSDL = "))
    S4 = S3 * DNN
    Accept.extend([S0,S1,S2,S3,S4])
    print(Accept)
    return Accept

print("________ACCEPTANCE________")
tthh_inc = Acceptance(tthh_inc, "tthh_inc")
tthh_had = Acceptance(tthh_had, "tthh_had")
tthh_semi = Acceptance(tthh_semi, "tthh_semi")
tthh_di = Acceptance(tthh_di, "tthh_di")

Acc = {     # Normalize Acceptance to the cross section * L.
    "tthh_inc" : [tthh_inc, x_tthh_inc/tthh_inc[0]],
    "tthh_had" : [tthh_had, x_tthh_had/tthh_had[0]], 
    "tthh_semi" : [tthh_semi, x_tthh_semi/tthh_semi[0]], 
    "tthh_di" : [tthh_di, x_tthh_di/tthh_di[0]], 
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

'''
# Significance # Modify! # 
for i in range(0,5): # 5 Cuts + 1 No cut
    print("Significance of ES :", i)
    Significance = CF["tthh"][0][i]/np.sqrt(CF["tth"][0][i] + CF["ttbbh"][0][i] + CF["ttzh"][0][i] + CF["ttvv"][0][i] + CF["ttbbv"][0][i] +CF["ttbb"][0][i] + CF["ttbbbb"][0][i] + CF["tttt"][0][i])
    print("Significance: {:.2f}".format(Significance))
    print("----Done----")
'''
