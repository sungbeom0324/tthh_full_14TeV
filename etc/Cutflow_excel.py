import ROOT
import numpy as np
import pandas as pd

# Criteria
PRE = "CF_DI"
Tree = "Delphes"
print(PRE)

# Input
tthh = "samples1/" + PRE + "_tthh.root"
tth = "samples1/" + PRE + "_tth.root"
ttbbh = "samples1/" + PRE + "_ttbbh.root"
ttzh = "samples1/" + PRE + "_ttzh.root"
ttvv = "samples1/" + PRE + "_ttvv.root"
ttbbv = "samples1/" + PRE + "_ttbbv.root"
ttbb = "samples1/" + PRE + "_ttbb.root"
ttbbbb = "samples1/" + PRE + "_ttbbbb.root"
tttt = "samples1/" + PRE + "_tttt.root"

TreeName = "Delphes"

# Luminosity [fb^-1]
L = 3000
#BR = 1.0 # Inclusive
#BR = 0.457 # Hadronic
#BR = 0.438 # Semileptonic
BR = 0.105 # Dileptonic

# [HLLHC : Inclusive, fb]
crossx = {
    "x_tthh": 0.948 * L * BR,
    "x_tth": 612 * L * BR,
    "x_ttbbh": 2.919 * L * BR,
    "x_ttzh": 1.543 * L * BR,
    "x_ttvv": 14.62 * L * BR,
    "x_ttbbv": 5.011 * L * BR,
    "x_ttbb": 2661 * 0.9027 * L * BR,
    "x_ttbbbb": 8.918 * L * BR,
    "x_tttt": 11.81 * L
}

for key, val in crossx.items():
    print(key + " : " + str(round(val / 3000., 2)))

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

# Acceptance Function
def Acceptance(df, df_name):
    Accept = []
    S0 = float(df.Count().GetValue())
    df = df.Filter("Lep_size == 2 && SS_OS_DL == -1")
    S1 = float(df.Count().GetValue())
    df = df.Filter("bJet_size >= 5")
    S2 = float(df.Count().GetValue())
    df = df.Filter("j_ht > 300")
    S3 = float(df.Count().GetValue())
    DNN = float(input(f"{df_name} DNN Efficiency SSDL = "))
    S4 = S3 * DNN
    Accept.extend([S0, S1, S2, S3, S4])
    print(Accept)
    return Accept

# Calculate Acceptance
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

# Create Dictionary for Acceptance
Acc = {  
    "tthh": tthh,
    "tth": tth,
    "ttbbh": ttbbh,
    "ttzh": ttzh,
    "ttvv": ttvv,
    "ttbbv": ttbbv,
    "ttbb": ttbb,
    "ttbbbb": ttbbbb,
    "tttt": tttt
}

# Create DataFrame for Acceptance
df_acceptance = pd.DataFrame.from_dict(Acc, orient='index', columns=['S0', 'S1', 'S2', 'S3', 'S4'])

# Normalize Acceptance and Prepare Cutflow
def Cutflow(Acc):
    cutflow_dict = {}
    for key, value in Acc.items():
        weight = crossx["x_" + key] / value[0]
        value = [element * weight for element in value]
        rounded = [round(num, 2) for num in value]
        cutflow_dict[key] = rounded
        print(key, rounded)
    return cutflow_dict

print("__________CUTFLOW__________")
CF = Cutflow(Acc)

# Create DataFrame for Cutflow
df_cutflow = pd.DataFrame.from_dict(CF, orient='index', columns=['S0', 'S1', 'S2', 'S3', 'S4'])

# Calculate Significance
significance = []
for i in range(5):  # 5 steps in the cutflow
    sig = CF["tthh"][i] / np.sqrt(
        CF["tth"][i] + CF["ttbbh"][i] + CF["ttzh"][i] + CF["ttvv"][i] + CF["ttbbv"][i] + CF["ttbb"][i] + CF["ttbbbb"][i] + CF["tttt"][i]
    )
    significance.append(round(sig, 2))

# Append Significance to Cutflow DataFrame
df_cutflow.loc["Significance"] = significance

# Save to Excel with multiple sheets
with pd.ExcelWriter("output_acceptance_cutflow.xlsx") as writer:
    df_cutflow.to_excel(writer, sheet_name="Cutflow")
    df_acceptance.to_excel(writer, sheet_name="Acceptance")

print("Data has been written to output_acceptance_cutflow.xlsx with two sheets: Cutflow and Acceptance")

