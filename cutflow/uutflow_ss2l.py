import ROOT
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
from array import array

# Input input_open
input_cf = [
    "SS_OS_DL", "j_ht", "Lep_size"
]

input_dnn = [
     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m",

     "bJet_size",
     "b1b2_dr", "b1b3_dr",
     "b2b3_dr",

     "Lep1_pt", "Lep1_eta", "Lep1_phi", 
     "Lep2_pt", "Lep2_eta", "Lep2_phi", 
     "MET_E",

   # Kinematic
#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_m", "bb_twist",
]

input_open = input_cf + input_dnn

# Criteria
PRE = "CF_INC"
print("PRE")
Tree = "Delphes"

# Input files
tthh = "samples1/" + PRE + "_tthh.root"
tth = "samples1/" + PRE + "_tth.root"
ttbbh = "samples1/" + PRE + "_ttbbh.root"
ttzh = "samples1/" + PRE + "_ttzh.root"
ttvv = "samples1/" + PRE + "_ttvv.root"
ttbbv = "samples1/" + PRE + "_ttbbv.root"
ttbb = "samples1/" + PRE + "_ttbb.root"
ttbbbb = "samples1/" + PRE + "_ttbbbb.root"
tttt = "samples1/" + PRE + "_tttt.root"

# Load data using uproot
tthh = uproot.open(tthh)[Tree].arrays(input_open, library="pd")
tth = uproot.open(tth)[Tree].arrays(input_open, library="pd")
ttbbh = uproot.open(ttbbh)[Tree].arrays(input_open, library="pd")
ttzh = uproot.open(ttzh)[Tree].arrays(input_open, library="pd")
ttvv = uproot.open(ttvv)[Tree].arrays(input_open, library="pd")
ttbbv = uproot.open(ttbbv)[Tree].arrays(input_open, library="pd")
ttbb = uproot.open(ttbb)[Tree].arrays(input_open, library="pd")
ttbbbb = uproot.open(ttbbbb)[Tree].arrays(input_open, library="pd")
tttt = uproot.open(tttt)[Tree].arrays(input_open, library="pd")

# Luminosity [fb^-1]
L = 3000
BR = 1.0  # Inclusive

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

# Acceptance Function with Basic event selection
def Acceptance(df, df_name):
    Accept = []
    S0 = len(df)
    df = df[(df['Lep_size'] == 2) & (df['SS_OS_DL'] == 1)]
    S1 = len(df)
    df = df[df['bJet_size'] >= 3]
    S2 = len(df)
    df = df[df['j_ht'] > 300]
    S3 = len(df)
#    DNN = float(input(f"{df_name} DNN Efficiency SSDL = "))
#    S4 = S3 * DNN
    Accept.extend([S0, S1, S2, S3])
#    print(Accept)
    return Accept, df

# Calculate Acceptance
print("________Basic selection acceptance________")
tthh_acc, tthh = Acceptance(tthh, "tthh")
tth_acc, tth = Acceptance(tth, "tth")
ttbbh_acc, ttbbh = Acceptance(ttbbh, "ttbbh")
ttzh_acc, ttzh = Acceptance(ttzh, "ttzh")
ttvv_acc, ttvv = Acceptance(ttvv, "ttvv")
ttbbv_acc, ttbbv = Acceptance(ttbbv, "ttbbv")
ttbb_acc, ttbb = Acceptance(ttbb, "ttbb")
ttbbbb_acc, ttbbbb = Acceptance(ttbbbb, "ttbbbb")
tttt_acc, tttt = Acceptance(tttt, "tttt")

# Import pre-trained DNN model
dnn_dir = "/home/stiger97/github/tthh_full_14TeV/DNN_result/SSDL_SEMI_0817/LetsFind_tthh/bCat_higgs5_2Mat/best_model.h5"
dnn_model = tf.keras.models.load_model(dnn_dir)
dnn_model.summary()

# label and concat the filtered data
tthh["category"]   = 0 ; tthh["process"]   = 0
tth["category"]    = 1 ; tth["process"]    = 1
ttzh["category"]   = 1 ; ttzh["process"]   = 2
ttbbh["category"]  = 1 ; ttbbh["process"]  = 3
ttvv["category"]   = 2 ; ttvv["process"]   = 4
ttbbv["category"]  = 2 ; ttbbv["process"]  = 5
ttbb["category"]   = 2 ; ttbb["process"]   = 6
ttbbbb["category"] = 2 ; ttbbbb["process"] = 7
tttt["category"]   = 3 ; tttt["process"]   = 8
df_test = pd.concat([tthh, tth, ttzh, ttbbh, ttvv, ttbbv, ttbb, ttbbbb, tttt]).sample(frac=1).reset_index(drop=True)
x_test = np.array(df_test.filter(items = input_dnn))
y_test = np.array(df_test.filter(items = ["category"]))
pred_test = dnn_model.predict(x_test); print("pred_test :", pred_test); pred_test_arg = np.argmax(pred_test, axis=1)
print("Is it similar?")
print("Answer :     " , y_test.T)
print("Prediction : " , pred_test_arg)

results_test = dnn_model.evaluate(x_test, y_test)
loss_test = results_test[0]
acc_test = results_test[1]
print(f"Test accuracy: {acc_test * 100:.2f}%")

df_test["G1"] = np.array(pred_test.T[0])
df_test["G2"] = np.array(pred_test.T[1])
df_test["G3"] = np.array(pred_test.T[2])
df_test["G4"] = np.array(pred_test.T[3])
df_test["DNN"] = df_test["G1"]/(df_test["G2"]+df_test["G3"]+df_test["G4"])

print("________ACCEPTANCE________")
def Acceptance2(df_test, acc, process):
    S3 = len(df_test)
    df = df_test[ (df_test["process"] == process) & (df_test["DNN"]>0.25) ]
    S4 = len(df)
    acc.append(S4)
    print(acc)
    return acc, df
 
tthh_acc, tthh = Acceptance2(df_test, tthh_acc, 0)
tth_acc, tth = Acceptance2(df_test, tth_acc, 1)
ttbbh_acc, ttbbh = Acceptance2(df_test, ttbbh_acc, 2)
ttzh_acc, ttzh = Acceptance2(df_test, ttzh_acc, 3)
ttvv_acc, ttvv = Acceptance2(df_test, ttvv_acc, 4)
ttbbv_acc, ttbbv = Acceptance2(df_test, ttbbv_acc, 5)
ttbb_acc, ttbb = Acceptance2(df_test, ttbb_acc, 6)
ttbbbb_acc, ttbbbb = Acceptance2(df_test, ttbbbb_acc, 7)
tttt_acc, tttt = Acceptance2(df_test, tttt_acc, 8)

###################################################
#                  Write  TTRee                   #
###################################################
# Prepare dataframe to be written
df_score = pd.DataFrame()
df_score["category"] = y_test.T[0]
df_score["process"] = df_test["process"]
df_score["G1"] = np.array(pred_test.T[0])
df_score["G2"] = np.array(pred_test.T[1])
df_score["G3"] = np.array(pred_test.T[2])
df_score["G4"] = np.array(pred_test.T[3])
df_score["DNN"] = df_score["G1"]/(df_score["G2"]+df_score["G3"]+df_score["G4"])

f = ROOT.TFile("SSDL_score.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")

# Empty array
category_array = array('f', [0])
process_array = array('f', [0])
G1_array = array('f', [0])
G2_array = array('f', [0])
G3_array = array('f', [0])
G4_array = array('f', [0])
DNN_array = array('f', [0])

# Attatch array on TTree
tree.Branch('category', category_array, 'category/F')
tree.Branch('process', process_array, 'process/F')
tree.Branch('G1', G1_array, 'G1/F')
tree.Branch('G2', G2_array, 'G2/F')
tree.Branch('G3', G3_array, 'G3/F')
tree.Branch('G4', G4_array, 'G4/F')
tree.Branch('DNN', DNN_array, 'DNN/F')

# DataFrame의 데이터를 TTree에 채우기
for _, row in df_score.iterrows():
    category_array[0] = row['category']
    process_array[0] = row['process']
    G1_array[0] = row['G1']
    G2_array[0] = row['G2']
    G3_array[0] = row['G3']
    G4_array[0] = row['G4']
    DNN_array[0] = row['DNN']
    tree.Fill()

# 파일 저장 및 종료
f.Write()
f.Close()
#print(outdir)
print("DNN score written :")
print("SSDL_dnn_score.root")
#infile = "./DNN_result/SEMI_SSDL/LetsFind_tthh/bCat_higgs5_2Mat//dnn2_result.root"
#drawHistoStack_Single(infile, tree, "DNN Discriminator", "DNN", "Normalized Entries", "DNN", 20, 0, 1, "SEMI_SSDL", lepo=1)


# Create Dictionary for Acceptance
Acc = {
    "tthh": tthh_acc,
    "tth": tth_acc,
    "ttbbh": ttbbh_acc,
    "ttzh": ttzh_acc,
    "ttvv": ttvv_acc,
    "ttbbv": ttbbv_acc,
    "ttbb": ttbb_acc,
    "ttbbbb": ttbbbb_acc,
    "tttt": tttt_acc
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
print("Sig :", significance)

# Save to Excel with multiple sheets (optional)
with pd.ExcelWriter("output_acceptance_cutflow.xlsx") as writer:
    df_cutflow.to_excel(writer, sheet_name="Cutflow")
    df_acceptance.to_excel(writer, sheet_name="Acceptance")

print("Data has been written to output_acceptance_cutflow.xlsx with two sheets: Cutflow and Acceptance")

