import ROOT
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
from array import array
from utils.var_functions import *
import json

# Input input_dnn
with open('./dnn/dnn_input.json', 'r') as file:
    data = json.load(file)
input_0 = ["Lep_size", "SS_OS_DL", "j_ht"] # event selection
input_1 = data["input_1"] # bfh, bft.
input_2 = data["input_2"]
input_3 = data["input_3"]

input_open = input_0 + input_1 + input_2
input_dnn = input_1 + input_2


# Criteria
PRE = "inc"
print("PRE")
Tree = "Delphes"

# Input files
tthh = "skimmed/" + PRE + "_tthh.root"
tth = "skimmed/" + PRE + "_tth.root"
ttbbh = "skimmed/" + PRE + "_ttbbh.root"
ttzh = "skimmed/" + PRE + "_ttzh.root"
ttvv = "skimmed/" + PRE + "_ttvv.root"
ttbbv = "skimmed/" + PRE + "_ttbbv.root"
ttbb = "skimmed/" + PRE + "_ttbb.root"
ttbbbb = "skimmed/" + PRE + "_ttbbbb.root"
tttt = "skimmed/" + PRE + "_tttt.root"

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
    df = df[(df['Lep_size'] == 2) & (df['SS_OS_DL'] == -1)]
    S1 = len(df)
    df = df[df['bJet_size'] >= 5]
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
df_total = pd.concat([tthh, tth, ttzh, ttbbh, ttvv, ttbbv, ttbb, ttbbbb, tttt]).sample(frac=1).reset_index(drop=True)
x_bCat  = np.array(df_total.filter(items = input_1))

###################################################
#               bJet Classification               #
###################################################
# Load b-tagger model.
bfh_dir = "dnn_result/bfh/best_model.h5"
bft_dir = "dnn_result/bft/best_model.h5"
bfh_model = tf.keras.models.load_model(bfh_dir)
bft_model = tf.keras.models.load_model(bft_dir)

# Predict bJet origin.
_pred_bfh = bfh_model.predict(x_bCat); print("bJet from Higgs score : ", _pred_bfh)
_pred_bft = bft_model.predict(x_bCat); print("bJet from Top quark score : ", _pred_bft)
pred_bfh = np.argmax(_pred_bfh, axis=1); print("bJet from Higgs : ", pred_bfh)

# Define new variables
df_total["pred_bfh"] = pred_bfh
for i in range(10):
    column_name = f"2bfh_{i + 1}" # 동적으로 변수명 할당.
    df_total[column_name] = _pred_bfh[:, i]
for i in range(5):
    column_name = f"bft_{i + 1}"
    df_total[column_name] = _pred_bft[:, i]

df_total["higgs_mass_list"] = df_total.apply(higgs_5_2, axis = 1) # MODIFY #
df_total["higgs_mass"] = df_total["higgs_mass_list"].apply(lambda x: x[0])
df_total["higgs_mass_sub"] = df_total["higgs_mass_list"].apply(lambda x: x[1])
df_total["higgs_mass_sum"] = df_total["higgs_mass_list"].apply(lambda x: x[2])
df_total["X_higgs"] = df_total.apply(X_higgs, axis = 1) # MODIFY #

df_total["bfh_Vars"] = df_total.apply(bfh_Vars, axis = 1)
df_total["bfh_dr"] = df_total["bfh_Vars"].apply(lambda x: x[0])
df_total["bfh_Ht"] = df_total["bfh_Vars"].apply(lambda x: x[1])
df_total["bfh_dEta"] = df_total["bfh_Vars"].apply(lambda x: x[2])
df_total["bfh_dPhi"] = df_total["bfh_Vars"].apply(lambda x: x[3])
df_total["bfh_mbmb"] = df_total["bfh_Vars"].apply(lambda x: x[4])

######################################################
#                Event classification                #
######################################################
# Import pre-trained model
dnn_dir = "dnn_result/os2l/best_model.h5"
dnn_model = tf.keras.models.load_model(dnn_dir)
dnn_model.summary()
    
x_test = np.array(df_total.filter(items = input_dnn + input_3))
y_test = np.array(df_total.filter(items = ["category"]))
pred_test = dnn_model.predict(x_test); print("pred_test :", pred_test); pred_test_arg = np.argmax(pred_test, axis=1)
print("Is it similar?")
print("Answer :     " , y_test.T)
print("Prediction : " , pred_test_arg)

results_test = dnn_model.evaluate(x_test, y_test)
loss_test = results_test[0]
acc_test = results_test[1]
print(f"Test accuracy: {acc_test * 100:.2f}%")

df_total["G1"] = np.array(pred_test.T[0])
df_total["G2"] = np.array(pred_test.T[1])
df_total["G3"] = np.array(pred_test.T[2])
df_total["G4"] = np.array(pred_test.T[3])
df_total["DNN"] = df_total["G1"]/(df_total["G2"]+df_total["G3"]+df_total["G4"])

print("________ACCEPTANCE________")
def Acceptance2(df_total, acc, process):
    S3 = len(df_total)
    df = df_total[ (df_total["process"] == process) & (df_total["DNN"]>0.3) ]
    S4 = len(df)
    acc.append(S4)
    print(acc)
    return acc, df
 
tthh_acc, tthh = Acceptance2(df_total, tthh_acc, 0)
tth_acc, tth = Acceptance2(df_total, tth_acc, 1)
ttbbh_acc, ttbbh = Acceptance2(df_total, ttbbh_acc, 2)
ttzh_acc, ttzh = Acceptance2(df_total, ttzh_acc, 3)
ttvv_acc, ttvv = Acceptance2(df_total, ttvv_acc, 4)
ttbbv_acc, ttbbv = Acceptance2(df_total, ttbbv_acc, 5)
ttbb_acc, ttbb = Acceptance2(df_total, ttbb_acc, 6)
ttbbbb_acc, ttbbbb = Acceptance2(df_total, ttbbbb_acc, 7)
tttt_acc, tttt = Acceptance2(df_total, tttt_acc, 8)

###################################################
#                  Write  TTRee                   #
###################################################
# Prepare dataframe to be written
df_score = pd.DataFrame()
df_score["category"] = y_test.T[0]
df_score["process"] = df_total["process"]
df_score["G1"] = np.array(pred_test.T[0])
df_score["G2"] = np.array(pred_test.T[1])
df_score["G3"] = np.array(pred_test.T[2])
df_score["G4"] = np.array(pred_test.T[3])
df_score["DNN"] = df_score["G1"]/(df_score["G2"]+df_score["G3"]+df_score["G4"])

f = ROOT.TFile("os2l_score.root", "RECREATE")
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
print("os2l_score.root")


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
with pd.ExcelWriter("cutflow_os2l.xlsx") as writer:
    df_cutflow.to_excel(writer, sheet_name="Cutflow")
    df_acceptance.to_excel(writer, sheet_name="Acceptance")

print("Data has been written to cutflow_os2l.xlsx with two sheets: Cutflow and Acceptance")

