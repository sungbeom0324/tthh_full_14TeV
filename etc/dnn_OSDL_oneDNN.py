# ssh gpu-0-X
# conda activate py36
import os
import sys
import time
#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import tensorflow as tf
from tensorflow.python.eager import backprop
import shap
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
from matplotlib.backends.backend_pdf import PdfPages
from utils.var_functions import *

import ROOT
from array import array

start_time = time.time()

###################################################
#                     I/O                         #
###################################################
indir = "./samples1/"; PRE = "OSDL_DI_0817"
outdir = "./DNN_result/" + PRE + "/LetsFind_tthh/bCat_higgs5_2Mat/"    # modify #
os.makedirs(outdir, exist_ok=True)
process_names = ["G1", "G2", "G3", "G4"]

input_1 = [ # Used for analytic reason, like bjet categorization. 
     "bJet_size",

     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m",
     "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m",

     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr",
     "b2b3_dr", "b2b4_dr", "b2b5_dr",
     "b3b4_dr", "b3b5_dr",
     "b4b5_dr",
            ]

input_2 = [
     "Lep1_pt", "Lep1_eta", "Lep1_phi",
     "Lep2_pt", "Lep2_eta", "Lep2_phi",
     "MET_E",

    # Kinematic
#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_m", "bb_twist",
           ]
 
dnn_input = input_1 + input_2

###################################################
#                 PreProcessing                   #
###################################################
# Load samples and label
df_tthh   = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(dnn_input,library="pd")
df_tth   = uproot.open(indir+PRE+"_tth.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttzh   = uproot.open(indir+PRE+"_ttzh.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttbbh   = uproot.open(indir+PRE+"_ttbbh.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttvv   = uproot.open(indir+PRE+"_ttvv.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttbbv   = uproot.open(indir+PRE+"_ttvv.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttbb   = uproot.open(indir+PRE+"_ttbb.root")["Delphes"].arrays(dnn_input,library="pd")
df_ttbbbb = uproot.open(indir+PRE+"_ttbbbb.root")["Delphes"].arrays(dnn_input,library="pd")
df_tttt   = uproot.open(indir+PRE+"_tttt.root")["Delphes"].arrays(dnn_input,library="pd")
df_tthh["category"]   = 0
df_tth["category"]    = 1
df_ttzh["category"]   = 1
df_ttbbh["category"]  = 1
df_ttvv["category"]   = 2
df_ttbbv["category"]  = 2
df_ttbb["category"]   = 2
df_ttbbbb["category"] = 2
df_tttt["category"]   = 3

# Concat into groups
df_1 = pd.concat([df_tthh], ignore_index=True)
df_2 = pd.concat([df_tth, df_ttzh, df_ttbbh], ignore_index=True)
df_3 = pd.concat([df_ttvv, df_ttbbv, df_ttbb, df_ttbbbb], ignore_index=True)
df_4 = pd.concat([df_tttt], ignore_index=True)

# Undersampling for balance statistic
n1 = len(df_1)
n2 = len(df_2)
n3 = len(df_3)
n4 = len(df_4)
print("n1 = ", n1)
print("n2 = ", n2)
print("n3 = ", n3)
print("n4 = ", n4)
ntrain = min(n1, n2, n3, n4)
print("ntrain = ", ntrain)
df_1 = df_1.sample(n=ntrain).reset_index(drop=True)
df_2 = df_2.sample(n=ntrain).reset_index(drop=True)
df_3 = df_3.sample(n=ntrain).reset_index(drop=True)
df_4 = df_4.sample(n=ntrain).reset_index(drop=True)

# Horizontal merging
df_total = pd.concat([df_1, df_2, df_3, df_4])
df_total = df_total.sample(frac=1).reset_index(drop=True)

# Casting to ndarray & train_validation partitioning
x_total = np.array(df_total.filter(items = dnn_input))
y_total = np.array(df_total.filter(items = ["category"]))
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)

###################################################
#               Correlation Matrix                #
###################################################
# (Optional) Sample Features
#print("Plotting corr_matrix total"); plot_corrMatrix(df_total, outdir,"total")
#print("Plotting corr_matrix tthh"); plot_corrMatrix(df_tthh, outdir,"tthh")
#print("Plotting corr_matrix tthbb"); plot_corrMatrix(df_tthbb, outdir,"tthbb")
#print("Plotting corr_matrix ttbb"); plot_corrMatrix(df_ttbb, outdir,"ttbb")
#print("Plotting corr_matrix ttbbbb"); plot_corrMatrix(df_ttbbbb, outdir,"ttbbbb")

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 20; batch_size = 256; print("batch size :", batch_size)
activation_function='relu'
weight_initializer = 'random_normal'
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)

model = tf.keras.models.Sequential()
###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
###############    Hidden Layer     ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(10, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(15, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(10, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(len(process_names), activation="softmax"))
###################################################

###############    Compile Model    ###############    
model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), 
              loss="sparse_categorical_crossentropy", 
              metrics = ["accuracy", "sparse_categorical_accuracy"])
model.summary()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); print("pred_train :", pred_train); pred_train_arg = np.argmax(pred_train, axis=1)
pred_val   = model.predict(x_val);   print(pred_val); pred_val_arg = np.argmax(pred_val, axis=1)
print("Is it similar?")
print("Answer :     " , y_val.T)
print("Prediction : " , pred_val_arg)

train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[0]]).T, columns=["True", "Pred"]) # True0~4,Pred0~1<"tthh" 
val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[0]]).T, columns=["True", "Pred"])
# result is prediction score respect to [G1, G2, G3, G4], add up to 1.

###################################################
#                Confusion Matrix                 #
###################################################
print("#           CONFUSION MATRIX             #")
plot_confusion_matrix(y_val, pred_val_arg, classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val_arg, classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train_arg, classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train_arg, classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

print("train_result : ")
print(train_result)

#plot_output_dist(train_result, val_result, sig="tthh", savedir=outdir)
plot_output_dist2(train_result, val_result, sig="tthh", savedir=outdir)
plot_performance(hist=hist, savedir=outdir)

###################################################
#                    Accuracy                     #
###################################################
print("#               ACCURACY                  #")
train_results = model.evaluate(x_train, y_train) # Cause you set two : "accuracy", "sparse_categorical_accuracy"
train_loss = train_results[0]
train_acc = train_results[1]
print(f"Train accuracy: {train_acc * 100:.2f}%")
val_results = model.evaluate(x_val, y_val)
val_loss = val_results[0]
val_acc = val_results[1]
print(f"Validation accuracy: {val_acc * 100:.2f}%")

###################################################
#                Feature Importance               #
###################################################
# SHAP
print("#           SHAP Feature importane            ")
explainer = shap.KernelExplainer(model.predict, x_train[:100])
shap_values = explainer.shap_values(x_train[:100])
shap.summary_plot(shap_values, x_train, plot_type='bar', max_display=50, feature_names=dnn_input, show=False)
plt.savefig(outdir+'/shap_summary_plot.pdf')

###################################################
#                  Write  TTRee                   #
###################################################
'''
# ROOT file & Tree object.
f = ROOT.TFile(outdir+"/dnn2_output.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")
# Empty array.
category_array = array('f', [0])
higgs_mass_array = array('f', [0])
higgs_mass_sub_array = array('f', [0])
higgs_mass_sum_array = array('f', [0])
bfh_dr_array = array('f', [0])
bfh_Ht_array = array('f', [0])
bfh_dEta_array = array('f', [0])
bfh_mbmb_array = array('f', [0])

# Attatch tree branches with empty array.
tree.Branch('category', category_array, 'category/F')
tree.Branch('higgs_mass', higgs_mass_array, 'higgs_mass/F')
tree.Branch('higgs_mass_sub', higgs_mass_sub_array, 'higgs_mass_sub/F')
tree.Branch('higgs_mass_sum', higgs_mass_sum_array, 'higgs_mass_sum/F')
tree.Branch('bfh_dr', bfh_dr_array, 'bfh_dr/F')
tree.Branch('bfh_Ht', bfh_Ht_array, 'bfh_Ht/F')
tree.Branch('bfh_dEta', bfh_dEta_array, 'bfh_dEta/F')
tree.Branch('bfh_mbmb', bfh_mbmb_array, 'bfh_mbmb/F')

# Fill array with data.
for _, row in df_total.iterrows():
    category_array[0] = row['category']
    higgs_mass_array[0] = row['higgs_mass']
    higgs_mass_sub_array[0] = row['higgs_mass_sub']
    higgs_mass_sum_array[0] = row['higgs_mass_sum']
    bfh_dr_array[0] = row['bfh_dr']
    bfh_Ht_array[0] = row['bfh_Ht']
    bfh_dEta_array[0] = row['bfh_dEta']
    bfh_mbmb_array[0] = row['bfh_mbmb']
    tree.Fill()

f.Write()
f.Close()
print(outdir)
print("root file is written.")
'''
# Prepare dataframe to be written
df_val = pd.DataFrame()
df_val["category"] = y_val.T[0]
df_val["G1"] = np.array(pred_val.T[0])
df_val["G2"] = np.array(pred_val.T[1])
df_val["G3"] = np.array(pred_val.T[2])
df_val["G4"] = np.array(pred_val.T[3])
df_val["DNN"] = df_val["G1"]/(df_val["G2"]+df_val["G3"]+df_val["G4"])

f = ROOT.TFile(outdir+"/dnn2_result.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")

# Empty array
category_array = array('f', [0])
G1_array = array('f', [0])
G2_array = array('f', [0])
G3_array = array('f', [0])
G4_array = array('f', [0])
DNN_array = array('f', [0])

# Attatch array on TTree
tree.Branch('category', category_array, 'category/F')
tree.Branch('G1', G1_array, 'G1/F')
tree.Branch('G2', G2_array, 'G2/F')
tree.Branch('G3', G3_array, 'G3/F')
tree.Branch('G4', G4_array, 'G4/F')
tree.Branch('DNN', DNN_array, 'DNN/F')

# DataFrame의 데이터를 TTree에 채우기
for _, row in df_val.iterrows():
    category_array[0] = row['category']
    G1_array[0] = row['G1']
    G2_array[0] = row['G2']
    G3_array[0] = row['G3']
    G4_array[0] = row['G4']
    DNN_array[0] = row['DNN']
    tree.Fill()

# 파일 저장 및 종료
f.Write()
f.Close()
print(outdir)
print("DNN score written :")
print(outdir+"/dnn2_result.root")

###################################################
#                     Time                        #
###################################################
print("Number of full data (Group x 4): ", ntrain*4)
end_time = time.time()
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
print("---Far Done---")
