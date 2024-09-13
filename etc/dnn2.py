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

###################################################
#                     I/O                         #
###################################################
indir = "./samples1/"; PRE = "FULL_0620"
outdir = "./DNN_result/" + PRE + "/LetsFind_tthh/bCat_higgs5_2Mat/"    # modify #
os.makedirs(outdir, exist_ok=True)
process_names = ["tthh", "tthbb", "ttbb", "ttbbbb"]

dnn1vars = [
     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_m",
     "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_m",

     # bb_dr
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr",
     "b2b3_dr", "b2b4_dr", "b2b5_dr",
     "b3b4_dr", "b3b5_dr",
     "b4b5_dr",
            ]


addvars = [

     "bJet_size",

     # Lepton
     "Lep_size",
     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",
     "MET_E",


    # Defined Kinematic vars
     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_m", "bb_twist",
#     "chi2_Higgs_pt", "chi2_Higgs_eta", "chi2_Higgs_phi", "chi2_Higgs_m"
           ]
 
inputvars = dnn1vars + addvars
newvars = [
#        "pred_bfh", "2bfh_1", "2bfh_2", "2bfh_3", "2bfh_4", "2bfh_5", "2bfh_6", "2bfh_7", "2bfh_8", "2bfh_9", "2bfh_10",
#        "bft_1", "bft_2", "bft_3", "bft_4", "bft_5",        
        "higgs_mass", "higgs_mass_sub", "X_higgs", "bfh_dr", "bfh_Ht", "bfh_dEta", "bfh_Phi", "bfh_mbmb"
        ]
inputvars_2 = inputvars #+ newvars
openvars = inputvars

###################################################
#                PreProcessing_1                  #
###################################################
df_tthh   = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(openvars,library="pd")
df_tthbb  = uproot.open(indir+PRE+"_tthbb.root")["Delphes"].arrays(openvars,library="pd")
df_ttbbbb = uproot.open(indir+PRE+"_ttbbbb.root")["Delphes"].arrays(openvars,library="pd")
df_ttbb   = uproot.open(indir+PRE+"_ttbb.root")["Delphes"].arrays(openvars,library="pd")
df_tthh["category"]   = 0
df_tthbb["category"]  = 1
df_ttbb["category"]   = 2
df_ttbbbb["category"] = 3
print("Columns", df_tthh.columns)

ntthh   = len(df_tthh)
ntthbb  = len(df_tthbb)
nttbb   = len(df_ttbb) 
nttbbbb = len(df_ttbbbb)
ntrain  = min(ntthh, ntthbb, nttbb, nttbbbb)
print("ntthh = ", ntthh)
print("ntthbb = ", ntthbb)
print("nttbb = ", nttbb)
print("nttbbbb = ", nttbbbb)
print("Minimum = ", ntrain)

df_tthh   = df_tthh.sample(n=ntrain).reset_index(drop=True)
df_tthbb  = df_tthbb.sample(n=ntrain).reset_index(drop=True)
df_ttbb   = df_ttbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbbb = df_ttbbbb.sample(n=ntrain).reset_index(drop=True)

# X-Y Partition
df_total = pd.concat([df_tthh, df_tthbb, df_ttbb, df_ttbbbb])
df_total = df_total.sample(frac=1).reset_index(drop=True)
x_bCat  = np.array(df_total.filter(items = dnn1vars))
x_total = np.array(df_total.filter(items = inputvars))
y_total = np.array(df_total.filter(items = ["category"]))

###################################################
#               bJet Classification               #
###################################################
bfh_dir = "DNN_result/B_0620/bJetCassification/bCat_higgs5_2Mat/best_model.h5"  ## modify ##
bft_dir = "DNN_result/B_0620/bJetCassification/bCat_top_1/best_model.h5"
bfh_model = tf.keras.models.load_model(bfh_dir)
bft_model = tf.keras.models.load_model(bft_dir)
bfh_model.summary()
bft_model.summary()
print("x_total: ", x_total); print("x_total_shape: ", x_total.shape)
_pred_bfh = bfh_model.predict(x_bCat); print("pred_bfh : ", _pred_bfh)
_pred_bft = bft_model.predict(x_bCat); print("_pred_bft : ", _pred_bft.shape)
pred_bfh = np.argmax(_pred_bfh, axis=1) # (arg 0~10 that with the biggest prob)

### NEW VARIABLES ###
df_total["pred_bfh"] = pred_bfh
for i in range(10):
    column_name = f"2bfh_{i + 1}"
    df_total[column_name] = _pred_bfh[:, i]
for i in range(5):
    column_name = f"bft_{i + 1}"
    df_total[column_name] = _pred_bft[:, i]
print(df_total)

df_total["higgs_mass_list"] = df_total.apply(higgs_5_2, axis = 1) # MODIFY #
df_total["higgs_mass"] = df_total["higgs_mass_list"].apply(lambda x: x[0])
df_total["higgs_mass_sub"] = df_total["higgs_mass_list"].apply(lambda x: x[1])
df_total["higgs_mass_sum"] = df_total["higgs_mass_list"].apply(lambda x: x[2])
df_total["X_higgs"] = df_total.apply(X_higgs, axis = 1) # MODIFY #
print(df_total)

df_total["bfh_Vars"] = df_total.apply(bfh_Vars, axis = 1)
df_total["bfh_dr"] = df_total["bfh_Vars"].apply(lambda x: x[0])
df_total["bfh_Ht"] = df_total["bfh_Vars"].apply(lambda x: x[1])
df_total["bfh_dEta"] = df_total["bfh_Vars"].apply(lambda x: x[2])
df_total["bfh_dPhi"] = df_total["bfh_Vars"].apply(lambda x: x[3])
df_total["bfh_mbmb"] = df_total["bfh_Vars"].apply(lambda x: x[4])

df_plot_tthh = df_total[(df_total["category"] == 0)]
df_plot_tthbb = df_total[(df_total["category"] == 1)]
df_plot_ttbb = df_total[(df_total["category"] == 2)]
df_plot_ttbbbb = df_total[(df_total["category"] == 3)]


###################################################
#                PreProcessing_2                  #
###################################################
_x_total = df_total.filter(items = inputvars_2)
#_x_total = _x_total.drop('pred_bfh', axis=1) # modify. w/ dnnvars un#
x_total = np.array(_x_total)
y_total = np.array(df_total.filter(items = ["category"]))

print("Final x = ", x_total)
print("Final y = ", y_total)

# Data Set Partioning
ntotal = len(y_total)
train_len = int(0.7*ntotal)
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

start_time = time.time()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])
end_time = time.time()

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); print("pres_train :");print(pred_train); pred_train_arg = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val); print(pred_val); pred_val_arg = np.argmax(pred_val, axis=1)
print("Is it similar?")
print("Prediction for validation set: ", pred_val_arg)
print("Answer for train set:         ", y_val.T)

train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[0]]).T, columns=["True", "Pred"]) # True0~4,Pred0~1<"tthh" 
val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[0]]).T, columns=["True", "Pred"])
# 0123 -> tthh ; tthbb ; ttbb ; ttbbbb : So upper ones are tthh prediction score. 

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
test_results = model.evaluate(x_val, y_val)
test_loss = test_results[0]
test_acc = test_results[1]
print(f"Test accuracy: {test_acc * 100:.2f}%")
###################################################
#                Feature Importance               #
###################################################
'''
print("#          FEATURE IMPORTANCE             #")
model_dir = outdir + '/best_model.h5'
plot_feature_importance(model_dir, x_val, _x_total.columns, outdir)
print("#        FEATURE IMPORTANCE DONE         #")

# GradientTape를 사용하여 그래디언트 계산
x_sample = x_train[0:10]  # 입력 데이터 중 하나의 샘플 선택
# BatchNormalization 레이어 생성
batch_norm_layer = tf.keras.layers.BatchNormalization()
# x_sample에 대한 BatchNormalization 적용
x_sample = batch_norm_layer(x_sample, training=False)
x_sample_tensor = tf.convert_to_tensor(x_sample, dtype=tf.float32)  # NumPy 배열을 TensorFlow 텐서로 변환
with tf.GradientTape() as tape:
    tape.watch(x_sample_tensor)  # 입력 데이터를 감시 대상으로 설정
    predictions = model(x_sample_tensor)  # 모델의 예측 계산

# 그래디언트 계산
gradients = tape.gradient(predictions, x_sample_tensor)

# 결과 출력
print("Input Data:", x_sample)
print("Model Predictions:", predictions)
print("Gradients:", gradients)

###### SHAP #####
print("Try SHAP Feature importance test")
# Calculate SHAP
explainer = shap.KernelExplainer(model.predict, x_train)
shap_values = explainer.shap_values(x_val[:100])
plt.figure()
shap.summary_plot(shap_values, x_val[:100], show=False)
plt.savefig('shap_summary_plot.pdf')
'''

###################################################
#                     Time                        #
###################################################
print("Number of full data: ", ntrain*5)
colnames = _x_total.columns; print("Columns :",colnames)
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
print("---Far Done---")

###################################################
#                  Write  TTRee                   #
###################################################
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
