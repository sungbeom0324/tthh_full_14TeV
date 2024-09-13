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
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
from matplotlib.backends.backend_pdf import PdfPages
from utils.var_functions import *

###################################################
#                     I/O                         #
###################################################
indir = "./samples2/"; PRE = "FULL_1204_14TeV"
outdir = "./CNN_result/" + PRE + "/LetsFind_tthh/bCat_higgs5_2Mat/"    # MODIFY  #
os.makedirs(outdir, exist_ok=True)
process_names = ["tthh", "tthbb", "ttbb", "ttbbbb", "ttbbcc"]

var_jet = [
# 32 = x_train.reshape(-1_nEvents, 8, 4, 1(Black/White))

     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass",
     "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass",

     # bb_dr
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr",
     "b2b3_dr", "b2b4_dr", "b2b5_dr",
     "b3b4_dr", "b3b5_dr",
     "b4b5_dr",
     "bJet_size", "bb_avg_dr",
]

var_event = [
     # Lepton
     "Lep_size",
     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",
     "MET_E", # why decrease..


     # Defined Kinematic vars
     "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
     "close_Higgs_pt", "close_Higgs_phi", "close_Higgs_mass", "close_Higgs_eta",
        ]

addvars = []
 
inputvars = var_jet + var_event + addvars
newvars = ["pred_bcat", "higgs_mass", "X_higgs", "bfh_dr", "bfh_Ht", "bfh_dEta", "bfh_Phi", "bfh_mbmb"]
inputvars_2 = inputvars + newvars
openvars = inputvars

###################################################
#                PreProcessing_1                  #
###################################################
df_tthh   = uproot.open(indir+PRE+"_tthh_di.root")["Delphes"].arrays(openvars,library="pd")
df_tthbb  = uproot.open(indir+PRE+"_tthbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbbbb = uproot.open(indir+PRE+"_ttbbbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbb   = uproot.open(indir+PRE+"_ttbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbbcc = uproot.open(indir+PRE+"_ttbbcc_di.root")["Delphes"].arrays(openvars,library="pd")
df_tthh["category"]   = 0
df_tthbb["category"]  = 1
df_ttbb["category"]   = 2
df_ttbbbb["category"] = 3
df_ttbbcc["category"] = 4
print("Columns", df_tthh.columns)

ntthh   = len(df_tthh)
ntthbb  = len(df_tthbb)
nttbb   = len(df_ttbb) 
nttbbbb = len(df_ttbbbb)
nttbbcc = len(df_ttbbcc)
ntrain  = min(ntthh, ntthbb, nttbb, nttbbbb, nttbbcc)

df_tthh   = df_tthh.sample(n=ntrain).reset_index(drop=True)
df_tthbb  = df_tthbb.sample(n=ntrain).reset_index(drop=True)
df_ttbb   = df_ttbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbbb = df_ttbbbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbcc = df_ttbbcc.sample(n=ntrain).reset_index(drop=True)

# X-Y Partition
df_total = pd.concat([df_tthh, df_tthbb, df_ttbb, df_ttbbbb, df_ttbbcc])
df_total = df_total.sample(frac=1).reset_index(drop=True)
x_jet   = np.array(df_total.filter(items = var_jet))
x_event = np.array(df_total.filter(items = var_event))
x_total = np.array(df_total.filter(items = inputvars))
y_total = np.array(df_total.filter(items = ["category"]))

###################################################
#               bJet Classification               #
###################################################
model_dir = "CNN_result/B_1217_14TeV/bJetCassification/bCat_higgs5_2Mat/best_model.h5"  # MODIFY #
bft_model = tf.keras.models.load_model(model_dir)
bft_model.summary()
pred_bcat = bft_model.predict([x_jet.reshape(-1,8,4,1), x_event]); print("pred_bcat : ", pred_bcat) # MODIFY #
pred_bcat = np.argmax(pred_bcat, axis=1) # (arg 0~10 that with the biggest prob)

### NEW VARIABLES ###
df_total["pred_bcat"] = pred_bcat
df_total["higgs_mass"] = df_total.apply(higgs_5_2, axis = 1) # MODIFY #
df_total["X_higgs"] = df_total.apply(X_higgs, axis = 1) # MODIFY #

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
df_plot_ttbbcc = df_total[(df_total["category"] == 4)]

### PLOTTING ###
#plot_histogram(df_plot['bfh_dr'], "bfh_dr", "bfh_dr", "AU", outdir)
#plot_histogram(df_plot['bfh_Ht'], "bfh_Ht", "bfh_Ht", "AU", outdir)
#plot_histogram(df_plot['bfh_dEta'], "bfh_dEta", "bfh_dEta", "AU", outdir)
#plot_histogram(df_plot['bfh_dPhi'], "bfh_dPhi", "bfh_dPhi", "AU", outdir)
#plot_histogram(df_plot['bfh_mbmb'], "bfh_mbmb", "bfh_mbmb", "AU", outdir)
#plot_2histogram(df_plot, 'bfh_dr', 'b1b2_dr', "bfh_dr vs b1b2_dr", "bfh_dr", "AU", outdir)
plot_multi_histograms(df_plot_tthh, "bfh_dr vs All_dr_tthh", "dr", "Events", outdir,
                     ('bfh_dr', 'dfh_dr', "red", "stepfilled"),
                     ('b1b2_dr', 'b1b2_dr', "blue", "step"),
                     ('b1b3_dr', 'b1b3_dr', "green", "step"),
                     ('b1b4_dr', 'b1b4_dr', "yellow", "step"),
                     ('b1b5_dr', 'b1b5_dr', "c", "step"),
                     ('b2b3_dr', 'b2b3_dr', "m", "step"),
                     ('b2b4_dr', 'b2b4_dr', "k", "step"),
                     ('b2b5_dr', 'b2b5_dr', "purple", "step"),
                     ('b3b4_dr', 'b3b4_dr', "green", "step"),
                     ('b3b5_dr', 'b3b5_dr', "blue", "step"),
                     ('b4b5_dr', 'b4b5_dr', "yellow", "step"),
)
plot_multi_histograms_dfs("bfh_dEta", "bfh_dEta_Comparison", "dEta", "Normalized Entries", outdir,
                          (df_plot_tthh, 'tthh', 'red', 'stepfilled'),
                          (df_plot_ttbbbb, 'ttbbbb', 'blue', 'step'),
                          (df_plot_ttbbcc, 'ttbbcc', 'purple', 'step'),
                          (df_plot_ttbb, 'ttbb', 'yellow', 'step'),
                          (df_plot_tthbb, 'tthbb', 'orange', 'stepfilled')
        )
plot_multi_histograms_dfs("bfh_dPhi", "bfh_dPhi_Comparison", "dPhi", "Normalized Entries", outdir,
                          (df_plot_tthh, 'tthh', 'red', 'stepfilled'),
                          (df_plot_ttbbbb, 'ttbbbb', 'blue', 'step'),
                          (df_plot_ttbbcc, 'ttbbcc', 'purple', 'step'),
                          (df_plot_ttbb, 'ttbb', 'yellow', 'step'),
                          (df_plot_tthbb, 'tthbb', 'orange', 'stepfilled')
        )
plot_multi_histograms_dfs("bfh_mbmb", "bfh_mbmb_Comparison", "mass [Gev]", "Normalized Entries", outdir,
                          (df_plot_tthh, 'tthh', 'red', 'stepfilled'),
                          (df_plot_ttbbbb, 'ttbbbb', 'blue', 'step'),
                          (df_plot_ttbbcc, 'ttbbcc', 'purple', 'step'),
                          (df_plot_ttbb, 'ttbb', 'yellow', 'step'),
                          (df_plot_tthbb, 'tthbb', 'orange', 'stepfilled')
        )
plot_multi_histograms_dfs("higgs_mass", "higgs_mass_Comparison", "mass [GeV]", "Normalized Entries", outdir,
                          (df_plot_tthh, 'tthh', 'red', 'stepfilled'),
                          (df_plot_ttbbbb, 'ttbbbb', 'blue', 'step'),
                          (df_plot_ttbbcc, 'ttbbcc', 'purple', 'step'),
                          (df_plot_ttbb, 'ttbb', 'yellow', 'step'),
                          (df_plot_tthbb, 'tthbb', 'orange', 'stepfilled')
        )
plot_multi_histograms_dfs("X_higgs", "X_higgs_Comparison", "dMass [GeV]", "Normalized Entries", outdir,
                          (df_plot_tthh, 'tthh', 'red', 'stepfilled'),
                          (df_plot_ttbbbb, 'ttbbbb', 'blue', 'step'),
                          (df_plot_ttbbcc, 'ttbbcc', 'purple', 'step'),
                          (df_plot_ttbb, 'ttbb', 'yellow', 'step'),
                          (df_plot_tthbb, 'tthbb', 'orange', 'stepfilled')
        )



### PLOT HIGGS ###
higgs_mass = np.array(df_plot_tthh.filter(items = ["higgs_mass"]))
hist, bin_edges = np.histogram(higgs_mass, bins=39, range = (10, 400), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
initial_guess = [125, 30, 1] # mu, sigma, A
popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess)

mu, sigma, A = popt
lower_bound = mu - sigma
upper_bound = mu + sigma

# Assuming df_new is your DataFrame containing 'close_Higgs_mass'
close_higgs_mass = np.array(df_plot_tthh.filter(items=["close_Higgs_mass"]))
hist_close, bin_edges_close = np.histogram(close_higgs_mass, bins=39, range=(10, 400), density=True) 
bin_centers_close = (bin_edges_close[:-1] + bin_edges_close[1:]) / 2

# Define Gaussian function
def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Initial guess for curve fitting
initial_guess = [125, 30, 1]  # mu, sigma, A

# Perform curve fitting
popt_close, _ = curve_fit(gaussian, bin_centers_close, hist_close, p0=initial_guess)

# Extract parameters from the fitting
mu_close, sigma_close, A_close = popt_close
lower_bound_close = mu_close - sigma_close
upper_bound_close = mu_close + sigma_close

# Plotting
plt.hist(higgs_mass, bins=39, range=(10, 400), density=True, alpha=0.2, color='g', label='Higgs_mass Histogram')
plt.plot(bin_centers, gaussian(bin_centers, *popt), 'r-', label='Higgs_mass Fit')

plt.hist(close_higgs_mass, bins=39, range=(10, 400), density=True, alpha=0.2, color='b', label='close_Higgs_mass Histogram')
plt.plot(bin_centers_close, gaussian(bin_centers_close, *popt_close), 'b-', label='close_Higgs_mass Fit')

#plt.axvline(x=lower_bound, color='orange', linestyle='--', label='Higgs_mass 1 Sigma')
#plt.axvline(x=upper_bound, color='orange', linestyle='--')
#plt.axvline(x=lower_bound_close, color='cyan', linestyle='--', label='close_Higgs_mass 1 Sigma')
#plt.axvline(x=upper_bound_close, color='cyan', linestyle='--')

plt.axvline(x=125, color='magenta', linestyle='--', label='x=125')

plt.fill_betweenx([0, A], lower_bound, upper_bound, color='orange', alpha=0.3, label='Higgs_mass 1 Sigma Range')
plt.fill_betweenx([0, A_close], lower_bound_close, upper_bound_close, color='cyan', alpha=0.3, label='close_Higgs_mass 1 Sigma Range')

plt.xlabel('Higgs Mass [GeV]')
plt.ylabel('Probability Density')
plt.title(f'Higgs_mass: mu = {popt[0]:.2f}, sigma = {popt[1]:.2f}\n close_Higgs_mass: mu = {popt_close[0]:.2f}, sigma = {popt_close[1]:.2f}')
plt.legend()
plt.savefig(outdir+'/higgs_mass_and_close_higgs_mass_dnn.pdf')
plt.show()

'''
plt.hist(higgs_mass, bins=39, range = (10,400), density=True, alpha=0.2, color='g', label='Histogram')
plt.plot(bin_centers, gaussian(bin_centers, *popt), 'r-', label='Fit')
plt.axvline(x=lower_bound, color='b', linestyle='--', label='1 Sigma')
plt.axvline(x=upper_bound, color='b', linestyle='--')
plt.axvline(x=125, color='magenta', linestyle='--', label='x=125')
plt.fill_betweenx([0, A], lower_bound, upper_bound, color='blue', alpha=0.3, label='1 Sigma Range')
plt.xlabel('Higgs Mass [GeV]')
plt.ylabel('Probability Density')
plt.title(f'Higgs_mass: mu = {popt[0]:.2f}, sigma = {popt[1]:.2f}')
plt.legend()
plt.savefig(outdir+'/higgs_mass_dnn.pdf')
plt.show()
'''

###################################################
#                PreProcessing_2                  #
###################################################
x_total = np.array(df_total.filter(items = inputvars_2))
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
#print("Plotting corr_matrix ttbbcc"); plot_corrMatrix(df_ttbbcc, outdir,"ttbbcc")

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 50; batch_size = 1024; print("batch size :", batch_size)
activation_function='relu'
weight_initializer = 'random_normal'
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)

model = tf.keras.models.Sequential()
###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
###############    Hidden Layer     ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dropout(0.2))    
model.add(tf.keras.layers.Dense(50, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(len(process_names), activation="softmax"))
###################################################

###############    Compile Model    ###############    
model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), 
              loss="sparse_categorical_crossentropy", 
              metrics = ["accuracy", "sparse_categorical_accuracy"])
model.summary()
###################################################
start_time = time.time()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])
end_time = time.time()

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); print(pred_train); pred_train_arg = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val); print(pred_val); pred_val_arg = np.argmax(pred_val, axis=1)
print("Is it similar?")
print("Prediction for validation set: ", pred_val_arg)
print("Answer for train set:         ", y_val.T)

train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[0]]).T, columns=["True", "Pred"]) # True0~4,Pred0~1<"tthh" 
val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[0]]).T, columns=["True", "Pred"])

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

plot_performance(hist=hist, savedir=outdir)
plot_output_dist(train_result, val_result,sig="tthh", savedir=outdir)

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

###################################################
#                     Time                        #
###################################################
print("Number of full data: ", ntrain*5)
colnames = df_total.columns; print("Columns :",colnames)
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
print("---Done---")
