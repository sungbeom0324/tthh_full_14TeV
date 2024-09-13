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
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
from matplotlib.backends.backend_pdf import PdfPages
from utils.var_functions import *

start_time = time.time()
###################################################
#                     I/O                         #
###################################################    
indir = "./samples2/"; OS = "Full1106"
outdir = "./DNN_result/" + OS + "/LetsFind_tthh/binary/higgs_2" ####MODIFY!!#####
os.makedirs(outdir, exist_ok=True)
dnn1vars = [ 
     "bJet_size",
     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass",
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b2b3_dr", "b2b4_dr", "b3b4_dr",

#     "Lep_size",
#     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
#     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",

#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
#     "close_Higgs_pt", "close_Higgs_eta", "close_Higgs_phi", "close_Higgs_mass"   
            ]

addvars = []
inputvars = dnn1vars + addvars
newvars = ["pred_bcat", "higgs_mass"]
inputvars_2 = inputvars + newvars
process_names = ["tthh_signal", "backgrounds"]
openvars = inputvars

###################################################
#                PreProcessing_1                  #
###################################################
df_tthh   = uproot.open(indir+OS+"_tthh_di.root")["Delphes"].arrays(openvars,library="pd")
df_tthbb  = uproot.open(indir+OS+"_tthbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbbbb = uproot.open(indir+OS+"_ttbbbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbb   = uproot.open(indir+OS+"_ttbb_di.root")["Delphes"].arrays(openvars,library="pd")
df_ttbbcc = uproot.open(indir+OS+"_ttbbcc_di.root")["Delphes"].arrays(openvars,library="pd")
print(df_tthh.columns)

# Taylor Samples
ntthh   = len(df_tthh)
ntthbb  = len(df_tthbb)
nttbb   = len(df_ttbb) 
nttbbbb = len(df_ttbbbb)
nttbbcc = len(df_ttbbcc)
ntrain = min(ntthh, ntthbb, nttbb, nttbbbb, nttbbcc)

df_tthh   = df_tthh.sample(n=ntrain).reset_index(drop=True)
df_tthbb  = df_tthbb.sample(n=ntrain).reset_index(drop=True)
df_ttbb   = df_ttbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbbb = df_ttbbbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbcc = df_ttbbcc.sample(n=ntrain).reset_index(drop=True)
df_bkg = pd.concat([df_tthbb, df_ttbb, df_ttbbbb, df_ttbbcc]); nbkg = len(df_bkg)
ntrain_ = min(ntrain, nbkg); print("ntrain_", ntrain_)
df_tthh   = df_tthh.sample(n=ntrain_).reset_index(drop=True)
df_bkg   = df_bkg.sample(n=ntrain_).reset_index(drop=True)

df_tthh["category"]   = 0
df_bkg["category"] = 1

# X-Y Partioning
df_total = pd.concat([df_tthh, df_bkg])
df_total = df_total.sample(frac=1).reset_index(drop=True)
x_total = np.array(df_total.filter(items = inputvars))
y_total = np.array(df_total.filter(items = ["category"]))

###################################################
#               bJet Classification               #
###################################################
model_dir = "DNN_result/Test1105_bJet/bJetCassification/bCat_higgs_2/best_model.h5"  ###### MODIFY!!! #######
bft_model = tf.keras.models.load_model(model_dir)
bft_model.summary()
pred_bcat = bft_model.predict(x_total); print("pred_bcat : ", pred_bcat) # (ndarray; nEvt * 11; probability)
pred_bcat = np.argmax(pred_bcat, axis=1) # (a value 0~10 that with the biggest prob)
df_total["pred_bcat"] = pred_bcat
######## Define New Variables ################
df_total["higgs_mass"] = df_total.apply(higgs_higgs_2, axis = 1) ##### MODIFY!!! #######

df_tthh_new = df_total[df_total["category"]==0]
higgs_mass = np.array(df_tthh_new.filter(items = ["higgs_mass"]))
hist, bin_edges = np.histogram(higgs_mass, bins=39, range = (10, 400), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
initial_guess = [125, 30, 1] # mu, sigma, A
popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess)

mu, sigma, A = popt
lower_bound = mu - sigma
upper_bound = mu + sigma

#plt.hist(higgs_mass, bins=39, range = (10,400), density=True, alpha=0.2, color='g', label='Histogram')
#plt.plot(bin_centers, gaussian(bin_centers, *popt), 'r-', label='Fit')
#plt.axvline(x=lower_bound, color='b', linestyle='--', label='1 Sigma')
#plt.axvline(x=upper_bound, color='b', linestyle='--')
#plt.axvline(x=125, color='magenta', linestyle='--', label='x=125')
#plt.fill_betweenx([0, A], lower_bound, upper_bound, color='blue', alpha=0.3, label='1 Sigma Range')
#plt.title(f'Higgs_mass: mu = {popt[0]:.2f}, sigma = {popt[1]:.2f}')
#plt.legend()
#plt.savefig(outdir+'/histogram_with_fit.pdf')
#plt.show()

###################################################
#                PreProcessing_2                  #
###################################################
x_total = np.array(df_total.filter(items = inputvars_2))
y_total = np.array(df_total.filter(items = ["category"]))
print("New x_total = ", x_total)
print("New y_total = ", y_total)

# Data Set Partioning
ntotal = len(y_total)
train_len = int(0.7*ntotal)
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)
    
# (Optional) Sample Features
colnames = df_total.columns; print("Columns :",colnames)
#print("Plotting corr_matrix total"); plot_corrMatrix(df_total, outdir,"total")
#print("Plotting corr_matrix tthh"); plot_corrMatrix(df_tthh, outdir,"tthh")
#print("Plotting corr_matrix tthbb"); plot_corrMatrix(df_tthbb, outdir,"tthbb")
#print("Plotting corr_matrix ttbb"); plot_corrMatrix(df_ttbb, outdir,"ttbb")
#print("Plotting corr_matrix ttbbbb"); plot_corrMatrix(df_ttbbbb, outdir,"ttbbbb")
#print("Plotting corr_matrix ttbbcc"); plot_corrMatrix(df_ttbbcc, outdir,"ttbbcc")

 
###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 10; batch_size = 1024; print("batch size :", batch_size)
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
model.add(tf.keras.layers.Dense(50, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(len(process_names), activation="softmax"))
###################################################

###############    Compile Model    ###############    
model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), loss="sparse_categorical_crossentropy", metrics = ["accuracy"])
model.summary()
###################################################

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); print(pred_train); pred_train = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val); print(pred_val); pred_val = np.argmax(pred_val, axis=1)
pred_val_2 = model.predict(x_val); print(pred_val_2); pred_val_2 = np.argmax(pred_val_2, axis=1)

###################################################
#                Confusion Matrix                 #
###################################################
plot_confusion_matrix(y_val, pred_val, classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

###################################################
#                    Accuracy                     #
###################################################
print("#               ACCURACY                  #")
train_loss, train_acc = model.evaluate(x_train, y_train)
print(f"Train accuracy: {train_acc * 100:.2f}%")
test_loss, test_acc = model.evaluate(x_val, y_val)
print(f"Test accuracy: {test_acc * 100:.2f}%")

end_time = time.time()
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
###################################################
#                      ROC                        #
###################################################
print("#              ROC-Curve                 #")
fpr, tpr, thresholds = roc_curve(y_val, pred_val)
roc_auc = roc_auc_score(y_val, pred_val)
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc='lower right')
pdf_pages = PdfPages(outdir+'/roc_curve.pdf')  # 파일 이름을 원하는 이름으로 변경
pdf_pages.savefig()
pdf_pages.close()

###################################################
#                Feature Importance               #
###################################################
'''
model_2 = tf.keras.models.load_model(outdir+"/best_model.h5")
model_2.summary()

print("feature importance starts")
input_data = tf.convert_to_tensor(x_val, dtype=tf.float32)
print(input_data)
name_inputvar = inputvars_2
n_evts = len(x_val)
n_var = len(name_inputvar)
mean_grads = n_var*[0.0]
all_grads = []
#mean_jacobian = np.zeros(len(name_inputvar))
#jacobian_matrix = np.zeros((len(name_inputvar),len(name_inputvar)))
for i, event in enumerate(x_val):
#    print(i,"/",n_evts)
    with tf.GradientTape() as tape:
        inputs = tf.Variable([event])
        tape.watch(inputs)
        prediction = model(inputs)[:, 1]
    first_order_gradients = tape.gradient(prediction, inputs)
    gradiants = tf.abs(first_order_gradients)
    numpy_array = gradiants.numpy()[0]
    all_grads.append(numpy_array)
    #print(i,numpy_array , len(numpy_array))
    for n in range(len(mean_grads)):
        mean_grads[n] += abs(numpy_array[n])/n_evts

print(mean_grads) # This is just final iteration (final event), not important yet.
df = pd.DataFrame(all_grads)
print(df)
df.to_csv('data.csv', index=False)

gradiants = tf.abs(first_order_gradients)
numpy_array = gradiants.numpy()
df = pd.DataFrame(numpy_array)
print(df)
df.to_csv(outdir+'/data.csv', index=False)
feature_importance_first_order  = tf.reduce_mean(tf.abs(first_order_gradients), axis=0)
feature_importance_dict = dict(zip(name_inputvar, feature_importance_first_order.numpy())) # Mapping
#feature_importance_Secondorder  = tf.reduce_mean(tf.abs(second_order_gradients), axis=0)
feature_importance_series = pd.Series(feature_importance_dict)

print("Feature importance series?")
print(feature_importance_series)
max_importance_score = feature_importance_series.max()
for i, importance_score in enumerate(feature_importance_first_order):
    print(f"Feature {i+1} , {name_inputvar[i]} Importance: {max_importance_score-importance_score:.5f}")

print(feature_importance_series.index, feature_importance_series.values)

# Order the Feature Importance
sorted_indices = np.argsort(feature_importance_series.values)[::-1]
sorted_importance = feature_importance_series.values[sorted_indices]
sorted_feature_names = feature_importance_series.index[sorted_indices]

plt.figure(figsize=(10, 10))
plt.barh(sorted_feature_names, max_importance_score - sorted_importance)
plt.xlabel('Feature Importance')
plt.ylabel('Feature Names')
plt.savefig(outdir+'/first_order_gradient_importance.png')
plt.title('First-Order Gradient Feature Importance')
plt.show()
'''
###################################################

print("---Done---")
