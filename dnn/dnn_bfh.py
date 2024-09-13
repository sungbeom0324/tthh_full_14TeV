# ssh gpu-0-X ; conda activate py36
import os
import sys
import time
#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import shap
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
import ROOT
from array import array
import json

###################################################
#                     I/O                         #
###################################################
indir = "./skimmed/"; PRE = "b" # Should be apart from data for event classification.
outdir = "./dnn_result/" + PRE + "fh/" # modify #
os.makedirs(outdir, exist_ok=True)
class_names = ["Cat1", "Cat2", "Cat3", "Cat4", "Cat5", "Cat6", "Cat7", "Cat8", "Cat9", "Cat10", "NoCat"]

with open('./dnn/dnn_input.json', 'r') as file:
    data = json.load(file)

input_1 = data["input_1"]

input_cat = ["bCat_higgs5_2Mat", "bCat_higgs5_2Mat_1", "bCat_higgs5_2Mat_2", "isMatchable"] # MODIFY #
openvars = input_1 + input_cat

###################################################
#                 PreProcessing                   #
###################################################
df_tthh = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(openvars,library="pd")
df_cat1 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 0]
df_cat2 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 1]
df_cat3 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 2]
df_cat4 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 3]
df_cat5 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 4]
df_cat6 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 5]
df_cat7 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 6]
df_cat8 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 7]
df_cat9 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 8]
df_cat10 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 9]
df_cat11 = df_tthh.loc[df_tthh["bCat_higgs5_2Mat"] == 10]

nCat1 = len(df_cat1) 
nCat2 = len(df_cat2)
nCat3 = len(df_cat3)
nCat4 = len(df_cat4)
nCat5 = len(df_cat5)
nCat6 = len(df_cat6)
nCat7 = len(df_cat7)
nCat8 = len(df_cat8)
nCat9 = len(df_cat9)
nCat10 = len(df_cat10)
nCat11 = len(df_cat11)
ntrain = min(nCat1, nCat2, nCat3, nCat4, nCat5, nCat6, nCat7, nCat8, nCat9, nCat10)
print("ntrain = ", ntrain)

df_cat1 = df_cat1.sample(n=ntrain).reset_index(drop=True)
df_cat2 = df_cat2.sample(n=ntrain).reset_index(drop=True)
df_cat3 = df_cat3.sample(n=ntrain).reset_index(drop=True)
df_cat4 = df_cat4.sample(n=ntrain).reset_index(drop=True)
df_cat5 = df_cat5.sample(n=ntrain).reset_index(drop=True)
df_cat6 = df_cat6.sample(n=ntrain).reset_index(drop=True)
df_cat7 = df_cat7.sample(n=ntrain).reset_index(drop=True)
df_cat8 = df_cat8.sample(n=ntrain).reset_index(drop=True)
df_cat9 = df_cat9.sample(n=ntrain).reset_index(drop=True)
df_cat10 = df_cat10.sample(n=ntrain).reset_index(drop=True)
df_cat11 = df_cat11.sample(n=ntrain).reset_index(drop=True)
print("df_cat1", df_cat1)

# Horiznotal merging
df_total = pd.concat([df_cat1, df_cat2, df_cat3, df_cat4, df_cat5, df_cat6, df_cat7, df_cat8, df_cat9, df_cat10, df_cat11])
df_total = df_total.sample(frac=1).reset_index(drop=True)

##################################################
#             Train-Val Partitioning             #
##################################################
# Casting to ndarray & train_validation partitioning
x_total = np.array(df_total.filter(items = input_1))
#y_total = np.array(df_total.filter(items = ["bCat_higgs5_2Mat"]))
y_total = np.array(df_total.filter(items = input_cat))
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)

y_val_1 = y_val.T[1]; print("y_val_1\n", y_val_1)
y_val_2 = y_val.T[2]; print("y_val_2\n", y_val_2)

df_val = pd.DataFrame(x_val, columns=input_1)
df_val['bCat_higgs5_2Mat'] = np.array(y_val.T[0])
df_val['bCat_higgs5_2Mat_1'] = np.array(y_val.T[1])
df_val['bCat_higgs5_2Mat_2'] = np.array(y_val.T[2])
df_val['isMatchable'] = np.array(y_val.T[3])

y_train = y_train.T[0].T
y_val = y_val.T[0].T

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 10; batch_size = 512; print("batch size :", batch_size)
activation_function='relu'
weight_initializer = 'random_normal'
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)

model = tf.keras.models.Sequential()
###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
###############    Hidden Layer     ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(30, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
print("class_names : ", len(class_names))
model.add(tf.keras.layers.Dense(len(class_names), activation="softmax"))
###################################################
start_time = time.time()

model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), 
                   loss="sparse_categorical_crossentropy", 
                   metrics = ["accuracy", "sparse_categorical_accuracy"])
model.summary()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])
end_time = time.time()

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); pred_train = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val) ; pred_val = np.argmax(pred_val, axis=1)
print("Is it similar?")
print("Prediction for validation set: ", pred_val)
print("Answer for train set:         ", y_val.T)

###################################################
#         Confusion Matrix, Acc Curve             #
###################################################
print("#           CONFUSION MATRIX             #")
plot_confusion_matrix(y_val, pred_val, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

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
#              Feature Importance                 #
###################################################
'''
# SHAP
print("#           SHAP Feature importane            ")
explainer = shap.KernelExplainer(model.predict, x_train[:100])
shap_values = explainer.shap_values(x_train[:100])
shap.summary_plot(shap_values, x_train, plot_type='bar', max_display=50, feature_names=input_dnn+input_3, show=False)
plt.savefig(outdir+'/shap_summary_plot.pdf')
'''
###################################################
#                     Time                        #
###################################################
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
print("---Done---")

###################################################
#                  Write  TTRee                   #
###################################################
# Prepare file.
f = ROOT.TFile(outdir+"/bfh.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")

# Prepare dataframe to be written
df_val["pred_val"] = np.array(pred_val)

# Empty array.
pred_val_array = array('f', [0])
bCat_higgs5_2Mat_array = array('f', [0])
bCat_higgs5_2Mat_1_array = array('f', [0])
bCat_higgs5_2Mat_2_array = array('f', [0])
#isMatchable_array = array('f', [0])

# Attach tree branches with empty array.
tree.Branch('pred_val', pred_val_array, 'pred_val/F')
tree.Branch('bCat_higgs5_2Mat', bCat_higgs5_2Mat_array, 'bCat_higgs5_2Mat/F')
tree.Branch('bCat_higgs5_2Mat_1', bCat_higgs5_2Mat_1_array, 'bCat_higgs5_2Mat_1/F')
tree.Branch('bCat_higgs5_2Mat_2', bCat_higgs5_2Mat_2_array, 'bCat_higgs5_2Mat_2/F')
#tree.Branch('isMatchable', isMatchable_array, 'isMatchable/F')

# Fill array with validation data.
for _, row in df_val.iterrows():
    pred_val_array[0] = row['pred_val']    
    bCat_higgs5_2Mat_array[0] = row['bCat_higgs5_2Mat']
    bCat_higgs5_2Mat_1_array[0] = row['bCat_higgs5_2Mat_1']
    bCat_higgs5_2Mat_2_array[0] = row['bCat_higgs5_2Mat_2']
    tree.Fill()

f.Write()
f.Close()
print("DNN score written :")
print(outdir+"bfh.root")
