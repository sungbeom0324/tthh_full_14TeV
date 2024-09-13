# ssh gpu-0-X
# conda activate py36
import os
import sys

#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split

# Criteria #
OS = "J30L20_0.2dpT_NodR"
class_names = ["tthh","ttbb", "ttbbbb", "ttbbcc"]

epochs = 1000
inputvars = [
            "bJet_size", "nMuon",
            "bJet1_pt","bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt",
            "bJet1_eta","bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta",
            "bJet1_phi","bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi",
            "Muon1_pt", "Muon2_pt",
            "Muon1_eta", "Muon2_eta",
            "Muon1_phi", "Muon2_phi"
            ]


df_sig_tt_list = []


# Input #
indir = "./samples2/"
df_tthh   = uproot.open(indir+OS+"_tthh_di.root")["Delphes"].arrays(inputvars,library="pd")
df_ttbb   = uproot.open(indir+OS+"_ttbb_di.root")["Delphes"].arrays(inputvars,library="pd")
df_ttbbbb = uproot.open(indir+OS+"_ttbbbb_di.root")["Delphes"].arrays(inputvars,library="pd")
df_ttbbcc = uproot.open(indir+OS+"_ttbbcc_di.root")["Delphes"].arrays(inputvars,library="pd")
print(type(df_tthh))


# Output #
outdir = "./DNN_result/" + OS + "/"
os.makedirs(outdir, exist_ok=True)


########## Start ############
print("Start multi Training")

ntthh   = len(df_tthh)
nttbb   = len(df_ttbb)
nttbbbb = len(df_ttbbbb)
nttbbcc = len(df_ttbbcc)

ntrain = min(ntthh, nttbb, nttbbbb, nttbbcc)
print(ntrain)

df_tthh   = df_tthh.sample(n=ntrain).reset_index(drop=True)
df_ttbb   = df_ttbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbbb = df_ttbbbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbcc = df_ttbbcc.sample(n=ntrain).reset_index(drop=True)

df_tthh["category"]   = 0
df_ttbb["category"]   = 1
df_ttbbbb["category"] = 2
df_ttbbcc["category"] = 3

pd_data = pd.concat([df_tthh, df_ttbb, df_ttbbbb, df_ttbbcc])
colnames = pd_data.columns
print(pd_data.head())
print("Col names:",colnames)

print("Plotting corr_matrix total")
plot_corrMatrix(pd_data,outdir,"total")
print("Plotting corr_matrix tthh")
plot_corrMatrix(df_tthh,outdir,"tthh")
print("Plotting corr_matrix ttbb")
plot_corrMatrix(df_ttbb,outdir,"ttbb")
print("Plotting corr_matrix ttbbbb")
plot_corrMatrix(df_ttbbbb,outdir,"ttbbbb")
print("Plotting corr_matrix ttbbcc")
plot_corrMatrix(df_ttbbcc,outdir,"ttbbcc")

pd_data = pd_data.sample(frac=1).reset_index(drop=True)

print(pd_data.head())

x_total = np.array(pd_data.filter(items = inputvars))
y_total = np.array(pd_data.filter(items = ['category']))

# Splitting between training set and cross-validation set
numbertr = len(y_total)
trainlen = int(0.7*numbertr) # Fraction used for training

x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)
print(len(x_train),len(x_val),len(y_train),len(y_val))

patience_epoch = 10
# Early Stopping with Validation Loss for Best Model
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)
print("xtrain shape:",x_train.shape)



###################################################
#                      Model                      #
###################################################
activation_function='relu'
weight_initializer = 'random_normal'

model = tf.keras.models.Sequential()
###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
###############    Hidden Layer     ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(len(class_names), activation="softmax"))
###################################################


batch_size = 1024
print("batch size :", batch_size)
model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), loss="sparse_categorical_crossentropy", metrics = ["accuracy"])

model.summary()
hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])


pred_train = model.predict_classes(x_train)
print("pred_train", pred_train)
print("orig train", y_train.T)
#train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[1]]).T, columns=["True", "Pred"])
pred_val = model.predict_classes(x_val)
print("pred_val", pred_val)
print("orig train", y_val.T)
print("conf matrix on train set ")
print(confusion_matrix(y_train, pred_train))
print("conf matrix on val set ")
print(confusion_matrix(y_val, pred_val))

plot_confusion_matrix(y_val, pred_val, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

pred_val = model.predict(x_val)
pred_train = model.predict(x_train)

exit()
print(pred_val)
print(pred_val.T)
print(y_val)
print(y_val.T)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[1]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[1]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="tt",savedir=outdir)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[2]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[2]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="st",savedir=outdir)
