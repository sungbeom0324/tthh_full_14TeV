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
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
import json

###################################################
#                     I/O                         #
###################################################
indir = "./skimmed/"; PRE = "b" # Should be apart from data for event classification.
outdir = "./dnn_result/" + PRE + "ft/" # MODIFY #
os.makedirs(outdir, exist_ok=True)
class_names = ["Cat1", "Cat2", "Cat3", "Cat4", "Cat5", "NoCat"]

with open('./dnn/dnn_input.json', 'r') as file:
    data = json.load(file)

input_1 = data["input_1"]

input_cat = ["bCat_top_1"] # MODIFY #
openvars = input_1 + input_cat

###################################################
#                 PreProcessing                   #
###################################################
pd_data = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(openvars,library="pd")
pd_cat1 = pd_data.loc[pd_data["bCat_top_1"] == 0]
pd_cat2 = pd_data.loc[pd_data["bCat_top_1"] == 1]
pd_cat3 = pd_data.loc[pd_data["bCat_top_1"] == 2]
pd_cat4 = pd_data.loc[pd_data["bCat_top_1"] == 3]
pd_cat5 = pd_data.loc[pd_data["bCat_top_1"] == 4]
pd_cat6 = pd_data.loc[pd_data["bCat_top_1"] == 5]

nCat1 = len(pd_cat1) 
nCat2 = len(pd_cat2)
nCat3 = len(pd_cat3)
nCat4 = len(pd_cat4)
nCat5 = len(pd_cat5)
nCat6 = len(pd_cat6)
ntrain = min(nCat1, nCat2, nCat3, nCat4, nCat5, nCat6)
print("ntrain = ", ntrain)

pd_cat1 = pd_cat1.sample(n=ntrain).reset_index(drop=True)
pd_cat2 = pd_cat2.sample(n=ntrain).reset_index(drop=True)
pd_cat3 = pd_cat3.sample(n=ntrain).reset_index(drop=True)
pd_cat4 = pd_cat4.sample(n=ntrain).reset_index(drop=True)
pd_cat5 = pd_cat5.sample(n=ntrain).reset_index(drop=True)
pd_cat6 = pd_cat6.sample(n=ntrain).reset_index(drop=True)
print("pd_cat1", pd_cat1)

pd_data = pd.concat([pd_cat1, pd_cat2, pd_cat3, pd_cat4, pd_cat5, pd_cat6])
pd_data = pd_data.sample(frac=1).reset_index(drop=True)
print("pd_data", pd_data)
x_total = np.array(pd_data.filter(items = input_1))
y_total = np.array(pd_data.filter(items = ['bCat_top_1'])) # MODIFY #

# Training and Cross-Validation Set
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)
print("x_train: ",len(x_train),"x_val: ", len(x_val),"y_train: ", len(y_train),"y_val", len(y_val))

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 30; batch_size = 512; print("batch size :", batch_size)
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
model.add(tf.keras.layers.Dense(20, activation=activation_function))
#model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dropout(0.2))    
model.add(tf.keras.layers.Dense(20, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
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
print("#          FEATURE IMPORTANCE             #")
model_dir = outdir + '/best_model.h5'
plot_feature_importance(model_dir, x_val, input_1, outdir)
 '''
###################################################
#                     Time                        #
###################################################
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
print("---Done---")
