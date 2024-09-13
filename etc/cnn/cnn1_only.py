import os
import time
os.environ["CUDA_VISIBLE_DEVICES"] = "3"

import uproot
import numpy as np
import tensorflow as tf
from keras.layers import Conv2D, Flatten, Dense, BatchNormalization
from keras.callbacks import EarlyStopping, ModelCheckpoint
import pandas as pd
import matplotlib.pyplot as plt
from utils.plots import *
from sklearn.model_selection import train_test_split

start_time = time.time()
###################################################
#                     I/O                         #
###################################################
indir = "./samples1/"; PRE = "B_1204_14TeV" # Should be apart from data for event classification.
outdir = "./CNN_result/" + PRE + "/bJetCassification/bCat_higgs5_2Mat/" # MODIFY #
os.makedirs(outdir, exist_ok=True)
class_names = ["Cat1", "Cat2", "Cat3", "Cat4", "Cat5", "Cat6", "Cat7", "Cat8", "Cat9", "Cat10", "NoCat"]

inputvars = [
# 30 = reshape(-1_nEvents, 5, *6, 1(Black/White))


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
     "bJet_size",

     # Lepton
#     "Lep_size",
#     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
#     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",
     "MET_E", # why decrease..


     # Defined Kinematic vars
#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
#     "close_Higgs_pt", "close_Higgs_phi", "close_Higgs_mass", "close_Higgs_eta",
        ]

catvar = ["bCat_higgs5_2Mat"] # MODIFY #
openvars = inputvars + catvar

###################################################
#                 PreProcessing                   #
###################################################
pd_data = uproot.open(indir+PRE+"_tthh_di.root")["Delphes"].arrays(openvars,library="pd")
pd_cat1 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 0]
pd_cat2 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 1]
pd_cat3 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 2]
pd_cat4 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 3]
pd_cat5 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 4]
pd_cat6 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 5]
pd_cat7 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 6]
pd_cat8 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 7]
pd_cat9 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 8]
pd_cat10 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 9]
pd_cat11 = pd_data.loc[pd_data["bCat_higgs5_2Mat"] == 10]

nCat1 = len(pd_cat1)
nCat2 = len(pd_cat2)
nCat3 = len(pd_cat3)
nCat4 = len(pd_cat4)
nCat5 = len(pd_cat5)
nCat6 = len(pd_cat6)
nCat7 = len(pd_cat7)
nCat8 = len(pd_cat8)
nCat9 = len(pd_cat9)
nCat10 = len(pd_cat10)
nCat11 = len(pd_cat11)
ntrain = min(nCat1, nCat2, nCat3, nCat4, nCat5, nCat6, nCat7, nCat8, nCat9, nCat10)#, nCat11)
print("ntrain = ", ntrain)

pd_cat1 = pd_cat1.sample(n=ntrain).reset_index(drop=True)
pd_cat2 = pd_cat2.sample(n=ntrain).reset_index(drop=True)
pd_cat3 = pd_cat3.sample(n=ntrain).reset_index(drop=True)
pd_cat4 = pd_cat4.sample(n=ntrain).reset_index(drop=True)
pd_cat5 = pd_cat5.sample(n=ntrain).reset_index(drop=True)
pd_cat6 = pd_cat6.sample(n=ntrain).reset_index(drop=True)
pd_cat7 = pd_cat7.sample(n=ntrain).reset_index(drop=True)
pd_cat8 = pd_cat8.sample(n=ntrain).reset_index(drop=True)
pd_cat9 = pd_cat9.sample(n=ntrain).reset_index(drop=True)
pd_cat10 = pd_cat10.sample(n=ntrain).reset_index(drop=True)
pd_cat11 = pd_cat11.sample(n=ntrain).reset_index(drop=True)
print("pd_cat1", pd_cat1)

pd_data = pd.concat([pd_cat1, pd_cat2, pd_cat3, pd_cat4, pd_cat5, pd_cat6, pd_cat7, pd_cat8, pd_cat9, pd_cat10, pd_cat11])
pd_data = pd_data.sample(frac=1).reset_index(drop=True)
print("pd_data", pd_data)
x_total = np.array(pd_data.filter(items = inputvars))
y_total = np.array(pd_data.filter(items = ['bCat_higgs5_2Mat'])) # MODIFY #

# Training and Cross-Validation Set
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)
x_train = x_train.reshape(-1, 8, 4, 1); x_val = x_val.reshape(-1, 8, 4, 1) # MODIFY as input_shape #
print("x_train: ",len(x_train),"x_val: ", len(x_val),"y_train: ", len(y_train),"y_val", len(y_val))

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch=100; batch_size = 1024; print("batch_size :", batch_size) 
activation_function = tf.nn.relu
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)

model = tf.keras.models.Sequential()
model.add(Conv2D(filters=16, kernel_size=2, strides=1, activation=activation_function, use_bias=False, input_shape=(8, 4, 1)))
model.add(BatchNormalization())    
model.add(Conv2D(filters=32, kernel_size=2, strides=1, activation=activation_function, use_bias=False))
model.add(BatchNormalization())
model.add(Conv2D(filters=64, kernel_size=2, strides=1, activation=activation_function, use_bias=False))
model.add(BatchNormalization())    
model.add(Flatten())
model.add(Dense(50, activation=activation_function))
model.add(BatchNormalization())    
model.add(Dense(11, activation=tf.nn.softmax))
###################################################
model.compile(loss='sparse_categorical_crossentropy',
              optimizer = 'adam',
              metrics=['accuracy', 'sparse_categorical_accuracy'])
model.summary()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, validation_data=(x_val, y_val),
        callbacks=[es,mc])

end_time = time.time()
###################################################
#                  Prediction                     #
###################################################
pred_train = model.predict(x_train); print(pred_train); pred_train_arg = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val); print(pred_val); pred_val_arg = np.argmax(pred_val, axis=1)
train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[0]]).T, columns=["True", "Pred"]) # True0~4,Pred0~1<"tthh"
val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[0]]).T, columns=["True", "Pred"])

###################################################
#                Confusion Matrix                 #
###################################################
print("#           CONFUSION MATRIX             #")
plot_confusion_matrix(y_val, pred_val_arg, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix.pdf")

plot_confusion_matrix(y_val, pred_val_arg, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix.pdf")

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

'''
###################################################
#              Feature Importance                 #
###################################################
print("#          FEATURE IMPORTANCE             #")
model_dir = outdir + '/best_model.h5'
plot_feature_importance(model_dir, x_val.reshape(-1), inputvars, outdir)
'''
###################################################
#                     Time                        #
###################################################
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")
