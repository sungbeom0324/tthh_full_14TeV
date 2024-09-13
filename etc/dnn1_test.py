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
indir = "./samples1/"; PRE = "TEST_0625"
outdir = "./DNN_result/" + PRE + "/LetsFind_tthh/bCat_higgs5_2Mat"    # modify #
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

catvar = ["bCat_higgs5_2Mat", "bCat_higgs5_2Mat_multi", "bCat_higgs5_2Mat_1", "bCat_higgs5_2Mat_2", "isMatchable"] # modify #
openvars = dnn1vars + catvar

###################################################
#                PreProcessing_1                  #
###################################################
df_tthh = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(openvars,library="pd")
x_bCat  = np.array(df_tthh.filter(items = dnn1vars))

###################################################
#               bJet Classification               #
###################################################
bfh_dir = "DNN_result/B_0620/bJetCassification/bCat_higgs5_2Mat/best_model.h5"  ## modify ##
bfh_model = tf.keras.models.load_model(bfh_dir)
bfh_model.summary()
_pred_bfh = bfh_model.predict(x_bCat); print("pred_bfh : ", _pred_bfh)
pred_bfh = np.argmax(_pred_bfh, axis=1) # (arg 0~10 that with the biggest prob)
print("pred_bfh argmax : ", pred_bfh)

isMatchable_np = np.array(df_tthh[["isMatchable"]])
y_np = np.array(df_tthh[["bCat_higgs5_2Mat"]])
y_np_1 = np.array(df_tthh[["bCat_higgs5_2Mat_1"]])
y_np_2 = np.array(df_tthh[["bCat_higgs5_2Mat_2"]])

###################################################
#                  Write  TTRee                   #
###################################################
# ROOT file & Tree object.
f = ROOT.TFile(outdir+"/bfh_test.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")

# Empty array.
isMatchable_array = array('f', [0])
bCat_higgs5_2Mat_array = array('f', [0])
bCat_higgs5_2Mat_1_array = array('f', [0])
bCat_higgs5_2Mat_2_array = array('f', [0])
#bCat_higgs5_2Mat_array_multi = array('f', [0])
pred_val_array = array('f', [0])

# Attach tree branches with empty array.
tree.Branch('isMatchable', isMatchable_array, 'isMatchable/F')
tree.Branch('bCat_higgs5_2Mat', bCat_higgs5_2Mat_array, 'bCat_higgs5_2Mat/F')
tree.Branch('bCat_higgs5_2Mat_1', bCat_higgs5_2Mat_1_array, 'bCat_higgs5_2Mat_1/F')
tree.Branch('bCat_higgs5_2Mat_2', bCat_higgs5_2Mat_2_array, 'bCat_higgs5_2Mat_2/F')
#tree.Branch('bCat_higgs5_2Mat_multi', bCat_higgs5_2Mat_array_multi, 'bCat_higgs5_2Mat_multi/F')
tree.Branch('pred_val', pred_val_array, 'pred_val/F')

# Fill array with validation data.
for i in range(len(x_bCat)):
    isMatchable_array[0] = isMatchable_np[i]
    bCat_higgs5_2Mat_array[0] = y_np[i]
    bCat_higgs5_2Mat_1_array[0] = y_np_1[i]
    bCat_higgs5_2Mat_2_array[0] = y_np_2[i]
#    bCat_higgs5_2Mat_array_multi[0] = y_val_multi[i]
    pred_val_array[0] = pred_bfh[i]
    tree.Fill()

f.Write()
f.Close()

print(outdir+"/bfh_test.root")




