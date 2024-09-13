# AE reconstructing tthh, ttbbbb, ttbbcc, tthbb.  
# It is unsupervised though, extract dixcriminative regions.
import os
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization
from sklearn.preprocessing import MinMaxScaler

# Data
indir = "./samples1/"; PRE = "CF_1223_14TeV"
outdir = "./DNN_result/" + PRE + "/LetsFind_tthh/ae/"    # MODIFY  #
os.makedirs(outdir, exist_ok=True)
process_names = ["tthh", "tthbb", "ttbb", "ttbbbb", "ttbbcc"]
openvars = [
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
#     "MET_E", # why decrease..


    # Defined Kinematic vars
#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
#     "close_Higgs_pt", "close_Higgs_eta", "close_Higgs_phi", "close_Higgs_mass"
            ]
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
ntrain = min(ntthh, ntthbb, nttbb, nttbbbb, nttbbcc)

df_tthh   = df_tthh.sample(n=ntrain).reset_index(drop=True)
df_tthbb  = df_tthbb.sample(n=ntrain).reset_index(drop=True)
df_ttbb   = df_ttbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbbb = df_ttbbbb.sample(n=ntrain).reset_index(drop=True)
df_ttbbcc = df_ttbbcc.sample(n=ntrain).reset_index(drop=True)

# X-Y Partition
scaler = MinMaxScaler()    
df_total = pd.concat([df_tthh, df_tthbb, df_ttbb, df_ttbbbb, df_ttbbcc])
df_total = df_total.sample(frac=1).reset_index(drop=True)
x_total = np.array(df_total.filter(items = openvars))
x_total = scaler.fit_transform(x_total)
y_total = np.array(df_total.filter(items = ["category"]))
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)

print("x_train : ", x_train)
print("y_train : ", y_train)

# Autoencoder Model
input_img = Input(shape=(31,))
encoded = Dense(16, activation='relu')(input_img)
#encoded = Dense(8, activation='relu')(encoded)
encoded = Dense(3, activation='relu')(encoded)  # 3차원 latent space

#decoded = Dense(8, activation='relu')(encoded)
decoded = Dense(16, activation='relu')(encoded)
decoded = Dense(31, activation='sigmoid')(decoded)

autoencoder = Model(input_img, decoded)

# Compile
autoencoder.compile(optimizer='adam', loss='mean_squared_error')

# Fit
autoencoder.fit(x_train.reshape((len(x_train), np.prod(x_train.shape[1:]))),
                x_train.reshape((len(x_train), np.prod(x_train.shape[1:]))),
                epochs=100,
                batch_size=1024,
                shuffle=True,
                validation_data=(x_val.reshape((len(x_val), np.prod(x_val.shape[1:]))),
                                 x_val.reshape((len(x_val), np.prod(x_val.shape[1:])))))

# Encoder
encoder = Model(input_img, encoded)

# Color
colors = ['b', 'g', 'r', 'c', 'brown']

# Visualize in 3D
encoded_imgs = encoder.predict(x_val.reshape((len(x_val), np.prod(x_val.shape[1:]))))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(3):
    indices = np.where(y_val == i)
    ax.scatter(encoded_imgs[indices, 0], encoded_imgs[indices, 1], encoded_imgs[indices, 2], c=colors[i], label=str(i), s=20, alpha=0.9)

ax.legend()
plt.savefig("fig_0110_3D.pdf")

