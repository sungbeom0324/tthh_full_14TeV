import os
import sys
import time
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from utils.plots import *
from matplotlib.backends.backend_pdf import PdfPages
from utils.var_functions import *
from utils.drawHistoModules import *
import ROOT
from array import array
import shap

start_time = time.time()

# Setting up directories and process names
indir = "./samples1/"; PRE = "SSDL_SEMI_0817"
outdir = "./Transformer_result/" + PRE + "/LetsFind_tthh/bCat_higgs5_2Mat/"
os.makedirs(outdir, exist_ok=True)
process_names = ["G1", "G2", "G3", "G4"]

# Define the inputs
input_1 = ["bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_m",
           "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_m",
           "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_m",
           "b1b2_dr", "b1b3_dr", "b2b3_dr"]

input_2 = ["bJet_size", "Lep1_pt", "Lep1_eta", "Lep1_phi",
           "Lep2_pt", "Lep2_eta", "Lep2_phi", "MET_E"]

input_transformer = input_1 + input_2

# Load data
df_tthh = uproot.open(indir+PRE+"_tthh.root")["Delphes"].arrays(input_transformer,library="pd")
df_tth = uproot.open(indir+PRE+"_tth.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttzh = uproot.open(indir+PRE+"_ttzh.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttbbh = uproot.open(indir+PRE+"_ttbbh.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttvv = uproot.open(indir+PRE+"_ttvv.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttbbv = uproot.open(indir+PRE+"_ttvv.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttbb = uproot.open(indir+PRE+"_ttbb.root")["Delphes"].arrays(input_transformer,library="pd")
df_ttbbbb = uproot.open(indir+PRE+"_ttbbbb.root")["Delphes"].arrays(input_transformer,library="pd")
df_tttt = uproot.open(indir+PRE+"_tttt.root")["Delphes"].arrays(input_transformer,library="pd")

df_tthh["category"] = 0; df_tthh["process"] = 0
df_tth["category"] = 1; df_tth["process"] = 1
df_ttzh["category"] = 1; df_ttzh["process"] = 2
df_ttbbh["category"] = 1; df_ttbbh["process"] = 3
df_ttvv["category"] = 2; df_ttvv["process"] = 4
df_ttbbv["category"] = 2; df_ttbbv["process"] = 5
df_ttbb["category"] = 2; df_ttbb["process"] = 6
df_ttbbbb["category"] = 2; df_ttbbbb["process"] = 7
df_tttt["category"] = 3; df_tttt["process"] = 8

# Concatenate into groups and balance
df_1 = pd.concat([df_tthh], ignore_index=True)
df_2 = pd.concat([df_tth, df_ttzh, df_ttbbh], ignore_index=True)
df_3 = pd.concat([df_ttvv, df_ttbbv, df_ttbb, df_ttbbbb], ignore_index=True)
df_4 = pd.concat([df_tttt], ignore_index=True)

ntrain = min(len(df_1), len(df_2), len(df_3), len(df_4))
df_1 = df_1.sample(n=ntrain).reset_index(drop=True)
df_2 = df_2.sample(n=ntrain).reset_index(drop=True)
df_3 = df_3.sample(n=ntrain).reset_index(drop=True)
df_4 = df_4.sample(n=ntrain).reset_index(drop=True)

df_total = pd.concat([df_1, df_2, df_3, df_4])
df_total = df_total.sample(frac=1).reset_index(drop=True)

x_total = np.array(df_total.filter(items=input_transformer))
y_total = np.array(df_total.filter(items=["category"]))
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)

# Define Transformer model
def transformer_encoder(inputs, head_size, num_heads, ff_dim, dropout=0):
    # Normalization and Attention
    x = tf.keras.layers.LayerNormalization(epsilon=1e-6)(inputs)
    x = tf.keras.layers.MultiHeadAttention(
        key_dim=head_size, num_heads=num_heads, dropout=dropout
    )(x, x)
    x = tf.keras.layers.Dropout(dropout)(x)
    res = x + inputs

    # Feed Forward Part
    x = tf.keras.layers.LayerNormalization(epsilon=1e-6)(res)
    x = tf.keras.layers.Conv1D(filters=ff_dim, kernel_size=1, activation="relu")(x)
    x = tf.keras.layers.Dropout(dropout)(x)
    x = tf.keras.layers.Conv1D(filters=inputs.shape[-1], kernel_size=1)(x)
    return x + res

def build_model(input_shape, head_size, num_heads, ff_dim, num_transformer_blocks, mlp_units, dropout=0, mlp_dropout=0):
    inputs = tf.keras.Input(shape=input_shape)
    x = inputs
    for _ in range(num_transformer_blocks):
        x = transformer_encoder(x, head_size, num_heads, ff_dim, dropout)

    x = tf.keras.layers.GlobalAveragePooling1D(data_format="channels_first")(x)
    for dim in mlp_units:
        x = tf.keras.layers.Dense(dim, activation="relu")(x)
        x = tf.keras.layers.Dropout(mlp_dropout)(x)
    outputs = tf.keras.layers.Dense(len(process_names), activation="softmax")(x)
    return tf.keras.Model(inputs, outputs)

input_shape = x_train.shape[1:]
model = build_model(
    input_shape,
    head_size=256,
    num_heads=4,
    ff_dim=4,
    num_transformer_blocks=4,
    mlp_units=[128],
    dropout=0.25,
    mlp_dropout=0.4,
)

model.compile(
    loss="sparse_categorical_crossentropy",
    optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4),
    metrics=["sparse_categorical_accuracy"],
)

model.summary()

# Train model
epochs = 100
es = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=20)
mc = tf.keras.callbacks.ModelCheckpoint(outdir+'/best_transformer_model.h5', monitor='val_loss', mode='min', save_best_only=True)

history = model.fit(
    x_train, y_train,
    validation_data=(x_val, y_val),
    epochs=epochs,
    batch_size=256,
    callbacks=[es, mc]
)

# Prediction and evaluation
pred_train = model.predict(x_train)
pred_val = model.predict(x_val)

train_result = pd.DataFrame(np.array([y_train.T[0], np.argmax(pred_train, axis=1)]).T, columns=["True", "Pred"])
val_result = pd.DataFrame(np.array([y_val.T[0], np.argmax(pred_val, axis=1)]).T, columns=["True", "Pred"])

# Confusion Matrix
plot_confusion_matrix(y_val, np.argmax(pred_val, axis=1), classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, np.argmax(pred_val, axis=1), classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, np.argmax(pred_train, axis=1), classes=process_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, np.argmax(pred_train, axis=1), classes=process_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

# Accuracy
train_results = model.evaluate(x_train, y_train)
train_acc = train_results[1]
print(f"Train accuracy: {train_acc * 100:.2f}%")
val_results = model.evaluate(x_val, y_val)
val_acc = val_results[1]
print(f"Validation accuracy: {val_acc * 100:.2f}%")

# SHAP Feature Importance
print("# SHAP Feature Importance #")
explainer = shap.KernelExplainer(model.predict, x_train[:100])
shap_values = explainer.shap_values(x_train[:100])
shap.summary_plot(shap_values, x_train, plot_type='bar', max_display=50, feature_names=input_transformer, show=False)
plt.savefig(outdir+'/shap_summary_plot_transformer.pdf')

# Write TTRee
df_val = pd.DataFrame()
df_val["category"] = y_val.T[0]
df_val["G1"] = pred_val.T[0]
df_val["G2"] = pred_val.T[1]
df_val["G3"] = pred_val.T[2]
df_val["G4"] = pred_val.T[3]
df_val["DNN"] = df_val["G1"] / (df_val["G2"] + df_val["G3"] + df_val["G4"])

f = ROOT.TFile(outdir+"/transformer_result.root", "RECREATE")
tree = ROOT.TTree("Delphes", "Example Tree")

# Empty array
category_array = array('f', [0])
G1_array = array('f', [0])
G2_array = array('f', [0])
G3_array = array('f', [0])
G4_array = array('f', [0])
DNN_array = array('f', [0])

# Attach array to TTree
tree.Branch('category', category_array, 'category/F')
tree.Branch('G1', G1_array, 'G1/F')
tree.Branch('G2', G2_array, 'G2/F')
tree.Branch('G3', G3_array, 'G3/F')
tree.Branch('G4', G4_array, 'G4/F')
tree.Branch('DNN', DNN_array, 'DNN/F')

for _, row in df_val.iterrows():
    category_array[0] = row['category']
    G1_array[0] = row['G1']
    G2_array[0] = row['G2']
    G3_array[0] = row['G3']
    G4_array[0] = row['G4']
    DNN_array[0] = row['DNN']
    tree.Fill()

f.Write()
f.Close()

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
print("---Complete---")

