# Build functions make additional variables with DNN result, "pred_bfh"
import os
import numpy as np
import matplotlib.pyplot as plt

def plot_histogram(data, title, x_label, y_label, save_dir):
    hist_values, bin_edges = np.histogram(data, bins=20)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]
    errors = np.sqrt(hist_values)
    
    plt.bar(bin_centers, hist_values, width=bin_width, color='blue', alpha=0.7)
    plt.errorbar(bin_centers, hist_values, yerr=errors, fmt='d', color='red', markersize=6, capsize=5)
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_name = title.replace(" ", "_").lower() + ".pdf"
    save_path = os.path.join(save_dir, file_name)
    plt.savefig(save_path, format='pdf')
    plt.close()

def plot_2histogram(df, col1, col2, title, x_label, y_label, save_dir):
    data1 = df[col1]
    data2 = df[col2]

    hist_values1, bin_edges1 = np.histogram(data1, bins=20)
    bin_centers1 = 0.5 * (bin_edges1[:-1] + bin_edges1[1:])
    bin_width1 = bin_edges1[1] - bin_edges1[0]
    errors1 = np.sqrt(hist_values1)

    hist_values2, bin_edges2 = np.histogram(data2, bins=20)
    bin_centers2 = 0.5 * (bin_edges2[:-1] + bin_edges2[1:])
    bin_width2 = bin_edges2[1] - bin_edges2[0]
    errors2 = np.sqrt(hist_values2)

    plt.bar(bin_centers1, hist_values1, width=bin_width1, color='blue', alpha=0.5, label=col1)
    plt.bar(bin_centers2, hist_values2, width=bin_width2, color='green', alpha=0.5, label=col2)
    plt.errorbar(bin_centers1, hist_values1, yerr=errors1, fmt='none', color='blue', markersize=6, capsize=5)
    plt.errorbar(bin_centers2, hist_values2, yerr=errors2, fmt='none', color='green', markersize=6, capsize=5)

    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    file_name = title.replace(" ", "_").lower() + ".pdf"
    save_path = os.path.join(save_dir, file_name)
    plt.savefig(save_path, format='pdf')
    plt.close()

def plot_multi_histograms(dataframe, histogram_name, x_label, y_label, outdir, *args):
    plt.figure(figsize=(10, 6))

    for arg in args:
        column_name, label, color, histtype = arg
        data = dataframe[column_name]
        data = data[data > 0]

        n, bins, patches = plt.hist(data, bins=30, density=True, alpha=0.5, label=label, color=color, histtype=histtype)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_width = bins[1] - bins[0]
        bin_heights = n / (len(data) * bin_width)
        errors = np.sqrt(n) / (len(data) * bin_width)
        plt.errorbar(bin_centers, bin_heights, yerr=errors, fmt='none', color=color)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(histogram_name)
    plt.legend()
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    save_filename = f'{histogram_name.replace(" ", "_")}.pdf'
    save_path = os.path.join(outdir, save_filename)
    plt.savefig(save_path)
    plt.show()

def plot_multi_histograms_dfs(column_name, histogram_name, x_label, y_label, outdir, *args):
    plt.figure(figsize=(10, 6))

    for dataframe, label, color, histtype in args:
        data = dataframe[column_name]
        data = data[data > 0]

        n, bins, patches = plt.hist(data, bins=30, density=True, alpha=0.5, label=label, color=color, histtype=histtype)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_width = bins[1] - bins[0]
        bin_heights = n / (len(data) * bin_width)
        errors = np.sqrt(n) / (len(data) * bin_width)
        plt.errorbar(bin_centers, bin_heights, yerr=errors, fmt='none', color=color)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(histogram_name)
    plt.legend()

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    save_filename = f'{histogram_name.replace(" ", "_")}.pdf'
    save_path = os.path.join(outdir, save_filename)
    plt.savefig(save_path)
    plt.show()

def higgs5_3(row):
    cat = row["pred_bfh"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m", "bJet5_m"]
    arrow = {0:[0,1,2], 1:[0,1,3], 2:[0,1,4], 3:[0,2,3], 4:[0,2,4], 5:[0,3,4], 6:[1,2,3], 7:[1,2,4], 8:[1,3,4], 9:[2,3,4]}
    if cat == 10: return 0;

    else: 
        pt1, pt2, pt3 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]], row[pts[arrow[cat][2]]]
        eta1, eta2, eta3 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]], row[etas[arrow[cat][2]]]
        phi1, phi2, phi3 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]], row[phis[arrow[cat][2]]]
        m1, m2, m3 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]], row[ms[arrow[cat][2]]]
   
        px1, px2, px3 = pt1 * np.cos(phi1), pt2 * np.cos(phi2), pt3 * np.cos(phi3)
        py1, py2, py3 = pt1 * np.sin(phi1), pt2 * np.sin(phi2), pt3 * np.sin(phi3)
        pz1, pz2, pz3 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2), pt3 * np.sinh(eta3)
        E1,E2,E3=np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2),np.sqrt(pt3**2+m3**2)*np.cosh(eta3)
        pxs = [px1, px2, px3]; pys = [py1, py2, py3]; pzs = [pz1, pz2, pz3]; Es = [E1, E2, E3]

    close_mass = 999
    for i in range(len(pxs)):
        px_1, py_1, pz_1, E_1 = pxs[i], pys[i], pzs[i], Es[i]
        for j in range(i+1, len(pxs)):
            px_2, py_2, pz_2, E_2 = pxs[j], pys[j], pzs[j], Es[j]

            px = px_1 + px_2 
            py = py_1 + py_2 
            pz = pz_1 + pz_2 
            E = E_1 + E_2

            mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
            if abs(mass-125) < abs(close_mass-125): 
                close_mass = mass

    if close_mass != 999:
        return close_mass
    else: return 0


def calculate_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1, E2 =np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
    pxs = [px1, px2]; pys = [py1, py2]; pzs = [pz1, pz2]; Es = [E1, E2]

    px_1, py_1, pz_1, E_1 = pxs[0], pys[0], pzs[0], Es[0]
    px_2, py_2, pz_2, E_2 = pxs[1], pys[1], pzs[1], Es[1]

    px = px_1 + px_2
    py = py_1 + py_2
    pz = pz_1 + pz_2
    E = E_1 + E_2

    mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    return mass


def higgs_5_2(row):
    cat = row["pred_bfh"]
    bfts = ["bft_1", "bft_2", "bft_3", "bft_4", "bft_5"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m", "bJet5_m"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[0,4], 4:[1,2], 5:[1,3], 6:[1,4], 7:[2,3], 8:[2,4], 9:[3,4]}
    arrow2 = {0:[2,3,4], 1:[1,3,4], 2:[1,2,4], 3:[1,2,3], 4:[0,3,4], 5:[0,2,4], 6:[0,2,3], 7:[0,1,4], 8:[0,1,3], 9:[0,1,2]}
    if cat == 10: return [0, 0, 0]; # modify output shape #

    # higgs_1
    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
    higgs_mass_1 = calculate_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)

    # higgs_2
    bfts = np.array([ row[bfts[arrow2[cat][0]]], row[bfts[arrow2[cat][1]]], row[bfts[arrow2[cat][2]]] ])   
    bft_idx = np.argmax(bfts)
    cand = arrow2[cat]
    del cand[bft_idx]
     
    pt1, pt2 = row[pts[cand[0]]], row[pts[cand[1]]]
    eta1, eta2 = row[etas[cand[0]]], row[etas[cand[1]]]
    phi1, phi2 = row[phis[cand[0]]], row[phis[cand[1]]]
    m1, m2 = row[ms[cand[0]]], row[ms[cand[1]]]
    higgs_mass_2 = calculate_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)

    # higgs 12
    higgs_mass_12 = higgs_mass_1 + higgs_mass_2 
 
    out = [higgs_mass_1, higgs_mass_2, higgs_mass_12]
    return out   


def Second_higgs_5_2(row):
    cat = row["pred_bfh"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m", "bJet5_m"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[0,4], 4:[1,2], 5:[1,3], 6:[1,4], 7:[2,3], 8:[2,4], 9:[3,4]}
    if cat == 10: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1, E2 =np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
    pxs = [px1, px2]; pys = [py1, py2]; pzs = [pz1, pz2]; Es = [E1, E2]

    close_mass = 999
    px_1, py_1, pz_1, E_1 = pxs[0], pys[0], pzs[0], Es[0]
    px_2, py_2, pz_2, E_2 = pxs[1], pys[1], pzs[1], Es[1]

    px = px_1 + px_2 
    py = py_1 + py_2 
    pz = pz_1 + pz_2 
    E = E_1 + E_2

    mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    if abs(mass-125) < abs(close_mass-125): 
        close_mass = mass

    if close_mass != 999:
        return close_mass
    else: return 0

def X_higgs(row):
    HiggsMass = 125.0; HiggsWidth = 4.07;
    mass = row["higgs_mass"]
    X_Higgs = ((HiggsMass - mass)/HiggsWidth)**2
    return X_Higgs

def bb_dr(row):
    cat = row["pred_bfh"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m", "bJet5_m"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[1,2], 4:[1,3], 5:[2,3], 6:[]}
    if cat == 6: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
    
    out = np.sqrt( (eta1-eta2)**2 + (phi1-phi2)**2 )
    return out

def bfh_Vars(row):
    cat = row["pred_bfh"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m", "bJet5_m"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[0,4], 4:[1,2], 5:[1,3], 6:[1,4], 7:[2,3], 8:[2,4], 9:[3,4]}
    if cat == 10: return [0, 0, 0, 0, 0];

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    bb_dr = np.sqrt( (eta1-eta2)**2 + (phi1-phi2)**2 )
    bb_Ht = pt1 + pt2
    bb_dEta = abs(eta1-eta2)
    bb_dPhi = abs(phi1-phi2)
    bb_mbmb = m1+m2
    out = [bb_dr, bb_Ht, bb_dEta, bb_dPhi, bb_mbmb]
    return out


def higgs_4(row):
    cat = row["pred_bfh"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_m", "bJet2_m", "bJet3_m", "bJet4_m"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[1,2], 4:[1,3], 5:[2,3], 6:[]}
    if cat == 6: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1,E2 = np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
    px = px1+px2; py = py1+py2; pz = pz1+pz2; E = E1+E2
    mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    return mass
