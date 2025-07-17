import os.path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_modification(mods, dir_name="Default", save_path="./", mod_type="Modifications", seq_name="Sequence"):
    ### deprecated ####
    return
    ##### get positions and mods ############
    modifications = np.array(mods, dtype=int)
    positions = modifications[:, 0].tolist()
    mods = modifications[:, 1].tolist()
    moz = [i for i, n in enumerate(mods) if n > 1]
    mods = [mods[i] for i in moz]
    positions = [positions[i] for i in moz]
    positions_size = 50
    num_plots = math.floor(len(positions)/positions_size)


    #### set color based on mod type #######
    bar_color = "b"
    if mod_type == "Deletions":
        bar_color = "r"
    elif mod_type == "Insertions":
        bar_color = "g"

    #### set figure length #####
    figlen = len(mods) * .3

    if figlen < 6.4:
        figlen = 6.4
    elif figlen > 100:
        figlen = 100

    #### plot modifications #####
    #fig, ax = plt.figure(figsize=(figlen, 4.8)) #6.4, 4.8 default size
    #fig, ax = plt.subplots(num_plots)
    #fig.set_size_inches(10, 50)
    #fig.set_figheight(figlen)
    #fig.suptitle(seq_name +' '+ dir_name +' ' + mod_type, fontsize="xx-large")

    #fig.set_figwidth(10)

    dir_name = save_path + dir_name + '_Modification_Plots'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    curr_pos = 0
    for i in range(num_plots):
        mod_pos = np.array(positions[curr_pos:curr_pos+positions_size])
        mod_num = np.array(mods[curr_pos:curr_pos+positions_size])
        df = pd.DataFrame({'Reference Nucleotide Position': mod_pos,
                           'Number of Modifications': mod_num})
        ax = df.plot.bar(x='Reference Nucleotide Position', y='Number of Modifications', rot=90, figsize=(figlen, 10),
                    color=bar_color, xlabel='Reference Nucleotide Position', ylabel='Number of ' + mod_type,
                    title=seq_name +' ' + mod_type, legend=False)
        # pps = ax.bar(X, np.array(mods[curr_pos:curr_pos+positions_size]),
        #                 color=bar_color, align="center", edgecolor="black", width=.4)
        #ax.bar_label(ax, label_type='edge')
        curr_pos = curr_pos + positions_size
        fig = plt.gcf()
        #plt.showblock=False)
        figname = os.path.join(dir_name + '/' + seq_name + '_' + mod_type + '_' + str(i) + '.png')
        fig.savefig(figname)


def plot_modification_summary(mods, dir_name="Default", save_path="./", \
                              mod_type="Modification Summary", seq_name="Sequence"):
    ##### deprecated for now #########
    return
    ##### get positions and mods ############
    modifications = np.array(mods, dtype=int)
    positions = modifications[:, 0].tolist()
    dels = (modifications[:, 1]/np.linalg.norm(modifications[:, 1])).tolist()
    ins = (modifications[:, 2]/np.linalg.norm(modifications[:,2])).tolist()
    mismatch = (modifications[:, 3]/np.linalg.norm(modifications[:,3])).tolist()
    positions_size = 50
    num_plots = math.floor(len(positions) / positions_size)

    #### set figure length #####
    figlen = len(mods) * .3
    if figlen < 6.4:
        figlen = 6.4
    elif figlen > 6.4:
        figlen = 30

    dir_name = save_path + dir_name + '_Modification_Plots'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    curr_pos = 0
    for i in range(num_plots):
        mod_pos = np.array(positions[curr_pos:curr_pos + positions_size])
        mod_del = np.array(dels[curr_pos:curr_pos + positions_size])
        mod_ins = np.array(ins[curr_pos:curr_pos + positions_size])
        mod_mis = np.array(mismatch[curr_pos:curr_pos + positions_size])

        df = pd.DataFrame({'Reference Nucleotide Position': mod_pos,
                           'Deletions': mod_del, 'Insertions': mod_ins,
                           'Mismatches': mod_mis})
        ax = df.plot.bar(x='Reference Nucleotide Position', rot=90, figsize=(figlen, 10),
                         xlabel='Reference Nucleotide Position', ylabel='Number of ' + mod_type,
                         title=seq_name + ' ' + mod_type, legend=True, ylim=(0,.01))
        # pps = ax.bar(X, np.array(mods[curr_pos:curr_pos+positions_size]),
        #                 color=bar_color, align="center", edgecolor="black", width=.4)
        # ax.bar_label(ax, label_type='edge')
        curr_pos = curr_pos + positions_size
        fig = plt.gcf()
        #plt.showblock=False)
        figname = os.path.join(dir_name + '/' + seq_name + '_' + mod_type + '_' + str(i) + '.png')
        fig.savefig(figname)


    # #### plot modifications #####
    # fig = plt.figure(figsize=(figlen, 4.8)) #6.4, 4.8 default size
    # ax = fig.add_axes([.1, .2, .8, .7])
    # ax.margins(.005, .5)
    # ax.set_xticks(range(0, len(positions)))
    # ax.set_xticklabels(positions)
    # ax.tick_params(direction='out', length=6, rotation=90)
    # ax.set_xlabel('Reference Nucleotide Position', fontsize="xx-large")
    # ax.set_ylabel('Number of ' + mod_type, fontsize="xx-large")
    # ax.set_title(seq_name +' '+ dir_name +' ' + mod_type, fontsize="xx-large")
    # ax.grid(axis='y')
    # # proper ticks
    # X = np.arange(len(mods))
    # ax.plot(X, dels, color='r', label="Deletions")
    # #ax.bar_label(pps, label_type='edge')
    # ax.plot(X, ins, color='g', label="Insertions")
    # ax.plot(X, mismatch, color='b', label="Mismatches")
    # ax.legend()
    # fig = plt.gcf()
    # dir_name = save_path + dir_name + '_Modification_Plots'
    # if not os.path.exists(dir_name):
    #     os.makedirs(dir_name)
    # figname = os.path.join(dir_name + '/' + seq_name + '_' + mod_type + '.png')
    # fig.savefig(figname)
    # plt.tight_layout()
    # #plt.show)

def plot_mismatch(mismatches, seq_name="Sequence", dir_name="Default", save_path="./",  mod_type="Mismatched Bases"):
    #### deprecated for now ####
    return
    ##### get positions and mods ############
    data = np.array(mismatches)
    positions = data[:,0].tolist()
    num_mismatch = data[:,1].tolist()
    nuc_mismatch = data[:,2].tolist()
    ref_seq_label = data[:,3].tolist()

    ref_label = []
    for i in range(len(ref_seq_label)):
        ref_label.append(str(positions[i]) +" "+ str(ref_seq_label[i]))

    print(seq_name)

    #### set figure length #####
    figlen = len(positions) * .3
    if figlen < 6.4:
        figlen = 6.4
    elif figlen > 100:
        figlen = 100

    #### plot modifications #####
    fig = plt.figure(figsize=(figlen, 6))  # 6.4, 4.8 default size
    ax = fig.add_axes([.1, .2, .8, .7])
    ax.margins(.005, .2)
    ax.set_xticks(range(0, len(positions)))
    ax.set_xticklabels(ref_label)
    ax.tick_params(direction='out', length=6, rotation=90)
    ax.set_xlabel('Reference Nucleotide Position', fontsize="xx-large")
    ax.set_ylabel('Number of ' + mod_type, fontsize="xx-large")
    ax.set_title(seq_name + ' ' + mod_type, fontsize="xx-large")
    ax.grid(axis='y')
    # proper ticks
    X = np.arange(len(num_mismatch))
    pps = ax.bar(X, np.array(num_mismatch).astype(float), color="b", align="center", edgecolor="black", width=.2)
    ax.bar_label(pps, nuc_mismatch, label_type='edge')
    #### save plot #####
    fig = plt.gcf()
    dir_name = save_path + dir_name + '_Modification_Plots'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    figname = os.path.join(dir_name + '/' + seq_name + '_Mismatched_Bases'  + '.png')
    fig.savefig(figname)
    #plt.showblock=False)

def plot_average_mod_rate(df, dir_name="Default", save_path="./"):
    plt.rcParams.update({'font.size': 20})
    df = df.drop("Condition", axis=1)
    #x = df.values  # returns a numpy array
    #min_max_scaler = preprocessing.MinMaxScaler()
    #x_scaled = min_max_scaler.fit_transform(x)
    #df = pd.DataFrame(x_scaled)
    df = df.iloc[:,:] * 100
    fig = df.plot(kind="bar", rot=30, title="Overall Modification Rates", figsize=(12, 12))
    #### save plot #####
    fig = plt.gcf()
    dir_name = save_path + dir_name + '_Modification_Plots'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    figname = os.path.join(dir_name + '/Average_Modification_Rate' + '.png')
    fig.savefig(figname)
    #plt.showblock=False)

def plot_average_mod_by_pos_rate(df, dir_name="Default", save_path="./"):
    plt.rcParams.update({'font.size': 20})
    #x = df.values  # returns a numpy array
    #min_max_scaler = preprocessing.MinMaxScaler()
    #x_scaled = min_max_scaler.fit_transform(x)
    #df = pd.DataFrame(x_scaled)
    df = df.iloc[:, :] * 100
    fig = df.plot(kind="line", rot=45, title="Modification Rates by Position", subplots=True, figsize=(12,8),
                  ylim=(0,100), color='purple')
    #fig = df.plot(kind="bar", rot=90, title="Modification Rates by Position", figsize=(12,8), stacked=True)
    #### save plot #####
    fig = plt.gcf()
    dir_name = save_path + dir_name + '_Modification_Plots'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    figname = os.path.join(dir_name + '/Position_Modification_Rate' + '.png')
    fig.savefig(figname)
    #plt.showblock=False)
