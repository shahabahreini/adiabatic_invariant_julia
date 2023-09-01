import pandas as pd
import re
import math
import matplotlib.pyplot as plt
import os
from lib import *
import lib 

# Import argreletextrema
import numpy as np
from scipy.signal import argrelextrema, find_peaks_cwt, find_peaks
from scipy.misc import electrocardiogram
from findpeaks import findpeaks


# ---------------------------------- Config ---------------------------------- #
save_file_name = r"2D_omega"
save_file_extension = ".svg"

is_multi_files = False
target_folder_multi_files = "/csv/1deg/"
plots_folder = "plots"

fpath = "/csv/exact/"

# doesn't need to update the parameters if individual file is used
parameter_dict = {
    "eps": 0.0177,
    "epsphi": 0.0,
    "kappa": 1.0,
    "deltas": 1.0,
    "beta": 90.0,
    "alpha": 30.0,
    "theta": 30.0,
    "time": 100,
}

simulation_time = 100
method = "Feagin14 Method"
# ------------------------------------ --- ----------------------------------- #


def my_plot(results):
    df = results['df']
    plt.plot(df['x'], df['y'], label=r'$V_X$')
    plt.plot(df['x'][df['peak']], df['y'][df['peak']], 'rx', label=r'peak')
    plt.xlabel('Steps (DataPoint Index)')
    plt.ylabel(r'$V_X$')
    plt.legend()
    plt.show()

def peakfinder_(X):
    # Initialize
    fp = findpeaks(method="peakdetect", lookahead=1)
    results = fp.fit(X)
    my_plot(results)

    peak_idx = peak_couter(results)

    # Plot
    fp.plot(xlabel='Steps (DataPoint Index)', ylabel=r'$V_X$')

    fp = findpeaks(method="topology", lookahead=1)
    results = fp.fit(X)
    #fp.plot()
    #fp.plot_persistence()

    return peak_idx


def peak_couter(results):
    df_result = results["df"]

    count_valleys = df_result[df_result["valley"] == True].shape[0]
    count_peaks = df_result[df_result["peak"] == True].shape[0]

    chosen_key = "peak" if count_peaks < count_valleys else "valley"
    peak_indexes = df_result[df_result[chosen_key] == True]["x"]
    peak_indexes = peak_indexes.tolist()

    return peak_indexes

def plotter(path_, fname_):
    global parameter_dict

    df = lib.read_exported_csv_2Dsimulation(path_, fname_ + ".csv")
    varibale_to_find_peaks_with = df["drho"]
    variable_to_show_peaks_on = df["rho"]


    peak_idxx = peakfinder_(varibale_to_find_peaks_with)

    lib.adiabtic_calculator(df["drho"], df["rho"], peak_idxx)

    x_axis_data = [df["timestamp"].tolist()[i] for i in peak_idxx]
    y_axis_data = [i for i in range(1, len(peak_idxx) + 1)]

    # append the end of the simulation to the list
    # x_axis_data.append(x_.tolist()[-1])
    # y_axis_data.append(y_axis_data[-1])

    # get the list of all csv files
    if is_multi_files:
        path_ = os.path.dirname(__file__) + target_folder_multi_files
        filelst = os.listdir(path_)
    else:
        path_ = ""
        filelst = indivisual_file_names_to_read

    parameter_dict = extract_parameters_by_file_name(filelst[0])

    eps = parameter_dict["eps"]
    epsphi = parameter_dict["epsphi"]
    kappa = parameter_dict["kappa"]
    deltas = parameter_dict["deltas"]
    beta = parameter_dict["beta"]
    alpha = parameter_dict["alpha"]
    theta = parameter_dict["theta"]
    simulation_time = parameter_dict["time"]
    
    plt.plot(
        x_axis_data,
        y_axis_data,
        marker="o",
        markerfacecolor="#344e41",
        markersize=3,
        label=f"$\epsilon_\phi = {parameter_dict['epsphi']}$",
    )
    plt.rcParams["figure.dpi"] = 150
    plt.ylabel(r"$\qquad \sum Cycles$")
    plt.xlabel(r"$\tau$")
    plt.suptitle(r"Accumulative Number Of Cycles per Dimenssionless Time ($\tau$)", fontsize=12)
    plt.title(
        #f"$\\theta_0 = {theta}^{{\circ}}$ , $\\alpha_0={alpha}^{{\circ}}$ , $\\beta_0 = {beta}^{{\circ}}$, $\\phi_0 = 0.0^{{\circ}}$, $\\kappa = {kappa}$, $\\delta_* = {deltas}$, $\\epsilon_ = {eps}$",
        f"$\\theta_0 = {theta}^{{\circ}}$ , $\\alpha_0={alpha}^{{\circ}}$ , $\\beta_0 = {beta}^{{\circ}}$, $\\phi_0 = 0.0^{{\circ}}$, $\\delta_* = {deltas}$, $\\epsilon_\phi = {epsphi}$, $\\epsilon = {eps}$",
        loc="right",
        fontsize=8,
        color="grey",
        style="italic",
    )
    plt.title(
        f"{method}, $\\tau = {simulation_time}$",
        loc="left",
        fontsize=8,
        color="grey",
        style="italic",
    )
    
    if is_multi_files:
        plt.legend()

    plt.tight_layout()
    path_to_save = os.path.join(plots_folder, str(save_file_name + save_file_extension))
    plt.savefig(path_to_save, dpi=600)
    path_to_save = os.path.join(plots_folder, str(save_file_name + ".png"))
    plt.savefig(path_to_save, dpi=600)
    plt.show()


if not is_multi_files:
    chosen_csv = search_for_export_csv()
    indivisual_file_names_to_read = [chosen_csv]
else:
    chosen_csv = "multi_plot"

plotter(fpath, chosen_csv.replace(".csv", ""))