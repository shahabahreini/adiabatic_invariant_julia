import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def search_for_export_csv():
    # get all files in current directory
    files = os.listdir()

    # filter out files with csv extension
    csv_files = [file for file in files if file.endswith(".csv")]

    # print the list of csv files
    print("CSV files in current directory:")
    for i, file in enumerate(csv_files):
        print(f"{i+1}. {file}")

    # ask user to choose a file
    while True:
        choice = input("Choose a file (enter a number from the list): ")
        try:
            choice = int(choice)
            if choice > 0 and choice <= len(csv_files):
                break
        except ValueError:
            pass
        print("Invalid choice, please try again.")

    # print the selected file name
    selected_file = csv_files[choice - 1]
    print(f"Selected file: {selected_file}")

    return selected_file


def extract_parameters_by_file_name(fname):
    numbers = {}

    for match in re.finditer(
        r"(eps|epsphi|kappa|deltas|beta|alpha|theta|time)(\d+\.\d+)", fname
    ):
        key = match.group(1)
        value = float(match.group(2))
        numbers[key] = value

    return numbers


def read_exported_csv_simulation(path_, fname_):
    """Gets the folder path and desired file name and load the data into Pandas DataFrame"""

    data = pd.read_csv(path_ + fname_)
    df = pd.DataFrame(
        data,
        columns=[
            "timestamp",
            "value1",
            "value2",
            "value3",
            "value4",
            "value5",
            "value6",
        ],
    )
    df.rename(
        columns={
            "value1": "dR",
            "value2": "dphi",
            "value3": "dZ",
            "value4": "R",
            "value5": "phi",
            "value6": "Z",
        },
        inplace=True,
    )
    return df


def read_exported_csv_2Dsimulation(path_, fname_):
    """Gets the folder path and desired file name and load the data into Pandas DataFrame"""

    data = pd.read_csv(fname_)
    df = pd.DataFrame(
        data, columns=["timestamp", "omega_rho", "omega_z", "rho", "z", "drho", "dz"]
    )

    return df


def adiabtic_calculator(v_x, x, extremum_idx):
    velocity = v_x
    position = x
    
    # Compute the changes in the components of X
    delta_X = position.diff()

    # Compute the cumulative sum of velocity times delta_X
    adiabatic = np.cumsum(velocity * delta_X)

    # Compute the integral of V.dX between sequential extremum indexes
    integral_VdX = []
    for i in range(len(extremum_idx) - 1):
        start_idx = extremum_idx[i]
        end_idx = extremum_idx[i + 1]
        integral = np.sum(velocity[start_idx:end_idx] * delta_X[start_idx:end_idx])
        integral_VdX.append(integral)
    
    # Plot the integral versus cycles
    plt.plot(range(len(integral_VdX)), integral_VdX)
    plt.xlabel('Cycles')
    plt.ylabel(r'$\oint\, V.\, dX$')
    plt.title('Closed Path Integral Of Radial Velocity per Cycles')
    plt.show()

    return adiabatic
