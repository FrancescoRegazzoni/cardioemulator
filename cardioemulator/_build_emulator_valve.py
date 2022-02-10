import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

try:
    import matplotlib.pyplot as plt
except:
    pass

from ._emulator_valve import EmulatorValve

def build_emulator_valve(
        file_PV_loops,
        label_pressure_jump = 'p_jump_aortic_valve',
        label_flow_rate = 'QAV',
        label_opening_coefficient = 'opening_coeff_aortic_valve',
        output_file = None,
        verbose = True,
        show_fig = False):
    """
    Build an emulator for a cardiac valve.

    Parameters
    ----------
    file_PV_loops : str
        CSV file containing pressure jumps and flow rates across the valve.
    label_pressure_jump : str
        Label of the variable in the CSV file containing the pressure jump
        across the valve.
    label_flow_rate : str
        Label for the variable in the CSV file containing the flow rate through
        the valve.
    label_opening_coefficient : str
        Label for the variable in the CSV file containing the opening
        coefficient of the valve. A value of 1 indicates a fully open valve, and
        a value of 0 indicates a fully closed valve. Intermediate values
        represent intermediate configurations, and are ignored.
    output_file : str, optional
        Path to a json file to store the dictionary encoding the emulator.
    verbose : bool, optional
        Flag to toggle verbosity.
    show_fig : bool, optional
        Flag to toggle plots visualization.

    Returns
    -------
    emulator : EmulatorValve
        Emulator.

    """

    data_loop = pd.read_csv(file_PV_loops)

    line = lambda x, R : R * x

    # Fit the "valve open" data points to get R_min.
    open_points = data_loop[data_loop[label_opening_coefficient] == 1]
    R_min, _    = curve_fit(line, np.array(open_points[label_flow_rate]), \
        np.array(open_points[label_pressure_jump]), p0 = 0.0075)
    R_min = float(R_min)
    print("R_min = %0.4f mmHg / (ml/s)" % R_min)

    # We set the maximum resistance always to 75000, regardless of data. This
    # prevents spurious effects (due to e.g. bad valve resolution in the RIIS
    # model) from being incorporated in the emulator.
    R_max = 75000

    # Plot.
    if show_fig:
        plt.plot(data_loop[label_pressure_jump],
                 data_loop[label_flow_rate],
                 label = 'sample PQ')
        plt.plot(open_points[label_pressure_jump],
                 open_points[label_pressure_jump] / R_min,
                 label = 'fitted, open')
        plt.plot(closed_points[label_pressure_jump],
                 closed_points[label_pressure_jump] / R_max,
                 label = 'fitted, closed')
        plt.legend();

    # Build the emulator instance and save it to file.
    emulator_valve = dict()
    emulator_valve['R_min'] = R_min
    emulator_valve['R_max'] = R_max
    emulator_instance = EmulatorValve(emulator_valve)

    if output_file is not None:
        emulator_instance.save_file(output_file, verbose = verbose)
