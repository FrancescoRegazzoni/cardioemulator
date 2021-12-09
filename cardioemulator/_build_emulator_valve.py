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

    data_loop = pd.read_csv(file_PV_loops)

    line = lambda x, R : R * x

    # Fit the "valve open" data points to get R_min.
    open_points = data_loop[data_loop[label_opening_coefficient] == 1]
    R_min, _    = curve_fit(line, np.array(open_points[label_flow_rate]), \
        np.array(open_points[label_pressure_jump]), p0 = 0.0075)
    R_min = float(R_min)
    print("R_min = %0.4f mmHg / (ml/s)" % R_min)

    # Fit the "valve closed" data points to get R_max.
    closed_points = data_loop[data_loop[label_opening_coefficient] == 0]
    # R_max, _      = curve_fit(line, np.array(closed_points[label_flow_rate]), \
    #     np.array(closed_points[label_pressure_jump]), p0 = 75000)
    # R_max = float(R_max)
    # if R_max < 0:
    #     print("Warning: fitted R_max < 0, setting R_max = 75000 mmHg / (ml/s)")
    #     R_max = 75000

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
