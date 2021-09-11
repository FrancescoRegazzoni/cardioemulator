import pandas as pd
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit, root_scalar
import json
try:
    import matplotlib.pyplot as plt
except:
    pass
from ._emulator import Emulator

def build_emulator(
        file_PV_loops,
        period,
        file_EDPV = None,
        cycles_to_drop = 0,
        num_cycles = None,
        PV_loops_sampling = 1,
        E_base_min = 2.0,
        E_base_max = 3.0,
        ESPV_E_min = 0.0,
        ESPV_E_max = 10.0,
        pressure_unit_PV_loops = 'mmHg', volume_unit_PV_loops = 'mL',
        pressure_unit_EDPV     = 'mmHg', volume_unit_EDPV     = 'mL',
        label_time = 'time',
        label_pressure_PV_loops = 'pressure', label_volume_PV_loops = 'volume',
        label_pressure_EDPV     = 'pressure', label_volume_EDPV     = 'volume',
        label_activation = None,
        threshold_activation = 1e-2,
        ED_ini = None,
        ED_end = None,
        klotz_An = 28.2,
        klotz_Bn = 2.79,
        activation_times_reinterpolation = False,
        output_file = None,
        output_file_fig = None,
        verbose = True,
        show_fig = False):
    """
    Build an emulator.

    Parameters
    ----------
    file_PV_loops : str
        CSV file containing pressure-volume transients.
    period : float [s]
        Heartbeat period.
    file_EDPV : str, optional
        CSV file containing the end-diastolic pressure volume relationship.
        If this input is provided, the `simulated EDPVR approach` is used;
        otherwise, the `fitted EDPVR approach` is used.
    cycles_to_drop : int, optional
        If greater than 0, the first ``cycles_to_drop`` heartbeats are discarded
        when constructing the emulator.
    num_cycles : int, optional
        If this input is provided, ``num_cycles`` heartbeats are considered to construct
        the emulator. Otherwise, all the provided heartbeat are used.
        Please notice that the total number of heartbeats cannot exceed ``cycles_to_drop``
        + ``num_cycles``.
    PV_loops_sampling : int, optional
        Sampling rate of the CSV file. If greater than 1, only one row every
        ``PV_loops_sampling`` rows is employed.
    E_base_min : float [mmHg/mL], optional
        Lower bound of elastance interval used to detect the end systolic points.
    E_base_max : float [mmHg/mL], optional
        Upper bound of elastance interval used to detect the end systolic points.
    ESPV_E_min : float [mmHg/mL], optional
        Lower bound of end systolic elastance. If this limit is exceeded, a warning is raised.
    ESPV_E_max : float [mmHg/mL], optional
        Upper bound of end systolic elastance. If this limit is exceeded, a warning is raised.
    pressure_unit_PV_loops : str, optional
        Measure unit for pressure in the CSV file ``file_PV_loops``. Allowed values are
        ``"mmHg"``, ``"Pa"``, ``"kPa"``.
    volume_unit_PV_loops : str, optional
        Measure unit for volume in the CSV file ``file_PV_loops``. Allowed values are ``"L"``,
        ``"mL"``, ``"m^3"``.
    pressure_unit_EDPV : str, optional
        Measure unit for pressure in the CSV file ``file_EDPV``. Allowed values are ``"mmHg"``,
        ``"Pa"``, ``"kPa"``.
    volume_unit_EDPV : str, optional
        Measure unit for volume in the CSV file ``file_EDPV``. Allowed values are ``"L"``,
        ``"mL"``, ``"m^3"``.
    label_time : str, optional
        Column label indicating time in the CSV file ``file_PV_loops``.
    label_pressure_PV_loops: str, optional
        Column label indicating pressure in the CSV file ``file_PV_loops``.
    label_volume_PV_loops : str, optional
        Column label indicating volume in the CSV file ``file_PV_loops``.
    label_pressure_EDPV : str, optional
        Column label indicating pressure in the CSV file ``file_EDPV``.
    label_volume_EDPV : str, optional
        Column label indicating volume in the CSV file ``file_EDPV``.
    label_activation : str, optional
        Column label indicating activation in the CSV file ``file_PV_loops``.
        If this input is provided, then the end diastolic phase is detected as those
        time instants for which activation is within the lower ``threshold_activation``
        fraction of the range.
    threshold_activation : float, optional
        Threshold on activation (relative w.r.t. range) defining the end diastolic phase
        (meaningful only if ``label_activation`` input is provided).
    ED_ini : float [s], optional
        Time from the beginning of each heartbeat defining the beginning of the diastolic
        phase (meaningful only if ``label_activation`` input is not provided).
    ED_end : float [s], optional
        Time from the beginning of each heartbeat defining the end of the diastolic
        phase (meaningful only if ``label_activation`` input is not provided).
    klotz_An : float [mmHg], optional
        Klotz's curve An parameter.
    klotz_Bn : float [-], optional
        Klotz's curve Bn parameter.
    activation_times_reinterpolation : bool, optional
        Flag to toggle the reinterpolation of activation function on a different time grid.
    output_file : str, optional
        Path to a json file to store the dictionary encoding the emulator.
    verbose : bool, optional
        Flag to toggle verbosity.
    show_fig : bool, optional
        Flag to toggle plots visualization.

    Returns
    -------
    emulator : Emulator
        Emulator.

    """

    if verbose:
        print('Building 0D emulator:')
        print('   file PV loops:         %s' % file_PV_loops)
        if file_EDPV is None:
            print('   EDPV approach:         fitted')
        else:
            print('   EDPV approach:         simulated')
            print('   file EDPV:             %s' % file_EDPV)
        print('   period:                %0.3f s' % period)

    make_plot = show_fig or output_file_fig is not None
    if make_plot:
        if label_activation is None:
            fig_out, axs_out = plt.subplots(1, 2, figsize = (9, 4))
            axs_PV  = axs_out[0]
            axs_act = axs_out[1]
        else:
            fig_out, axs_out = plt.subplots(2, 2, figsize = (9, 8))
            axs_PV  = axs_out[0,0]
            axs_act = axs_out[0,1]
            axs_T_loops = axs_out[1,0]
            axs_T_time  = axs_out[1,1]

    data_loop = pd.read_csv(file_PV_loops)
    data_loop = data_loop[::PV_loops_sampling]
    convert_pressure(data_loop, label_pressure_PV_loops, pressure_unit_PV_loops)
    convert_volume(data_loop, label_volume_PV_loops, volume_unit_PV_loops)

    if file_EDPV is not None:
        data_EDPV = pd.read_csv(file_EDPV)
        convert_pressure(data_EDPV, label_pressure_EDPV, pressure_unit_EDPV)
        convert_volume(data_EDPV, label_volume_EDPV, volume_unit_EDPV)

    if num_cycles is None:
        num_cycles = int(np.floor(data_loop[label_time].max() / period) - cycles_to_drop)
    if verbose:
        print('   num cycles:            %d' % num_cycles)

    data_loop_used = data_loop[np.logical_and(data_loop[label_time] >= cycles_to_drop * period,
                                              data_loop[label_time] <= (cycles_to_drop + num_cycles) * period)]

    if make_plot:
        axs_PV.plot(data_loop[label_volume_PV_loops],
                    data_loop[label_pressure_PV_loops],
                    label = 'sample PV loops')
        axs_PV.plot(data_loop_used[label_volume_PV_loops],
                    data_loop_used[label_pressure_PV_loops],
                    label = 'sample PV loops (used)')
        if file_EDPV is not None:
            axs_PV.plot(data_EDPV[label_volume_EDPV],
                        data_EDPV[label_pressure_EDPV],
                        label = 'simulated EDPVR')

    ######## EDPV (nonlinear)

    if label_activation is None:
        if ED_ini is None: ED_ini = period * 7.0/8.0
        if ED_end is None: ED_end = period
        idx_deactivated = np.logical_or.reduce([np.logical_and(
            data_loop_used[label_time] >= (cycles_to_drop + n_cycle) * period + ED_ini,
            data_loop_used[label_time] <= (cycles_to_drop + n_cycle) * period + ED_end)
            for n_cycle in range(num_cycles)])
    else:
        Ta = data_loop_used[label_activation]
        Ta_min = Ta.min()
        Ta_max = Ta.max()
        idx_deactivated = Ta < Ta_min + threshold_activation * (Ta_max - Ta_min)

    EDPV_V_points = data_loop_used[label_volume_PV_loops][idx_deactivated]
    EDPV_p_points = data_loop_used[label_pressure_PV_loops][idx_deactivated]

    if file_EDPV is not None:
        p_EDPV = np.array(data_EDPV[label_pressure_EDPV])
        V_EDPV = np.array(data_EDPV[label_volume_EDPV])

        EDPV = interpolate.interp1d(V_EDPV, p_EDPV, fill_value = 'extrapolate')
    else:
        klotz = lambda V, V0, V30: klotz_An * abs((V - V0) / (V30 - V0))**klotz_Bn

        p_opt, _ = curve_fit(klotz, EDPV_V_points, EDPV_p_points, p0 = (30, 200))
        EDPV_V0, EDPV_V30 = p_opt
        EDPV = lambda V: klotz(V, EDPV_V0, EDPV_V30)

        if verbose:
            print('   EDPV V_0:              %0.2f mL' % EDPV_V0)
            print('   EDPV V_30:             %0.2f mL' % EDPV_V30)

        if make_plot:
            axs_PV.plot(data_loop_used[label_volume_PV_loops][idx_deactivated],
                        data_loop_used[label_pressure_PV_loops][idx_deactivated],
                        '+', label = 'ED points')
            if label_activation is not None:
                axs_T_loops.plot(data_loop_used[label_time], Ta, label = 'activation')
                axs_T_loops.plot(data_loop_used[label_time][idx_deactivated], Ta[idx_deactivated],
                                '+', label = 'ED points')
                axs_T_loops.set_title(label_activation)
                axs_T_loops.set_xlabel('time [s]')
                axs_T_loops.legend(loc = 'upper left', frameon = False, fontsize = 7, labelspacing = 0.2)

    ######## ESPV

    p_list = []
    V_list = []
    if abs(E_base_min - E_base_max) < 1e-12:
        E_base_list = [E_base_min]
    else:
        E_base_list = np.linspace(E_base_min, E_base_max, 10)
    for i in range(num_cycles):
        data_cycle = data_loop_used[np.logical_and(data_loop_used[label_time] >  (cycles_to_drop + i    ) * period,
                                                   data_loop_used[label_time] <= (cycles_to_drop + i + 1) * period)]
        p_cycle = np.array(data_cycle[label_pressure_PV_loops])
        V_cycle = np.array(data_cycle[label_volume_PV_loops])
        for E_base in E_base_list:
            idx = np.argmax(p_cycle / E_base - V_cycle)
            p_list.append(p_cycle[idx])
            V_list.append(V_cycle[idx])

    if make_plot:
        axs_PV.plot(V_list, p_list, 'o', label = 'ES points')

    coeffs = np.polyfit(V_list, p_list, 1)

    ESPV_V0 = -coeffs[1]/coeffs[0]
    ESPV_E  = coeffs[0]

    if verbose:
        print('   ESPV V_0:              %0.2f mL' % ESPV_V0)
        print('   ESPV elastance:        %0.3f mmHg/mL' % ESPV_E)

    if ESPV_E < ESPV_E_min or ESPV_E > ESPV_E_max:
        if verbose:
            print('*** ESPV elastance out of bounds! Resorting to EDPV V_0...')
        if file_EDPV is not None:
            root_result = root_scalar(EDPV, x0 = 10.0, x1 = 30.0)
            if root_result.converged:
                EDPV_V0 = root_result.root
                if verbose:
                    print('*** EDPV V_0 search did not converge.')
            else:
                EDPV_V0 = 0.0
        ESPV_V0 = EDPV_V0
        ESPV_E = np.mean(p_list / (V_list - ESPV_V0))

    ESPV = lambda V : ESPV_E * (V - ESPV_V0)

    ######## EDPV (linear)

    EDPV_E = np.mean(EDPV_p_points / (EDPV_V_points - ESPV_V0))
    EDPV_linear = lambda V : EDPV_E * (V - ESPV_V0)

    if verbose:
        print('   EDPV elastance:        %0.3f mmHg/mL' % EDPV_E)

    if label_activation is None:
        time_C = 0.0
        duration_C = 0.0
        duration_R = 0.0
    else:
        data_cycle = data_loop_used[np.logical_and(data_loop_used[label_time] >  (cycles_to_drop    ) * period,
                                                data_loop_used[label_time] <= (cycles_to_drop + 1) * period)]

        t_peak = np.array(data_cycle[label_time])[np.argmax(data_cycle[label_activation])]

        data_C = data_cycle[data_cycle[label_time] <= t_peak]
        data_R = data_cycle[data_cycle[label_time] > t_peak]

        Ta_min_C = data_C[label_activation].min()
        Ta_min_R = data_R[label_activation].min()
        t_C = np.array(data_C[label_time])[np.argmax(data_C[label_activation] > Ta_min_C + threshold_activation * (Ta_max - Ta_min_C))]
        t_R = np.array(data_R[label_time])[np.argmax(data_R[label_activation] < Ta_min_R + threshold_activation * (Ta_max - Ta_min_R))]

        time_C = t_C - cycles_to_drop * period
        duration_C = t_peak - t_C
        duration_R = t_R - t_peak
        if verbose:
            print('   time contraction:      %0.3f mmHg/mL' % time_C)
            print('   duration contraction:  %0.3f mmHg/mL' % duration_C)
            print('   duration relaxation:   %0.3f mmHg/mL' % duration_R)

    if make_plot:
        V_plot = np.linspace(0, 200, 1000)
        axs_PV.plot(V_plot, EDPV(V_plot), label = 'EDPVR')
        axs_PV.plot(V_plot, ESPV(V_plot), label = 'ESPVR')

        if label_activation is not None:
            axs_T_time.plot(data_cycle[label_time], data_cycle[label_activation])
            axs_T_time.axvline(t_C   , color = 'k', linestyle = '--')
            axs_T_time.axvline(t_peak, color = 'k', linestyle = '--')
            axs_T_time.axvline(t_R   , color = 'k', linestyle = '--')
            axs_T_time.set_title(label_activation)
            axs_T_time.set_xlabel('time [s]')

    ######## Activation kinetics

    PV = lambda V, phi: (1-phi) * EDPV(V) + phi * ESPV(V)

    if not activation_times_reinterpolation:
        tt = data_loop[label_time]
        ff = (data_loop[label_pressure_PV_loops] - EDPV(data_loop[label_volume_PV_loops])) / (ESPV(data_loop[label_volume_PV_loops]) - EDPV(data_loop[label_volume_PV_loops]))
    else:
        tt = period * np.linspace(cycles_to_drop, cycles_to_drop + num_cycles, num_cycles * 500)
        p_curr_func = interpolate.interp1d(data_loop[label_time], data_loop[label_pressure_PV_loops])
        V_curr_func = interpolate.interp1d(data_loop[label_time], data_loop[label_volume_PV_loops])
        phi_vals = np.linspace(0, 1, 1000)

        ff = list()
        for t in tt:
            p_curr = p_curr_func(t)
            V_curr = V_curr_func(t)
            p_theor = PV(V_curr, phi_vals)
            if p_curr < np.min(p_theor):
                phi = 0.0
            elif p_curr > np.max(p_theor):
                phi = 1.0
            else:
                phi = interpolate.interp1d(p_theor, phi_vals)(p_curr)
            ff.append(phi)
        ff = np.array(ff)

    activation_t = np.linspace(0, period, 1001)
    activation_v = interpolate.interp1d(tt, ff)(np.clip(activation_t + (cycles_to_drop + num_cycles - 1) * period, tt.min(), tt.max()))

    activation_t[-1] = period + 1e-6

    if make_plot:
        for i in range(num_cycles):
            axs_act.plot(tt - (i + cycles_to_drop) * period, ff, label = 'cycle %d' % (i + 1))
        axs_act.plot(activation_t, activation_v, 'k--', label = 'final')
        axs_act.set_title('Activation function')
        axs_act.set_xlim([0, period])
        axs_act.set_xlabel('time [s]')
        axs_act.legend(loc = 'upper left', frameon = False, fontsize = 7, labelspacing = 0.2)

        axs_PV.set_xlim((0, 1.2 * np.max(data_loop[label_volume_PV_loops])))
        axs_PV.set_ylim((0, 1.2 * np.max(data_loop[label_pressure_PV_loops])))
        axs_PV.set_xlabel('volume [mL]')
        axs_PV.set_ylabel('pressure [mmHg]')
        axs_PV.set_title('pressure-volume relationships')
        axs_PV.legend(loc = 'upper left', frameon = False, fontsize = 7, labelspacing = 0.2)

        fig_out.tight_layout()

    ######## Output
    emulator = dict()

    emulator['period'] = period

    emulator['EDPV'] = dict()
    if file_EDPV is not None:
        emulator['EDPV']['klotz'] = False
        emulator['EDPV']['P'] = list(p_EDPV)
        emulator['EDPV']['V'] = list(V_EDPV)
    else:
        emulator['EDPV']['klotz'] = True
        emulator['EDPV']['An'] = klotz_An
        emulator['EDPV']['Bn'] = klotz_Bn
        emulator['EDPV']['V0'] = EDPV_V0
        emulator['EDPV']['V30'] = EDPV_V30

    emulator['ESPV'] = dict()
    emulator['ESPV']['V_0'] = ESPV_V0
    emulator['ESPV']['E'] = ESPV_E

    emulator['activation'] = dict()
    emulator['activation']['t'] = list(activation_t)
    emulator['activation']['v'] = list(activation_v)

    emulator['linear'] = dict()
    emulator['linear']['EA'] = ESPV_E
    emulator['linear']['EB'] = EDPV_E
    emulator['linear']['TC'] = duration_C / period
    emulator['linear']['TR'] = duration_R / period
    emulator['linear']['tC'] = time_C / period
    emulator['linear']['V0'] = ESPV_V0

    emulator_instance = Emulator(emulator)
    if output_file is not None:
        emulator_instance.save_file(output_file, verbose = verbose)

    if output_file_fig is not None:
        fig_out.savefig(output_file_fig)

    return emulator_instance

def convert_pressure(dataset, field_name, measure_unit):
    if measure_unit == 'Pa':
        dataset[field_name] *= 1.0 / 133.322
    elif measure_unit == 'kPa':
        dataset[field_name] *= 1e-3 / 133.322
    elif measure_unit == 'mmHg':
        pass
    else:
        raise Exception('Unknown pressure unit %s' % measure_unit)

def convert_volume(dataset, field_name, measure_unit):
    if measure_unit == 'm^3':
        dataset[field_name] *= 1e6
    elif measure_unit == 'L':
        dataset[field_name] *= 1e-3
    elif measure_unit == 'mL':
        pass
    else:
        raise Exception('Unknown volume unit %s' % measure_unit)