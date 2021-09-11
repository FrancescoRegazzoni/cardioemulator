import json
from scipy import interpolate
import copy
from ._emulator import Emulator

def transform_ES_elastance(emulator_data, factor):
    """
    Transform an emulator by applying a multiplicative factor to end systolic elastance.

    Parameters
    ----------
    emulator_data : str, dict or Emulator
        Emulator to be transformed.
    factor : float
        Multiplicative factor for end systolic elastance.

    Returns
    -------
    emulator : Emulator
        Transformed emulator.

    """

    if isinstance(emulator_data, str):
        with open(emulator_data) as json_file:
            emulator_data = json.load(json_file)
    if isinstance(emulator_data, Emulator):
        emulator_data = emulator_data.data

    emulator_data_new = copy.deepcopy(emulator_data)
    emulator_data_new['ESPV']['E'] *= factor

    return Emulator(emulator_data_new)

def transform_time_shift(emulator_data, time_lag):
    """
    Transform an emulator by shifting in time the activation pattern.

    Parameters
    ----------
    emulator_data : str, dict or Emulator
        Emulator to be transformed.
    time_lag : float [s]
        Time lag.

    Returns
    -------
    emulator : Emulator
        Transformed emulator.

    """

    if isinstance(emulator_data, str):
        with open(emulator_data) as json_file:
            emulator_data = json.load(json_file)
    if isinstance(emulator_data, Emulator):
        emulator_data = emulator_data.data

    activation_base = interpolate.interp1d(emulator_data['activation']['t'], emulator_data['activation']['v'])
    activation = lambda t: activation_base(t % emulator_data['period'])

    emulator_data_new = copy.deepcopy(emulator_data)
    emulator_data_new['activation']['v'] = activation(emulator_data['activation']['t'] - time_lag)

    return Emulator(emulator_data_new)