import json
from ._emulator import Emulator_parametric

def build_emulator_parametric(
        emulator_A,
        emulator_B,
        param_A,
        param_B,
        output_file = None,
        verbose = True):
    """
    Build a parametric emulator, by interpolation of two emulators.

    Parameters
    ----------
    emulator_A : dict or str
        Dictionary encoding emulator A.
        This input can be provided as a dictionary or as the path to json
        file containing the dictionary.
    emulator_B : dict or str
        Dictionary encoding emulator B.
        This input can be provided as a dictionary or as the path to json
        file containing the dictionary.
    param_A : float
        Parameter associated with emulator A.
    param_B : float
        Parameter associated with emulator B.
    output_file : str, optional
        Path to a json file to store the dictionary encoding the parametric emulator.
    verbose : bool, optional
        Flag to toggle verbosity.

    Returns
    -------
    emulator : Emulator_parametric
        Parametric emulator.

    """

    if isinstance(emulator_A, str):
        with open(emulator_A) as json_file:
            emulator_A = json.load(json_file)
    if isinstance(emulator_B, str):
        with open(emulator_B) as json_file:
            emulator_B = json.load(json_file)

    if emulator_A['period'] != emulator_B['period']:
        raise Exception('The two emulators should have the same period')

    emulator = dict()
    emulator['period'] = emulator_A['period']
    emulator['parametric'] = True
    emulator['paramA'] = param_A
    emulator['paramB'] = param_B
    emulator['emulatorA'] = emulator_A
    emulator['emulatorB'] = emulator_B

    if output_file is not None:
        with open(output_file, 'w') as f:
            json.dump(emulator, f, indent=4)
        if verbose:
            print('saved file %s' % output_file)

    return Emulator_parametric(emulator)