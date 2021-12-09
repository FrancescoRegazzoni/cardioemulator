import json
import numpy as np

class EmulatorValve:
    """
    0D emulator of a valve.

    Represents the valve as a function Q = delta_p / R(delta_p; R_min, R_max),
    finding the values of Rmin and Rmax that fit given data.
    """

    def __init__(self, emulator_data):
        """
        Parameters
        ----------
        emulator_data : dict or str
            Dictionary encoding the emulator.
            This input can be provided as a dictionary or as the path to json
            file containing the dictionary.
        """

        if isinstance(emulator_data, str):
            with open(emulator_data) as json_file:
                emulator_data = json.load(json_file)
        self.data = emulator_data

        self.R_min = self.data["R_min"]
        self.R_max = self.data["R_max"]

    def R(self, p_up, p_down):
        """
        Valve resistance as a function of the pressure jump.
        """

        heavisideMY = lambda x: np.arctan( np.pi / 2 * x * 200 ) * 1 / np.pi + 0.5
        return 10.**( np.log10( self.R_min ) + ( np.log10( self.R_max ) - np.log10( self.R_min ) ) * heavisideMY( p_down - p_up ) )

    def save_file(self, output_file, verbose = False):
        """
        Save the emulator in json format.

        Parameters
        ----------
        output_file : str
            File path.
        """
        with open(output_file, 'w') as f:
            json.dump(self.data, f, indent = 4)
        if verbose:
            print('Saved file %s' % output_file)

def get_valve_R_function(emulator_data):
    if isinstance(emulator_data, str):
        with open(emulator_data) as json_file:
            emulator_data = json.load(json_file)

    return EmulatorValve(emulator_data).R
