from scipy import interpolate
import json

class Emulator:
    """
    0D emulator of a cardiac chamber.
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

        if self.data['EDPV']['klotz']:
            self.EDPV = lambda V: self.data['EDPV']['An'] * ((V - self.data['EDPV']['V0']) / (self.data['EDPV']['V30'] - self.data['EDPV']['V0']))**self.data['EDPV']['Bn']
        else:
            self.EDPV = interpolate.interp1d(self.data['EDPV']['V'], self.data['EDPV']['P'], fill_value = 'extrapolate')

        self.activation_base = interpolate.interp1d(self.data['activation']['t'], self.data['activation']['v'])

    def ESPV(self, V):
        """
        End-systolic pressure-volume relationship.

        Parameters
        ----------
        V : float [mL]
            Chamber volume.

        Returns
        -------
        p : float [mmHg]
            Chamber pressure.
        """
        return self.data['ESPV']['E'] * ( V - self.data['ESPV']['V_0'] )

    def EDPV(self, V):
        """
        End-distolic pressure-volume relationship.

        Parameters
        ----------
        V : float [mL]
            Chamber volume.

        Returns
        -------
        p : float [mmHg]
            Chamber pressure.
        """
        raise Exception('EDPV function should be overwritten in the class constructor!')

    def activation(self, t):
        """
        Time-dependent activation function.

        Parameters
        ----------
        t : float [s]
            time.

        Returns
        -------
        act : float [-]
            Activation (0 = fully relaxed, 1 = fully activated).
        """
        return self.activation_base(t % self.data['period'])

    def PV(self, V, t):
        """
        Time-dependent pressure-volume relationship.

        Parameters
        ----------
        V : float [mL]
            Chamber volume.
        t : float [s]
            time.

        Returns
        -------
        p : float [mmHg]
            Chamber pressure.
        """
        return (1 - self.activation(t)) * self.EDPV(V) + self.activation(t) * self.ESPV(V)

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

class Emulator_parametric:
    """
    Parametric 0D emulator of a cardiac chamber.
    """
    def __init__(self, emulator_data):
        """
        Parameters
        ----------
        emulator_data : dict or str
            Dictionary encoding the parametric emulator.
            This input can be provided as a dictionary or as the path to json
            file containing the dictionary.
        """

        if isinstance(emulator_data, str):
            with open(emulator_data) as json_file:
                emulator_data = json.load(json_file)
        self.data = emulator_data

        self.emulatorA = Emulator(emulator_data['emulatorA'])
        self.emulatorB = Emulator(emulator_data['emulatorB'])
        self.paramA = emulator_data['paramA']
        self.paramB = emulator_data['paramB']

    def PV(self, V, t, param):
        """
        Time-dependent pressure-volume relationship.

        Parameters
        ----------
        V : float [mL]
            Chamber volume.
        t : float [s]
            time.
        param : float
            Emulator parameter.

        Returns
        -------
        p : float [mmHg]
            Chamber pressure.
        """
        return (param - self.paramA) / (self.paramB - self.paramA) * self.emulatorB.PV(V,t) \
             + (param - self.paramB) / (self.paramA - self.paramB) * self.emulatorA.PV(V,t)
             
    def save_file(self, output_file, verbose = False):
        """
        Save the parametric emulator in json format.

        Parameters
        ----------
        output_file : str
            File path.
        """
        with open(output_file, 'w') as f:
            json.dump(self.data, f, indent = 4)
        if verbose:
            print('Saved file %s' % output_file)

def get_PV_relationship(emulator_data, emulator_param = None):
    """
    Get the time-dependent pressure-volume relationship associated
    with an emulator.

    Parameters
    ----------
    emulator_data : dict or str
        Dictionary encoding the emulator or parametric emulator.
        This input can be provided as a dictionary or as the path to json
        file containing the dictionary.
    emulator_param : float, optional
        In case the ``emulator_data`` encodes a parametric emulator, this input
        represents the parameter associated with the concrete emulator.

    Returns
    -------
    fun : callable
        Time-dependent pressure-volume relationship

            ``fun(volume, time) -> pressure``

        with input volume [mL], time [s] and output pressure [mmHg].

    """

    if isinstance(emulator_data, str):
        with open(emulator_data) as json_file:
            emulator_data = json.load(json_file)

    if emulator_data.get('parametric', False):
        return lambda V, t: Emulator_parametric(emulator_data).PV(V, t, emulator_param)
    else:
        return Emulator(emulator_data).PV
