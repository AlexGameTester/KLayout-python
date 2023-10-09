from collections import OrderedDict
from itertools import product

import numpy as np
from datetime import datetime

from classLib.chipDesign import ChipDesign
from .sonnetLab import SonnetLab, SimulationBox

class SimulatedDesign(ChipDesign):
    def __init__(self, cell_name):
        super().__init__(cell_name)

        self.SL = None # matlab interface for simulating here

        # structure is {"sweep_par_name":sweep_par_values_list}
        # simulation is intended to happen across tensor product
        # of all parameters
        self._swept_pars = OrderedDict()

        # data storage references definition
        self.sMatrices = None  # np.array

        # fixed parameters definition
        self.freqs = None  # np.linspace(freq_start, freq_end, freqs_N)
        self.simulated_layer = self.layer_ph
        self.simulation_type = None

        # fixed or variable parameters definition
        # can be provided in set_fixed_params as well as
        # key: values_array pair in set_swept_params
        self.simBox = SimulationBox()  # box size and grid resolution

        # variable that depend on swept_pars
        self.ports = []  # list of SonnetPort class instances

        # additional control and visualizing variables
        self._start_time = None

        # provided by programmer (OPTIONAL)
        self._name = "default"  # here is optional measurement is stored

    def __reopen_socket(self):
        del self.SL
        self.SL = SonnetLab()  # new socket for every iteration possible memory leakage

    def calculate_ports(self, design_params):
        """
        virtual method
        """
        raise NotImplementedError

    def draw_simulation(selfs, iter_params_dict):
        raise NotImplementedError

    def set_fixed_parameters(self, freqs, simBox=None,  simulated_layer=None,
                             simulation_type="ABS"):
        self.simulation_type = simulation_type
        self.simBox = simBox
        self.freqs = freqs
        if simulated_layer is not None:
            self.simulated_layer = simulated_layer  # by default it is self.layer_ph

    def set_swept_parameters(self, sweep_parameters):
        """
        e.g. "simBox": list of SimulationBox() instances \n
        with, maybe, different amount of cells
        """
        self._swept_pars = sweep_parameters

    def allocate_sMatrices(self, freqs_n):
        self.post_freqs = np.zeros(tuple(
            len(swept_par_list) for swept_par_list in self._swept_pars.values()
        )+(freqs_n, ), dtype=np.complex128)
        self.sMatrices = np.zeros(tuple(
            len(swept_par_list) for swept_par_list in self._swept_pars.values()
        )+(freqs_n, len(self.ports), len(self.ports)),
                                  dtype=np.complex128)

    def get_Sij(self, i, j):
        return self.sMatrices[..., i+1, j+1]

    def set_measurement_name(self, name):
        self._name = name

    def simulate_sweep(self):
        self._start_time = datetime.now()

        vals_prod = product(*self._swept_pars.values())
        vals_length_list = list(map(lambda x: len(x), list(self._swept_pars.values())))
        _idxs_iterables = [range(vals_length_list[i]) for i in range(len(self._swept_pars))]
        idxs_prod = product(*_idxs_iterables)

        iter_i = 0
        for idxs, values in zip(idxs_prod, vals_prod):
            iter_params_dict = OrderedDict([(key, val) for key, val in zip(self._swept_pars.keys(), values)])
            self.draw_simulation(iter_params_dict)
            if( iter_i == 0 ):
                freqs, sMatrices = self.simulate_design(iter_params_dict)
                self.allocate_sMatrices(len(freqs))
            else:
                freqs, sMatrices = self.simulate_design(iter_params_dict)

            self.post_freqs[idxs] = freqs
            self.sMatrices[idxs] = sMatrices

    def simulate_design(self, iter_params_dict):
        self.__reopen_socket()

        ### parameters that can be both fixed or swept START ###
        if "simBox" in iter_params_dict:
            self.simBox = iter_params_dict["simBox"]
        elif self.simBox is None:
            print("simulate_design has no boxProps property")
        ### parameters that can be both fixed or swept END ###

        self.SL.clear()
        self.SL.set_boxProps(self.simBox)
        if self.simulation_type == "LINEAR":
            self.SL.set_linspace_sweep(self.freqs[0]/1e9, self.freqs[-1]/1e9, len(self.freqs))
        elif self.simulation_type == "ABS":
            self.SL.set_ABS_sweep(self.freqs[0]/1e9, self.freqs[-1]/1e9)
        else:
            self.SL.set_ABS_sweep(self.freqs[0]/1e9, self.freqs[-1]/1e9)

        self.calculate_ports(self.design_pars)
        self.SL.set_ports(self.ports)
        reg2sim = self._reg_from_layer(self.simulated_layer)
        self.SL.send_polygons(reg2sim)  # only 1 cell is supported
        # print("starting simulation")
        self.SL.start_simulation(wait=True)
        self.SL.release()
        return self.SL.get_s_params()

    def get_save_path(self, path=None):
        import os

        if not os.path.exists("data"):
            os.makedirs("data")

        sample_directory = os.path.join('data', self.__class__.__name__)  # both for python 2.x and 3.x
        if not os.path.exists(sample_directory):
            os.makedirs(sample_directory)

        date_directory = os.path.join(sample_directory,
                                      self._start_time.strftime("%b %d %Y"))
        if not os.path.exists(date_directory):
            os.makedirs(date_directory)

        time_directory = os.path.join(date_directory,
                                      self._start_time.strftime("%H-%M-%S")
                                      + " - " + self._name)

        if not os.path.exists(time_directory):
            os.makedirs(time_directory)

        return time_directory

    def save(self):
        import os
        import pickle as pkl
        with open(os.path.join(self.get_save_path(), self._name + '.pkl'), 'w+b') as f:
            pkl.dump(self, f)

