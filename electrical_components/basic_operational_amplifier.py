from PySpice.Spice.Netlist import SubCircuitFactory
from PySpice.Unit import *


class BasicOperationalAmplifier(SubCircuitFactory):
    Name = 'BasicOperationalAmplifier'
    NODES = ('non_inverting_input', 'inverting_input', 'output')

    ##############################################

    def __init__(self):
        super().__init__()

        # Input impedance
        self.R('input', 'non_inverting_input', 'inverting_input', mega(10))

        # dc gain=100k and pole1=100hz
        # unity gain = dcgain x pole1 = 10MHZ
        self.VCVS('gain', 'non_inverting_input', 'inverting_input', 1, self.gnd, voltage_gain=kilo(100))
        self.R('P1', 1, 2, kilo(1))
        self.C('P1', 2, self.gnd, micro(1.5915))

        # Output buffer and resistance
        self.VCVS('buffer', 2, self.gnd, 3, self.gnd, 1)
        self.R('out', 3, 'output', 10)
