from PySpice.Spice.Netlist import SubCircuitFactory
from PySpice.Unit import *
from basic_operational_amplifier import BasicOperationalAmplifier


class InvertingOpAmp(SubCircuitFactory):
    """
    Inverting OpAmp SubCir
    Termanals:
        Vin
        Vout
    Parms:
        R1[Ohms]
        R2[Ohms]
    """

    Name = 'InvertingOpAmp'
    Nodes = ('Vin', 'Vout')

    def __init__(self, R1=1, R2=1):
        super().__init__()
        self.R1 = R1;
        self.R2 = R2

        self.R('1', 'Vin', '2', R1 @ u_Ω)
        self.R('2', '2', 'Vout', R2 @ u_Ω)

        self.subcircuit(BasicOperationalAmplifier())
        self.X('op', 'BasicOperationalAmplifier', self.gnd, '2', 'Vout')
        self.Theory()
