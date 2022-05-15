from PySpice.Spice.Netlist import SubCircuitFactory
from PySpice.Unit import *


class BasicComparator(SubCircuitFactory):
    __name__ = 'BasicComparator'
    __nodes__ = ('non_inverting_input', 'inverting_input',
                 'voltage_plus', 'voltage_minus',
                 'output')

    ##############################################

    def __init__(self, ):
        super().__init__()

        # Fixme: ngspice is buggy with such subcircuit

        # Fixme: how to pass voltage_plus, voltage_minus ?
        # output_voltage_minus, output_voltage_plus = 0, 15

        # to plug the voltage source
        self.R(1, 'voltage_plus', 'voltage_minus', 1 @ u_MΩ)
        self.NonLinearVoltageSource(1, 'output', 'voltage_minus',
                                    expression='V(non_inverting_input, inverting_input)',
                                    # table=((-micro(1), output_voltage_minus),
                                    #       (micro(1), output_voltage_plus))
                                    table=(('-1uV', '0V'), ('1uV', '15V'))
                                    )
