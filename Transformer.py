# Transformer class
import numpy


class Transformer:

    # Transformer has base parameters name, bus1, bus2, apparentpower, v1rated, v2rated, impedance, xrratio, Sbase, Zt_connection1, Zt_grounding1, Zt_value1, Zt_connection2, Zt_grounding2, Zt_value2
    def __init__(self, name, bus1, bus2, apparentpowerrating, v1rated, v2rated, impedance, xrratio, Sbase, Zt_connection1, Zt_grounding1, Zt_value1, Zt_connection2, Zt_grounding2, Zt_value2):

        # Set base values
        self.name = name
        self.bus1 = bus1
        self.bus2 = bus2
        self.apparentpowerrating = apparentpowerrating
        self.v1rated = v1rated
        self.v2rated = v2rated
        self.impedance = impedance
        self.xrratio = xrratio
        self.powerloss = 0
        self.Zt_connection1 = Zt_connection1
        self.Zt_grounding1 = Zt_grounding1
        self.Zt_value1 = Zt_value1
        self.Zt_connection2 = Zt_connection2
        self.Zt_grounding2 = Zt_grounding2
        self.Zt_value2 = Zt_value2

        # Establish Sbase and Vbase
        self.Sbase = Sbase  # MVA
        Vbase = v2rated  # kV

        # Calculate Z Real and Z imaginary for the transformer
        self.Rpu = impedance * (v2rated * v2rated / apparentpowerrating)/(Vbase * Vbase/Sbase) * numpy.cos(numpy.arctan(xrratio))
        self.Xpu = impedance * (v2rated * v2rated / apparentpowerrating)/(Vbase * Vbase/Sbase) * numpy.sin(numpy.arctan(xrratio))

    def store_power_loss(self, powerloss):
        self.powerloss = powerloss