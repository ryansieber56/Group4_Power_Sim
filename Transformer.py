# Transformer class
import numpy


class Transformer:

    # Transformer has base parameters name, bus1, bus2, apparentpower, v1rated, v2rated, impedance, and xrratio
    def __init__(self, name, bus1, bus2, apparentpowerrating, v1rated, v2rated, impedance, xrratio):

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
        # Establish Sbase and Vbase
        Sbase = 100  # MVA
        Vbase = v2rated  # kV

        # Calculate Z Real and Z imaginary for the transformer
        self.Rpu = impedance * (v2rated * v2rated / apparentpowerrating)/(Vbase * Vbase/Sbase) * numpy.cos(numpy.arctan(xrratio))
        self.Xpu = impedance * (v2rated * v2rated / apparentpowerrating)/(Vbase * Vbase/Sbase) * numpy.sin(numpy.arctan(xrratio))
        #print(self.name, " Rpu ", self.Rpu, " XPU ", self.Xpu)

    def store_power_loss(self, powerloss):
        self.powerloss = powerloss