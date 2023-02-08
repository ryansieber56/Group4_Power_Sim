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

        # Establish Sbase
        Sbase = 100  # MVA

        # Calculate Z Real and Z imaginary for the transformer
        self.Zpur = impedance * (v1rated * v1rated / apparentpowerrating)/(v1rated * v1rated/Sbase) * numpy.cos(numpy.arctan(xrratio))
        self.Zpui = impedance * (v1rated * v1rated / apparentpowerrating)/(v1rated * v1rated/Sbase) * numpy.sin(numpy.arctan(xrratio))