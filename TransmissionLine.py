# Transmission Line class

import math
import numpy

# TransmissionLineBundles contains codeword, DSC, DSL, and resistance per feet information
from TransmissionLineBundles import TransmissionLineBundles


class TransmissionLine:

    # By default, transmission line has a name, bus to bus location, length  (miles), coordinates for each phase,
    # A codeword for the conductor type, number of bundles per phase, seperation of bundles per phase
    def __init__(self, name: str, bus1: str, bus2: str, lengthmi: float,
                 axaxis: float, ayaxis: float, bxaxis: float, byaxis: float, cxaxis: float, cyaxis: float,
                 codeword: str, numberofbundles: int, seperationdistance: float):

        # Set name, buses, and length of line
        self.name = name
        self.bus1 = bus1
        self.bus2 = bus2
        self.lengthmi = lengthmi

        # Set S base, Vbase, frequency, and calculate Z base
        Sbase = 100
        Vbase = 230
        Zbase = Vbase**2/Sbase
        frequency = 60

        # Calculate Equivalent Distance between lines (in feet)
        Dab = ((axaxis - bxaxis) ** 2 - (ayaxis - byaxis) ** 2) ** (1/2)
        Dbc = ((bxaxis - cxaxis) ** 2 - (byaxis - cyaxis) ** 2) ** (1/2)
        Dca = ((cxaxis - axaxis) ** 2 - (cyaxis - ayaxis) ** 2) ** (1/2)
        dequivalent = (Dab * Dbc * Dca) ** (1/3)

        # Create a Bundles object that contains the desired values from the TranmissionLineBundles class
        Bundles = TransmissionLineBundles(numberofbundles, seperationdistance, codeword)

        # Store the variables from that class
        DSL = Bundles.DSL
        DSC = Bundles.DSC
        R = Bundles.resistanceperft * self.lengthmi * 5280 / numberofbundles # 5280 converts miles to ft

        # Calculate Capacitance values
        CFpermi = ((2 * numpy.pi * 8.8541878128 * (10 ** (-12))) / (math.log(dequivalent * 12 / DSC))) * 1609.344 # F/mi
        Ctotal = CFpermi * lengthmi # Farads

        # Calculate Inductance Values
        LFpermi = (2 * (10 ** (-7)) * math.log(dequivalent * 12 / DSL)) * 1609.344 # F/mi
        Ltotal = LFpermi * lengthmi # Farads

        # Use R to get R per unit
        self.Rpu = R / Zbase

        # Use L to get X per unit
        self.Xpu = Ltotal * 2 * numpy.pi * frequency / Zbase

        # Use C to get B per unit
        self.Bpu = Ctotal * 2 * numpy.pi * frequency * Zbase
