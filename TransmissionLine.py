# Transmission Line class

import math
import numpy

# TransmissionLineBundles contains codeword, DSC, DSL, and resistance per feet information
from TransmissionLineBundles import TransmissionLineBundles


class TransmissionLine:

    # By default, transmission line has a name, bus to bus location, length  (miles), coordinates for each phase,
    # A codeword for the conductor type, number of bundles per phase, separation of bundles per phase, Vbase, and Sbase
    def __init__(self, name: str, bus1: str, bus2: str, lengthmi: float,
                 axaxis: float, ayaxis: float, bxaxis: float, byaxis: float, cxaxis: float, cyaxis: float,
                 codeword: str, numberofbundles: int, seperationdistance: float, Vbase: float, Sbase: float):

        # Set name, buses, and length of line
        self.name = name
        self.bus1 = bus1
        self.bus2 = bus2
        self.lengthmi = lengthmi
        self.numberofbundles = numberofbundles
        self.powerloss = None

        # Set S base, frequency, and calculate Z base
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
        Rpermi = Bundles.R

        # Calculate Capacitance values
        CFpermi = ((2 * numpy.pi * 8.8541878128 * (10 ** (-12))) / (math.log(dequivalent * 12 / DSC))) * 1609.344 # F/mi
        Ctotal = CFpermi * lengthmi # Farads

        # Calculate Inductance Values
        LFpermi = (2 * (10 ** (-7)) * math.log(dequivalent * 12 / DSL)) * 1609.344 # F/mi
        Ltotal = LFpermi * lengthmi # Farads

        # Calculate total resistance of line
        self.Rtotal = Rpermi * lengthmi

        # Use R to get R per unit
        self.Rpu = self.Rtotal / Zbase

        # Use L to get X per unit
        self.Xpu = Ltotal * 2 * numpy.pi * frequency / Zbase

        # Use C to get B per unit
        self.Bpu = Ctotal * 2 * numpy.pi * frequency * Zbase

    def store_power_loss(self, powerloss):
        self.powerloss = powerloss
