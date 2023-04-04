# import Dictionaries, List, and numpy
from typing import Dict, List
import numpy
import sys

# import files the Grid contains
from Bus import Bus
from Generator import Generator
from Transformer import Transformer
from TransmissionLine import TransmissionLine

# Create grid class to contain all buses, generators, transformers, and transmission lines
class Grid:

    # Pass the name of the Grid upon initialization
    def __init__(self, name: str):

        # Set name of Grid
        self.name: str = name

        # Establish bus orders and a dictionary for the buses
        self.buses_order: List[str] = list()
        self.buses: Dict[str, Bus] = dict()

        # Create a dictionary for generators, transformers, and transmissionlines, with the key being their name
        self.generators: Dict[str, Generator] = dict()
        self.transformers: Dict[str, Transformer] = dict()
        self.transmissionline: Dict[str, TransmissionLine] = dict()

        # Parameters for all Files
        self.Sbase = 100 # In units of MVA
        self.Vbase = 230 # In units of kV
        self.convergencevalue = 0.0001 #pu
        self.Q_k_limit = 175000000 #VA

    # Function to add a bus, by first making sure the bus does not already exist
    def __add_bus(self, bus):
        if bus not in self.buses.keys():
            self.buses[bus] = Bus(bus)
            self.buses_order.append(bus)

    # Function to add a transformer to the grid. It takes a name, initial bus, final bus, apparentpower (MVA),
    # low side voltage rating (kV), high side voltage rating (kV), impedance (pu), and xrratio
    def add_transformer(self, name: str, bus1: str, bus2: str, apparentpower: float,
                        v1rated: float, v2rated: float, impedance: float, xrratio: float):

        # Check for errors before adding transformer
        self.error_check_transformer(bus1, bus2, apparentpower, v1rated, v2rated, impedance, xrratio)

        # Add transformer to dictionary with all of its values
        self.transformers[name] = Transformer(name, bus1, bus2, apparentpower, v1rated, v2rated, impedance, xrratio, self.Sbase)

        # Add the buses it is connected to
        self.__add_bus(bus1)
        self.__add_bus(bus2)

    # Function to add a transmission line to the grid.
    # It takes parameters name, initial bus, final bus, length (miles), line coordinates (ft.(inches/12)),
    # Codeword of wire, number of bundles, and spacing between bundles (ft.(inches/12))
    # This function will throw an error if two phases are in the same exact location,
    # the number of bundles is not an integer 1-4, or if the codeword entered does not have data stored
    def add_transmissionline(self, name: str, bus1: str, bus2: str, lengthmi: float,
                             axaxis: float, ayaxis: float, bxaxis: float, byaxis: float, cxaxis: float, cyaxis: float,
                             codeword: str, numberofbundles: int, seperationdistance: float):

        # Check for errors before adding line
        self.error_check_transmission_line(bus1, bus2, lengthmi, axaxis, ayaxis, bxaxis, byaxis, cxaxis, cyaxis, codeword, numberofbundles, seperationdistance)

        # Add Vbase as the high voltage from the first transformer entered
        Vbase = self.transformers[list(self.transformers.keys())[0]].v2rated

        # Add transmission line to dictionary
        self.transmissionline[name] = TransmissionLine(name, bus1, bus2, lengthmi, axaxis, ayaxis, bxaxis, byaxis, cxaxis, cyaxis, codeword, numberofbundles, seperationdistance, Vbase, self.Sbase)

        # Add the buses it is connected to
        self.__add_bus(bus1)
        self.__add_bus(bus2)

    # Function to add a generator. It takes parameters name, initial bus, and nominal power (MVA).
    def add_generator(self, name, bus1, nominalpower):

        # Check for errors before adding generator
        self.error_check_generator(nominalpower)

        # Add generator to dictionary
        self.generators[name] = Generator(name, bus1, nominalpower)
        self.__add_bus(bus1)

    # Function to calculate the Ybus. There are no input parameters, and it outputs the Ybus matrix.
    # This function only works for a 7 bus network with the specified grid given in class.
    def calculate_Ybus(self):

        # Create a matrix with all zeros depending on how many buses are in the system of complex variables
        self.Ybus = numpy.zeros((len(self.buses_order), len(self.buses_order)), dtype=complex)

        # Assign Values to each component

        # Set Non-diagonals just using -1/Z
        self.Ybus[0][1] = -1 / (self.transformers["T1"].Rpu + 1j * self.transformers["T1"].Xpu) # T1
        self.Ybus[1][0] = self.Ybus[0][1] # T1
        self.Ybus[1][2] = -1 / (self.transmissionline["L2"].Rpu + 1j * self.transmissionline["L2"].Xpu) # L2
        self.Ybus[2][1] = self.Ybus[1][2] # L2
        self.Ybus[3][1] = -1 / (self.transmissionline["L1"].Rpu + 1j * self.transmissionline["L1"].Xpu) # L1
        self.Ybus[1][3] = self.Ybus[3][1] # L1
        self.Ybus[4][2] = -1 / (self.transmissionline["L3"].Rpu + 1j * self.transmissionline["L3"].Xpu) # L3
        self.Ybus[2][4] = self.Ybus[4][2] # L3
        self.Ybus[4][3] = -1 / (self.transmissionline["L6"].Rpu + 1j * self.transmissionline["L6"].Xpu) # L6
        self.Ybus[3][4] = self.Ybus[4][3] # L6
        self.Ybus[5][3] = -1 / (self.transmissionline["L4"].Rpu + 1j * self.transmissionline["L4"].Xpu) # L4
        self.Ybus[3][5] = self.Ybus[5][3] # L4
        self.Ybus[5][4] = -1 / (self.transmissionline["L5"].Rpu + 1j * self.transmissionline["L5"].Xpu) # L5
        self.Ybus[4][5] = self.Ybus[5][4] # L5
        self.Ybus[6][5] = -1 / (self.transformers["T2"].Rpu + 1j * self.transformers["T2"].Xpu) # T2
        self.Ybus[5][6] = self.Ybus[6][5] # T2

        # Set diagonals, do not include generators, include capacitances here
        # Use previous matrix values plus shunt charging values for necessary transmission lines
        self.Ybus[0][0] = -self.Ybus[0][1]  # G1, T1
        self.Ybus[1][1] = -self.Ybus[0][1] - self.Ybus[1][3] - self.Ybus[1][2] + (1j * self.transmissionline["L1"].Bpu / 2) + 1j * self.transmissionline["L2"].Bpu / 2# T1, L1, L2
        self.Ybus[2][2] = -self.Ybus[1][2] - self.Ybus[2][4] + (1j * self.transmissionline["L2"].Bpu / 2) + 1j * self.transmissionline["L3"].Bpu / 2           # L2, L3
        self.Ybus[3][3] = -self.Ybus[1][3] - self.Ybus[3][5] - self.Ybus[3][4] + (1j * self.transmissionline["L1"].Bpu / 2) + 1j * self.transmissionline["L4"].Bpu / 2 + 1j * self.transmissionline["L6"].Bpu / 2  # L1, L4, L6
        self.Ybus[4][4] = -self.Ybus[2][4] - self.Ybus[4][5] - self.Ybus[3][4] + (1j * self.transmissionline["L3"].Bpu / 2) + 1j * self.transmissionline["L5"].Bpu / 2 + 1j * self.transmissionline["L6"].Bpu / 2  # L3, L5, L6
        self.Ybus[5][5] = -self.Ybus[3][5] - self.Ybus[4][5] - self.Ybus[5][6] + (1j * self.transmissionline["L4"].Bpu / 2) + 1j * self.transmissionline["L5"].Bpu / 2  # L4, L5, T2
        self.Ybus[6][6] = -self.Ybus[5][6]  # G2, T2

        # Print the Y-bus matrix -> rechecked, correct
        #print("Y-bus matrix:")
        #i = 0
        #while i < 7:
        #    j = 0
        #    print("\nRow " + str(i + 1))
        #    while j < 7:
        #      print(self.Ybus[i][j])
        #      j += 1
        #    i = i + 1

    # Store bus data from main
    def setBusData(self, bus: str, bustype: str,  real_P: float, Q_or_V: float):
        self.buses[bus].setbusdata(bustype, real_P, Q_or_V)

    # Store power loss calculated from Power Flow class
    def store_power_loss(self, name: str, powerloss):
        self.transmissionline[name].store_power_loss(powerloss)

    def store_power_loss_transformer(self, name:str, powerloss):
        self.transformers[name].store_power_loss(powerloss)

    # Function to check errors in transmission line
    def error_check_transmission_line(self, bus1, bus2, lengthmi, axaxis, ayaxis, bxaxis, byaxis, cxaxis, cyaxis, codeword, numberofbundles, seperationdistance):

        # List of acceptable codewords
        codewordlist = ["Partridge"]

        # If two phases are in the same location, throw an error
        if (axaxis == bxaxis and ayaxis == byaxis) or (bxaxis == cxaxis and byaxis == cyaxis) or (
                axaxis == cxaxis and ayaxis == cyaxis):
            sys.exit("Error. Phases can not be in the same location. Enter phases at different coordinates.")

        # If number of bundles was less than 1 or more than 4, throw an error
        if numberofbundles < 1 or numberofbundles > 4:
            sys.exit("Error. Number of Bundles not accepted. Enter a value 1-4")

        # If codeword does not have data stored, throw an error
        if codeword not in codewordlist:
            sys.exit("Error: Codeword not accepted. Enter a different conductor type.")

        # If transmission line begins and ends at same bus, throw error
        if bus1 == bus2:
            sys.exit("Error: Cannot end a transmission line to the same bus it begins. Enter different buses.")

        # If transmission line length is negative or 0, throw an error
        if lengthmi <= 0:
            sys.exit("Error: Cannot have a negative or zero transmission line length. Enter a positive value.")

        # If separation distance is negative, throw an error
        if seperationdistance < 0:
            sys.exit("Error: Cannot have a negative bundle separation distance. Enter a positive value.")

    # Function to check errors in generator
    def error_check_generator(self, nominalpower):

        # If nominal power is negative, throw an error
        if nominalpower < 0:
            sys.exit("Error: Cannot have a negative generator power. Enter a positive value.")

    # Function to check errors in generator
    def error_check_transformer(self, bus1, bus2, apparentpower, v1rated, v2rated, impedance, xrratio):

        # If transmission line begins and ends at same bus, throw error
        if bus1 == bus2:
            sys.exit("Error: Cannot end a transformer at the same bus it begins. Enter different buses.")

        # If apparent power is less than zero, throw error
        if apparentpower < 0:
            sys.exit("Error: Cannot have a negative apparent power in the transformer. Enter a positive value.")

        # If v1rated is negative, throw an error
        if v1rated < 0:
            sys.exit("Error: Cannot have a negative v1 rated voltage. Enter a positive value.")

        # If v2rated is negative, throw an error
        if v2rated < 0:
            sys.exit("Error: Cannot have a negative v2 rated voltage. Enter a positive value.")

        # If impedance is negative, throw an error
        if impedance < 0:
            sys.exit("Error: Cannot have a negative transformer impedance. Enter a positive value.")

        # If xrratio is negative, throw an error
        if xrratio < 0:
            sys.exit("Error: Cannot have a negative transformer xrratio. Enter a positive value.")


