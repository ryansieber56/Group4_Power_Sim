# Main File

from Grid import Grid
from Newton_Raphson_Power_Flow import NewtonRhapson
from DC_Power_Flow_Solver import DCPowerFlow
from Fast_Decoupled_Solver import FastDecoupled
from Sequence_Networks import SequenceNet
from Fault_Calculation import FaultCalculation

# Create Power Grid
MainGrid = Grid("MainGrid")

# Add Generators, Parameters: name, bus, and nominalpower (MVA)
MainGrid.add_generator("G1", "Bus1", 100)
MainGrid.add_generator("G2", "Bus7", 100)

# Add Transformers, Paramters: name, Starting Bus, Ending Bus, Apparentpower (MVAR), v1rated, v2rated, impedance, and
# xrratio
MainGrid.add_transformer("T1", "Bus1", "Bus2", 125, 20, 230, 0.085, 10)
MainGrid.add_transformer("T2", "Bus7", "Bus6", 200, 18, 230, 0.105, 12)

# Add Transmission Line Connections
# Parameters: name, bus to bus location, length  (miles), coordinates for each phase, codeword for the conductor type,
# number of bundles per phase, and separation of bundles per phase
MainGrid.add_transmissionline("L1", "Bus2", "Bus4", 10, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L2", "Bus2", "Bus3", 25, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L3", "Bus3", "Bus5", 20, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L4", "Bus4", "Bus6", 20, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L5", "Bus5", "Bus6", 10, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L6", "Bus4", "Bus5", 35, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)

# Calculate Y_bus Matrix
MainGrid.calculate_Ybus()

# Set bus types, Parameters: name, Bus type, real power, Q or V depending on Bus type
MainGrid.setBusData("Bus1", "Slack Bus", 0, 0)
MainGrid.setBusData("Bus2", "Load Bus", 0, 0)
MainGrid.setBusData("Bus3", "Load Bus", 110, 50)
MainGrid.setBusData("Bus4", "Load Bus", 100, 70)
MainGrid.setBusData("Bus5", "Load Bus", 100, 65)
MainGrid.setBusData("Bus6", "Load Bus", 0, 0)
MainGrid.setBusData("Bus7", "Voltage Controlled Bus", 200, 1)

# Calculate Power Flow, Parameters: Grid, capacitor_bank
# Set integer to 0 if you want to solve the system by altering the non-slack generator out of PV mode
# Set integer to 1 if a capacitor bank should be added instead
#NewtonRhapson(MainGrid, 1)
FastDecoupled(MainGrid, 1)
#DCPowerFlow(MainGrid)

# Solve Sequence Network
#SeqNet = SequenceNet(MainGrid, 0.12, 0.14, 0.05, "Solid ground", 0, "Resistor", 1, "Delta", "N/A", 0, "Grounded Wye", "Resistor", 1, "Delta", "N/A", 0, "Grounded Wye", "Resistor", 1)
#FaultCalculation(MainGrid, SeqNet, "Symmetrical Fault", 7, 0)
