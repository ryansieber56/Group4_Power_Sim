# Main File
# Need help with:
#   DC Power Flow Solver - Can't get DC Power Flow Correct  -> Ybus exact same, Bbus exact same for what it should be
# To Do:


from Grid import Grid
from Newton_Raphson_Power_Flow import NewtonRhapson
from DC_Power_Flow_Solver import DCPowerFlow
from Fast_Decoupled_Solver import FastDecoupled
from Sequence_Networks import SequenceNet
from Fault_Calculation import FaultCalculation
# Create Power Grid
MainGrid = Grid("MainGrid")

# Add Generators
MainGrid.add_generator("G1", "Bus1", 100)
MainGrid.add_generator("G2", "Bus7", 100)

# Add Transformers
MainGrid.add_transformer("T1", "Bus1", "Bus2", 125, 20, 230, 0.085, 10)
MainGrid.add_transformer("T2", "Bus7", "Bus6", 200, 18, 230, 0.105, 12)

# Add Transmission Line Connections
MainGrid.add_transmissionline("L1", "Bus2", "Bus4", 10, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L2", "Bus2", "Bus3", 25, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L3", "Bus3", "Bus5", 20, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L4", "Bus4", "Bus6", 20, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L5", "Bus5", "Bus6", 10, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L6", "Bus4", "Bus5", 35, 0, 0, 19.5, 0, 39, 0, "Partridge", 2, 1.5)

# Calculate Y_bus Matrix
MainGrid.calculate_Ybus()

# Set bus types
MainGrid.setBusData("Bus1", "Slack Bus", 0, 0)
#MainGrid.setBusData("Bus7", "Slack Bus", 0, 0)

MainGrid.setBusData("Bus2", "Load Bus", 0, 0)
MainGrid.setBusData("Bus3", "Load Bus", 110, 50)
#MainGrid.setBusData("Bus4", "Load Bus", 100, 70)
MainGrid.setBusData("Bus4", "Load Bus", 100, 70)

MainGrid.setBusData("Bus5", "Load Bus", 100, 65)
MainGrid.setBusData("Bus6", "Load Bus", 0, 0)
MainGrid.setBusData("Bus7", "Voltage Controlled Bus", 200, 1)
#MainGrid.setBusData("Bus1", "Voltage Controlled Bus", 200, 1)

# Calculate Power Flow
NewtonRhapson(MainGrid)
#DCPowerFlow(MainGrid)  # should have no losses
#FastDecoupled(MainGrid)

# Solve Sequence Network
SeqNet = SequenceNet(MainGrid, 0.12, 0.14, 0.05, "Solid ground", 0, "Resistor", 1, "Delta", "N/A", 0, "Grounded Wye", "Resistor", 1, "Delta", "N/A", 0, "Grounded Wye", "Resistor", 1)
FaultCalculation(MainGrid, SeqNet, "Symmetrical Fault", 7, 0)