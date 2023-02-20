# Main File
from Grid import Grid
from Power_Flow import PowerFlow

# Create Power Grid
MainGrid = Grid("MainGrid")

# Add Generators
MainGrid.add_generator("G1", "Bus1", 100)
MainGrid.add_generator("G2", "Bus7", 100)

# Add Transformers
MainGrid.add_transformer("T1", "Bus1", "Bus2", 125, 20, 230, 0.085, 10)
MainGrid.add_transformer("T2", "Bus7", "Bus6", 200, 18, 230, 0.105, 12)

# Add Transmission Line Connections
MainGrid.add_transmissionline("L1", "Bus2", "Bus4", 10, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L2", "Bus2", "Bus3", 25, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L3", "Bus3", "Bus5", 20, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L4", "Bus4", "Bus6", 20, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L5", "Bus5", "Bus6", 10, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)
MainGrid.add_transmissionline("L6", "Bus4", "Bus5", 35, 0, 0, 19.5, 0, 40, 0, "Partridge", 2, 1.5)

# Calculate and Display Y_bus Matrix
MainGrid.calculate_Ybus()

# Calculate Power Flow
PowerFlow(MainGrid)

