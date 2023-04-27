# Main File

from Grid import Grid
from Newton_Raphson_Power_Flow import NewtonRhapson
from DC_Power_Flow_Solver import DCPowerFlow
from Fast_Decoupled_Solver import FastDecoupled
from Sequence_Networks import SequenceNet
from Fault_Calculation import FaultCalculation

# Create Power Grid
MainGrid = Grid("MainGrid")

# Add Generators, Parameters: name, bus, and nominalpower (MVA), x1, x2, x0 per unit impedances, grounding type, grounding value
MainGrid.add_generator("G1", "Bus1", 100, 0.12, 0.14, 0.05, "Solid ground", 0)
MainGrid.add_generator("G2", "Bus7", 100, 0.12, 0.14, 0.05, "Resistor", 1)

# Add Transformers, Paramters: name, Starting Bus, Ending Bus, Apparentpower (MVAR), v1rated, v2rated, impedance,
# xrratio, Connection 1 type, Grounding type, Value of grounding, Connection 2 type, Grounding type, Value of grounding
MainGrid.add_transformer("T1", "Bus1", "Bus2", 125, 20, 230, 0.085, 10, "Delta", "N/A", 0, "Grounded Wye", "Resistor", 1)
MainGrid.add_transformer("T2", "Bus7", "Bus6", 200, 18, 230, 0.105, 12, "Delta", "N/A", 0, "Ungrounded Wye", "Ungrounded", 0)

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
# Capacitor bank parameter will only be used if VAR limit exceeded on non-slack generator
print("Which method of Power Flow would you like to solve? Type the number associated with your desired method.")
method = input("1) Newton-Rhapson\n2) Fast Decoupled\n3) DC Power Flow\n")
if method == "1":
    print("\nSelected Newton-Rhapson Power Flow")
elif method == "2":
    print("\nSelected Fast Decoupled Power Flow")
elif method == "3":
    print("\nSelected DC Power Flow")
else:
    print("\nInvalid Selection. Program exiting.")
    exit(-1)

faulttype = None
buslocation = None
faulting_impedance = None

# Solve Sequence Network
selection = input("Would you like a Fault Calculation? Type: YES or NO\n")
if selection == "YES":
    print("\nWhat type of Fault analysis would you like? Type the integer of the desired analysis.")
    faulttype = input("1)Symmetrical Fault\n2)Single Line to Ground Fault\n3)Line to Line Fault\n4)Double Line to Ground Fault\n")
    if faulttype == "1":
        faulttype = "Symmetrical Fault"
    elif faulttype == "2":
        faulttype = "Single Line to Ground Fault"
    elif faulttype == "3":
        faulttype = "Line to Line Fault"
    elif faulttype == "4":
        faulttype = "Double Line to Ground Fault"
    else:
        print("\nInvalid Selection. Program Ending.")
        exit(-1)

    buslocation = input("\nWhat bus do you want the fault to occur at? Type the integer of the desired bus between 1 and 7.\n")
    if(int(buslocation) < 1 or int(buslocation) > 7):
        print("Bus location out of allowed range. Default value of Bus 4 being used to continue.")
        buslocation = "4"
    faulting_impedance = input("\nIf there is a faulting impedance, enter the value of that impedance. If not enter 0.\n")

elif selection == "NO":
    print("\nNo Fault analysis will be conducted.")
else:
    print("\nInvalid selection. No Fault analysis will be conducted.")

print("\nPrinting Power Flow Information:")
if method == "1":
    NewtonRhapson(MainGrid, 0)
elif method == "2":
    FastDecoupled(MainGrid, 0)
elif method == "3":
    DCPowerFlow(MainGrid)
if selection == "YES":
    print("\nPrinting Relevant Fault Calculation Information:")
    SeqNet = SequenceNet(MainGrid)
    FaultCalculation(MainGrid, SeqNet, faulttype, int(buslocation), float(faulting_impedance))
