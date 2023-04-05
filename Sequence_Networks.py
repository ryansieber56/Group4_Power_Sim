# File to develop the positive, negative, and zero sequence bus impedance matrices
import numpy as np
class SequenceNet:

    # Power flow
    def __init__(self, Grid, Zg1_grounding: str, Zg1_value:float, Zg2_grounding: str, Zg2_value:float, Zt1_grounding1: str, Zt1_value1:float, Zt1_grounding2: str, Zt1_value2:float, Zt2_grounding1: str, Zt2_value1:float, Zt2_grounding2: str, Zt2_value2:float):
        # Per unit values, update later to make more generic
        self.x1generators_oldpu = 0.12
        self.x2generators_oldpu = 0.14
        self.x0generators_oldpu = 0.05

        # Establish Which Sequence Networks you want to create
        Generate_Zbus_0 = 1
        Generate_Zbus_1 = 1
        Generate_Zbus_2 = 1

        # Establish Base Values
        self.Sbase = Grid.Sbase * 1000000 # In VA
        self.Vbase_main = Grid.Vbase * 1000
        self.Vbase_bus1 = Grid.transformers[list(Grid.transformers.keys())[0]].v1rated * 1000
        self.Vbase_bus7 = Grid.transformers[list(Grid.transformers.keys())[1]].v1rated * 1000
        nominalpower1 = Grid.generators[list(Grid.generators.keys())[0]].nominalpower * 1000000
        nominalpower2 = Grid.generators[list(Grid.generators.keys())[1]].nominalpower * 1000000

        # Establish Other Parameters
        self.length = len(Grid.buses)
        self.Zbus0 = np.zeros(self.length, self.length, type=complex)
        self.Zbus1 = np.zeros(self.length, self.length, type=complex)
        self.Zbus2 = np.zeros(self.length, self.length, type=complex)
        self.Ybus0 = np.zeros(self.length, self.length, type=complex)
        self.Ybus1 = np.zeros(self.length, self.length, type=complex)
        self.Ybus2 = np.zeros(self.length, self.length, type=complex)

        # Transformer Grounding ->transformer1 through 1 Ohm on Y side, transformer2 Y ungrounded

        # Generate Zbus1 -> Positive Sequence
        if Generate_Zbus_1 == 1:

            # Switch Generators to System Per Unit Values instead of individual
            self.x1generators_newpu1 = self.x1generators_oldpu * Grid.Sbase / nominalpower1
            self.x1generators_newpu2 = self.x1generators_oldpu * Grid.Sbase / nominalpower2

            # Update Bus 1 with generator information

            # Update Bus 7 with generator information



        # Generate Zbus2 -> Negative Sequence
        if Generate_Zbus_2 == 1:
            # Switch Generators to System Per Unit Values instead of individual
            self.x2generators_newpu1 = self.x2generators_oldpu * Grid.Sbase / nominalpower1
            self.x2generators_newpu2 = self.x2generators_oldpu * Grid.Sbase / nominalpower2

            # Update Bus 1 with generator information

            # Update Bus 7 with generator information



        # Generate Zbus0 -> 0 Sequence
        if Generate_Zbus_0 == 1:

            # This is Generator Information
            if Zg1_grounding == "Solid ground":
                self.Zg1 = 0
            if Zg1_grounding == "Ungrounded":
                self.Zg1 = float('inf')
            if Zg1_grounding == "Resistor":
                self.Zg1 = Zg1_value

            if Zg2_grounding == "Solid ground":
                self.Zg2 = 0
            if Zg2_grounding == "Ungrounded":
                self.Zg2 = float('inf')
            if Zg2_grounding == "Resistor":
                self.Zg2 = Zg2_value

            # This is Transformer Information -> Depends on grounding, initial remains here for reference
            self.Ybus0[6][5] = -1 / (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu)  # T2
            self.Ybus0[5][6] = self.Ybus0[6][5]  # T2
            self.Ybus0[0][1] = -1 / (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu)  # T1
            self.Ybus0[1][0] = self.Ybus0[0][1]  # T1

            # Grounded Wye can go to Grounded Wye or Delta to not be 0
            if Zt1_grounding1 == "Grounded Wye":
                if Zt1_grounding2 == "Grounded Wye":
                    self.Zt11 = Zt1_value1 * 3 * (self.Vbase_bus1**2 / nominalpower1) * (Grid.Sbase / nominalpower2)
                    self.Zt12 = Zt1_value2 * 3 * (self.Vbase_main**2 / Grid.Sbase)
                    self.Ybus0[0][1] = -1 / (Grid.transformers["T1"].Rpu + self.Zt11 + 1j * Grid.transformers["T1"].Xpu)  # T1
                    self.Ybus0[1][0] = self.Ybus0[0][1]  # T1
                if Zt1_grounding2 == "Delta":
                    self.Zt11 = Zt1_value1 * 3 * (self.Vbase_bus1**2 / nominalpower1) * (Grid.Sbase / nominalpower2)
                    self.Ybus0[0][1] = -1 / (Grid.transformers["T1"].Rpu + self.Zt11 + 1j * Grid.transformers["T1"].Xpu)  # T1 Number
                    self.Ybus0[1][0] = 0  # T1
                if Zt1_grounding2 == "Ungrounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1

            # Delta to Grounded Wye Only Non-Zero Combination
            if Zt1_grounding1 == "Delta":
                if Zt1_grounding2 == "Grounded Wye":
                    self.Zt12 = Zt1_value2 * 3 * (self.Vbase_main**2 / Grid.Sbase)
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = -1 / (Grid.transformers["T1"].Rpu + self.Zt12 + 1j * Grid.transformers["T1"].Xpu)  # T1  # T1
                if Zt1_grounding2 == "Delta":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1
                if Zt1_grounding2 == "Ungrounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1

            # All Ungrounded Wye have 0
            if Zt1_grounding1 == "Ungrounded Wye":
                if Zt1_grounding2 == "Grounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1
                if Zt1_grounding2 == "Delta":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1
                if Zt1_grounding2 == "Ungrounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1

            # Zt2 Needs done
            if Zt2_grounding1 == "Grounded Wye to Grounded Wye":
                self.Ybus0[6][5] = -1 / (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu)  # T2
                self.Ybus0[5][6] = self.Ybus0[6][5]  # T2
            if Zt2_grounding1 == "Delta to Grounded Wye":
                self.Ybus0[6][5] = -1 / (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu)  # T2
                self.Ybus0[5][6] = self.Ybus0[6][5]  # T2
            if Zt2_grounding1 == "Ungrounded Wye to Grounded Wye":
                self.Ybus0[6][5] = -1 / (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu)  # T2
                self.Ybus0[5][6] = self.Ybus0[6][5]  # T2

            # Establish Generators, takes into account 3Zg as well -> Include in Ybus0 calculation, done, not sure about - signs
            self.totalx0generators_newpu1 = (self.x0generators_oldpu + 3 * self.Zg1 * (self.Vbase_bus1**2)/nominalpower1) * (Grid.Sbase / nominalpower1)
            self.totalx0generators_newpu2 = (self.x0generators_oldpu + 3 * self.Zg2 * (self.Vbase_bus7**2)/nominalpower2) * (Grid.Sbase / nominalpower2)

            # Establish Lines -> Z0 = 3Z1 -> Included in Ybus0 calculation, done
            # Set Non-diagonals just using -1/Z
            self.Ybus0[1][2] = -1 / (3 * (Grid.transmissionline["L2"].Rpu + 1j * Grid.transmissionline["L2"].Xpu))  # L2
            self.Ybus0[2][1] = self.Ybus0[1][2]  # L2
            self.Ybus0[3][1] = -1 / (3 * (Grid.transmissionline["L1"].Rpu + 1j * Grid.transmissionline["L1"].Xpu))  # L1
            self.Ybus0[1][3] = self.Ybus0[3][1]  # L1
            self.Ybus0[4][2] = -1 / (3 * (Grid.transmissionline["L3"].Rpu + 1j * Grid.transmissionline["L3"].Xpu))  # L3
            self.Ybus0[2][4] = self.Ybus0[4][2]  # L3
            self.Ybus0[4][3] = -1 / (3 * (Grid.transmissionline["L6"].Rpu + 1j * Grid.transmissionline["L6"].Xpu)) # L6
            self.Ybus0[3][4] = self.Ybus0[4][3]  # L6
            self.Ybus0[5][3] = -1 / (3 * (Grid.transmissionline["L4"].Rpu + 1j * Grid.transmissionline["L4"].Xpu))  # L4
            self.Ybus0[3][5] = self.Ybus0[5][3]  # L4
            self.Ybus0[5][4] = -1 / (3 * (Grid.transmissionline["L5"].Rpu + 1j * Grid.transmissionline["L5"].Xpu))  # L5
            self.Ybus0[4][5] = self.Ybus0[5][4]  # L5

            # Set diagonals, do not include generators, include capacitances here
            # Use previous matrix values plus shunt charging values for necessary transmission lines
            self.Ybus0[0][0] = -self.Ybus0[0][1] + self.totalx0generators_newpu1  # G1, T1
            self.Ybus0[1][1] = -self.Ybus0[0][1] - self.Ybus0[1][3] - self.Ybus0[1][2] + (1j * Grid.transmissionline["L1"].Bpu / 2) + 1j * Grid.transmissionline["L2"].Bpu / 2  # T1, L1, L2
            self.Ybus0[2][2] = -self.Ybus0[1][2] - self.Ybus0[2][4] + (1j * Grid.transmissionline["L2"].Bpu / 2) + 1j * Grid.transmissionline["L3"].Bpu / 2  # L2, L3
            self.Ybus0[3][3] = -self.Ybus0[1][3] - self.Ybus0[3][5] - self.Ybus0[3][4] + (1j * Grid.transmissionline["L1"].Bpu / 2) + 1j * Grid.transmissionline["L4"].Bpu / 2 + 1j * Grid.transmissionline["L6"].Bpu / 2  # L1, L4, L6
            self.Ybus0[4][4] = -self.Ybus0[2][4] - self.Ybus0[4][5] - self.Ybus0[3][4] + (1j * Grid.transmissionline["L3"].Bpu / 2) + 1j * Grid.transmissionline["L5"].Bpu / 2 + 1j * Grid.transmissionline["L6"].Bpu / 2  # L3, L5, L6
            self.Ybus0[5][5] = -self.Ybus0[3][5] - self.Ybus0[4][5] - self.Ybus0[5][6] + (1j * Grid.transmissionline["L4"].Bpu / 2) + 1j * Grid.transmissionline["L5"].Bpu / 2  # L4, L5, T2
            self.Ybus0[6][6] = -self.Ybus0[5][6] + self.totalx0generators_newpu2  # G2, T2

            # Update Bus 1 with generator information, depends on grounding
            # Update Bus 1 with transformer information, depends on grounding

            # Update Bus 7 with generator information, depends on grounding
            # Update Bus 7 with transformer information, depends on grounding
