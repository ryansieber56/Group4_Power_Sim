# File to develop the positive, negative, and zero sequence bus impedance matrices
import numpy as np
class SequenceNet:

    # Power flow
    def __init__(self, Grid):
        # Values for generator
        self.x1generators_oldpu1 = Grid.generators["G1"].x1gen
        self.x2generators_oldpu1 = Grid.generators["G1"].x2gen
        self.x0generators_oldpu1 = Grid.generators["G1"].x0gen
        self.x1generators_oldpu2 = Grid.generators["G2"].x1gen
        self.x2generators_oldpu2 = Grid.generators["G2"].x2gen
        self.x0generators_oldpu2 = Grid.generators["G2"].x0gen
        Zg1_grounding = Grid.generators["G1"].grounding_type
        Zg1_value = Grid.generators["G1"].grounding_value
        Zg2_grounding = Grid.generators["G2"].grounding_type
        Zg2_value = Grid.generators["G2"].grounding_value

        # Transformer values
        Zt1_connection1 = Grid.transformers["T1"].Zt_connection1
        Zt1_grounding1 = Grid.transformers["T1"].Zt_grounding1
        Zt1_value1 = Grid.transformers["T1"].Zt_value1
        Zt1_connection2 = Grid.transformers["T1"].Zt_connection2
        Zt1_grounding2 = Grid.transformers["T1"].Zt_grounding2
        Zt1_value2 = Grid.transformers["T1"].Zt_value2
        Zt2_connection1 = Grid.transformers["T2"].Zt_connection1
        Zt2_grounding1 = Grid.transformers["T2"].Zt_grounding1
        Zt2_value1 = Grid.transformers["T2"].Zt_value1
        Zt2_connection2 = Grid.transformers["T2"].Zt_connection2
        Zt2_grounding2 = Grid.transformers["T2"].Zt_grounding2
        Zt2_value2 = Grid.transformers["T2"].Zt_value2

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
        self.Zbasemain = self.Vbase_main ** 2 / self.Sbase
        self.Zbase1 = self.Vbase_bus1 ** 2 / nominalpower1
        self.Zbase7 = self.Vbase_bus7 ** 2 / nominalpower2

        # Establish Other Parameters -> [Row][Column]
        self.length = len(Grid.buses)
        self.Zbus0 = np.zeros((self.length, self.length), dtype=complex)
        self.Zbus1 = np.zeros((self.length, self.length), dtype=complex)
        self.Zbus2 = np.zeros((self.length, self.length), dtype=complex)
        self.Ybus0 = np.zeros((self.length, self.length), dtype=complex)
        self.Ybus1 = np.zeros((self.length, self.length), dtype=complex)
        self.Ybus2 = np.zeros((self.length, self.length), dtype=complex)

        # Transformer Grounding ->transformer1 through 1 Ohm on Y side, transformer2 Y ungrounded

        # Generate Zbus1 -> Positive Sequence
        if Generate_Zbus_1 == 1:

            # Switch Generators to System Per Unit Values instead of individual
            self.x1generators_newpu1 = self.x1generators_oldpu1 * self.Sbase / nominalpower1
            self.x1generators_newpu2 = self.x1generators_oldpu2 * self.Sbase / nominalpower2

            # Ybu1 is slightly modified Ybus from before
            self.Ybus1 = Grid.Ybus

            # Update Bus 1 with generator information
            self.Ybus1[0][0] = 1 / (Grid.transformers["T1"].Rpu + 1j * self.x1generators_newpu1 + 1j * Grid.transformers["T1"].Xpu)

            # Update Bus 7 with generator information
            self.Ybus1[6][6] = 1 / (Grid.transformers["T2"].Rpu + 1j * self.x1generators_newpu2 + 1j * Grid.transformers["T2"].Xpu)

            # Zbus1 is inverse of Ybus1
            self.Zbus1 = np.linalg.inv(self.Ybus1)
            #print("Ybus1")
            #self.printmatrix(self.Ybus1)
            #print("Zbus1")
            #self.printmatrix(self.Zbus1)

        # Generate Zbus2 -> Negative Sequence
        if Generate_Zbus_2 == 1:
            # Switch Generators to System Per Unit Values instead of individual
            self.x2generators_newpu1 = self.x2generators_oldpu1 * self.Sbase / nominalpower1
            self.x2generators_newpu2 = self.x2generators_oldpu2 * self.Sbase / nominalpower2

            # Ybu2 is slightly modified Ybus from before
            self.Ybus2 = Grid.Ybus

            # Update Bus 1 with generator information
            self.Ybus2[0][0] = 1 / (Grid.transformers["T1"].Rpu + 1j * self.x2generators_newpu1 + 1j * Grid.transformers["T1"].Xpu)

            # Update Bus 7 with generator information
            self.Ybus2[6][6] = 1 / (Grid.transformers["T2"].Rpu + 1j * self.x2generators_newpu2 + 1j * Grid.transformers["T2"].Xpu)

            # Zbus2 is inverse of Ybus2
            self.Zbus2 = np.linalg.inv(self.Ybus2)
            #print("Ybus2")
            #self.printmatrix(self.Ybus2)
            #print("Zbus2")
            #self.printmatrix(self.Zbus2)

        # Generate Zbus0 -> 0 Sequence, Need Help
        if Generate_Zbus_0 == 1:

            # Generator 1 Grounding Information
            if Zg1_grounding == "Solid ground":
                self.Zg1 = 0
            elif Zg1_grounding == "Ungrounded":
                self.Zg1 = float('inf')
            elif Zg1_grounding == "Resistor":
                self.Zg1 = Zg1_value
            else:
                print("Invalid Grounding Information For Generator 1")
                exit(-1)

            # Generator 2 Grounding Information
            if Zg2_grounding == "Solid ground":
                self.Zg2 = 0
            elif Zg2_grounding == "Ungrounded":
                self.Zg2 = float('inf')
            elif Zg2_grounding == "Resistor":
                self.Zg2 = Zg2_value
            else:
                print("Invalid Grounding Information For Generator 2")
                exit(-1)

            # Establish Generators, takes into account 3Zg as well -> Include in Ybus0 calculation
            self.totalx0generators_newpu1 = (1j * self.x0generators_oldpu1 + 3 * self.Zg1 / self.Zbase1) * (self.Sbase / nominalpower1)
            self.totalx0generators_newpu2 = (1j * self.x0generators_oldpu2 + 3 * self.Zg2 / self.Zbase7) * (self.Sbase / nominalpower2)

            # Transformer 1 Grounding Information, Side 1
            if Zt1_grounding1 == "Solid ground":
                self.Zt1_value1 = 0
            elif Zt1_grounding1 == "Ungrounded":
                self.Zt1_value1 = float('inf')
            elif Zt1_grounding1 == "Resistor":
                self.Zt1_value1 = Zt1_value1
            elif Zt1_grounding1 == "N/A":
                self.Zt1_value1 = None
            else:
                print("Invalid Grounding Information For Transformer 1 Side 1")
                exit(-1)

            # Transformer 1 Grounding Information, Side 2
            if Zt1_grounding2 == "Solid ground":
                self.Zt1_value2 = 0
            elif Zt1_grounding2 == "Ungrounded":
                self.Zt1_value2 = float('inf')
            elif Zt1_grounding2 == "Resistor":
                self.Zt1_value2 = Zt1_value2
            elif Zt1_grounding2 == "N/A":
                self.Zt1_value2 = None
            else:
                print("Invalid Grounding Information For Transformer 1 Side 2")
                exit(-1)

            # Transformer 2 Grounding Information, Side 1
            if Zt2_grounding1 == "Solid ground":
                self.Zt2_value1 = 0
            elif Zt2_grounding1 == "Ungrounded":
                self.Zt2_value1 = float('inf')
            elif Zt2_grounding1 == "Resistor":
                self.Zt2_value1 = Zt2_value1
            elif Zt2_grounding1 == "N/A":
                self.Zt2_value1 = None
            else:
                print("Invalid Grounding Information For Transformer 2 Side 1")
                exit(-1)

            # Transformer 2 Grounding Information, Side 2
            if Zt2_grounding2 == "Solid ground":
                self.Zt2_value2 = 0
            elif Zt2_grounding2 == "Ungrounded":
                self.Zt2_value2 = float('inf')
            elif Zt2_grounding2 == "Resistor":
                self.Zt2_value2 = Zt2_value2
            elif Zt2_grounding2 == "N/A":
                self.Zt2_value2 = None
            else:
                print("Invalid Grounding Information For Transformer 2 Side 2")
                exit(-1)

            # This is Transformer Information -> Depends on grounding, initial remains here for reference only while I develop the proper code
            #self.Ybus0[6][5] = -1 / (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu)  # T2
            #self.Ybus0[5][6] = self.Ybus0[6][5]  # T2
            #self.Ybus0[0][1] = -1 / (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu)  # T1
            #self.Ybus0[1][0] = self.Ybus0[0][1]  # T1

            # Updating Ybus Depending on Transformer 1 Connection
            self.Ybus0 = Grid.Ybus

            # Grounded Wye can go to Grounded Wye or Delta to not be 0
            # Grounded Wye to Grounded Wye on Sequence Networks-> Set Ohms for each side using each side's base
            if Zt1_connection1 == "Grounded Wye":
                if Zt1_connection2 == "Grounded Wye":
                    self.Zt11 = Zt1_value1 * 3 / self.Zbase1
                    self.Zt12 = Zt1_value2 * 3 / self.Zbasemain
                    self.Ybus0[0][1] = -1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt11 + self.Zt12)  # T1
                    self.Ybus0[1][0] = self.Ybus0[0][1]  # T1
                    self.Ybus0[0][0] = 1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt11 + self.Zt12 + 1j * self.totalx0generators_newpu1)

                elif Zt1_connection2 == "Delta":
                    self.Zt11 = Zt1_value1 * 3 / self.Zbase1
                    self.Ybus0[0][1] = -1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt11)  # T1 Number
                    self.Ybus0[1][0] = 0  # T1
                    self.Ybus0[0][0] = 1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt11 + 1j * self.totalx0generators_newpu1)

                elif Zt1_connection2 == "Ungrounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1
                    self.Ybus0[0][0] = 1 / (1j * self.totalx0generators_newpu1)

                else:
                    print("Invalid Information For Transformer 1 Connection 2")
                    exit(-1)

            # Delta to Grounded Wye Only Non-Zero Combination
            elif Zt1_connection1 == "Delta":
                if Zt1_connection2 == "Grounded Wye":
                    #print("Transformer 1: Delta to Grounded Wye")
                    self.Zt12 = Zt1_value2 * 3 / self.Zbasemain
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = -1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt12)  # T1  # T1
                    self.Ybus0[0][0] = 1 / (1j * self.totalx0generators_newpu1)
                    #print("self.totalx0generators_newpu1", self.totalx0generators_newpu1)
                elif Zt1_connection2 == "Delta" or Zt1_connection2 == "Ungrounded Wye":
                    self.Ybus0[0][1] = 0  # T1
                    self.Ybus0[1][0] = 0  # T1
                    self.Ybus0[0][0] = 1 / (1j * self.totalx0generators_newpu1)

                else:
                    print("Invalid Information For Transformer 1 Connection 2")
                    exit(-1)
            # All Ungrounded Wye have 0
            elif Zt1_connection1 == "Ungrounded Wye":
                self.Ybus0[0][1] = 0  # T1
                self.Ybus0[1][0] = 0  # T1
                self.Ybus0[0][0] = 1 / (1j * self.totalx0generators_newpu1)

            else:
                print("Invalid Information For Transformer 1 Connection 1")
                exit(-1)

            # Updating Ybus Depending on Transformer 1 Connection
            # Grounded Wye can go to Grounded Wye or Delta to not be 0
            if Zt2_connection1 == "Grounded Wye":
                if Zt2_connection2 == "Grounded Wye":
                    self.Zt21 = Zt2_value1 * 3 / self.Zbase7
                    self.Zt22 = Zt2_value2 * 3 / self.Zbasemain
                    self.Ybus0[6][5] = -1 / (3 * (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu) + self.Zt21 + self.Zt22)  # T2
                    self.Ybus0[5][6] = self.Ybus0[6][5]
                    self.Ybus0[6][6] = 1 / (3 * (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu) + self.Zt21 + self.Zt22 + 1j * self.totalx0generators_newpu2)

                elif Zt2_connection2 == "Delta":
                    self.Zt21 = Zt2_value1 * 3 / self.Zbase7
                    self.Ybus0[6][5] = -1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt21)  # T2
                    self.Ybus0[5][6] = 0  # T2
                    self.Ybus0[6][6] = 1 / (3 * (Grid.transformers["T1"].Rpu + 1j * Grid.transformers["T1"].Xpu) + self.Zt21 + 1j * self.totalx0generators_newpu2)

                elif Zt2_connection2 == "Ungrounded Wye":
                    self.Ybus0[6][5] = 0  # T2
                    self.Ybus0[5][6] = 0  # T2
                    self.Ybus0[6][6] = 1 / (1j * self.totalx0generators_newpu2)

                else:
                    print("Invalid Information For Transformer 2 Connection 2")
                    exit(-1)

            # Delta to Grounded Wye Only Non-Zero Combination
            elif Zt2_connection1 == "Delta":
                if Zt2_connection2 == "Grounded Wye":
                    #print("Transformer 2: Delta to Grounded Wye")
                    self.Zt22 = Zt2_value2 * 3 / self.Zbasemain
                    self.Ybus0[6][5] = 0  # T2
                    self.Ybus0[5][6] = -1 / (3 * (Grid.transformers["T2"].Rpu + 1j * Grid.transformers["T2"].Xpu) + self.Zt22)  # T2
                    self.Ybus0[6][6] = 1 / (1j * self.totalx0generators_newpu2)

                elif Zt2_connection2 == "Delta" or Zt2_connection2 == "Ungrounded Wye":
                    self.Ybus0[6][5] = 0  # T2
                    self.Ybus0[5][6] = 0  # T2
                    self.Ybus0[6][6] = 1 / (1j * self.totalx0generators_newpu2)

                else:
                    print("Invalid Information For Transformer 2 Connection 2")
                    exit(-1)

            # All Ungrounded Wye have 0
            elif Zt2_connection1 == "Ungrounded Wye":
                self.Ybus0[6][5] = 0  # T2
                self.Ybus0[5][6] = 0  # T2
                self.Ybus0[6][6] = 1 / (1j * self.totalx0generators_newpu2)

            else:
                print("Invalid Information For Transformer 2 Connection 1")
                exit(-1)

            # Establish Lines -> Z0 = 3Z1
            # Set Non-diagonals just using -1/Z, Diagonals are same as original Ybus, or updated for generators in if segment

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

            # Set Diagonals, That do not have a generator attached
            self.Ybus0[1][1] /= 3
            self.Ybus0[2][2] /= 3
            self.Ybus0[3][3] /= 3
            self.Ybus0[4][4] /= 3
            self.Ybus0[5][5] /= 3

            # Print the Z0-bus matrix
            self.Zbus0 = np.linalg.inv(self.Ybus0)
            #print("Ybus0 Matrix")
            #self.printmatrix(self.Ybus0)
            #print("Z0-bus Matrix")
            #self.printmatrix(self.Zbus0)

    # Function to easily print Matrices
    def printmatrix(self, Matrix):
        i = 0
        while i < 7:
            j = 0
            print("\nRow " + str(i + 1))
            while j < 7:
                print(Matrix[i][j])
                j += 1
            i = i + 1
        print("\n")
