# File to Perform Fault Calculations
import numpy as np
class FaultCalculation:
    # Draw Zero Sequence for submission
    # Fault Calculation
    def __init__(self, Grid, SeqNet, faulttype: str, faultlocation: int, faulting_impedance: float):
        self.length = len(Grid.buses)
        self.V_0 = np.zeros(self.length, dtype=complex)
        self.V_1 = np.zeros(self.length, dtype=complex)
        self.V_2 = np.zeros(self.length, dtype=complex)
        self.V_012 = np.zeros((self.length, 3), dtype=complex)
        self.V_abc = np.zeros((self.length, 3), dtype=complex)
        self.In_1 = None
        self.In_2 = None
        self.In_0 = None
        self.I_012 = None
        self.I_abc = None

        self.Znn_1 = SeqNet.Zbus1[faultlocation-1][faultlocation-1]
        self.Znn_2 = SeqNet.Zbus2[faultlocation-1][faultlocation-1]
        self.Znn_0 = SeqNet.Zbus0[faultlocation-1][faultlocation-1]
        self.a1 = -0.5 + 1j * 0.8660254
        self.a2 = -0.5 - 1j * 0.8660254
        self.A = [[1, 1, 1], [1, self.a2, self.a1], [1, self.a1, self.a2]]
        # Vf will always be 1 per unit because we are neglecting pre fault currents
        self.VF = 1

        # A bolted fault will have a Zf of 0
        self.Zf = faulting_impedance

        # Using the Fault Type Information, Calculate the Currents
        if faulttype == "Symmetrical Fault":
            self.In_1 = self.VF/self.Znn_1
            self.In_2 = 0
            self.In_0 = 0

        elif faulttype == "Single Line to Ground Fault":
            self.In_0 = self.VF/((self.Znn_0 +self.Znn_1 + self.Znn_2) + 3 * self.Zf)
            self.In_1 = self.In_0
            self.In_2 = self.In_0

        elif faulttype == "Line to Line Fault":
            self.In_0 = 0
            self.In_1 = self.VF/(self.Znn_1 + self.Znn_2 + self.Zf)
            self.In_2 = -self.In_1

        elif faulttype == "Double Line to Ground Fault":
            self.In_1 = self.VF/(self.Znn_1+(self.Znn_2 * (self.Znn_0 + 3 * self.Zf))/(self.Znn_2 + self.Znn_0+3 * self.Zf))
            self.In_2 = -self.In_1 * (self.Znn_0 + 3 * self.Zf) / (self.Znn_0 + 3 * self.Zf + self.Znn_2)
            self.In_0 = -self.In_1 * self.Znn_2 / (self.Znn_0 + 3 * self.Zf + self.Znn_2)

        else:
            print("Invalid Fault Type")
            exit(-1)

        # Set I_012
        self.I_012 = np.array([[self.In_0], [self.In_1], [self.In_2]])

        # Solve for all V_012 Values and V-012
        for k in range(self.length):
            self.V_0[k] = -(SeqNet.Zbus0[k][faultlocation-1] * self.In_0)
            self.V_1[k] = self.VF - (SeqNet.Zbus1[k][faultlocation-1] * self.In_1)
            self.V_2[k] = -(SeqNet.Zbus2[k][faultlocation-1] * self.In_2)
            self.V_012[k] = np.array([[self.V_0[k]], [self.V_1[k]], [self.V_2[k]]]).flatten()

        # Calculate I_abc and V_abc
        self.I_abc = np.dot(self.A, self.I_012)
        for k in range(self.length):
            self.V_abc[k] = np.dot(self.A, self.V_012[k])
        print("\nI_abc")
        for i in range(len(self.I_abc)):
            print(self.I_abc[i])
        print("\nAll V_abc")
        for k in range(self.length):
            print("V_abc[", k, "] = ", self.V_abc[k])
