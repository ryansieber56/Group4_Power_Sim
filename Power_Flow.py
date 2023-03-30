import numpy as np
from DC_Power_Flow_Solver import DCPowerFlow
from Newton_Raphson_Power_Flow import NewtonRhapson
from Fast_Decoupled_Solver import FastDecoupled
class PowerFlow:

    # Power flow
    def __init__(self, result):
        if isinstance(result, DCPowerFlow):
            self.obj = DCPowerFlow
            print("DC Power Flow selected.")
        elif isinstance(result, NewtonRhapson):
            self.obj = NewtonRhapson
            print("Newton Rhapson Power Flow selected.")
        elif isinstance(result, FastDecoupled):
            self.obj = FastDecoupled
            print("Fast Decoupled Power Flow selected.")

        self.V_complex = self.obj.V_complex






# Milestone 3: Solve Power System Power Flow
# Due Wednesday, March 17th

# Students will use the previous input power flow data to solve the Newton-Raphson algorithm. This should be done with a tolerance of maximum mismatch of 0.0001 pu.
# Students will calculate the current flowing through each line and determine if that current exceeds line ampacity (this is on datasheet, remember that the conductors are bundled!)
# Students will calculate the power flowing through the lines on both the sending a receiving end
# Students will calculate the power loss in each line, and the total power loss for the system
# Students will solve for the real and reactive power injected into the slack bus and the reactive power injected into the PV bus





length = len(self.V_complex)
        # Calculate per unit current value
        I_values_per_unit = np.zeros((length, length), dtype=complex)

        for i in range(length):
            for j in range(length):
                if j <= i:
                    continue
                # double checked Ybus -> correct
                #I_values_per_unit[i, j] = (self.V_complex[j] - self.V_complex[i]) * (Grid.Ybus[i, j])
                I_values_per_unit[i, j] = (self.V_complex[i] - self.V_complex[j]) * (Grid.Ybus[i, j]) + (self.V_complex[i] ** 2 * Grid.transmissionline["L" + str(i + 1)].Bpu) / 2 - (self.V_complex[j] ** 2 * Grid.transmissionline["L" + str(i + 1)].Bpu) / 2

                # I_values_per_unit[i, j] = (V_complex[j]-V_complex[i]) * (Grid.Ybus[i, j]+ Grid.transmissionline["L"+str(i+1)].Bpu)
                print("YBUS", i, j, Grid.Ybus[i, j])
                # if abs(Grid.Ybus[i, j]) == 0:
                #    continue
                # print("V_complex[",j+1,"]-V_complex[",i+1,"]", V_complex[j]-V_complex[i])

        # Calculate actual current values
        self.I_values = (self.Sbase / (self.Vbase * np.sqrt(3))) * I_values_per_unit  # / number_of_bundles
        if self.slack_bus == 0:
            self.I_values[0, 1] *= self.Vbase / 20000
        if self.slack_bus == 6:
            self.I_values[5, 6] *= self.Vbase / 18000
        # BUS4 to BUS5 is wrong -> Future code is also incorrect

        # print("I_values_per_unit")
        # print(I_values_per_unit)
        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V
        self.P_values = np.zeros((length, length), dtype=complex)
        self.Q_values = np.zeros((length, length), dtype=complex)
        self.S_values = np.zeros((length, length), dtype=complex)
        I_conjugate = np.conjugate(self.I_values)

        # Take voltage values out of per unit
        self.V_complex *= self.Vbase

        # Calculate S, P, and Q values
        for i in range(length):
            for j in range(length):
                self.S_values[i,j] = I_conjugate[i,j]*self.V_complex[i]
        if self.slack_bus == 0:
            self.S_values[0, 1] *= self.Vbase / 20000
        if self.slack_bus == 6:
            self.S_values[5, 6] *= self.Vbase / 18000

        self.P_values = np.real(self.S_values)
        self.Q_values = np.imag(self.S_values)

        # Power and System Loss
        # exit("Rest of code incomplete")
        # Calculate power loss in each line
        P_loss = np.zeros((length, length), dtype=complex)

        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        for i1 in range(length):
            i2 = 0
            for i, value in Grid.transmissionline.items():
                P_loss[i1, i2] = abs(self.I_values[i1, i2]) * abs(self.I_values[i1, i2]) * value.Rtotal
                # Store the power loss in that transmission line
                Grid.store_power_loss(i, P_loss[i1, i2])
                i2 += 1

        # Find Total System loss
        self.system_loss = 0

        for i in range(length):
            self.system_loss += P_loss[i]
