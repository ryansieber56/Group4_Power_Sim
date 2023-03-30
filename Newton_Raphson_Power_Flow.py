# Newton_Rhapson
import numpy as np
import pandas as pd

from Grid import Grid


class NewtonRhapson:

    # Power flow
    def __init__(self, Grid):
        self.P_loss = None
        self.system_loss = None
        self.S_values = np.zeros((7, 7), dtype=complex)
        self.Q_values = np.zeros((7, 7))
        self.P_values = np.zeros((7, 7))
        self.I_values = np.zeros((7, 7), dtype=complex)
        self.Parr = []
        self.Qarr = []
        self.convergencemet = 0
        convergencevalue = 0.0001
        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.number_of_bundles = Grid.transmissionline[
            list(Grid.transmissionline.keys())[0]].numberofbundles  # Obtain number of bundles
        self.Sbase = 100000000  # VA
        self.Vbase_slack1 = 20000
        self.Vbase_slack2 = 18000
        # set blank array for given values
        P_given = np.zeros(len(Grid.buses))
        Q_given = np.zeros(len(Grid.buses))

        # find slack bus
        self.slack_bus = None
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i
                break

        # set given values
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = -Grid.buses["Bus" + str(i + 1)].Q / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = Grid.buses["Bus" + str(i + 1)].Q / 100

        # set intial guess
        # delta will be in radians because cos and sin take radians as arguments, and np.angle will return the Ybus angle in radians
        V = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        delta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # establish array for calculated value, without slack bus
        Parr = np.zeros(len(Grid.buses))
        Qarr = np.zeros(len(Grid.buses))
        iteration = 1
        while self.convergencemet == 0 and iteration < 30:
            # print("\n\nIteration = ", iteration)
            for i in range(len(Grid.buses)):
                Parr[i] = 0
                Qarr[i] = 0
            # calculate mismatch, ignoring the slack bus
            for i in range(len(Grid.buses)):
                # if slack bus, skip
                if i == self.slack_bus:
                    continue

                for j in range(len(Grid.buses)):
                    Parr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.cos(
                        delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    if i == self.voltage_controlled_bus:
                        continue
                    Qarr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.sin(
                        delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # P does not include slack bus
            P_mismatch = P_given - Parr
            if self.slack_bus == 0:
                P_mismatch = P_mismatch[1:7]
            if self.slack_bus == 6:
                P_mismatch = P_mismatch[0:6]

            # Q not include slack bus or VCB
            Q_mismatch = Q_given - Qarr
            Q_mismatch = Q_mismatch[1:6]

            # print("P_mismatch")
            # print(P_mismatch)
            # print("Q_mismatch")
            # print(Q_mismatch)

            self.convergencemet = 1

            # Check Power Mismatch for Convergence
            mismatchPQ = np.concatenate((P_mismatch, Q_mismatch))
            for i in range(len(mismatchPQ)):
                if mismatchPQ[i] > convergencevalue:
                    self.convergencemet = 0
                    break

            if self.convergencemet == 1:
                print("It converged after iteration ", iteration - 1)
                break

            # Calculate Jacobian Matrix
            J11 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
            J12 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
            J21 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
            J22 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))

            skipterm = 0

            for i in range(len(Grid.buses)):
                # if slack bus skip
                if i == self.slack_bus:
                    skipterm = 1
                    continue
                for j in range(len(Grid.buses)):
                    if i == j:
                        for z in range(len(Grid.buses)):
                            J12[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))
                            J22[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                            if z == i:
                                continue

                            J11[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))
                            J21[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                        # Non-Summation Modification
                        J11[i - skipterm, j - skipterm] = -J11[i - skipterm, j - skipterm] * V[i]
                        J12[i - skipterm, j - skipterm] = J12[i - skipterm, j - skipterm] + (
                                    V[i] * abs(Grid.Ybus[i, j]) * np.cos(np.angle(Grid.Ybus[i, j])))
                        J21[i - skipterm, j - skipterm] = (J21[i - skipterm, j - skipterm] * V[i])
                        J22[i - skipterm, j - skipterm] = J22[i - skipterm, j - skipterm] - (
                                    V[i] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))

                    else:
                        J11[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J12[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.cos(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J21[i - skipterm, j - skipterm] = -V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.cos(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J22[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # Combine Jacobian
            J = np.block([[J11, J12], [J21, J22]])

            # Delete slack bus 1 or 7, here if it is 1
            J_temp = np.delete(J, 11, 1)
            J = J_temp
            J_temp = np.delete(J, 11, 0)
            J = J_temp

            # print("Jacobian Matrix", iteration-1)
            # i = 0
            # while i < len(J):
            #    j = 0
            #    print("\nRow " + str(i + 1))
            #    while j < len(J):
            #        print(J[i][j])
            #        j += 1
            #    i = i + 1

            # calculate change in voltage and phase angle
            J_inv = np.linalg.inv(J)

            # combine mismatches
            mismatch = np.concatenate((P_mismatch, Q_mismatch))
            correction = np.dot(J_inv, mismatch)
            # print("correction")
            # print(correction)
            delta_correction = correction[:6]
            V_correction = correction[6:]

            # do not update the angle of the slack bus
            delta_correction = np.concatenate((delta_correction[:self.slack_bus], [0], delta_correction[self.slack_bus:]), axis=0)

            # do not update the voltage of the slack bus or voltage controlled bus
            V_correction = np.concatenate((V_correction[:6], [0], V_correction[6:]), axis=0)
            V_correction = np.concatenate((V_correction[:0], [0], V_correction[0:]), axis=0)

            # print("Delta Correction Data")
            # print(delta_correction)
            # print("V Correction Data")
            # print(V_correction)

            delta += delta_correction
            V += V_correction

            # print("Delta Updated Data")
            # print(delta)
            # print("V Updated Data")
            # print(V)
            iteration += 1

            # Calculate if VAR Limit has been exceeded
            #if Q_k >
        # Calculate Current Flowing
        # Set up array of zeros
        self.V_complex = np.zeros(len(V), dtype=complex)

        # Calculate the complex voltage
        for i in range(len(V)):
            self.V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])
            # print("V_complex for bus", i+1, ":", V_complex[i])
            # Print Values used -> double-checked and all voltages and angles are correct
            # print(i, " V ", V[i], "delta[i]", delta[i])
        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):
        length = len(self.V_complex)
        # Calculate per unit current value
        I_values_per_unit = np.zeros((length, length), dtype=complex)
        for i in range(length):
            for j in range(length):
                if j <= i:
                    continue
                I_values_per_unit[i, j] = (self.V_complex[j] - self.V_complex[i]) * (Grid.Ybus[i, j])

        # Calculate actual current values
        self.I_values = (self.Sbase / (self.Vbase * np.sqrt(3))) * I_values_per_unit
        if self.slack_bus == 0:
            self.I_values[0, 1] *= self.Vbase / self.Vbase_slack1 # -> Because different bases for slack bus
        if self.slack_bus == 6:
            self.I_values[5, 6] *= self.Vbase / self.Vbase_slack2 # -> Because different bases for slack bus

        # Check ampacity limit
        self.check_ampacity()

        # Take voltage values out of per unit
        self.V_complex *= self.Vbase
        I_conjugate = np.conjugate(self.I_values)
        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V * sqrt(3)
        for i in range(length):
            for j in range(length):
                self.S_values[i,j] = self.V_complex[i] * I_conjugate[i,j] * np.sqrt(3) # because 3 phase
                if self.slack_bus == 0 and i == 0:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack1
                if self.slack_bus == 6 and i == 6:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack2

        self.P_values = np.real(self.S_values)
        self.Q_values = np.imag(self.S_values)

        self.P_loss = np.zeros((length, length))

        self.system_loss = 0
        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        for i1 in range(length):
            for i2 in range(length):
                for i, value in Grid.transmissionline.items():
                    if value.bus1 == ("Bus" + str(i1+1)) and value.bus2 == ("Bus" + str(i2+1)):
                        powerloss = abs(self.I_values[i1, i2]) ** 2 * value.Rtotal * 3 # because 3 phase
                        Grid.store_power_loss(i, powerloss)
                        self.system_loss += powerloss

        # Same for transformers
        for i1 in range(length):
            for i2 in range(length):
                for i, value in Grid.transformers.items():
                    if value.bus1 == ("Bus" + str(i1+1)) and value.bus2 == ("Bus" + str(i2+1)) or value.bus1 == ("Bus" + str(i2+1)) and value.bus2 == ("Bus" + str(i1+1)):
                        if abs(self.I_values[i1, i2]) == 0:
                            continue
                        powerloss = abs(self.I_values[i1, i2]) ** 2 * value.Rpu * self.Vbase ** 2/self.Sbase * 3 # because 3 phase
                        if(self.slack_bus == 0 and value.name == "T1"):
                            powerloss = powerloss / (self.Vbase **2) * self.Vbase_slack1 ** 2
                        if(self.slack_bus == 6 and value.name == "T2"):
                            powerloss = powerloss / (self.Vbase **2) * self.Vbase_slack2 ** 2
                        Grid.store_power_loss_transformer(i, powerloss)
                        # print("Power Loss", i, powerloss)
                        self.system_loss += powerloss

        # NEED HELP WITH POWER AND Q INJECTIONS FOR SLACK + PV
        if self.slack_bus == 0:
            self.P_inj_slack = abs(self.P_values[0, 1])
            self.Q_inj_slack = abs(self.Q_values[0, 1])
            print("self.P_inj_slack:", self.P_inj_slack)
            print("self.Q_inj_slack:", self.Q_inj_slack)
            self.Q_inj_PV = 0
            for i in range(length):
                self.Q_inj_PV += Grid.Ybus[6, i] * abs(self.V_complex[i]) * np.sin(np.angle(self.V_complex[6])-np.angle(self.V_complex[i] - np.angle(Grid.Ybus[6,i]))) / self.Vbase
            self.Q_inj_PV *= self.Sbase * self.Vbase_slack2 **2
            print("self.Q_inj_PV:", abs(self.Q_inj_PV))

        if self.slack_bus == 6:
            self.P_inj_slack = abs(self.P_values[5, 6])
            self.Q_inj_slack = abs(self.Q_values[5, 6])
            print("self.P_inj_slack:", self.P_inj_slack)
            print("self.Q_inj_slack:", self.Q_inj_slack)
            self.Q_inj_PV = 0
            for i in range(length):
                self.Q_inj_PV = Grid.Ybus[0, i] * abs(self.V_complex[i]) * np.sin(
                    -np.angle(self.V_complex[i] - np.angle(Grid.Ybus[0, i]))) / self.Vbase
            self.Q_inj_PV *= self.Vbase_slack1 ** 2 / 3
            print("self.Q_inj_PV:", abs(self.Q_inj_PV))

        # Write to Excel file/ Print
        # Print Values
        print("\nI_Values")
        for i in range(length):
            for j in range(length):
                if abs(self.I_values[i, j]) == 0:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", self.I_values[i, j], abs(self.I_values[i, j]))

        print("\nS_values")
        for i in range(length):
            for j in range(length):
                if abs(self.S_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.S_values[i, j] / 1000000, "MVA", abs(self.S_values[i,j])/1000000, "MVA")

        print("\nP_values")
        for i in range(length):
            for j in range(length):
                if abs(self.P_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.P_values[i, j] / 1000000, "MW")

        print("\nQ_values")
        for i in range(length):
            for j in range(length):
                if abs(self.Q_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.Q_values[i, j] / 1000000, "MVAR")

        print("\nPower Loss Per Line")
        for i, value in Grid.transmissionline.items():
            print(i, value.powerloss/1000000, "MW")

        print("\nPower Loss in Transformers")
        for i, value in Grid.transformers.items():
            print(i, value.powerloss/1000000, "MW")

        print("\nTotal System Power Loss:", self.system_loss/1000000, "MW")

    def check_ampacity(self):
        i = 1
        j = 1
        while i < 6:
            while j < 6:
                if abs(self.I_values[i, j]) == 0:
                    j += 1
                    continue
                if self.I_values[i][j] >= 475:
                    print("Error: Ampacity limit exceeded")
                    exit(-1)
                j += 1
            j = 0
            i += 1

# Students will solve for the real and reactive power injected into the slack bus and the reactive power injected into the PV bus
