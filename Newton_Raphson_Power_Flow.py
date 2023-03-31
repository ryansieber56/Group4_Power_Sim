# Newton_Rhapson
import numpy as np
import pandas as pd

from Grid import Grid


class NewtonRhapson:

    # Power flow
    def __init__(self, Grid):
        self.length = len(Grid.buses)
        self.P_loss = np.zeros((self.length, self.length))
        self.system_loss = 0
        self.P_inj_slack = 0
        self.Q_inj_slack = 0
        self.Q_inj_PV = 0
        self.Q_k_limit = 175000000
        self.Q_limit_passed = 0
        self.S_values = np.zeros((self.length, self.length), dtype=complex)
        self.Q_values = np.zeros((self.length, self.length))
        self.P_values = np.zeros((self.length, self.length))
        self.I_values = np.zeros((self.length, self.length), dtype=complex)
        self.V_complex = np.zeros(self.length, dtype=complex)
        self.Parr = []
        self.Qarr = []
        self.convergencemet = 0
        self.convergencevalue = 0.0001
        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.number_of_bundles = Grid.transmissionline[
            list(Grid.transmissionline.keys())[0]].numberofbundles  # Obtain number of bundles
        self.Sbase = 100000000  # VA
        self.Vbase_slack1 = 20000
        self.Vbase_slack2 = 18000
        # set blank array for given values
        P_given = np.zeros(self.length)
        Q_given = np.zeros(self.length)

        # find slack bus
        self.slack_bus = None
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i

        # set given values
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = -Grid.buses["Bus" + str(i + 1)].Q / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = Grid.buses["Bus" + str(i + 1)].Q / 100

        # set intial guess
        V = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        delta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # establish array for calculated value, without slack bus
        Parr = np.zeros(self.length)
        Qarr = np.zeros(self.length)
        iteration = 1
        while self.convergencemet == 0 and iteration < 30:

            for i in range(self.length):
                Parr[i] = 0
                Qarr[i] = 0

            # calculate mismatch, ignoring the slack bus
            for i in range(self.length):
                # if slack bus, skip
                if i == self.slack_bus:
                    continue

                for j in range(self.length):
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

            self.convergencemet = 1

            # Check Power Mismatch for Convergence
            mismatchPQ = np.concatenate((P_mismatch, Q_mismatch))
            for i in range(len(mismatchPQ)):
                if mismatchPQ[i] > self.convergencevalue:
                    self.convergencemet = 0
                    break

            if self.convergencemet == 1:
                print("It converged after iteration ", iteration - 1)
                break

            # Calculate Jacobian Matrix
            J11 = np.zeros((self.length - 1, self.length - 1))
            J12 = np.zeros((self.length - 1, self.length - 1))
            J21 = np.zeros((self.length - 1, self.length - 1))
            J22 = np.zeros((self.length - 1, self.length - 1))

            skipterm = 0

            for i in range(self.length):
                # if slack bus skip
                if i == self.slack_bus:
                    skipterm = 1
                    continue
                for j in range(self.length):
                    if skipterm == 0 and j == self.slack_bus:
                        continue
                    if i == j:
                        for z in range(self.length):
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

            # Delete Voltage Controlled bus 1 or 7, here if it is 7
            if self.slack_bus == 0:
                J_temp = np.delete(J, 11, 1)
                J = J_temp
                J_temp = np.delete(J, 11, 0)
                J = J_temp

            # Here VCB is 1
            if self.slack_bus == 6:
                J_temp = np.delete(J, 5, 1)
                J = J_temp
                J_temp = np.delete(J, 5, 0)
                J = J_temp

            # calculate change in voltage and phase angle
            J_inv = np.linalg.inv(J)

            # combine mismatches
            mismatch = np.concatenate((P_mismatch, Q_mismatch))
            correction = np.dot(J_inv, mismatch)

            delta_correction = correction[:6]
            V_correction = correction[6:]

            # do not update the angle of the slack bus
            delta_correction = np.concatenate((delta_correction[:self.slack_bus], [0], delta_correction[self.slack_bus:]), axis=0)

            # do not update the voltage of the slack bus or voltage controlled bus
            V_correction = np.concatenate((V_correction[:6], [0], V_correction[6:]), axis=0)
            V_correction = np.concatenate((V_correction[:0], [0], V_correction[0:]), axis=0)

            delta += delta_correction
            V += V_correction

            iteration += 1

            self.Q_k = 0
            for i in range(self.length):
                self.Q_k += abs(Grid.Ybus[self.voltage_controlled_bus, i]) * V[i] * np.sin(delta[self.voltage_controlled_bus] - delta[i] - np.angle(Grid.Ybus[self.voltage_controlled_bus,i]))
            self.Q_k *= V[self.voltage_controlled_bus]
            self.Q_k *= self.Sbase
            if self.Q_k > self.Q_k_limit:
                print("EXCEEDED")
                self.solve_exceeded_var_power_flow(Grid)
                self.Q_limit_passed = 1
                break

        print("SETUP V_Complex")
        # Check to make sure V_complex was not set up by Var limit solver
        if self.Q_limit_passed == 0:
            # Calculate the complex voltage
            for i in range(self.length):
                self.V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])
            # print("V_complex for bus", i+1, ":", V_complex[i])
            # Print Values used -> double-checked and all voltages and angles are correct
            # print(i, " V ", V[i], "delta[i]", delta[i])
        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):

        # Calculate per unit current value
        I_values_per_unit = np.zeros((self.length, self.length), dtype=complex)
        for i in range(self.length):
            for j in range(self.length):
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
        for i in range(self.length):
            for j in range(self.length):
                self.S_values[i, j] = self.V_complex[i] * I_conjugate[i, j] * np.sqrt(3)  # because 3 phase
                if self.slack_bus == 0 and i == 0:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack1
                if self.slack_bus == 6 and i == 6:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack2

        self.P_values = np.real(self.S_values)
        self.Q_values = np.imag(self.S_values)

        self.system_loss = 0

        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        for i1 in range(self.length):
            for i2 in range(self.length):
                for i, value in Grid.transmissionline.items():
                    if value.bus1 == ("Bus" + str(i1+1)) and value.bus2 == ("Bus" + str(i2+1)):
                        powerloss = abs(self.I_values[i1, i2]) ** 2 * value.Rtotal * 3 # because 3 phase
                        Grid.store_power_loss(i, powerloss)
                        self.system_loss += powerloss

        # Same for transformers
        for i1 in range(self.length):
            for i2 in range(self.length):
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
                        self.system_loss += powerloss

        # POWER AND Q INJECTIONS FOR SLACK + PV
        for i in range(self.length):
            self.Q_inj_slack += abs(Grid.Ybus[self.slack_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.sin(np.angle(self.V_complex[self.slack_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.slack_bus, i]))
            self.P_inj_slack += abs(Grid.Ybus[self.slack_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.cos(np.angle(self.V_complex[self.slack_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.slack_bus, i]))
            self.Q_inj_PV += abs(Grid.Ybus[self.voltage_controlled_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.sin(np.angle(self.V_complex[self.voltage_controlled_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.voltage_controlled_bus, i]))
        self.P_inj_slack *= abs(self.V_complex[self.slack_bus])/self.Vbase * self.Sbase
        self.Q_inj_slack *= abs(self.V_complex[self.slack_bus])/self.Vbase * self.Sbase
        self.Q_inj_PV *= abs(self.V_complex[self.voltage_controlled_bus])/self.Vbase * self.Sbase
        print("self.P_inj_slack:", self.P_inj_slack/1000000, "MW")
        print("self.Q_inj_slack:", self.Q_inj_slack/1000000, "MVAR")
        print("self.Q_inj_PV:", self.Q_inj_PV/1000000, "MVAR")

        # Write to Excel file/ Print
        # Print Values
        print("\nI_Values")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.I_values[i, j]) == 0:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", abs(self.I_values[i, j]))

        #print("\nS_values")
        #for i in range(self.length):
        #    for j in range(self.length):
        #        if abs(self.S_values[i, j]) == 0:
        #            continue
        #        print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.S_values[i, j] / 1000000, "MVA", abs(self.S_values[i,j])/1000000, "MVA")

        print("\nP_values")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.P_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.P_values[i, j] / 1000000, "MW")

        print("\nQ_values")
        for i in range(self.length):
            for j in range(self.length):
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

    def solve_exceeded_var_power_flow(self, Grid):
        # set blank array for given values
        P_given = np.zeros(self.length)
        Q_given = np.zeros(self.length)

        # find slack bus
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i

        # set given values
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = -Grid.buses["Bus" + str(i + 1)].Q / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = Grid.buses["Bus" + str(i + 1)].Q / 100
        Q_given[self.voltage_controlled_bus] = 175000000 #self.Q_k_limit
        # Set

        # set intial guess
        V = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        delta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # establish array for calculated value, without slack bus
        Parr = np.zeros(self.length)
        Qarr = np.zeros(self.length)
        iteration = 1

        while self.convergencemet == 0 and iteration < 30:
            # print("\n\nIteration = ", iteration)
            for i in range(self.length):
                Parr[i] = 0
                Qarr[i] = 0
            # calculate mismatch, ignoring the slack bus
            for i in range(self.length):
                # if slack bus, skip
                if i == self.slack_bus:
                    continue

                for j in range(self.length):
                    Parr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.cos(
                        delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    if i == self.voltage_controlled_bus:
                        continue
                    Qarr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.sin(
                        delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # P or Q does not include slack bus
            P_mismatch = P_given - Parr
            Q_mismatch = Q_given - Qarr
            if self.slack_bus == 0:
                P_mismatch = P_mismatch[1:7]
                Q_mismatch = Q_mismatch[1:7]

            if self.slack_bus == 6:
                P_mismatch = P_mismatch[0:6]
                Q_mismatch = Q_mismatch[0:6]

            self.convergencemet = 1

            # Check Power Mismatch for Convergence
            mismatchPQ = np.concatenate((P_mismatch, Q_mismatch))
            for i in range(len(mismatchPQ)):
                if mismatchPQ[i] > self.convergencevalue:
                    self.convergencemet = 0
                    break

            if self.convergencemet == 1:
                print("It converged after iteration ", iteration - 1)
                break

            # Calculate Jacobian Matrix
            J11 = np.zeros((self.length - 1, self.length - 1))
            J12 = np.zeros((self.length - 1, self.length - 1))
            J21 = np.zeros((self.length - 1, self.length - 1))
            J22 = np.zeros((self.length - 1, self.length - 1))

            skipterm = 0

            for i in range(self.length):
                # if slack bus skip
                if i == self.slack_bus:
                    skipterm = 1
                    continue
                for j in range(self.length):
                    if i == j:
                        for z in range(self.length):
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


            delta_correction = correction[:6]
            V_correction = correction[6:]

            # do not update the angle of the slack bus
            delta_correction = np.concatenate((delta_correction[:self.slack_bus], [0], delta_correction[self.slack_bus:]), axis=0)

            # do not update the voltage of the slack bus
            V_correction = np.concatenate((V_correction[:self.slack_bus], [0], V_correction[self.slack_bus:]), axis=0)

            # Correction
            delta += delta_correction
            V += V_correction
            iteration += 1

        # Setup V_complex
        for i in range(self.length):
            self.V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])

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