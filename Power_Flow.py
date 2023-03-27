# Power Flow Class
import numpy as np
import pandas as pd

from Grid import Grid


class PowerFlow:

    # Power flow
    def __init__(self, Grid):
        self.Parr = []
        self.Qarr = []
        self.convergencemet = 0
        convergencevalue = 0.0001
        Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        number_of_bundles = Grid.transmissionline[
            list(Grid.transmissionline.keys())[0]].numberofbundles  # Obtain number of bundles
        Sbase = 100000000  # VA

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
            P_mismatch = P_mismatch[1:7]

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
            delta_correction = np.concatenate((delta_correction[:0], [0], delta_correction[0:]), axis=0)

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

        # Calculate Current Flowing
        # Set up array of zeros
        V_complex = np.zeros(len(V), dtype=complex)

        # Calculate the complex voltage
        for i in range(len(V)):
            V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])
            # print("V_complex for bus", i+1, ":", V_complex[i])
            # Print Values used -> double-checked and all voltages and angles are correct
            # print(i, " V ", V[i], "delta[i]", delta[i])

        # Calculate per unit current value
        I_values_per_unit = np.zeros((len(V), len(V)), dtype=complex)
        for i in range(len(V_complex)):
            for j in range(len(V_complex)):
                if j <= i:
                    continue
                # double checked Ybus -> correct
                I_values_per_unit[i, j] = (V_complex[j] - V_complex[i]) * (Grid.Ybus[i, j])

                # I_values_per_unit[i, j] = (V_complex[j]-V_complex[i]) * (Grid.Ybus[i, j]+ Grid.transmissionline["L"+str(i+1)].Bpu)
                print("YBUS", i, j, Grid.Ybus[i, j])
                # if abs(Grid.Ybus[i, j]) == 0:
                #    continue
                # print("V_complex[",j+1,"]-V_complex[",i+1,"]", V_complex[j]-V_complex[i])

        # Calculate actual current values
        self.I_values = (Sbase / (Vbase * np.sqrt(3))) * I_values_per_unit  # / number_of_bundles
        if self.slack_bus == 0:
            self.I_values[0, 1] *= Vbase / 20000
        if self.slack_bus == 6:
            self.I_values[5, 6] *= Vbase / 18000
        # BUS4 to BUS5 is wrong -> Future code is also incorrect








        # print("I_values_per_unit")
        # print(I_values_per_unit)
        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V
        self.P_values = np.zeros((len(V), len(V)), dtype=complex)
        self.Q_values = np.zeros((len(V), len(V)), dtype=complex)
        self.S_values = np.zeros((len(V), len(V)), dtype=complex)
        I_conjugate = np.conjugate(self.I_values)

        # Take voltage values out of per unit
        V_complex *= Vbase

        # Calculate S, P, and Q values
        self.S_values = np.dot(I_conjugate, V_complex)
        if self.slack_bus == 0:
            self.S_values[0, 1] *= Vbase / 20000
        if self.slack_bus == 6:
            self.S_values[5, 6] *= Vbase / 18000

        self.P_values = np.real(self.S_values)
        self.Q_values = np.imag(self.S_values)


        # Power and System Loss
        # exit("Rest of code incomplete")
        # Calculate power loss in each line
        P_loss = np.zeros((7, 7), dtype=complex)
        i1 = 0
        i2 = 0

        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        if i1 < 7:
            for i, value in Grid.transmissionline.items():
                P_loss[i1, i2] = abs(self.I_values[i1, i2]) * abs(self.I_values[i1, i2]) * value.Rtotal
                # Store the power loss in that transmission line
                Grid.store_power_loss(i, P_loss[i1, i2])
                i2 += 1
            i1 += 1

        # Find Total System loss
        self.system_loss = 0

        for i in range(len(P_loss)):
            self.system_loss += P_loss[i]

        # Write to Excel file/ Print
        # Print Values
        print("\nI_Values")
        for i in range(len(V_complex)):
            for j in range(len(V_complex)):
                if abs(self.I_values[i, j]) == 0:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", self.I_values[i, j], abs(self.I_values[i, j]))

        print("\nP_values")
        for i in range(len(V_complex)):
            for j in range(len(V_complex)):
                if abs(self.P_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.P_values[i, j] / 100000000, "MW")

        print("\nQ_values")
        for i in range(len(V_complex)):
            for j in range(len(V_complex)):
                if abs(self.Q_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.Q_values[i, j] / 100000000, "MVAR")

        print("\nPower Loss Per Line")
        for i, value in Grid.transmissionline.items():
            print(i, value.powerloss)

        print("\nTotal System Power Loss:", self.system_loss)

# Milestone 3: Solve Power System Power Flow
# Due Wednesday, March 17th

# Students will use the previous input power flow data to solve the Newton-Raphson algorithm. This should be done with a tolerance of maximum mismatch of 0.0001 pu.
# Students will calculate the current flowing through each line and determine if that current exceeds line ampacity (this is on datasheet, remember that the conductors are bundled!)
# Students will calculate the power flowing through the lines on both the sending a receiving end
# Students will calculate the power loss in each line, and the total power loss for the system
# Students will solve for the real and reactive power injected into the slack bus and the reactive power injected into the PV bus
