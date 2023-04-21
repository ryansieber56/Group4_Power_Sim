# Fast_Decoupled_Solver
import numpy as np
import pandas as pd

from Grid import Grid


class FastDecoupled:

    # Fast Decoupled
    def __init__(self, Grid):

        # Initial parameters that will be used later
        self.length = len(Grid.buses)
        self.system_loss = 0
        self.P_inj_slack = 0
        self.Q_inj_slack = 0
        self.Q_inj_PV = 0
        self.Q_limit_passed = 0
        self.convergencemet = 0
        iteration = 1

        # Set to 0, if becomes 1 means an adjustment was recently made with the capacitor
        # bank and the program should reiterate
        self.capacitor_bank_adjustment = 0

        # B is added capacitance if VAR limit is exceeded
        self.B = 0

        # 1 is to use capacitor bank correction on an exceeded Var limit, 0 is to get rid
        # of "PV" status of the generator's bus
        self.add_cap = 1

        # Take the Q limit and convergence value from the Grid file
        self.Q_k_limit = Grid.Q_k_limit # In units of VA
        self.convergencevalue = Grid.convergencevalue


        # Base Unit Parameters
        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.Sbase = Grid.Sbase * 1000000  # VA
        self.Vbase_slack1 = Grid.transformers[list(Grid.transformers.keys())[0]].v1rated * 1000  # In Volts
        self.Vbase_slack2 = Grid.transformers[list(Grid.transformers.keys())[1]].v1rated * 1000  # In Volts

        # Initialize Empty Arrays/ Matrices
        self.S_values = np.zeros((self.length, self.length), dtype=complex)
        self.Q_values = np.zeros((self.length, self.length))
        self.P_values = np.zeros((self.length, self.length))
        self.I_values = np.zeros((self.length, self.length), dtype=complex)
        self.V_complex = np.zeros(self.length, dtype=complex)
        self.Parr = []
        self.Qarr = []
        self.P_loss = np.zeros((self.length, self.length))
        P_given = np.zeros(self.length)
        Q_given = np.zeros(self.length)
        P_mismatch = np.zeros(self.length-1)
        Q_mismatch = np.zeros(self.length-2)

        # Array for calculated P, Q, and Jacobian values
        Parr = np.zeros(self.length)
        Qarr = np.zeros(self.length)
        J11 = np.zeros((self.length - 1, self.length - 1))
        J22 = np.zeros((self.length - 1, self.length - 1))

        # Find slack bus and voltage controlled bus by looping through Buses Dictionary in Grid Class
        self.slack_bus = None
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i

        # Set given values by looping through Buses Dictionary in Grid Class, Loads are negative
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = -Grid.buses["Bus" + str(i + 1)].Q / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
                Q_given[i] = Grid.buses["Bus" + str(i + 1)].Q / 100

        # Add a statement to turn Q_given into self. so that it can be used in other functions
        self.Q_given = Q_given

        # Set flat start
        V = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        delta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # Make sure convergence has not been met and that a capacitor bank was not just added, or that the max
        # Iterations has not been exceeded
        while (self.convergencemet == 0 and self.capacitor_bank_adjustment == 1) or iteration < 30:

            # Reset the calculated values every time through
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

            # P does not include slack bus, Q does not include Voltage controlled bus or slack bus
            if self.slack_bus == 0:
                P_mismatch = P_given[1:] - Parr[1:]
                Q_mismatch = Q_given[1:6] - Qarr[1:6]
            if self.slack_bus == 6:
                P_mismatch = P_given[0:6] - Parr[0:6]
                Q_mismatch = Q_given[1:6] - Qarr[1:6]

            # Set the convergence met variable to 1
            self.convergencemet = 1

            # Check Power Mismatch for Convergence
            mismatchPQ = np.concatenate((P_mismatch, Q_mismatch))
            for i in range(len(mismatchPQ)):
                # If any value in the mismatched is greater than the convergence requirement, set convergence parameter
                # to not met
                if mismatchPQ[i] > self.convergencevalue:
                    self.convergencemet = 0
                    break
            # If convergence has been met and there was no capacitor adjustment, notify user it converged
            if self.convergencemet == 1 and self.capacitor_bank_adjustment == 0:
                print("It converged after iteration ", iteration - 1)
                break

            # Calculate Jacobian Matrix
            skipterm = 0
            for i in range(self.length):
                # if slack bus skip
                if i == self.slack_bus:
                    skipterm = 1
                    continue
                for j in range(self.length):
                    # if slack bus skip
                    if j == self.slack_bus:
                        continue

                    # Diagonals have their own formula
                    if i == j:
                        for z in range(self.length):
                            J22[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                            if z == i:
                                continue

                            J11[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                        # Non-Summation Modification
                        J11[i - skipterm, j - skipterm] = -J11[i - skipterm, j - skipterm] * V[i]
                        J22[i - skipterm, j - skipterm] = J22[i - skipterm, j - skipterm] - (
                                    V[i] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))

                    # If not a diagonal of the jacobian
                    else:
                        J11[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J22[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # Delete Voltage Controlled bus 1 or 7, here if it is 7
            if self.slack_bus == 0:
                J_temp22 = np.delete(J22, 5, 1)
                J22 = J_temp22
                J_temp22 = np.delete(J22, 5, 0)
                J22 = J_temp22

            # Here VCB is 1
            if self.slack_bus == 6:
                J_temp22 = np.delete(J22, 0, 1)
                J22 = J_temp22
                J_temp22 = np.delete(J22, 0, 0)
                J22 = J_temp22

            # calculate change in voltage and phase angle
            delta_correction = np.dot(np.linalg.inv(J11), P_mismatch)
            V_correction = np.dot(np.linalg.inv(J22), Q_mismatch)

            # do not update the angle of the slack bus
            delta_correction = np.concatenate((delta_correction[:self.slack_bus], [0], delta_correction[self.slack_bus:]), axis=0)

            # do not update the voltage of the slack bus or voltage controlled bus
            V_correction = np.concatenate((V_correction[:6], [0], V_correction[6:]), axis=0)
            V_correction = np.concatenate((V_correction[:0], [0], V_correction[0:]), axis=0)

            # Update the voltage and angle and increase the iteration
            delta += delta_correction
            V += V_correction
            iteration += 1

            # Parameter to store Q from Voltage controlled bus
            self.Q_k = 0

            # Calculate Q_k
            for i in range(self.length):
                self.Q_k += abs(Grid.Ybus[self.voltage_controlled_bus, i]) * V[i] * np.sin(delta[self.voltage_controlled_bus] - delta[i] - np.angle(Grid.Ybus[self.voltage_controlled_bus,i]))
            self.Q_k *= V[self.voltage_controlled_bus]
            self.Q_k *= self.Sbase

            # If Q from VCB has exceeded the limit, and a capacitor bank is not wanted
            if self.Q_k > self.Q_k_limit and self.add_cap == 0:
                print("VAR LIMIT EXCEEDED: Generator bus will no longer be a Voltage Controlled Bus")
                self.solve_exceeded_var_power_flow(Grid)
                self.Q_limit_passed = 1
                break

            # Reset the term to 0 to see if another adjustment needs made
            self.capacitor_bank_adjustment = 0

            # If Q from VCB has exceeded the limit, and a capacitor bank is wanted
            if self.Q_k > self.Q_k_limit and self.add_cap == 1:
                # Adjustment needed
                print("LIMIT EXCEEDED, INCREASING CAPACITOR BANK")
                self.capacitor_bank_adjustment = 1
                # Add capacitor bank to highest MVAR load
                j = 0
                for i in range(self.length):
                    if self.Q_given[j] >= self.Q_given[i]:
                        j = i
                Grid.Ybus[j, j] += -1j * self.Q_given[j] * self.Sbase/ (self.Vbase ** 2)
                self.B += -1j * self.Q_given[j] * self.Sbase/ (self.Vbase ** 2)

                # Reset iteration so that it does not exceed the iteration limit while adding banks
                iteration = 0

        # Check to make sure V_complex was not set up by Var limit solver
        if self.Q_limit_passed == 0:
            # Calculate the complex voltage
            for i in range(self.length):
                self.V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])

        # Function to solve power flow
        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):

        # Calculate per unit current value
        I_values_per_unit = np.zeros((self.length, self.length), dtype=complex)
        for i in range(self.length):
            for j in range(self.length):
                I_values_per_unit[i, j] = (self.V_complex[j] - self.V_complex[i]) * (Grid.Ybus[i, j])

        # Calculate actual current values
        self.I_values = (self.Sbase / (self.Vbase * np.sqrt(3))) * I_values_per_unit
        self.I_values[0, 1] *= self.Vbase / self.Vbase_slack1 # -> Because from Bus 1 to 2, transformer
        self.I_values[6, 5] *= self.Vbase / self.Vbase_slack2 # -> Becasue from Bus 7 to 6, transformer

        # Check ampacity limit on lines
        self.check_ampacity()

        # Take voltage values out of per unit
        self.V_complex *= self.Vbase

        # Calcualte conjugate of current values
        I_conjugate = np.conjugate(self.I_values)

        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V * sqrt(3)
        # Find S first
        for i in range(self.length):
            for j in range(self.length):
                self.S_values[i, j] = self.V_complex[i] * I_conjugate[i, j] * np.sqrt(3)  # because 3 phase
                # Transformer has different base
                if i==0:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack1
                # Transformer has different base
                if i==6:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack2

        # P is the real part of S, Q is the imaginary part of S
        self.P_values = np.real(self.S_values)
        self.Q_values = np.imag(self.S_values)

        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        for i1 in range(self.length):
            for i2 in range(self.length):
                for i, value in Grid.transmissionline.items():
                    if value.bus1 == ("Bus" + str(i1+1)) and value.bus2 == ("Bus" + str(i2+1)):
                        powerloss = abs(self.I_values[i1, i2]) ** 2 * value.Rtotal * 3 # because 3 phase
                        # Store the power loss, so it can be accessed through the grid class
                        Grid.store_power_loss(i, powerloss)
                        # Update the system loss
                        self.system_loss += powerloss

        # Same for transformers
        for i1 in range(self.length):
            for i2 in range(self.length):
                for i, value in Grid.transformers.items():
                    if value.bus1 == ("Bus" + str(i1+1)) and value.bus2 == ("Bus" + str(i2+1)) or value.bus1 == ("Bus" + str(i2+1)) and value.bus2 == ("Bus" + str(i1+1)):
                        if abs(self.I_values[i1, i2]) == 0:
                            continue
                        powerloss = abs(self.I_values[i1, i2]) ** 2 * value.Rpu * self.Vbase ** 2/self.Sbase * 3 # because 3 phase
                        # Transformers have a different base
                        if value.name == "T1":
                            powerloss = powerloss / (self.Vbase ** 2) * self.Vbase_slack1 ** 2
                        if value.name == "T2":
                            powerloss = powerloss / (self.Vbase ** 2) * self.Vbase_slack2 ** 2
                        # Store the power loss in the transformer in the Grid class also
                        Grid.store_power_loss_transformer(i, powerloss)
                        # Update system power loss
                        self.system_loss += powerloss

        # POWER AND Q INJECTIONS FOR SLACK + PV
        for i in range(self.length):
            # Iterate through each Bus to sum up these values
            self.Q_inj_slack += abs(Grid.Ybus[self.slack_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.sin(np.angle(self.V_complex[self.slack_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.slack_bus, i]))
            self.P_inj_slack += abs(Grid.Ybus[self.slack_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.cos(np.angle(self.V_complex[self.slack_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.slack_bus, i]))
            self.Q_inj_PV += abs(Grid.Ybus[self.voltage_controlled_bus, i]) * abs(self.V_complex[i])/self.Vbase * np.sin(np.angle(self.V_complex[self.voltage_controlled_bus]) - np.angle(self.V_complex[i]) - np.angle(Grid.Ybus[self.voltage_controlled_bus, i]))

        # Convert to the values to their actual values instead of per unit
        self.P_inj_slack *= abs(self.V_complex[self.slack_bus])/self.Vbase * self.Sbase
        self.Q_inj_slack *= abs(self.V_complex[self.slack_bus])/self.Vbase * self.Sbase
        self.Q_inj_PV *= abs(self.V_complex[self.voltage_controlled_bus])/self.Vbase * self.Sbase

        # Print these values in "Mega" units
        print("self.P_inj_slack:", self.P_inj_slack/1000000, "MW")
        print("self.Q_inj_slack:", self.Q_inj_slack/1000000, "MVAR")
        print("self.Q_inj_PV:", self.Q_inj_PV/1000000, "MVAR")

        # Print I Values
        print("\nI_Values")
        for i in range(self.length):
            for j in range(self.length):
                # If the value is 0, do not print, if the value is negative, current is flowing the other way
                # So it will print out at a later iteration in the loop
                if abs(self.I_values[i, j]) == 0 or self.I_values[i, j] < 0:
                    continue
                # Only print Line currents, although transformer currents also matched with powerworld
                if i == 0 or j == 0 or i == 6 or j == 6:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", abs(self.I_values[i, j]))

        # Print the P values from the Sending End, meaning they are positive
        print("\nP_values Sending End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.P_values[i, j]) == 0 or self.P_values[i, j] < 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.P_values[i, j] / 1000000, "MW")

        # Print the P values from the Receiving End, meaning they are negative
        print("\nP_values Receiving End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.P_values[i, j]) == 0 or self.P_values[i, j] > 0:
                    continue
                print("Bus" + str(i + 1) + " from Bus" + str(j + 1) + " ", self.P_values[i, j] / 1000000, "MW")

        # Print the Q values from the Sending End, meaning they are positive
        print("\nQ_values Sending End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.Q_values[i, j]) == 0 or self.Q_values[i, j] < 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.Q_values[i, j] / 1000000, "MVAR")

        # Print the Q values from the Sending End, meaning they are negative
        print("\nQ_values Receiving End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.Q_values[i, j]) == 0 or self.Q_values[i, j] > 0:
                    continue
                print("Bus" + str(i + 1) + " from Bus" + str(j + 1) + " ", self.Q_values[i, j] / 1000000, "MVAR")

        # Print Power Loss Per Line
        print("\nPower Loss Per Line")
        for i, value in Grid.transmissionline.items():
            print(i, value.powerloss/1000000, "MW")

        # Print Power Loss in the Transformers
        print("\nPower Loss in Transformers")
        for i, value in Grid.transformers.items():
            print(i, value.powerloss/1000000, "MW")

        # Print the total system power loss
        print("\nTotal System Power Loss:", self.system_loss/1000000, "MW")

    # If the QVAR limit was exceeded previously, this function will be used if no capacitor bank is added.
    # It has the same functionality as before, however there is no longer a Voltage controlled bus
    def solve_exceeded_var_power_flow(self, Grid):
        # set blank array for given values
        P_given = np.zeros(self.length)
        Q_given = np.zeros(self.length)
        P_mismatch = np.zeros(self.length-1)
        Q_mismatch = np.zeros(self.length-1)
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
        Q_given[self.voltage_controlled_bus] = self.Q_k_limit / self.Sbase

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
                    Qarr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.sin(
                        delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # P does not include slack bus
            if self.slack_bus == 0:
                P_mismatch = P_given[1:7] - Parr[1:7]
                Q_mismatch = Q_given[1:7] - Qarr[1:7]
            if self.slack_bus == 6:
                P_mismatch = P_given[0:6] - Parr[0:6]
                Q_mismatch = Q_given[0:6] - Qarr[0:6]

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
            J22 = np.zeros((self.length - 1, self.length - 1))

            skipterm = 0

            for i in range(self.length):
                # if slack bus skip
                if i == self.slack_bus:
                    skipterm = 1
                    continue
                for j in range(self.length):
                    if j == self.slack_bus:
                        continue
                    if i == j:
                        for z in range(self.length):
                            J22[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))
                            if z == i:
                                continue
                            J11[i - skipterm, j - skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(
                                delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                        # Non-Summation Modification
                        J11[i - skipterm, j - skipterm] = -J11[i - skipterm, j - skipterm] * V[i]
                        J22[i - skipterm, j - skipterm] = J22[i - skipterm, j - skipterm] - (
                                V[i] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))

                    else:
                        J11[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J22[i - skipterm, j - skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.sin(
                            delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # Combine Jacobian
            J_inv11 = np.linalg.inv(J11)
            J_inv22 = np.linalg.inv(J22)

            # calculate change in voltage and phase angle
            delta_correction = np.dot(J_inv11, P_mismatch)
            V_correction = np.dot(J_inv22, Q_mismatch)

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
            print(self.V_complex[i])

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
