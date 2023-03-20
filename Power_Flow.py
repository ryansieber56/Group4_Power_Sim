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
        Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000 # In Volts
        number_of_bundles = Grid.transmissionline[list(Grid.transmissionline.keys())[0]].numberofbundles # Obtain number of bundles
        Sbase = 100000000 #VA

        # set blank array for given values
        P_given = np.zeros(len(Grid.buses))
        Q_given = np.zeros(len(Grid.buses))

        # find slack bus
        self.slack_bus = None
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i+1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i+1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i
                break

        # set given values
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i+1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i+1)].P / 100
                Q_given[i] = -Grid.buses["Bus" + str(i+1)].Q / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i+1)].P / 100
                Q_given[i] = Grid.buses["Bus" + str(i+1)].Q / 100

        # set intial guess
        # delta will be in radians because cos and sin take radians as arguments, and np.angle will return the Ybus angle in radians
        V = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        delta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # establish array for calculated value, without slack bus
        Parr = np.zeros(len(Grid.buses))
        Qarr = np.zeros(len(Grid.buses))
        iteration = 1
        while self.convergencemet == 0 and iteration < 30:
            #print("\n\nIteration = ", iteration)
            for i in range(len(Grid.buses)):
                Parr[i] = 0
                Qarr[i] = 0
            # calculate mismatch, ignoring the slack bus
            for i in range(len(Grid.buses)):
                # if slack bus, skip
                if i == self.slack_bus:
                    continue

                for j in range(len(Grid.buses)):
                    Parr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    if i == self.voltage_controlled_bus:
                        continue
                    Qarr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))

            # P does not include slack bus
            P_mismatch = P_given - Parr
            P_mismatch = P_mismatch[1:7]

            # Q not include slack bus or VCB
            Q_mismatch = Q_given - Qarr
            Q_mismatch = Q_mismatch[1:6]

            #print("P_mismatch")
            #print(P_mismatch)
            #print("Q_mismatch")
            #print(Q_mismatch)

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
                            J12[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))
                            J22[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                            if z == i:
                                continue

                            J11[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))
                            J21[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(delta[i] - delta[z] - np.angle(Grid.Ybus[i, z]))

                        # Non-Summation Modification
                        J11[i-skipterm, j-skipterm] = -J11[i-skipterm, j-skipterm] * V[i]
                        J12[i-skipterm, j-skipterm] = J12[i-skipterm, j-skipterm] + (V[i] * abs(Grid.Ybus[i, j]) * np.cos(np.angle(Grid.Ybus[i, j])))
                        J21[i-skipterm, j-skipterm] = (J21[i-skipterm, j-skipterm] * V[i])
                        J22[i-skipterm, j-skipterm] = J22[i-skipterm, j-skipterm] - (V[i] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))

                    else:
                        J11[i-skipterm, j-skipterm] = V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J12[i-skipterm, j-skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J21[i-skipterm, j-skipterm] = -V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J22[i-skipterm, j-skipterm] = V[i] * abs(Grid.Ybus[i, j]) * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))


            # Combine Jacobian
            J = np.block([[J11, J12], [J21, J22]])

            # Delete slack bus 1 or 7, here if it is 1
            J_temp = np.delete(J, 11, 1)
            J = J_temp
            J_temp = np.delete(J, 11, 0)
            J = J_temp

            #print("Jacobian Matrix", iteration-1)
            #i = 0
            #while i < len(J):
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
            #print("correction")
            #print(correction)
            delta_correction = correction[:6]
            V_correction = correction[6:]



            # do not update the angle of the slack bus
            delta_correction = np.concatenate((delta_correction[:0], [0], delta_correction[0:]), axis=0)

            # do not update the voltage of the slack bus or voltage controlled bus
            V_correction = np.concatenate((V_correction[:6], [0], V_correction[6:]), axis=0)
            V_correction = np.concatenate((V_correction[:0], [0], V_correction[0:]), axis=0)

            #print("Delta Correction Data")
            #print(delta_correction)
            #print("V Correction Data")
            #print(V_correction)

            delta += delta_correction
            V += V_correction

            #print("Delta Updated Data")
            #print(delta)
            #print("V Updated Data")
            #print(V)
            iteration += 1

        # Calculate Current Flowing
        # Set up array of zeros
        V_complex = np.zeros(len(V), dtype=complex)

        # Calculate the complex voltage
        for i in range(len(V)):
            V_complex[i] = V[i] * np.cos(delta[i]) + 1j * V[i] * np.sin(delta[i])

        # Calculate per unit current value: I = Ybus * V
        I_values_perunit = np.dot(Grid.Ybus, V_complex)

        # Calculate actual current values
        self.I_values = (Sbase/Vbase) * I_values_perunit / number_of_bundles

        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V
        P_values = np.zeros(len(self.I_values))
        Q_values = np.zeros(len(self.I_values))
        S_values = np.zeros(len(self.I_values), dtype=complex)
        I_conjugate = np.conjugate(self.I_values)

        # Take voltage values out of per unit
        V_complex *= Vbase

        # Calculate S, P, and Q values
        for i in range(len(self.I_values)):
            S_values[i] = I_conjugate[i] * V_complex[i]
            P_values[i] = np.real(S_values[i])
            Q_values[i] = np.imag(S_values[i])

        # Print out the current values
        print("\nI_values")
        for i in range(len(V)):
            print(self.I_values[i])

        print("\nP_values")
        for i in range(len(V)):
            print(P_values[i])

        print("\nQ_values")
        for i in range(len(V)):
            print(Q_values[i])

        # Calculate power loss in each line
        P_loss = np.zeros(len(Grid.transmissionline.items()))
        iteration_of_loop = 0
        # Loop over every item in the transmission line dictionary, i is the name and value is the object
        for i, value in Grid.transmissionline.items():
            P_loss[iteration_of_loop] = abs(self.I_values[iteration_of_loop]) * abs(self.I_values[iteration_of_loop]) * value.Rtotal

            # Store the power loss in that transmission line
            Grid.store_power_loss(i, P_loss[iteration_of_loop])
            iteration_of_loop += 1

        print("\nPower Loss Per Line")
        for i, value in Grid.transmissionline.items():
            print(i, value.powerloss)

        # Find Total System loss
        self.system_loss = 0

        for i in range(len(P_loss)):
            self.system_loss += P_loss[i]
        print("\nTotal System Power Loss")
        print(self.system_loss)
        # Write to Excel file

#Milestone 3: Solve Power System Power Flow
#Due Wednesday, March 17th

#Students will use the previous input power flow data to solve the Newton-Raphson algorithm. This should be done with a tolerance of maximum mismatch of 0.0001 pu.
#Students will calculate the current flowing through each line and determine if that current exceeds line ampacity (this is on datasheet, remember that the conductors are bundled!)
#Students will calculate the power flowing through the lines on both the sending a receiving end
#Students will calculate the power loss in each line, and the total power loss for the system
#Students will solve for the real and reactive power injected into the slack bus and the reactive power injected into the PV bus
