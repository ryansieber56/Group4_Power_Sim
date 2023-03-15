# Power Flow Class
import numpy as np
from Grid import Grid

class PowerFlow:

    # Power flow
    def __init__(self, Grid):
        self.Parr = []
        self.Qarr = []
        self.convergencemet = 0
        convergencevalue = 0.001
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
            print("\n\nIteration = ", iteration)
            iteration += 1
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
            print("P_mismatch")
            print(P_mismatch)
            print("Q_mismatch")
            print(Q_mismatch)
            self.convergencemet = 1
            # Check Power Mismatch for Convergence
            mismatchPQ = np.concatenate((P_mismatch, Q_mismatch))
            for i in range(len(mismatchPQ)):
                if mismatchPQ[i] > convergencevalue:
                    self.convergencemet = 0
                    break

            if self.convergencemet == 1:
                print("It converged after iteration ", iteration - 2)
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

            print("Delta Updated Data")
            print(delta)
            print("V Updated Data")
            print(V)

