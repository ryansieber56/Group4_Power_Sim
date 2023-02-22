# Power Flow Class
import numpy as np
from Grid import Grid

class PowerFlow:

    # Power flow
    def __init__(self, Grid):
        self.Parr = []
        self.Qarr = []

        # set given values
        P_given = np.zeros(len(Grid.buses)-1)
        Q_given = np.zeros(len(Grid.buses)-1)

        # TEMPORARY BUS DATA
        # Bus 1: Slack Bus -  𝑉=1.0 pu, 𝛿=0^∘
        # Bus 2: 𝑃_𝐿=0, 𝑄_𝐿=0
        # Bus 3:𝑃_𝐿=110 "MW", 𝑄_𝐿=50" Mvar"
        # Bus 4: 𝑃_𝐿=100" MW", 𝑄_𝐿=70" Mvar"
        # Bus 5: 𝑃_𝐿=100" MW",𝑄_𝐿=65" Mvar"
        # Bus 6: 𝑃_𝐿=0 "MW", 𝑄_𝐿=0" Mvar"
        # Bus 7:𝑃_𝐺=200 "MW", 𝑉=1.0,  𝑃_𝐿=0, 〖 𝑄〗_𝐿=0

        # do not include slack bus
        P_given[0] = -0
        P_given[1] = -110
        P_given[2] = -100
        P_given[3] = -100
        P_given[4] = -0
        P_given[5] = 200

        Q_given[0] = -0
        Q_given[1] = -50
        Q_given[2] = -70
        Q_given[3] = -65
        Q_given[4] = -0
        Q_given[5] = 0



        # set intial guess, length of 6 b/c do not include slack bus
        V = [1.0,1.0,1.0,1.0,1.0,1.0]
        delta = [0.0,0.0,0.0,0.0,0.0,0.0]

        # establish array for calculated value, without slack bus
        Parr = np.zeros(len(Grid.buses)-1)
        Qarr = np.zeros(len(Grid.buses)-1)

        skipterm = 0
        # calculate mismatch, ignoring the slack bus
        for i in range(len(Grid.buses)):
            # if slack bus, skip -> using example first bus is slack bus
            if i == 0:
                i += 1
                skipterm = 1
                continue
            for j in range(len(Grid.buses)):
                # if slack bus, skip -> using example first bus is slack bus
                if j == 0:
                    j += 1
                    continue
                Parr[i-skipterm] += V[i-skipterm] * V[j-skipterm] * abs(Grid.Ybus[i, j]) * np.cos(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                Qarr[i-skipterm] += V[i-skipterm] * V[j-skipterm] * abs(Grid.Ybus[i, j]) * np.sin(delta[j-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))

        P_mismatch = P_given - Parr
        Q_mismatch = Q_given - Qarr
        print("P_mismatch")
        print(P_mismatch)
        print("Q_mismatch")
        print(Q_mismatch)
        # Calculate Jacobian Matrix
        i = 0
        j = 0
        z = 0
        J11 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J12 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J21 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J22 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))

        skipterm = 0
        while i < len(Grid.buses):
            # if slack bus skip -> using example first bus is slack bus
            if i == 0:
                skipterm = 1
                i += 1
                continue
            while j < (len(Grid.buses)):
                # if slack bus skip -> using example first bus is slack bus
                if j == 0:
                    j += 1
                    continue
                if i == j:
                    while z < len(Grid.buses):
                        # if slack bus skip -> using example first bus is slack bus
                        if z == 0:
                            z += 1
                            continue

                        J12[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z-skipterm] * np.cos(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                        J22[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z-skipterm] * np.sin(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))

                        if z == i:
                            z += 1
                            continue

                        J11[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z-skipterm] * np.sin(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                        J21[i-skipterm, j-skipterm] += abs(Grid.Ybus[i, z]) * V[z-skipterm] * np.cos(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                        z += 1
                    z = 1
                    J11[i-skipterm, j-skipterm] = J11[i-skipterm, j-skipterm] * -V[i-skipterm]
                    J12[i-skipterm, j-skipterm] = J12[i-skipterm, j-skipterm] + (V[i-skipterm] * abs(Grid.Ybus[i, j]) * np.cos(np.angle(Grid.Ybus[i, j])))
                    J21[i-skipterm, j-skipterm] = J21[i-skipterm, j-skipterm] * V[i-skipterm]
                    J22[i-skipterm, j-skipterm] = J22[i-skipterm, j-skipterm] + -(V[i-skipterm] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))
                    i += 1
                    j += 1
                else:
                    J11[i-skipterm, j-skipterm] = V[i-skipterm] * abs(Grid.Ybus[i, j]) * V[j-skipterm] * np.sin(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                    J12[i-skipterm, j-skipterm] = V[i-skipterm] * abs(Grid.Ybus[i, j]) * np.cos(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                    J21[i-skipterm, j-skipterm] = -V[i-skipterm] * abs(Grid.Ybus[i, j]) * V[j-skipterm] * np.cos(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                    J22[i-skipterm, j-skipterm] = V[i-skipterm] * abs(Grid.Ybus[i, j]) * np.sin(delta[i-skipterm] - delta[j-skipterm] - np.angle(Grid.Ybus[i, j]))
                    i += 1
                    j += 1

        # get rid of row and column for voltage controlled bus, so it is then an 11*11 (-2 cuz of slack,-1 cuz of V.C.B.)
        # here getting rid of bus 7
        J = np.block([[J11, J12], [J21, J22]])

        print("Jacobian Matrix")
        print(J)

        J_temp = np.delete(J, 5, 1)
        J = J_temp

        J_temp = np.delete(J, 5, 0)
        J = J_temp

        # calculate change in voltage and phase angle
        J_inv = np.linalg.inv(J)

        # have J-1 which is 11X11, and the P/Q mismatch which is 1X12, so delete the row of voltage controlled bus
        Q_temp = np.delete(Q_mismatch, 5, 0)
        Q_mismatch = Q_temp

        # combine mismatches
        mismatch = np.concatenate((P_mismatch, Q_mismatch))
        correction = np.dot(J_inv, mismatch)

        delta_correction = correction[:6]
        V_correction = correction[6:]
        V_correction = np.concatenate((V_correction[:5], [0], V_correction[5:]), axis=0)

        # so last change of V is 0?
        delta += delta_correction
        V += V_correction

        print("Delta Updated Data")
        print(delta)
        print("V Updated Data")
        print(V)


#Milestone 2: Produce Power Flow Input Data, Jacobian, and Injection Equations
#Students will be given loading and generator setpoint data
#Students will write up all power flow input data which should include line data, transformer data, and bus data
#Students will write a computer program which calculates the first four steps of the Newton-Raphson power flow program for a flat start for one iteration
#Calculate power mismatch
#Calculate Jacobian matrix
#Calculate voltage and phase mismatches
#Compute voltage and phase angle for next iteration

#Bus Data:
#Bus 1: Slack Bus -  𝑉=1.0 pu, 𝛿=0^∘
#Bus 2: 𝑃_𝐿=0, 𝑄_𝐿=0
#Bus 3:𝑃_𝐿=110 "MW", 𝑄_𝐿=50" Mvar"
#Bus 4: 𝑃_𝐿=100" MW", 𝑄_𝐿=70" Mvar"
#Bus 5: 𝑃_𝐿=100" MW",𝑄_𝐿=65" Mvar"
#Bus 6: 𝑃_𝐿=0 "MW", 𝑄_𝐿=0" Mvar"
#Bus 7:𝑃_𝐺=200 "MW", 𝑉=1.0,  𝑃_𝐿=0, 〖 𝑄〗_𝐿=0



