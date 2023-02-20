# Power Flow Class
import numpy as np
from Grid import Grid

class PowerFlow:

    # Power flow
    def __init__(self, Grid):
        self.Parr = []
        self.Qarr = []

        # set given values
        P_given = np.zeros(len(Grid.buses))
        Q_given = np.zeros(len(Grid.buses))

        # TEMPORARY BUS DATA
        # Bus 1: Slack Bus -  𝑉=1.0 pu, 𝛿=0^∘
        # Bus 2: 𝑃_𝐿=0, 𝑄_𝐿=0
        # Bus 3:𝑃_𝐿=110 "MW", 𝑄_𝐿=50" Mvar"
        # Bus 4: 𝑃_𝐿=100" MW", 𝑄_𝐿=70" Mvar"
        # Bus 5: 𝑃_𝐿=100" MW",𝑄_𝐿=65" Mvar"
        # Bus 6: 𝑃_𝐿=0 "MW", 𝑄_𝐿=0" Mvar"
        # Bus 7:𝑃_𝐺=200 "MW", 𝑉=1.0,  𝑃_𝐿=0, 〖 𝑄〗_𝐿=0
        P_given[0] = 0
        P_given[1] = 0
        P_given[2] = 110
        P_given[3] = 100
        P_given[4] = 100
        P_given[5] = 0
        P_given[6] = 200
        Q_given[0] = 0
        Q_given[1] = 0
        Q_given[2] = 50
        Q_given[3] = 70
        Q_given[4] = 65
        Q_given[5] = 0
        Q_given[6] = 0



        # set intial guess
        V = [1.0,1.0,1.0,1.0,1.0,1.0,1.0]
        delta = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

        # estblish array for calculated value
        Parr = np.zeros(len(Grid.buses))
        Qarr = np.zeros(len(Grid.buses))

        # calculate mismatch
        for i in range(len(Grid.buses)):
            for j in range (len(Grid.buses)):
                Parr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                Qarr[i] += V[i] * V[j] * abs(Grid.Ybus[i, j]) * np.sin(delta[j] - delta[j] - np.angle(Grid.Ybus[i, j]))

        P_mismatch = P_given - Parr
        Q_mismatch = Q_given - Qarr
        P_mismatch = np.delete(P_mismatch, 0, axis = 0)
        Q_mismatch = np.delete(Q_mismatch, 0, axis = 0)

        # Calculate Jacobian Matrix
        i = 1
        j = 1
        z = 1
        J11 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J12 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J21 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))
        J22 = np.zeros((len(Grid.buses) - 1, len(Grid.buses) - 1))

        while i < len(Grid.buses):
            while j < (len(Grid.buses)):
                if i == j:
                    while z < len(Grid.buses):
                        J12[i-1, j-1] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J22[i-1, j-1] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        if z == i:
                            z += 1
                            continue

                        J11[i-1, j-1] += abs(Grid.Ybus[i, z]) * V[z] * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        J21[i-1, j-1] += abs(Grid.Ybus[i, z]) * V[z] * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                        z += 1
                    z = 1
                    J11[i-1, j-1] = J11[i-1, j-1] * -V[i]
                    J12[i-1, j-1] = J12[i-1, j-1] + (V[i] * abs(Grid.Ybus[i, j]) * np.cos(np.angle(Grid.Ybus[i, j])))
                    J21[i-1, j-1] = J21[i-1, j-1] * V[i]
                    J22[i-1, j-1] = J22[i-1, j-1] + -(V[i] * abs(Grid.Ybus[i, j]) * np.sin(np.angle(Grid.Ybus[i, j])))
                    i += 1
                    j += 1
                else:
                    J11[i-1, j-1] = V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    J12[i-1, j-1] = V[i] * abs(Grid.Ybus[i, j]) * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    J21[i-1, j-1] = -V[i] * abs(Grid.Ybus[i, j]) * V[j] * np.cos(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    J22[i-1, j-1] = V[i] * abs(Grid.Ybus[i, j]) * np.sin(delta[i] - delta[j] - np.angle(Grid.Ybus[i, j]))
                    i += 1
                    j += 1

        J = np.block([[J11, J12], [J21, J22]])

        # calculate change in voltage and phase angle
        J_inv = np.linalg.inv(J)
        mismatch = np.concatenate((P_mismatch, Q_mismatch))
        correction = np.dot(J_inv, mismatch)

        delta_correction = correction[:6]
        V_correction = correction[6:]

        delta_correction = np.insert(delta_correction, 0, 0)
        V_correction = np.insert(V_correction, 0, 0)

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



