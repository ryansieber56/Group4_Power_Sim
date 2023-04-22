import numpy as np
from Grid import Grid
class DCPowerFlow:

    # Power flow
    def __init__(self, Grid):
        # DC Power Flow Solver does not work properly
        # Initial parameters that will be used later
        self.Q_inj_PV = 0
        self.P_inj_slack = 0
        self.Q_inj_slack = 0
        self.length = len(Grid.buses)
        self.P_loss = 0
        self.system_loss = 0
        self.slack_bus = None
        self.voltage_controlled_bus = None

        # Establish Empty Arrays
        P_given = np.zeros(self.length)
        self.voltage_pu = np.ones(self.length)
        self.S_values = np.zeros((self.length, self.length), dtype=complex)
        self.Q_values = np.zeros((self.length, self.length))
        self.P_values = np.zeros((self.length, self.length))
        self.I_values = np.zeros((self.length, self.length), dtype=complex)
        self.V_complex = np.zeros(self.length-1, dtype=complex)
        self.Bbus = np.zeros((self.length, self.length))
        self.Bbus_inv = np.zeros((self.length, self.length))


        # Establish Base Parameters
        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.Vbase_slack1 = Grid.transformers[list(Grid.transformers.keys())[0]].v1rated * 1000  # In Volts
        self.Sbase = Grid.Sbase * 1000000  # VA

        # Find slack bus and voltage controlled bus by looping through Buses Dictionary in Grid Class
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i

        # Set given values by looping through Buses Dictionary in Grid Class, Loads are negative
        for i in range(self.length):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
        self.P_given = P_given

        # Calculate Bbus and delete slack row and column
        self.Bbus = np.imag(Grid.Ybus)
        self.Bbus_whole = self.Bbus
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 0)
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 1)

        # Calculate the Bbus inverse
        self.Bbus_inv = np.linalg.inv(self.Bbus)
        temp_P_given = P_given

        # Get rid of slack bus for the P_given array
        if self.slack_bus == 0:
            temp_P_given = temp_P_given[1:7]

        if self.slack_bus == 6:
            temp_P_given = temp_P_given[0:6]

        # Calculate the angle values by Bbus_inv * P_given
        self.delta = -1 * np.dot(self.Bbus_inv, temp_P_given)

        # Print angle and voltage information
        #for i in range(self.length-1):
        #    print("Bus:", i+2, "delta:", self.delta[i]*180/np.pi, "voltage:", self.voltage_pu[i], "P_given:", P_given[i+1])

        # Calculate the complex voltage
        for i in range(self.length-1):
            self.V_complex[i] = self.voltage_pu[i] * np.cos(self.delta[i]) + 1j * self.voltage_pu[i] * np.sin(self.delta[i])

        self.V_complex = np.insert(self.V_complex, self.slack_bus, 1)
        #for i in range(self.length):
        #    print("V_complex for bus", i+1, ":", self.V_complex[i])

        # Solve for Powers, Currents, and other important values and print them out
        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):

        # Calculate per unit current value
        I_values_per_unit = np.zeros((self.length, self.length), dtype=complex)

        # Going from the smaller bus to the larger bus
        for i in range(self.length):
            for j in range(self.length):
                I_values_per_unit[i, j] = -(self.V_complex[j] - self.V_complex[i]) * (self.Bbus_whole[i, j])

        # Calculate actual current values
        self.I_values = (self.Sbase / (self.Vbase * np.sqrt(3))) * I_values_per_unit
        self.I_values[0, 1] *= self.Vbase / self.Vbase_slack1  # -> Because from Bus1 to 2

        # Check ampacity limit
        self.check_ampacity()

        # Take voltage values out of per unit
        self.V_complex *= self.Vbase

        # Set up to find P
        self.P_inj_slack = 0
        for i in range(self.length):
            self.P_inj_slack += self.P_given[i]
        self.P_inj_slack = -self.P_inj_slack * 100
        print("self.P_inj_slack:", self.P_inj_slack, "MW")

        I_conjugate = np.conjugate(self.I_values)
        for i in range(self.length):
            for j in range(self.length):
                self.S_values[i, j] = self.V_complex[i] * I_conjugate[i, j] * np.sqrt(3)  # because 3 phase
                # Transformer has different base
                if i==0:
                    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack1

        # Print Values
        P = abs(self.S_values)
        P[0, 1] = self.P_inj_slack*1000000
        P[1, 0] = self.P_inj_slack*1000000

        print("\nI_Values")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.I_values[i, j]) == 0 or self.I_values[i, j] < 0:
                    continue
                # Only print Line currents
                if i == 0 or j == 0 or i == 6 or j == 6:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", abs(self.I_values[i, j]))

        # Print the P values from the Sending End, meaning they are positive
        print("\nP_values Sending End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(P[i, j]) == 0 or np.imag(self.S_values[i,j]) > 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", P[i, j] / 1000000, "MW")

        # Print the P values from the Receiving End, meaning they are negative
        print("\nP_values Receiving End")
        for i in range(self.length):
            for j in range(self.length):
                if abs(P[i, j]) == 0 or np.imag(self.S_values[i,j]) < 0:
                    continue
                print("Bus" + str(i + 1) + " from Bus" + str(j + 1) + " ", P[i, j] / 1000000, "MW")

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
