import numpy as np
from Grid import Grid
class DCPowerFlow:

    # Power flow
    def __init__(self, Grid):
        self.Q_inj_PV = 0
        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.Vbase_slack1 = Grid.transformers[list(Grid.transformers.keys())[0]].v1rated * 1000  # In Volts
        self.Sbase = Grid.Sbase * 1000000  # VA
        self.P_inj_slack = 0
        self.Q_inj_slack = 0
        self.length = len(Grid.buses)
        self.V_complex = np.zeros(self.length-1, dtype=complex)
        self.Bbus = np.zeros((self.length, self.length))
        self.Bbus_inv = np.zeros((self.length, self.length))

        # set blank array for given values
        P_given = np.zeros(self.length)

        self.voltage_pu = np.ones(self.length)
        self.P_loss = None
        self.system_loss = None
        self.S_values = np.zeros((self.length, self.length), dtype=complex)
        self.Q_values = np.zeros((self.length, self.length))
        self.P_values = np.zeros((self.length, self.length))
        self.I_values = np.zeros((self.length, self.length), dtype=complex)

        # find slack bus
        self.slack_bus = 0
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i + 1)].type == "Slack Bus":
                self.slack_bus = i
            if Grid.buses["Bus" + str(i + 1)].type == "Voltage Controlled Bus":
                self.voltage_controlled_bus = i

        # set given values
        for i in range(len(Grid.buses)):
            if Grid.buses["Bus" + str(i + 1)].type == "Load Bus":
                P_given[i] = -Grid.buses["Bus" + str(i + 1)].P / 100
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
        self.P_given = P_given
        # Calculate Bbus and it's inverse
        self.Bbus = np.imag(Grid.Ybus)
        self.Bbus_whole = self.Bbus
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 0)
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 1)

        # Calculate the Bbus inverse
        self.Bbus_inv = np.linalg.inv(self.Bbus)
        temp_P_given = P_given

        if self.slack_bus == 0:
            temp_P_given = temp_P_given[1:7]

        if self.slack_bus == 6:
            temp_P_given = temp_P_given[0:6]

        # Calculate the angle values by Bbus_inv * P_given
        self.delta = -1 * np.dot(self.Bbus_inv, temp_P_given)

        #self.delta = [-4.264455*np.pi/180,-5.294268*np.pi/180,-4.562926*np.pi/180,-4.697326*np.pi/180,-3.791440*np.pi/180,2.203836*np.pi/180,]

        # Print angle and voltage information -> Complex power is good
        for i in range(self.length-1):
            print("Bus:", i+2, "delta:", self.delta[i]*180/np.pi, "voltage:", self.voltage_pu[i], "P_given:", P_given[i+1])

        # Calculate the complex voltage
        for i in range(self.length-1):
            self.V_complex[i] = self.voltage_pu[i] * np.cos(self.delta[i]) + 1j * self.voltage_pu[i] * np.sin(self.delta[i])

        self.V_complex = np.insert(self.V_complex, self.slack_bus, 1)
        for i in range(self.length):
            print("V_complex for bus", i+1, ":", self.V_complex[i])

        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):

        # Calculate per unit current value
        I_values_per_unit = np.zeros((self.length, self.length), dtype=complex)
        for i in range(self.length):
            for j in range(self.length):
                if j <= i:
                    continue
                I_values_per_unit[i, j] = (self.V_complex[j] - self.V_complex[i]) * (self.Bbus_whole[i, j])

        # Calculate actual current values
        self.I_values = (self.Sbase / (self.Vbase * np.sqrt(3))) * I_values_per_unit
        self.I_values[0, 1] *= self.Vbase / self.Vbase_slack1  # -> Because from Bus1 to 2

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
                # if self.slack_bus == 0 and i == 0:
                if i == 0:
                    self.S_values[i, j] = self.S_values[i, j] / self.V_complex[i] * self.Vbase_slack1
                # if self.slack_bus == 6 and i == 6:
                #    self.S_values[i, j] = self.S_values[i, j]/self.V_complex[i] * self.Vbase_slack2

        self.P_values = self.S_values
        #self.P_values = abs(self.S_values)
        #self.P_values = np.zeros(self.length)
        #for i in range(self.length):
         #   for j in range(self.length):
          #      self.P_values[i] += self.Bbus_whole[i,j] * (self.V_complex[j] - self.V_complex[i])  # because 3 phase
                # if self.slack_bus == 0 and i == 0:
                #if i == 0:
                #    self.P_values[i, j] = self.P_values[i, j] / self.V_complex[i] * self.Vbase_slack1

        # POWER INJECTIONS FOR SLACK
        #for i in range(self.length):
        #    self.P_inj_slack += self.Bbus_whole[self.slack_bus, i] * self.V_complex[i] / self.Vbase
        #self.P_inj_slack *= abs(self.V_complex[self.slack_bus]) / self.Vbase * self.Sbase
        self.P_inj_slack = 0
        for i in range(self.length):
            self.P_inj_slack += self.P_given[i]
        self.P_inj_slack = self.P_inj_slack * 100
        print("self.P_inj_slack:", abs(self.P_inj_slack), "MW")

        # Write to Excel file/ Print
        # Print Values
        print("\nI_Values")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.I_values[i, j]) == 0:
                    continue
                print("B" + str(i + 1) + " to B" + str(j + 1) + " ", abs(self.I_values[i, j]))

        print("\nP_values")
        for i in range(self.length):
            for j in range(self.length):
                if abs(self.P_values[i, j]) == 0:
                    continue
                print("Bus" + str(i + 1) + " to Bus" + str(j + 1) + " ", self.P_values[i, j] / 1000000, "MW")

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

#Milestone 4: Var Limits and Different Solvers
#Due March 29th
#• Implement var limits in your program so that your second generator
#does not exceed 175 Mvar
#• Raise the vars on Bus 4 until the var control is used for the
#generator on Bus 7. Report the bus voltages and angles.
#• Implement a capacitor bank and place it until the generator is in
#voltage control

##• Implement a fast decoupled solver and compare the iterations
#versus those from the full Newton Raphson
#• Implement the DC power flow solver
#• Compare the voltages for all three solvers described above
