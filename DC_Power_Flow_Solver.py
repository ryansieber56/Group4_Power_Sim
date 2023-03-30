import numpy as np
from Grid import Grid
class DCPowerFlow:

    # Power flow
    def __init__(self, Grid):
        length = len(Grid.buses)
        self.V_complex = np.zeros(length-1, dtype=complex)
        self.Bbus = np.zeros((length, length))
        self.Bbus_inv = np.zeros((length, length))

        # set blank array for given values
        P_given = np.zeros(length)

        self.voltage_pu = np.ones(length)
        self.P_loss = None
        self.system_loss = None
        self.S_values = np.zeros((length, length), dtype=complex)
        self.Q_values = np.zeros((length, length))
        self.P_values = np.zeros((length, length))
        self.I_values = np.zeros((length, length), dtype=complex)

        self.Vbase = Grid.transformers[list(Grid.transformers.keys())[0]].v2rated * 1000  # In Volts
        self.Sbase = 100000000  # VA

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
            else:
                P_given[i] = Grid.buses["Bus" + str(i + 1)].P / 100
        # Calculate Bbus and it's inverse
        self.Bbus = np.imag(Grid.Ybus)
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 0)
        self.Bbus = np.delete(self.Bbus, self.slack_bus, 1)

        # Since inverse is not working, check to make sure determinant is not zero
        print("Determinant of Bbus matrix:", np.linalg.det(self.Bbus))

        # Calculate the Bbus inverse
        self.Bbus_inv = np.linalg.inv(self.Bbus)
        temp_P_given = P_given
        print(len(temp_P_given))
        if self.slack_bus == 0:
            temp_P_given = temp_P_given[1:7]

        if self.slack_bus == 6:
            temp_P_given = temp_P_given[0:6]

        # Calculate the angle values by Bbus_inv * P_given
        self.delta = -1 * np.dot(self.Bbus_inv, temp_P_given)



        # Print angle and voltage information -> Complex power is good
        for i in range(length-1):
            print(i, "delta:", self.delta[i]*180/np.pi, "voltage:", self.voltage_pu[i], "P_given:", P_given[i])


        # Calculate the slack Power
        #if self.slack_bus == 0:
        #if self.slack_bus == 6:

        # Calculate the complex voltage
        for i in range(length-1):
            self.V_complex[i] = self.voltage_pu[i] * np.cos(self.delta[i]) + 1j * self.voltage_pu[i] * np.sin(self.delta[i])

        self.V_complex = np.insert(self.V_complex, self.slack_bus, 1)
        for i in range(length):
            print("V_complex for bus", i+1, ":", self.V_complex[i])

        self.solve_power_flow(Grid)

    def solve_power_flow(self, Grid):
        length = len(self.V_complex)
        print(self.V_complex)
        print(length)
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
            self.I_values[0, 1] *= self.Vbase / 20000  # -> Because different bases for slack bus
        if self.slack_bus == 6:
            self.I_values[5, 6] *= self.Vbase / 18000  # -> Because different bases for slack bus

        # Take voltage values out of per unit
        self.V_complex *= self.Vbase
        I_conjugate = np.conjugate(self.I_values)
        # Set up to find P and Q
        # Find P, S = P + jQ, and S = I_conjugate * V * sqrt(3)
        for i in range(length):
            for j in range(length):
                self.S_values[i,j] = self.V_complex[i] * I_conjugate[i,j] * np.sqrt(3) # because 3 phase
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
                            powerloss = powerloss / (self.Vbase **2) * 20000 * 20000

                        if(self.slack_bus == 6 and (i1 or i2 == 6)):
                            powerloss = powerloss / (self.Vbase **2) * 18000 * 18000

                        Grid.store_power_loss_transformer(i, powerloss)
                        print("Power Loss", i, powerloss)
                        self.system_loss += powerloss

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
        print("This power loss does not include Transformers, fix this")





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
