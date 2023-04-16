# Bus Class

class Bus:

    def __init__(self, name: str):

        # Counter to track how many exist
        counter = 0

        # Create the bus with its name
        def __init__(self, name: str):
            self.name = name
            self.index = Bus.counter

            # Increase Counter
            Bus.counter = Bus.counter + 1

    # Power Flow Setting Bus Data
    def setbusdata(self, type: str, real_P: float, Q_or_V: float):
        # If the bus is a Slack Bus, set necessary parameters
        if type == "Slack Bus":
            self.type = "Slack Bus"
            self.V = 1.0
            self.delta = 0.0
            self.P = 0.0
            self.Q = 0.0

        # If the bus is a Load Bus, set necessary parameters
        elif type == "Load Bus":
            self.type = "Load Bus"
            self.V = 0.0
            self.delta = 0.0
            self.P = real_P
            self.Q = Q_or_V

        # If the bus is a Voltage Controlled Bus, set necessary parameters
        elif type == "Voltage Controlled Bus":
            self.type = "Voltage Controlled Bus"
            self.V = Q_or_V
            self.delta = 0.0
            self.P = real_P
            self.Q = 0.0

        # If the bus is typed wrong or invalid, print Error Message
        else:
            print("Type not accepted. Enter Slack Bus, Load Bus, or Voltage Controlled Bus.")
            exit(-1)