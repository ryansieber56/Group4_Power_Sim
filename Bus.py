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
        if type == "Swing Bus":
            self.V = 1.0
            self.delta = 0.0
            self.P = 0.0
            self.Q = 0.0

        elif type == "Load Bus":
            self.V = 0.0
            self.delta = 0.0
            self.P = real_P
            self.Q = Q_or_V

        elif type == "Voltage Controlled Bus":
            self.V = Q_or_V
            self.delta = 0.0
            self.P = real_P
            self.Q = 0.0

        else:
            print("Type not accepted. Enter Swing Bus, Load Bus, or Voltage Controlled Bus.")