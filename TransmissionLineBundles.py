# TransmissionLineBundles class
# Used to find DSL and DSC given the number of bundles, distance between bundles in feet, and codeword of conductor

class TransmissionLineBundles:
    def __init__(self, numberofbundles: int, distance: float, codeword: str):

        # Change distance in feet to inches
        distance = distance * 12

        # If statement to set GMR, r, and resistance values for that specific codeword
        if codeword == "Partridge":
            self.GMR = 0.2604  # in inches
            self.r = 0.321  # in inches
            self.resistanceperft = 0.0779 / 1000  # in Ohms per ft

        # Else when the codeword is not a stored type
        else:
            print("Error: Type not accepted. Try another conductor type.")

        # If there is 1 bundle, set values
        if numberofbundles == 1:
            self.DSL = self.GMR
            self.DSC = self.r

        # If there are 2 bundles, set values
        elif numberofbundles == 2:
            self.DSL = (self.GMR * distance) ** (1 / 2)
            self.DSC = (self.r * distance) ** (1 / 2)

        # If there are 3 bundles, set values
        elif numberofbundles == 3:
            self.DSL = (self.GMR * distance * distance) ** (1 / 3)
            self.DSC = (self.r * distance * distance) ** (1 / 3)

        # If there are 4 bundles, set values
        elif numberofbundles == 4:
            self.DSL = (self.GMR * distance ** 3) ** (1 / 4)
            self.DSC = 1.0941 * (self.r * distance ** 3) ** (1 / 4)
