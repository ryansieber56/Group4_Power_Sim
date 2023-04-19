# File to Perform Fault Calculations
import numpy as np
class FaultCalculation:

    # Power flow
    def __init__(self, Grid, SeqNet, faulttype: str, faultlocation: int):
        self.In_1 = None
        self.In_2 = None
        self.In_0 = None

        if faulttype == "Symmetrical Fault":
            self.In_1 = VF/SeqNet.Zbus1[faultlocation-1][faultlocation-1]
            self.In_2 = 0
            self.In_0 = 0

        elif faulttype == "Single Line to Ground Fault":
            self.In_0 = VF/
            self.In_1 = self.In_0
            self.In_2 = self.In_0
            pass
        elif faulttype == "Line to Line Fault":
            pass
        elif faulttype == "Double Line to Ground Fault":
            pass
        else:
            print("Invalid Fault Type")
            exit(-1)

