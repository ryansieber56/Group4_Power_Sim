# Generator Class
class Generator:

    # Generator Contains a name, bus, and nominalpower (MVA)
    def __init__(self, name, bus1, nominalpower):
        self.name = name
        self.bus1 = bus1
        self.nominalpower = nominalpower