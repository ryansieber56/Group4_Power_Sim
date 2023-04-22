# Generator Class
class Generator:

    # Generator Contains a name, bus, nominalpower (MVA), positive sequence impedance, negative sequence impedance, and the zero sequence impedance
    def __init__(self, name, bus1, nominalpower, x1gen, x2gen, x0gen, grounding_type, grounding_value):
        self.name = name
        self.bus1 = bus1
        self.nominalpower = nominalpower
        self.x1gen = x1gen
        self.x2gen = x2gen
        self.x0gen = x0gen
        self.grounding_type = grounding_type
        self.grounding_value = grounding_value
