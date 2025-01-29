# Class to track sector contacts. This class will have parameters
# holding dictionaries to track the amino acid composition of
# associated sectors.
class sectorContacts:
    # Simply constructor.
    def __init__(self):
        self.S1 = {};
        self.S2 = {};
        self.S3 = {};
        self.S4 = {};
        self.S5 = {};
        self.S6 = {};
        self.S7 = {};
        self.S8 = {};
        self.S9 = {};
        self.S10 = {};
        self.S11 = {};
        self.S12 = {};
        self.sectors = {"S1": self.S1, "S2": self.S2, "S3": self.S3, "S4": self.S4, "S5": self.S5, "S6": self.S6,
                        "S7": self.S7, "S8":  self.S8, "S9": self.S9, "S10": self.S10, "S11": self.S11, "S12": self.S12}
