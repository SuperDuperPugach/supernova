class Supernova:
    def __init__(self, ev):
        self.ev = ev                #neutrino energy
        self.e0  = self.ev - delt   #positron energy in null order
        self.ve0 = sqrt(self.e0*self.e0 - me*me)/self.e0  #positron velocity in null order
        self.p_probability = []
        self.p_energyBin   = []
        self.auxiliaryFunc() #init self.p_probability & self.p_energyBin
