from superSpectr import *
from firstorder import *
import math
import numpy

class SupernovaNeutrino:
    def __init__(self, eNu):
        s = SuperSpectr(eNu)
        self.enNu = s.shootEnu(100000)
        self.firstOrder = FirstOrder(10)

    def shootPos(self):
        posAngle = []
        posEn    = []
        kinN       = []
        for eNu in self.enNu:
            self.firstOrder = FirstOrder(eNu)
            cos = self.firstOrder.shootPosAngle(1)
            pE  = self.firstOrder.enPos(cos)
            kN  = self.firstOrder.kinNeutron(cos)
            posAngle.append(cos)
            posEn.append(pE)
            kinN.append(kN)
        return numpy.asarray(posAngle), numpy.asarray(posEn), numpy.asarray(kinN)
