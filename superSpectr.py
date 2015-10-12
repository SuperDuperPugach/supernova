import math
import numpy
#CONSTANTS
Dist             = 3.0856776*10**22 #cm
Enu_total        =6.2415096516*10**5 * 3.*10**53 /6. #MeV
Enu_average      = 15 #MeV
Delt             = 939.565378 - 938.272046 #MeV
Me               = 0.510999 #MeV
Nprotons         = 10**31 #numbers of proton's tagerts


class SuperSpectr :
    def __init__(self, eNu):
        self.eNu = eNu
        self.p_probability = []
        self.p_energyBin   = []
        self.auxiliary()

    def auxiliary(self):
        """
        init self.p_probability & self.p_energyBin
        """
        summ  = 0.
        norma = 0.
        try:
            init  = self.eNu[0]
            step  = self.eNu[1] - self.eNu[0]
        except IndexError:
            print "Nevernii format eNu"
        dN_dE = self.dNumb_dEnu(self.eNu)
        dN_dE[0] = 0.0
        for i in range(0, self.eNu.size):
            norma += dN_dE[i]
        for i in range(0, self.eNu.size):
            summ += dN_dE[i]/norma
            self.p_probability.append(summ)
            self.p_energyBin.append(init + i*step)

    def shootEnu(self, number):
        """
        :param number: number of request value eNu

        :return: array of eNu according dNumb_dEnu (MeV)
        """
        angleList = []
        val = numpy.random.ranf(number)
        deltaX, x, y = 0., 0., 0.
        for v in val :
            for i in range(0, self.eNu.size) :
             #   print "in second for"
                if(self.p_probability[i] >= v) :
                    deltaX = v - self.p_probability[i]
                    y = self.p_energyBin[i] - self.p_energyBin[i-1]
                    x = self.p_probability[i] - self.p_probability[i-1]
                    angleList.append(deltaX*y/x + self.p_energyBin[i])
                    break
            #    angleList.append(3.) # esli if ne srabotaet


        return numpy.asarray(angleList)

    def spectrMaxwBoltz(self, eNu):
        """
        :param eNu: neutrino energy in interval from eNu_min to infinity (MeV)

        :return: Maxwell-Boltzmann spectrum(supernova eNu distribution)
        """
        return 128./3. * eNu**3/Enu_average**4 * numpy.exp(-4*eNu/Enu_average)

    def dF_dEnu(self, eNu):
        """
        :param eNu: neutrino energy in interval from eNu_min to infinity (MeV)

        :return: time-integrated flux for a single neutrino flavor (MeV*cm**2)**(-1)
        """
        return 1./(4*math.pi*Dist**2) * Enu_total/Enu_average * self.spectrMaxwBoltz(eNu)

    def sigma_total(self, eNu):
        """
        :param eNu: neutrino energy in interval from eNu_min to infinity (MeV)

        :return: total cross section IBD (cm**2)
        """
        e0 = eNu - Delt
        return 0.0952 * 10.**(-42) * e0 * numpy.sqrt(e0**2 - Me**2)

    def dNumb_dEnu(self, eNu):
        """
        :param eNu: neutrino energy in interval from eNu_min to infinity (MeV)

        :return: number of events per eNu (1/MeV)
        """
        return Nprotons * self.sigma_total(eNu) * self.dF_dEnu(eNu)

