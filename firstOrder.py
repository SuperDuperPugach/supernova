# -*- coding: utf-8 -*-
from numpy import linspace, random, asarray
from math import sqrt, pi

#CONSTANTS:
f         = 1.
g         = 1.26
f2        = 3.706
delt      = 939.565378 - 938.272046
M         = (939.565378 + 938.272046)/2.
me        = 0.510999
mn        = 939.565378
mp        = 938.272046
cosUCab   = 0.97428
deltInRad = 0.024
Gfermi    = 1.166 * 10**(-11)
y         = (delt**2-me**2)/2.
sigma0    = ((Gfermi**2)*cosUCab**2)*(1+deltInRad)/pi

class FirstOrder :
    'Generate energy and angle in first order'
    cosTeta = linspace(-1,1,201)
    def __init__(self, ev):
        self.ev = ev                #neutrino energy
        self.e0  = self.ev - delt   #positron energy in null order
        self.ve0 = sqrt(self.e0*self.e0 - me*me)/self.e0  #positron velocity in null order
        self.p_probability = []
        self.p_energyBin   = []
        self.auxiliaryFunc() #init self.p_probability & self.p_energyBin

    def auxiliaryFunc(self):
        summ  = 0.
        norma = 0.
        init  = -1.
        step  = 0.01
        dS_dc = self.dSigma_dcos(self.bigGamma, self.enPos, FirstOrder.cosTeta) # значения сечения в угол для cosTeta = linspace(-1,1,201)
        dS_dc[0] = 0.0
        for i in range(0,FirstOrder.cosTeta.size):
            norma += dS_dc[i]

        for i in range(0,FirstOrder.cosTeta.size):
            summ += dS_dc[i]/norma
            self.p_probability.append(summ)
            self.p_energyBin.append(init + i*step)

    def dSigma_dcos(self, G, Ee1, cosTeta) :
        '''
        :param cosTeta:asdasdasd
        :param Ee1: sdfdsfsdfsdf
        '''
        temp = (sigma0/2.)*((f**2+3*g*g)+(f**2-g**2)*((((Ee1(cosTeta))**2 - me**2)**(1./2.)) / \
                                                      (Ee1(cosTeta)))*cosTeta)*Ee1(cosTeta)*((Ee1(cosTeta))**2 - me**2)**(1./2.) - \
                                                      (sigma0/2.)*(G(cosTeta)/M)*self.e0*(self.e0**2-me**2)**(1./2.)
        return temp

    def bigGamma(self, cosTeta) :
        #G = 2*(f+f2)*g*((2*self.e0 + delt)*(1-cosTeta))
        #+ (f**2 + g**2)*(delt*(1+cosTeta))
        #+ (f**2 + 3*g*g)*(3*(self.e0 + delt)*(1 - cosTeta) - delt)
        #+ (f**2 - g**2)*(3*(self.e0 + delt)*(1-cosTeta) - delt)*cosTeta

        G2 =  2*(f+f2)*g*((2*self.e0 + delt)*(1-self.ve0*cosTeta) - (me**2)/self.e0)+ \
             (f**2 + g**2)*(delt*(1+self.ve0*cosTeta) + (me**2)/self.e0)+ \
             (f**2 + 3*g*g)*((self.e0 + delt)*(1 - (1/self.ve0)*cosTeta) - delt) + \
             (f**2 - g**2)*((self.e0 + delt)*(1 - (1/self.ve0)*cosTeta) - delt)*self.ve0*cosTeta

        return G2

    def enPos(self, cosTeta) :
        temp = self.e0*(1-(self.ev/M)*(1-self.ve0*cosTeta))-y**2/M
        return temp
        

    def shootPosAngle(self, n):
        angleList = []
        val = random.ranf(n)
        deltaX, x, y = 0., 0., 0.
        j = 1
        for v in val : # цикл по множеству необходимых значений косинуса
            #print "in first for v = " , v
            
            for i in range(0, FirstOrder.cosTeta.size) :
             #   print "in second for"
                if(self.p_probability[i] >= v) :

                    deltaX = v - self.p_probability[i]
                    y = self.p_energyBin[i] - self.p_energyBin[i-1]
                    x = self.p_probability[i] - self.p_probability[i-1]
                    angleList.append(deltaX*y/x + self.p_energyBin[i])
                    if(deltaX*y/x + self.p_energyBin[i] < -0.999995):
                        print 'Yes, you got it!'
                    if(deltaX*y/x + self.p_energyBin[i] > 0.999995):
                        print 'Ups, you got it!'
                    break
            #    angleList.append(3.) # vot tuta xerznaet
    

        return asarray(angleList)

    def kinNeutron(self, cosTeta):
        temp = (self.ev*self.e0/M)*(1-self.ve0*cosTeta) + y**2/M
        return temp

    def neutronAngle(self, cosTeta, e_positr, t_neutr):
        #temp = (self.ev -((e_positr**2-me**2)**(1./2.))*cosTeta)/ \
               #(((t_neutr+mn)**2-mn**2)**(1./2.))
        temp = (mp**2 + me**2 - 2*mp*e_positr - mn**2 +2*self.ev*(t_neutr+mn))/ \
               (2*((t_neutr+mn)**2-mn**2)**(1./2.)*self.ev)
        return temp

    def neutronAngle2(self, cosTeta, e_positr, t_neutr):
        temp = (self.ev -((e_positr**2-me**2)**(1./2.))*cosTeta)/ \
               (((t_neutr+mn)**2-mn**2)**(1./2.))
        return temp

    def getAll(self, cosTeta):
        eP = self.enPos(cosTeta)
        kN = self.kinNeutron(cosTeta)
        nAngle = self.neutronAngle(cosTeta, eP,kN)
        nAngle2 = self.neutronAngle2(cosTeta, eP,kN)
        return eP, kN, nAngle, nAngle2
