import util

import matplotlib.pyplot as plt 
from math import exp, pi, sqrt, log
from collections import defaultdict

class Thermo:

    def __init__(self, species, completeLevels, tempRange, cutoffMethod, numberDensity=0):
        self.species = species
        self.completeLevels = completeLevels
        self.tempRange = tempRange
        self.numberDensity = numberDensity
        
        if type(cutoffMethod) in [float, int]:
            self.cutoffMethod = 'fixedLowering'
            self.cutoffEnergy = cutoffMethod
        else:
            self.cutoffMethod = cutoffMethod

        self.constants = util.Constants()

        self.epsilon = []
        self.J = []

        for config in self.completeLevels:
            for level in config.term.levels:
                if level.energy is not None:
                    self.epsilon.append(level.energy * self.constants.Cm_1ToJoules)
                    self.J.append(level.J)

    '''
    def getCpRange(self):
        CpRange = self.calcCpRange(self.tempRange)
        return CpRange
    '''

    def calcThermo(self, T):
        # function which evaluates the specific heat (J / mol K))
        # T - temperature at which to evaluate the specific heat (K)
        # cutoff - index of maximum energy level   

        self.cutoffIdx = self.calcCutoffIdx(T)
        #print(self.cutoffIdx)

        def Q(T):
            Q = 0.
            for m in range(self.cutoffIdx):
                Q += (2.*self.J[m]+1.) * exp(-self.epsilon[m] /(self.constants.kB * T))
                #print(Q)
            return Q
        
        def Qdot(T):
            Qdot = 0.
            for m in range(self.cutoffIdx):
                Qdot += (2.*self.J[m]+1.) * self.epsilon[m] \
                        * exp(-self.epsilon[m] /(self.constants.kB*T))

            return 1.0/(self.constants.kB * T**2.0) * Qdot

        def Qdot2(T):
            Qdot2 = 0.
            for m in range(self.cutoffIdx):
                Qdot2 += (2.*self.J[m]+1.) * self.epsilon[m] \
                        * (self.epsilon[m] - 2.*self.constants.kB*T) \
                        * exp(-self.epsilon[m] / (self.constants.kB*T))
                        
            return Qdot2 / (self.constants.kB**2.0 * T**4.0) 

        Q = Q(T)
        Qdot = Qdot(T)
        Qdot2 = Qdot2(T)

        Cp = self.constants.R*(T**2*Qdot2/Q - (T*Qdot/Q)**2 + 2.0*T*Qdot/Q + 5./2.)
        H = self.constants.R*(T**2.*Qdot/Q + 5./2.*T)

        T_1 = 1.0
        P_0 = 100.0E3



        sc = 5./2. + log(((2*pi*self.constants.atomicMassUnits * self.constants.kB * T_1) \
             / (self.constants.h**2.))**(3./2.) * (self.constants.kB*T_1)/P_0)
      
        #S = self.constants.R*(T/Q * Qdot + log(Q) + 3./2.*log(self.species.molarMass) + 5./2.*log(T) + sc)

        #print(self.species.mass)
        S = self.constants.R*(T/Q * Qdot + log(Q) + 3./2.*log(self.species.molarMass) + 5./2.*log(T) + sc)

        #print(S, S-sc)

        return Cp, H, S

    def calcCutoffIdx(self, T):
        
        if self.cutoffMethod == 'fixedLowering':
            # Capitelli style fixed cutoffEnergy method
            cutoffJoules = (self.species.Io - self.cutoffEnergy) * self.constants.Cm_1ToJoules
        elif self.cutoffMethod == 'CEATemper':
            # Gordon and McBride 'TEMPER' Method
            cutoffJoules = self.species.Io*self.constants.Cm_1ToJoules - T * self.constants.kB
        elif self.cutoffMethod == 'DebyeHuckel':
            # Debye-Huckel approximation of ionization lowering
            z = self.species.charge
            kB = self.constants.kB
            epsilon0 = self.constants.vacPermiativity
            qE = self.constants.electronCharge
            nE = self.numberDensity

            lambdaD = sqrt(epsilon0 * kB * T / (qE**2. * nE))
            loweringEnergy = qE**2. *  (z + 1)/(4*pi * epsilon0 * lambdaD)

            cutoffJoules = self.species.Io * self.constants.Cm_1ToJoules - loweringEnergy
        else:
            print('NO VALID IONIZATION LOWERING FOUND!')


        cutoff = len(self.epsilon)
        for e, energy in enumerate(self.epsilon):
            if energy > cutoffJoules:
                cutoff = e - 1
                break
        if cutoff < 1:
            cutoff = 1

        return cutoff


    def calcThermoPropsRange(self, temps):
        # calculates a range of Cp values (J / mol K)) for a range of temperatures
        # returns a 2D array of temperatures and values

        # temps - array of input temperatures (K)
        self.CpArray = []
        self.enthalpyArray = []
        self.entropyArray = []

        for T in temps:
            Cp, H, S = self.calcThermo(T)
            self.CpArray.append(Cp)
            self.enthalpyArray.append(H)
            self.entropyArray.append(S)
            
        return self.CpArray, self.enthalpyArray, self.entropyArray

            
if __name__ == '__main__':
    Capitelli2005oi = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                    15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                    30000, 34000, 40000, 44000, 50000],
                [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                    31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                    69.37, 51.52, 35.92, 30.6, 26.30]]

    Gordon1999 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000],
                [23.703, 22.734, 21.257, 21.040, 20.915, 20.827, 20.937, 21.799, 22.708, 23.150, 23.868, 24.722]]

    constants = util.Constants()

    #calcSpecificHeat(300)

    #oi = Species('oi', 1000)
    #oii = Species('oii', 1000)
    #ni = Species('ni', 1000)
    #nii = Species('nii', 1000)
    #ci = Species('ci', 1000)

    #oi.calcCpRange(range(1000, 51000, 1000))
    #oii.calcCpRange(range(1000, 51000, 1000))
    #ni.calcCpRange(range(1000, 51000, 1000))
    #nii.calcCpRange(range(1000, 51000, 1000))

    #he = Species('he', 1000)

    #for level in he.states:
    #    print(level)

    #plt.figure()

    #plt.plot(oi.CpArray[0], oi.CpArray[1], label='o-i')
    #plt.plot(oii.CpArray[0], oii.CpArray[1], label='o-ii')
    #plt.plot(ni.CpArray[0], ni.CpArray[1], label='n-i')
    #plt.plot(nii.CpArray[0], nii.CpArray[1], label='n-ii')

    #plt.plot(Gordon1999[0], Gordon1999[1], label='Gordon Coeffs')
    plt.plot(Capitelli2005oi[0], Capitelli2005oi[1], label='Capitelli o-i')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Cp (J/mol/K)')

    plt.legend()
    #plt.show()

    #O = self.species.
    