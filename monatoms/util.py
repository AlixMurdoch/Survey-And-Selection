class Constants:
    def __init__(self):
        self.kB = 1.38064852E-23                        #Boltzmann constant (NIST, 2014)
        self.electronCharge = 1.60217662E-19            #charge of electron (NIST, 2014)
        self.R = 8.3144598                              #specific gas constant (NIST, 2015)
        self.eVToCm_1 = 8065.544005                     #1 cm^-1 as measured in eV (NIST, 2015)
        self.IH = 13.59843449 * self.eVToCm_1           #Hydrogen ionization energy (NIST, 2018)
        self.Rydberg = 13.605693009 * self.eVToCm_1     #Rydberg constant (NIST, 2015)
        self.vacPermiativity = 8.854187817E-12          #vacuum permittivity (NIST, 2015)
        self.c = 2.99792458E8                           #speed of light (def)
        self.h = 6.626070040E-34                        #Planck's constant (Schlamminger et. al, 2014)
        self.atomicMassUnits = 1.660539040E-27          #atomic mass unit (NIST, 2018)
        self.Cm_1ToJoules = self.electronCharge * 1.23986E-4    #
        self.avagadro = 6.02214076E23                       #Mole (def)

        
subShellMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', \
               'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', \
               'w', 'x', 'y', 'z', 'ss', 'pp', 'dd', 'ff', \
               'gg', 'hh',' ii', 'jj', 'kk', 'll', 'mm', 'nn',\
               'oo', 'qq', 'rr', 'tt', 'uu', 'vv', 'ww', 'xx',\
               'yy', 'zz', 'sss', 'ppp', 'ddd', 'fff', \
               'ggg', 'hhh',' iii', 'jjj', 'kkk', 'lll', 'mmm', 'nnn',\
               'ooo', 'qqq', 'rrr', 'ttt', 'uuu', 'vvv', 'www', 'xxx',\
               'yyy', 'zzz']
               
termMap = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', \
           'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', \
           'W', 'X', 'Y', 'Z' 'SS', 'PP', 'DD', 'FF', 'GG',\
           'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO',\
           'QQ', 'RR', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY',\
           'ZZ', 'SSS', 'PPP', 'DDD', 'FFF', 'GGG',\
           'HHH', 'III', 'JJJ', 'KKK', 'LLL', 'MMM', 'NNN', 'OOO',\
           'QQQ', 'RRR', 'TTT', 'UUU', 'VVV', 'WWW', 'XXX', 'YYY',\
           'ZZZ']

def getSpectralDataDir():
    return 'spectraSource/'