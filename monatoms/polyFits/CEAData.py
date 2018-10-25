def getCEA():    
    'N+'
    CEA = {(200., 1000.): [
       5.237079210e+03,
       2.299958315e+00,
       2.487488821e+00,
       2.737490756e-05,
      -3.134447576e-08,
       1.850111332e-11,
      -4.447350984e-15,
       2.256284738e+05,
       5.076830786e+00,],

    (1000., 6000.): [
       2.904970374e+05,
      -8.557908610e+02,
       3.477389290e+00,
      -5.288267190e-04,
       1.352350307e-07,
      -1.389834122e-11,
       5.046166279e-16,
       2.310809984e+05,
      -1.994146545e+00,],

    (6000., 50000.): [
       1.646092148e+07,
      -1.113165218e+04,
       4.976986640e+00,
      -2.005393583e-04,
       1.022481356e-08,
      -2.691430863e-13,
       3.539931593e-18,
       3.136284696e+05,
      -1.706646380e+01,]}

    return CEA
