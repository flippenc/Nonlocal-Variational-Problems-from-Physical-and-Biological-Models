import numpy as np
from scipy.optimize import curve_fit

def fit_func(x,lambdaCoeff,constant):
    return lambdaCoeff*x+constant

# lambda vs a
points = np.array([(0.1,0.01567565),(0.1,0.0159005),(0.1,0.0159005),(0.2,0.0222320),(0.2,0.0225410),(0.2,0.0225410),(0.3,0.0272525),(0.3,0.027662),(0.3,0.027662),(0.4,0.031515),(0.4,0.0320295),(0.4,0.032028),(0.5,0.0352495),(0.5,0.0358465),(0.5,0.035845),(0.6,0.0386735),(0.6,0.0392855),(0.6,0.0392855),(0.7,0.0417355),(0.7,0.0424225),(0.7,0.0424235)])
x = points[:,0]
y = points[:,1]

params = curve_fit(fit_func,x,y)
[lambdaCoeff,constant] = params[0]
print(f'Lambda vs a: a = {lambdaCoeff}lambda + {constant}')

# lambda vs b
points = np.array([(0.1,0.06580900),(0.1,0.0660780),(0.1,0.0660790),(0.2,0.0955545),(0.2,0.0960010),(0.2,0.0960010),(0.3,0.1204050),(0.3,0.120810),(0.3,0.120810),(0.4,0.143345),(0.4,0.1438750),(0.4,0.143875),(0.5,0.1656700),(0.5,0.1664500),(0.5,0.166450),(0.6,0.1884350),(0.6,0.1892450),(0.6,0.1892450),(0.7,0.2124950),(0.7,0.2133950),(0.7,0.2133950)])
x = points[:,0]
y = points[:,1]

params = curve_fit(fit_func,x,y)
[lambdaCoeff,constant] = params[0]
print(f'Lambda vs b: b = {lambdaCoeff}lambda + {constant}')