import dynamics
import math
import poliastro.constants.general as constant

global mu
global Re
mu=constant.GM_earth.value
Re=constant.R_earth.value

def getV(r,a):
    return math.sqrt(2*mu/r-mu/a)

h1=500e3
h2=120e3
v1=getV(Re+h1,Re+h1)
v2=getV(Re+h1,Re+(h1+h2)/2)
print(v1)
print(v2)