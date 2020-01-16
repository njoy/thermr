import numpy as np
from random import random

Egrid = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
Egrid = [0.0, 1.0, 2.0, 3.0]

def f(x):
    return -x**2+5.0*x+2.0


def sufficientlyClose(x1,x2,tol=1e-3):
    return ( abs(f(x1)-f(x2)) <= tol )



braggs = sorted([random()*Egrid[-1] for i in range(1)])
braggs = [2.3,2.9,3.1,3.5,6.1]
braggs = [2.3]

def findLocation(Egrid, bragg, j=0):
    for i in range(j,len(Egrid)-1):
        if Egrid[i+1] > bragg:
            return i


finalE = []
i_old = -1

for bragg in braggs:
    i = findLocation(Egrid,bragg)
    if ( i > i_old ):
        finalE += Egrid[i_old+1:i+1]
    finalE.append(bragg)

    i_old = i


finalE += Egrid[i+1:]
print(finalE)
print( (([x in finalE for x in Egrid]==[True]*len(Egrid)) and \
        ([x in finalE for x in braggs]==[True]*len(braggs)) and \
        finalE == sorted(finalE)))



