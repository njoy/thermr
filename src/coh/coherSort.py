import numpy as np
from random import random

Egrid = list(np.linspace(0,10,100))

def f(x):
    return -x**2+5.0*x+2.0


def sufficientlyClose(x1,x2,tol=1e-1):
    return ( abs(f(x1)-f(x2)) <= tol )


braggs = sorted([random()*Egrid[-1] for i in range(1)])

def findLocation(Egrid, bragg, j=0):
    for i in range(j,len(Egrid)-1):
        if Egrid[i+1] > bragg:
            return i
    return len(Egrid)


def addPoint(x1,x2,vec):
    if (sufficientlyClose(x1,x2)):
        #print(x1,'and',x2,'are sufficiently close!')
        finalE.append(x2)
    else:
        #print("need to add a midpoint!")
        xm = (x1+x2)*0.5
        #print("midpoint = ",xm)
        addPoint(x1,xm,vec)
        #vec.append(xm)
        addPoint(xm,x2,vec)



finalE = []
i_old = 0
finalE = [Egrid[0]]
for bragg in braggs:
    i = findLocation(Egrid,bragg)
    if ( i > i_old ):
        for k in range(i_old+1,i+1):
            addPoint(finalE[-1],Egrid[k],Egrid)
    addPoint(finalE[-1],bragg,finalE)
    i_old = i

finalE += Egrid[i+1:]
print( (([x in finalE for x in Egrid]==[True]*len(Egrid)) and \
        ([x in finalE for x in braggs]==[True]*len(braggs)) and \
        finalE == sorted(finalE)))



