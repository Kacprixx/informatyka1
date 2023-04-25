#Poczatek projektu
from math import *
import numpy as np


    
def Np(f, self):
    N = self.a / np.sqrt(1 - self.ecc2 * np.sin(f)**2)
    return(N) 
    
def XYZ2flh(self,X,Y,Z):
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p * (1 - self.ecc2)))
    while True:
        N = Np(f,self)
        h = (p/np.cos(f)) - N
        fp = f    #f poprzednie 
        f = np.arctan(Z/(p*(1 - self.ecc2 * N/ (N +h))))
        if abs(fp - f) < (0.000001/206265):
            break
    l = np.arctan2(Y , X)
    return(f,l,h)    

def flh2XYZ(f, l, h, a, e2):
    N = Np(f, a, e2)
    X = (N + h) * np.cos(f) * np.cos(l)
    Y = (N + h) * np.cos(f) * np.sin(l)
    Z = (N * (1 - e2) + h) * np.sin(f)
    return(X, Y, Z)

def Rneu(f, l):
    R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                  [-np.sin(f) * np.sin(l), np.cos(l) , np.cos(f) * np.sin(l)],
                  [np.cos(f)             ,      0.   , np.sin(f)]])     
    return(R)

def XYZ2neu(dX, f, l):
    R = Rneu(f, l)
    return(R.T @ dX )

def GaussKruger(f,l,l0,a,e2):
    a2 = a**2
    b2 = a2 * (1 - e2)
    e_2 = (a2 - b2)/b2
    dl = l - l0
    dl2 = dl**2
    dl4 = dl**4
    t = tan(f)
    t2 = t**2
    t4 = t**4
    n2 = e_2 * (cos(f)**2)
    n4 = n2 ** 2
    N = Np(f,a,e2)
    e4 = e2**2
    e6 = e2**3
    A0 = 1 - (e2/4) - ((3*e4)/64) - ((5*e6)/256)
    A2 = (3/8) * (e2 + e4/4 + (15*e6)/128)
    A4 = (15/256) * (e4 + (3*e6)/4)
    A6 = (35*e6)/3072
    sigma = a * ((A0 * f) - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
    xgk = sigma + ((dl**2)/2) * N * sin(f) * cos(f) * (1 + ((dl**2)/12)*(cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
    ygk = dl * N * cos(f) * (1 + (dl2/6) * (cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
    return(xgk,ygk)

def GK2PL1992(xgk,ygk):
    x_92 = xgk * 0.9993 - 5300000
    y_92 = ygk * 0.9993 + 500000
    return(x_92, y_92)

def GK2PL2000(xgk,ygk,l0):
    strefa = int(l0 * 180/pi)/3
    x_00 = xgk * 0.999923
    y_00 = ygk * 0.999923 + strefa * 1000000 + 500000
    return(x_00, y_00)



o = object()

class Transformation:
    def init(self, model: str = "wgs84"):
    # Wykorzystywane Parametry elipsoid
    # a = dłuższa półos (promień rownikowy)
    # b = krótsza półos (promień południkowy)
    # flat = spłaszczenie
    # e2 = mimosród^2
    if model == "wgs84"
        self.a = 6378137
        self.b = 6356752.31424518
    elif model == "grs80":
        self.a = 6378137.0
        self.b = 6356752.31414036
    elif model == "mars":
        self.a = 3396900.0
        self.b = 3376097.80585952
    else:
        raise NotImplementedError(f"{model} model not implemented")
    self.flat = (self.a - self.b) / self.a    
    self.e = sqrt(2 * self.flat - self.flat ** 2) # mimosrod
    self.e2 = (2 * self.flat - self.flat ** 2)    #mimosrod^2