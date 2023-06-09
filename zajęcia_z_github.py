#Poczatek projektu
from math import *
from argparse import ArgumentParser


class Transformation:
    def __init__(self, model: str = "grs80"):
        '''       
        Wykorzystywane Parametry elipsoid
        a = dłuższa półos (promień rownikowy)
        b = krótsza półos (promień południkowy)
        flat = spłaszczenie
        e2 = mimosród^2
 
        Inicjuje obiekt elipsoidy na podstawie wybranego modelu.
        Dostępne modele to: WGS84, GRS80 i Mars.
        
        Argumenty:
        model (str): Łańcuch znaków określający model elipsoidy. 
                  Domyślnie ustawione na 'wgs84'.
        '''
        if model == "wgs84":
            self.a = 6378137 
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a        # splaszczenie
        self.e = sqrt(2 * self.flat - self.flat ** 2) # mimosrod
        self.e2 = (2 * self.flat - self.flat ** 2)    # mimosrod^2  
        
         
        
    def XYZ2flh(self,X,Y,Z):
        '''
        przeliczenie wspolrzednych geocentrycznych na 
        wspolrzedne geodezyjne (Algorytm Hirvonena)

        Parameters
        ----------
        X : float
            wspolrzedna geocentryczna [m]
        Y : float
            wspolrzedna geocentryczna [m]
        Z : float
            wspolrzedna geocentryczna [m]
        a : float
            dluzsza polos elipsoidy [m]
        e2: float
            mimosrod elipsoidy [niemianowane]

        Returns
        -------
        f : float
            szerokosc geodezyjna [radiany]
        l : float
            dlugosc geodezyjna [radiany]
        h : float
            wysokosc elipsoidalna [m]

        '''
        p = sqrt(X**2 + Y**2)
        f = atan(Z / (p * (1 - self.e2)))
        while True:
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)
            h = (p/cos(f)) - N
            fp = f    #f poprzednie 
            f = atan(Z/(p*(1 - self.e2 * N/ (N +h))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = atan2(Y , X)
        return(f,l,h)
    
   
    
    def XYZ_neu(self, x, y, z):
        """   
        Funkcja przelicza wspolrzedne geodezyjne  
        na wspolrzedne topograficzne.
    
        Parameters
        -------
        x : float
            wspolrzedna geocentryczna [m]
        y : float
            wspolrzedna geocentryczna [m]
        z : float
            wspolrzedna geocentryczna [m]
        a : float 
            dluzsza polos elipsoidy [m]
        e2: float 
            mimosrod elipsoidy [niemianowana]
      
        Returns
        -------
        N : float 
            wpolrzedna topocentryczna N (north) [m]
        E : float 
            wpolrzedna topocentryczna E (east) [m]
        U : float 
            wpolrzedna topocentryczna U (up) [m]
    
        """ 
        a = self.a
        e2 = self.e2
        f = self.XYZ2flh(x, y, z)[0]
        l = self.XYZ2flh(x, y, z)[1]
        #N = self.a / sqrt(1 - self.e2 * sin(f)**2)
      
        N = -sin(f) * cos(l) * x - sin(f) * sin(l) * y + cos(f) * z
        E = -sin(l) * x + cos(l) * y
        U = cos(f) * cos(l) * x + cos(f) * sin(l) * y  + sin(f) * z
        return (N, E, U)
    
    
    

    

    def XY_2000(self, x, y, z):
        '''
        

        Parameters
        ----------
        x : float
            wspolrzedna geocentryczna [m]
        y : float
            wspolrzedna geocentryczna [m]
        z : float
            wspolrzedna geocentryczna [m]

        Returns
        -------
        x_00 : float
            wspolrzedna geocentryczna w ukladzie PL-2000 [m]
        y_00 : float
            wspolrzedna geocentryczna w ukladzie PL-2000 [m]

        '''
        a = self.a
        e2 = self.e2
        f = self.XYZ2flh(x,y,z)[0]
        l = self.XYZ2flh(x,y,z)[1]
        if abs(degrees(l) - 15) <= 1.5:
            l0_deg = 15
        elif abs(degrees(l) - 18) < 1.5:
            l0_deg = 18
        elif abs(degrees(l) - 21) <= 1.5:
            l0_deg = 21
        else:
            l0_deg = 24
        l0 = radians(l0_deg)
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
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        e4 = e2**2
        e6 = e2**3
        A0 = 1 - (e2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (e2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = a * ((A0 * f) - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        xgk = sigma + ((dl**2)/2) * N * sin(f) * cos(f) * (1 + ((dl**2)/12)*(cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * cos(f) * (1 + (dl2/6) * (cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        strefa = int(l0 * 180/pi)/3
        x_00 = xgk * 0.999923
        y_00 = ygk * 0.999923 + strefa * 1000000 + 500000
        return x_00, y_00




    
    def XY_1992(self, x, y, z):
        """   
        przeliczenie wszpolrzednych geodezyjnych 
        na wspolrzedne ukladu 1992
            
        Parameters
        -------
        x : float
            wspolrzedna geocentryczna [m]
        y : float
            wspolrzedna geocentryczna [m]
        z : float
            wspolrzedna geocentryczna [m]
        a : [float] 
            dluzsza polos elipsoidy [m]
        e2: [float] 
            mimosrod elipsoidy [niemianowana]
            
        Returns
        -------
        x92 : [float] 
            wspolrzedna w ukladzie PL-1992 [m]
        y92 : [float]  
            wspolrzedna w ukladzie PL-1992 [m]  
        
        """ 
        a = self.a
        e2 = self.e2
        f = self.XYZ2flh(x,y,z)[0]
        l = self.XYZ2flh(x,y,z)[1]
        m = 0.9993
        
        N = self.a/sqrt(1-self.e2*sin(f)**2)
        e2p = self.e2/(1-self.e2)
        t = tan(f)
        n2 = e2p * cos(f)**2
        l0 = radians(19)
        lam = l - l0
        
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35 * (self.e2**3))/3072
        
        sig = self.a * ((A0*f) - (A2*sin(2*f)) + (A4*sin(4*f)) - (A6*sin(6*f)))
        x = sig + ((lam**2)/2) * (N*sin(f)*cos(f)) * (1 + ((lam**2)/12) * ((cos(f))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((lam**4)/360) * ((cos(f))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = lam * (N * cos(f)) * (1 + ((((lam**2)/6) * (cos(f))**2) * (1-(t**2) + n2)) +  (((lam**4)/(120)) * (cos(f)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
            
        x92 = x * m - 5300000
        y92 = y * m + 500000
            
        return(x92, y92)


      
X = []
Y = []
Z = []
F = []
L = []
H = []
X_92 = []
Y_92 = []
X_00 = []
Y_00 = []
N = []
E = []
U = []


parser = ArgumentParser()
parser.add_argument('-m', '--m', type = str, help = 'podaj elipsoide')
parser.add_argument('-t', '--t', type=str, help="Podaj nazwe pliku tekstowegoz danymi z rozszerzeniem txt")

args = parser.parse_args()
            
with open(args.t, 'r') as plik:
    lines = plik.readlines()
    t = 0
    for i in lines:
        t = t + 1
        if t > 0:
            x = i.split(',')
            X.append(float(x[0]))
            Y.append(float(x[1]))
            Z.append(float(x[2]))
        
        
if __name__ == "__main__":
#    parser = ArgumentParser()
#    parser.add_argument('-t', '--t', type=str, help="Podaj nazwe pliku tekstowegoz danymi z rozszerzeniem txt")
#    args = parser.parse_args()
    #tworze obiekt
    geo = Transformation(model = "wgs84")
    #wsp geocentryczne
#    x = 100; y = 120; z = 0
#    x00,y00 = geo.XY_2000(x,y,z)
#    x92,y92 = geo.XY_1992(x, y, z)
 #   f1, l1, h = geo.XYZ2flh(x,y,z)
 #   N, E, U = geo.XYZ_neu(x, y)
 #   f = f1 * 180 / pi
 #   l = l1 * 180 / pi
 #   print(x00,y00)
 #   print('')
 #   print(x92,y92)
#    print('')
 #   print(f, l, h)
 #   print('')
 #  print(N, E, U)
  #  print('')
#    X = []
#    Y = []
#    Z = []
#    F = []
#    L = []
#    H = []
#    X_92 = []
#    Y_92 = []
#    X_00 = []
#    Y_00 = []
#    N = []
#    E = []
#    U = []
    '''
    with open(args.t, 'r') as plik:
        lines = plik.readlines()
        t = 0
        for i in lines:
            t = t + 1
            if t > 0:
                x = i.split(',')
                X.append(float(x[0]))
                Y.append(float(x[1]))
                Z.append(float(x[2]))
      '''      
    for A,B,C in zip(X,Y,Z):
        f, l, h = geo.XYZ2flh(A, B, C)
        F.append(degrees(f))
        L.append(degrees(l))
        H.append(h)
        x92, y92 = geo.XY_1992(A, B, C)
        X_92.append(x92)
        Y_92.append(y92)
        x00, y00 = geo.XY_2000(A, B, C)
        X_00.append(x00)
        Y_00.append(y00)
        n, e, u = geo.XYZ_neu(A, B, C)
        N.append(n)
        E.append(e)
        U.append(u)
    



plik=open("wyniki.txt","w")
plik.write(f'Współrzędne flh, PL_1992, PL_2000, NEU stacji permanentnej GNSS \n')
plik.write(f'Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu \n')
plik.write(f'# ************************************* \n')
plik.write(f'# fLh **********************************\n')
plik.write(f'  f[d]         l[d]         h[m] \n')
plik.write(f'# ************************************* \n')
for A,B,C in zip(F,L,H):
    A = f'{A:7.4f}'
    B = f'{B:7.4f}'
    C = f'{C:7.4f}'
    plik.write(f'{A},      {B},      {C} \n')
    
  
plik.write(f'# ************************************* \n')
plik.write(f'# PL_2000 ************************************* \n')
plik.write(f'  X[m]         Y[m] \n')
plik.write(f'# ************************************* \n')
for A,B in zip(X_00,Y_00):
    A = f'{A:7.3f}'
    B = f'{B:7.3f}'
    plik.write(f'{A},   {B} \n')
    
plik.write(f'# ************************************* \n')
plik.write(f'# PL_1992 ************************************* \n')
plik.write(f'  X[m]         Y[m] \n')
plik.write(f'# ************************************* \n')
for A,B in zip(X_92,Y_92):
    A = f'{A:7.3f}'
    B = f'{B:7.3f}'
    plik.write(f'{A},   {B} \n')

plik.write(f'# ************************************* \n')
plik.write(f'# NEU ************************************* \n')
plik.write(f'  N[m]         E[m]         U[m] \n')
plik.write(f'# ************************************* \n')

for A,B,C in zip(N,E,U):
    A = f'{A:7.3f}'
    B = f'{B:7.3f}'
    C = f'{C:7.3f}'
    plik.write(f'{A},   {B},      {C} \n')
plik.close()





























