#Poczatek projektu
from math import *
import numpy as np


def Np(f, self):
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        return(N) 

   
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
        
        
    def Np(f, self):
        '''
        liczy promien krzywizny w I wertykale 

        Parameters
        ----------
        f : float
            szerokosc geodezyjna [radiany]
        

        Returns
        -------
        N : float
            promien krzywizny w I wertykale [m]

        '''
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        return(N) 
        
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
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z / (p * (1 - self.e2)))
        while True:
            N = Np(f,self)
            h = (p/np.cos(f)) - N
            fp = f    #f poprzednie 
            f = np.arctan(Z/(p*(1 - self.e2 * N/ (N +h))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y , X)
        return(f,l,h)
        
    
    
    def flh2XYZ(f, l, h, self):
        '''
        przeliczane wspolrzednych geodezyjnych 
        na wspolrzedne geocentryczne 

        Parameters
        ----------
        f : float
            szerokosc geodezyjna [radiany]
        l : float
            dlugosc geodezyjna [radiany]
        h : float
            wysokosc elipsoidalna [m]
        a : float
            dluzsza polos elipsoidy [m]
        e2: float
            mimosrod elipsoidy [niemianowane]

        Returns
        -------
        X : float
            wspolrzedna geocentryczna [m]
        Y : float
            wspolrzedna geocentryczna [m]
        Z : float
            wspolrzedna geocentryczna [m]

        '''
        N = Np(f, self)
        X = (N + h) * np.cos(f) * np.cos(l)
        Y = (N + h) * np.cos(f) * np.sin(l)
        Z = (N * (1 - self.e2) + h) * np.sin(f)
        return(X, Y, Z)
    
    
    def flh2neu(fl1, fl2, self):
        """   
        Funkcja przelicza wspolrzedne geodezyjne  
        na wspolrzedne topograficzne.
    
        Parameters
        -------
        fl1 : [list]  
             wpolrzedne geodezyjne punktu poczatkowego [radiany]
        fl2 : [list]  
             wpolrzedne geodezyjne punktu koncowego [radiany]
        a   : [float] 
             dluzsza polos elipsoidy [m]
        e2  : [float] 
             mimosrod elipsoidy [niemianowana]
      
        Returns
        -------
        N : [float] 
            wpolrzedna topocentryczna N (north) [m]
        E : [float] 
            wpolrzedna topocentryczna E (east) [m]
        U : [float] 
            wpolrzedna topocentryczna U (up) [m]
    
        """  
        X1, Y1, Z1 = self.flh2XYZ(fl1[0], fl1[1], fl1[2], self.a, self.e2)
        X2, Y2, Z2 = self.flh2XYZ(fl2[0], fl2[1], fl2[2], self.a, self.e2)
    
        dx = [X2 - X1, Y2 - Y1, Z2 - Z1]  
        R = np.array([[-np.sin(fl1[0]) * np.cos(fl1[1]), -np.sin(fl1[1]), np.cos(fl1[0]) * np.cos(fl1[1])],
                      [-np.sin(fl1[0]) * np.sin(fl1[1]), np.cos(fl1[1]) , np.cos(fl1[0]) * np.sin(fl1[1])],
                      [np.cos(fl1[0]) , 0 , np.sin(fl1[0])]])
    
        neu = R.T @ dx
        N = neu[0]
        E = neu[1]
        U = neu[2]
        return(N, E, U)
    
    
    

    

    def ukl2000(f, l, self):
        """   
        Funkcja przelicza współrzędne geodezyjne  
         na współrzędne układu 2000.
     
         Parameters
         -------
         fi  :  float] : szerokość geodezyjna [rad]
         lam : [float] : długość geodezyjna [rad]
         a   : [float] : dłuższa półoś elipsoidy [m]
         e2  : [float] : mimośrod elipsoidy [niemianowana]
     
         Returns
         -------
         x00 : [float] : współrzędna w układzie 2000 [m]
         y00 : [float] : współrzędna w układzie 2000 [m]
     
         """
        m = 0.999923
    
        N = self.a/np.sqrt(1-self.e2*np.sin(f)**2)
        e2p = self.e2/(1-self.e2)
        t = np.arctan(f)
        n2 = e2p * (np.cos(f))**2
        l = np.degrees(l)
    
        if (l.all() > 13.5 and l.all()) < 16.5:
            s = 5
            l0 = 15
        elif (l.all() > 16.5 and l.all() < 19.5):
            s = 6  
            l0 = 18
        elif (l.all() > 19.5 and l.all() < 22.5):
            s = 7
            l0 = 21
        elif (l.all() > 22.5 and l.all() < 25.5):
            s = 8
            l0 = 24
    
        l = np.radians(l)
        l0 = np.radians(l0)
        lam = l - l0
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35 * (self.e2**3))/3072
    
        sig = self.a * ((A0*f) - (A2*np.sin(2*f)) + (A4*np.sin(4*f)) - (A6*np.sin(6*f)))
        x = sig + ((lam**2)/2) * (N*np.sin(f)*np.cos(f)) * (1 + ((lam**2)/12) * ((np.cos(f))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((lam**4)/360) * ((np.cos(f))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 * (t**2))))
        y = l * (N*np.cos(f)) * (1 + ((((l**2)/6) * (np.cos(f))**2) * (1-(t**2) + n2)) + (((l**4)/(120)) * (np.cos(f)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 = m * x
        y00 = m * y + (s*1000000) + 500000
        return(x00, y00)





    def ukl1992(f, l, self):
            """   
            przeliczenie wszpolrzednych geodezyjnych 
            na wspolrzedne ukladu 1992
            
            Parameters
            -------
            f  : [float] 
                  szerokosc geodezyjna [rad]
            l : [float] 
                  długosc geodezyjna [rad]
            a : [float] 
                  dluzsza polos elipsoidy [m]
            e2: [float] 
                  mimosrod elipsoidy [niemianowana]
            
            Returns
            -------
            x92 : [float] 
                  wspolrzedna w ukladzie 1992 [m]
            y92 : [float]  
                  wspolrzedna w ukladzie 1992 [m]  
            
            """ 
            m = 0.9993
            
            N = self.a/math.sqrt(1-self.e2*math.sin(fi)**2)
            e2p = self.e2/(1-self.e2)
            t = math.tan(fi)
            n2 = e2p * math.cos(fi)**2
            l0 = math.radians(19)
            lam = l - l0
            
            A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
            A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
            A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
            A6 = (35 * (self.e2**3))/3072
            
            sig = self.a * ((A0*fi) - (A2*math.sin(2*fi)) + (A4*math.sin(4*fi)) - (A6*math.sin(6*fi)))
            x = sig + ((lam**2)/2) * (N*math.sin(fi)*math.cos(fi)) * (1 + ((lam**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((lam**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
            y = lam * (N * math.cos(fi)) * (1 + ((((lam**2)/6) * (math.cos(fi))**2) * (1-(t**2) + n2)) +  (((lam**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
            
            x92 = x * m - 5300000
            y92 = y * m + 500000
            
            return(x92, y92)






    
if __name__ == "__main__":
    #tworze obiekt
    geo = Transformation(model = "wgs84")
    #wsp geocentryczne 
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    f1, l1 , h = geo.XYZ2flh(X, Y, Z)
    f = f1 * 180 / pi 
    l = l1 * 180 / pi
    print(f, l, h)


