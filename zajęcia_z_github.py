#Poczatek projektu
from math import *
import numpy as np


def Np(f, self):
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        return(N) 

    
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
    
    
    

    

    def ukl2000(self, x, y, z):
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
        a = self.a
        e2 = self.e2
        f = self.XYZ2flh(x, y, z)[0]
        l = self.XYZ2flh(x, y, z)[1]
        
        
        def GaussKruger(fa, la, L0, a, e2):
            b = a*(sqrt(1-e2))
            e2 = ((a)**2-(b)**2)/((b)**2)
            t = np.tan(fa)
            eta2 = (e2)*((np.cos(fa))**2)
            dl = la-L0
            N = Np(fa, a, e2)
            sigma1 = sigma(e2, a, fa)
            Xgk = sigma1+((dl**2)/2)*N*(sin(fa))*(cos(fa))*(1+((dl**2)/12)*((cos(fa))**2)*(5-(t**2)+9*(eta2)+4*(eta2)**2)+((dl**4)/360)*((cos(fa))**4)*(61-58*(t**2)+(t**4)+270*(eta2)-330*(eta2)*(t**2)))
            Ygk = dl*N*(cos(fa))*(1+((dl**2)/6)*((cos(fa))**2)*(1-t**2+eta2)+((dl**4)/120)*((cos(fa))**4)*(5-18*(t**2)+(t**4)+14*eta2-58*(eta2)*(t**2)))
            return(Xgk, Ygk)




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


    def u92u00_2_GK(X, Y):
        """   
        Funkcja przelicza współrzędne układu 1992 lub układu 1992  
        na współrzędne Gaussa - Krugera.
        
        Parameters
        -------
        X    : [float] : współrzędna w układzie 1992/2000 [m]
        Y    : [float] : współrzędna w układzie 1992/2000 [m]
          
        Returns
        -------
        xGK  : [float] : współrzędna w układzie Gaussa - Krugera 
        yGK  : [float] : współrzędna w układzie Gaussa - Krugera 
        lam0 : [float] : południk osiowy [rad] 
        m    : [float] : elemntarna skala długości [niemianowana]
        
        """     
        if X < 1000000 and Y < 1000000:
            m92 = 0.9993
            xGK = (X + 5300000)/m92
            yGK = (Y - 500000)/m92
            l0 = math.radians(19)
            m = m92
            
        elif X > 1000000 and Y > 1000000:
            if Y > 5000000 and Y < 6000000:
                s = 5
                l0 = math.radians(15)
            elif Y > 6000000 and Y < 7000000:
                s = 6
                l0 =  math.radians(18)
            elif Y > 7000000 and Y < 8000000:
                s = 7
                l0 =  math.radians(21)
            elif Y > 8000000 and Y < 9000000:
                s = 8
                l0 =  math.radians(24)
            m00 = 0.999923
            xGK = X/m00
            yGK = (Y - (s * 1000000) - 500000)/m00
            m = m00
            
        return(xGK, yGK, l0, m)

    def GK_2_flh(xGK, yGK, m, l0, self):    
        """   
        przeliczenie wspolrzednych z ukladu Gaussa-Krugera
        na wspolrzedne geodezyjne 
        
        Parameters
        -------
        xGK  : [float] 
                wspolrzedna w ukladzie Gaussa - Krugera 
        yGK  : [float] 
                wspolrzedna w ukladzie Gaussa - Krugera 
        l0 : [float] 
                poludnik osiowy [radiany] 
        m    : [float] 
                elemntarna skala długości [niemianowana]
        a    : [float] 
                dluzsza polos elipsoidy [m]
        e2   : [float] 
                mimosrod elipsoidy [niemianowana]
          
        Returns
        -------
        f   : [float] 
                szerokosc geodezyjna [radiany]
        l  : [float] 
                dlugosc geodezyjna [radiany]
        
        """  
        A0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128))
        A4 = 15/256 * ((self.e2**2) + (3*(self.e2**3))/4)
        A6 = (35*(self.e2**3))/3072
        eps = 0.000001/3600 *math.pi/180
        b1 = xGK/(self.a * A0)
        
        while True:
            b0 = b1
            b1= (xGK/(self.a * A0 )) + ((A2/A0) * math.sin(2*b0)) - ((A4/A0) * math.sin(4*b0)) + ((A6/A0) * math.sin(6*b0))
            b = b1
            if math.fabs(b1 - b0) <= eps:
                break
            
        e_2 = self.e2/(1-self.e2)
        N = self.a/math.sqrt(1 - self.e2 * math.sin(b)**2)
        t = math.tan(b)
        n2 = e_2 * math.cos(b)**2
        
        fi = b - (t/2) * (((yGK/(N))**2) * (1 + n2) - (1/12) * ((yGK/( N))**4) * (5 + (3 * t**2) + (6*n2) - (6 * n2 * t**2) - (3 * n2**2) - (9 * t**2 * n2**2)) + (1/360) * ((yGK/(N))**6) * (61 + (90 * t**2) + (45 * t**4) + (107 * n2) - (162 * t**2 * n2) - (45 * (t**4) * n2)))
        lam = (1/math.cos(b)) * ((yGK/(N)) - ((1/6) * (yGK/(N))**3 * (1 + 2 * t**2 + n2)) + ((1/120) * (yGK/( N))**5 * (5 + (28 * t**2) + (24 * t**4) + (6 * n2) + (8 * n2 * t**2))))
        l = l0 + lam
        
        return(f, l)
    
    def u1992_u2000_2_flh(X, Y, self):
        """   
        Funkcja przelicza współrzędne układu 1992 lub układu 2000  
        na współrzędne geodezyjne.
        
        Parameters
        -------
        X   : [float] : współrzędna w układzie 1992/2000 [m]
        Y   : [float] : współrzędna w układzie 1992/2000 [m]
        a   : [float] : dłuższa półoś elipsoidy [m]
        e2  : [float] : mimośrod elipsoidy [niemianowana]
        
        Returns
        -------
        fi  : [float] : szerokość geodezyjna [rad]
        lam : [float] : długość geodezyjna [rad]
        
        """     
        xGK, yGK, lam0, m = self.u92u00_2_GK(X, Y)
        f, l = self.GK_2_flh(xGK, yGK, m, l0, self.e2, self.a)
        
        return(f, l)

    def rad2dms(rad):
        """   
        przeliczenie wartosci katow z radianow na stopnie 
        
        Parameters
        -------
        rad : [float] 
            kat w radianach [radiany]
       
        Returns
        -------
        dms : [list] 
            kat w stopniach, minutach, sekundach [d, m, s]
        
        """     
        dd = np.rad2deg(rad)
        dd = dd
        deg = int(np.trunc(dd))
        mnt = int(np.trunc((dd-deg) * 60))
        sec = ((dd-deg) * 60 - mnt) * 60
        dms = [deg, mnt, round(sec, 5)]
        
        return(dms)

    
if __name__ == "__main__":
    #tworze obiekt
    geo = Transformation(model = "wgs84")
    #wsp geocentryczne 
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    x, y = geo.ukl2000(X, Y, Z)
    #f = f1 * 180 / pi 
    #l = l1 * 180 / pi
    print(x, y)





























