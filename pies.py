import numpy as np

def Np(f,a,e2):
    return a / np.sqrt(1 - e2 * np.sin(f)**2)

def XYZ2flh(X,Y,Z,a,e2):
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p * (1 - e2)))
    while True:
        N = Np(f,a,e2)
        h = (p/np.cos(f)) - N
        fp = f    #f poprzednie 
        f = np.arctan(Z/(p*(1 - e2 * N/ (N +h))))
        if abs(fp - f) < (0.000001/206265):
            break
    l = np.arctan2(Y , X)
    return(f,l,h)

# Otwórz plik 'wsp_inp.txt' i wczytaj wartości X, Y, Z
with open('wsp_inp.txt', 'r') as f:
    lines = f.readlines()
    values = []
    for line in lines:
        if line.startswith('#'):
            continue
        values.append(list(map(float, line.strip().split(','))))

# Przelicz wartości za pomocą funkcji XYZ2flh()
a = 6378137
e2 = 0.00669438002290
results = []
for value in values:
    X, Y, Z = value
    f, l, h = XYZ2flh(X, Y, Z, a, e2)
    results.append((f,l,h))

# Zapisz wyniki do pliku 'wsp_out.txt'
with open('wsp_out.txt', 'w') as f:
    for result in results:
        f.write(f"Latitude (φ): {np.degrees(result[0]):.6f} degrees, ")
        f.write(f"Longitude (λ): {np.degrees(result[1]):.6f} degrees, ")
        f.write(f"Height (h): {result[2]:.2f} meters\n")
    
# Wyświetl informację o zapisie wyników
print("Wyniki zostały zapisane do pliku 'wsp_out.txt'")