# kalkulator współrzędnych 
# pełny skrypt projektu znajduję się w pliku "zajęcia_z_github.py"
# w pliku "pies.py" znajduje się skrypt, który pozwala wywołać nasz projekt za pomocą biblioteki arpgarse

# 1) DO CZEGO SŁUŻY PROGRAM I JAKĄ FUNKCJONALNOŚĆ OFERUJE:
    program został stworzony aby w łatwy sposób transformować współrzędne pomiędzy damymi układami  
    Danymi początkowymi są wsp. geocentryczne(X, Y, Z) 
    obsugiwane elipsoidy to : WGS84, GRS80
    współrzędne można transformować pomiędzy układami wspołrzędnych: geocentryczne, geodezyjne, topocentryczne, oraz geocentryczne w układzie PL-2000 i PL-1992 

# 2) WYMAGANIA DZIAŁANIA PROGRAMU:
    aby aplikacja działała w sposób prawidłowy należy użyć python w wersji 3.9 lub nowszej aktualizacji
    program jest dedykowany dla python w wersji 3.9, w innych wersjach mogą wystąpić niespodziewane w błędy 

# 3) WYMAGANY SYSTEM OPERACYJNY 
    program został napisany w systemie operacyjnym WINDOWS 

# 4) INSTRUKCJA KORZYSTANIA Z PROGRAMU 
    program działa na dwa sposoby: 
      - WCZYTANIE PLIKU.txt W COMMAND WINDOW
       plik.txt - jest to plik ze współrzędnymi geocentrycznymi(X,Y,Z)
       wspołrzędne należy wpisać w kolejności X,Y,Z odzielając je przecinkami
       przykład zapisu:

       3664940.500,1409153.590,5009571.170
       3664940.510,1409153.580,5009571.167
       3664940.520,1409153.570,5009571.167

       Dany plik ze współrzędnymi należy mieć w tym samym folderze co skrypt programu
       aby program wywołał, ten plik musi się nazywać "wsp_inp.txt" , w przeciwnym wypadku plik nie zostanie wczytany. 
       Plik z wynikami o nazwie wyniki.txt zostanie utworzony w folderze, którym znajduję sie skrypt programu 
       kolejność wyników: wsp. 1-geodezyjne, 2-geocentryczne PL-2000, 3-geocentryczne PL-1992, 4-topocentryczne 
       W otwartym command window należy wpisać "python zajęcia_z_github.py -m nazwa_modelu -t wsp_inp.txt (jako nazwa_modelu możemy przyjąć wgs84 lub grs80)
       - Z WYKORZYSTANIEM BIBLIOTEKI ARGPARSE (dla współrzędnych jednego punktu) 
       aby policzyć współrzędne należy otworzyć command window ( aby to zrobić należy otworzyć folder w którym znajduje się skrypt programu
       w ścieżce pliku należy kliknąc strzałkę, usunać tekts, który podświetli się na niebiesko, napisac "cmd" i nacisnąć enter
       w otwartym command window należy wpisać "python pies.py -m nazwa_modelu -x 100 -y 100 -z 100"     (jako nazwa_modelu możemy przyjąć wgs84 lub grs80, wsp są przykładowe )
       wartosći i argumenty odzielamy spacją
       kolejność wyników: wsp. 1-geodezyjne, 2-geocentryczne PL-2000, 3-geocentryczne PL-1992, 4-topocentryczny 
#    5) ZNANE BŁĘDY I NIETYPOWE ZACHOWANIA PROGRAMU
       nietypowym zachowaniem programu jest problem z wynikami topocentrycznymi(N,E,U) często wartości wyników są nie poprawne     
       Podczas poprawy programu aby można było wgrywać plik tekstowy w command window powstał konflikt przy próbie wywołania pojedyńczych współrzędnych. Dlatego          teraz gdy chce sie przeliczyć pojedyńcze współrzędne należy je wgrać do pliku textowego i postąpić zgodnie z instrukcją programu pierwszego sposobu. 
