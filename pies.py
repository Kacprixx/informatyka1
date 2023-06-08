
from math import *
from argparse import ArgumentParser
import zajęcia_z_github
from zajęcia_z_github import Transformation 

parser = ArgumentParser()
parser.add_argument('-m', '--m', type = str, help = 'podaj elipsoide')
parser.add_argument('-x', '--x', type=float)
parser.add_argument('-y', '--y', type=float)
parser.add_argument('-z', '--z', type=float)
args = parser.parse_args()


geo = Transformation(model = args.m)
f1, l1, h = geo.XYZ2flh(args.x, args.y, args.z) 
N, E, U = geo.XYZ_neu(args.x, args.y, args.y) 
x_92, y_92 = geo.XY_1992(args.x, args.y, args.z)
x_20, y_20 = geo.XY_2000(args.x, args.y, args.z)



f = f1 * 180 /pi 
l = l1 * 180 /pi 
print(round(f,5),round(l,5), round(h,5))
print(round(x_20,3), round(y_20,3))
print(round(x_92,3), round(y_92,3))
print(round(N,3), round(E,3), round(U,3))











       
 
