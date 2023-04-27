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
f, l, h = geo.XYZ2flh(args.x, args.y, args.z) 


#print(suma)
#print(get_product(args.x, args.y))
#print(module_12_14.__name__)













       
 
