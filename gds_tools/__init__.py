# Module definitions
import copy
from gds_tools import functions, classes, routing, geometry

# Handy constants
a = 'A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z'.split(',')
alphabet = copy.copy(a)
for i in range(1, 20):
    for c in a:
        alphabet.append(c + str(i))

# Aliases
GDSComponent = classes.GDSComponent
instruction_parse = functions.instruction_parse
add = functions.add
save = functions.save
flatten = functions.flatten
lattice = functions.lattice
lattice_cutter = functions.lattice_cutter

# Backwards compatibility
funcs = functions
transmission = routing
