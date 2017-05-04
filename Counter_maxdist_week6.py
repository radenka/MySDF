import argparse
import csv
import math
import pprint
import sys
from collections import Counter
import time

start_time = time.time()

parser = argparse.ArgumentParser(prog = 'SDF4ever', 
                                 description = 'SDF format processing', 
                                 epilog='End of help block. Now try it yourself.')
parser.add_argument('--foo', help='This thing serves for nothing.')
parser.add_argument('-v', '--verbose', action='store_true', help='Prints computed data in detail.')
args = parser.parse_args()
if args.verbose:
    print('Chatty output turned on.')
    
with open('ChI_BI_project.sdf') as file:     # sys.argv[1]

    final_stat = Counter()
    n = 1
    while True:
        try:
            for i in range(3):   # zahozeni prvnich tri radku souboru
                file.readline()
            info = file.readline()

            """ THE COUNTS LINE
                aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv

            Where::
            aaa = number of atoms
            bbb = number of bonds """

            atom_count = int(info[0:3])   # prvni pozice ctvrteho radku: pocet atomu
            bond_count = int(info[3:6])   # druha pozice: pocet vazeb

            atoms = list()
            coordinates_all = list()    # list of all atom coordinates created for each molecule in sdf file
            
            # block calculating maximal distance of atoms in molecule
            # Euclidean distance in three-dimensional space
            for i in range(atom_count):
                line = file.readline()                
                atoms.append(line[31:34].strip())  # atoms list is created (['N', 'O', ...])
                
                x, y, z = (float(line[l: r]) for l, r in [(3,10), (13,20), (23,30)])
                coordinates_all.append((x, y, z))
            
            maxdist = 0
            atom1 = 1
            atom2 = 1
            
            for i in range(len(coordinates_all)):    # creates matrix of distances of atoms
                for j in range(i+1, len(coordinates_all)):
                    distance = math.sqrt(sum( (coordinates_all[i][k] - coordinates_all[j][k])**2 for k in range(3) ))
                        
                    if distance > maxdist:
                        maxdist = distance
                        atom1 = i+1
                        atom2 = j+1
            
            if args.verbose:
                print(str(n) + ': Maxdist between atms ' + str(atom1) + ' ' + str(atom2) + '; ' + ('%.4f' % maxdist))                       
            #print(str(n) + ': Maxdist between atms ' + str(atom1) + ' ' + str(atom2) + '; ' + ('%.4f' % maxdist)) ####
                
            maxbonds = [0] * atom_count

            for i in range(bond_count):
                line = file.readline()
                
                first, second, bonds = (int(line[l: r]) for l, r in [(0,3), (3,6), (6,9)])

                # hodnota vazby se prepise pouze kdyz je aktualne ctena hodnota vazby je vyssi nez ulozena
                # je treba od pozice 'first' odecist 1, aby pozice atomu odpovidala pozici v 'maxbonds'
                if bonds > maxbonds[first-1]:
                    maxbonds[first-1] = bonds

                if bonds > maxbonds[second-1]:
                    maxbonds[second-1] = bonds

            pairs_tuples = list(zip(atoms, maxbonds))  # creates e.g. [('O', 2), ('N', 2), ('N', 1), ('H', 1)]
            mol_counter = Counter()
            
            for i in pairs_tuples:
                mol_counter[i] += 1  # mol_counter counts different atom types in molecule
            if args.verbose:
                PrintMol = pprint.PrettyPrinter(indent=2, width=50)
                PrintMol.pprint(mol_counter)
                
            for key in mol_counter:  # creates final statistics of atom types in whole sdf file
                final_stat[key] += mol_counter[key]

        except ValueError:
            break

        while True:
            line = file.readline()
            if '$$$$' in line:
                n += 1
                break

    print("\n### FINAL STATISTICS OF ATOM TYPES: ###")    
    PrintFinal = pprint.PrettyPrinter(indent=2)   # fcking awesome
    PrintFinal.pprint(final_stat)

print("--- %s seconds ---" % (time.time() - start_time))