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

with open('ChI_BI_project.sdf') as file, open('Periodic_Table_Of_Elements.csv') as csvfile:     # sys.argv[1]

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
            
            for i in range(atom_count):
                line = file.readline()                
                atoms.append(line[31:34].strip())  # vytvoreni seznamu prvku ['N', 'O', ...]
                x, y, z = (float(line[l: r]) for l, r in [(3,10), (13,20), (23,30)])
                coordinates_all.append((x, y, z))   # predavam jako tuple
            
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
                        
            #print(str(n) + ': Maxdist between atms ' + str(atom1) + ' ' + str(atom2) + '; ' + ('%.4f' % maxdist))
                
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

            pairs_tuples = list(zip(atoms, maxbonds))
            mol_counter = Counter()

            for i in pairs_tuples:
                mol_counter[i] += 1
            
            for key in mol_counter:
                final_stat[key] += mol_counter[key]

            # block calculating molecular mass (real weight of molecule in kgs)
            supertrouper = csv.DictReader(csvfile, delimiter=',')
            molecular_mass = 0
            
            for atom_type, count in mol_counter.items():
                csvfile.seek(0)
                for row in supertrouper:
                    if atom_type[0] == (row['Symbol']):
                        molecular_mass += float(row['Atomic Weight']) * count

            fin_molecular_mass = molecular_mass * 1.66053904e-27  
            #print("Weight = " + str("%.3e" % fin_molecular_mass) + " kg")

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