# zpracovani sdf souboru, extrahovani informaci o atomech a jejich vazbach
# na zaklade struktury sfd souboru vytvorit dvojice atom: max_vazba
# ulozit informace o vaznosti / napr. C 1 (jednovazny uhlik): 2, H 1: 12 apod.

import math
import pprint
import sys
from collections import Counter

with open('ChI_BI_project.sdf') as file:

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
                coordinates_all.append(list(float(line[l: r]) for l, r in [(3,10), (13,20), (23,30)]))
            
            maxdist = 0
            
            try:
                for i in range(len(coordinates_all)):    # creates matrix of distances of atoms
                    if len(coordinates_all) == 1:
                        print("Single atom. No distance calculated.")
                        
                    for j in range(len(coordinates_all)):                  
                        if i >= j:
                            distance = 0
                            
                        else:
                            distance = math.sqrt(sum((coordinates_all[i][k] - coordinates_all[j][k])**2 for k in range(3)))
                            
                        if distance > maxdist:
                            maxdist = distance
                            atom1 = i+1
                            atom2 = j+1
                            
                print(str(n) + ': Maxdist between atms ' + str(atom1) + ' ' + str(atom2) + '; ' + ('%.4f' % maxdist))
                
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