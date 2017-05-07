import argparse
import csv
import math
import pprint
from collections import Counter
import time

start_time = time.time()

parser = argparse.ArgumentParser(prog='SDF4ever',
                                 description='SDF format processing.\n',
                                 epilog='End of help block. Now try it yourself. Good luck!')
parser.add_argument('-v', '--verbose', action='store_true', help='Prints computed data in detail.')
parser.add_argument('filename')
args = parser.parse_args()
if args.verbose:
    print('Chatty output turned on.')


def open_read_file(filename):

    atom_counts = list()
    bond_counts = list()
    atoms_all = list()
    coordinates_all = list()
    bond_data = list()

    with open(filename) as file:
        while True:

            try:
                for i in range(3):
                    file.readline()                
                info = file.readline()

                atom_count = int(info[0:3])
                bond_count = int(info[3:6])

                atom_counts.append(atom_count)  # atom_count: position '012' in line
                bond_counts.append(bond_count)  # bond_count: position '345' in line

                coordinates = list()
                atoms = list()
                for i in range(atom_count):
                    line = file.readline()                
                    atoms.append(line[31:34].strip())
                    x, y, z = (float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                    coordinates.append((x, y, z))

                atoms_all.append(atoms)
                coordinates_all.append(coordinates)

                molecule_data = list()
                for i in range(bond_count):
                    line = file.readline()
                    first, second, bond = (int(line[l: r]) for l, r in [(0,3), (3,6), (6,9)])
                    molecule_data.append((first, second, bond))
                bond_data.append(molecule_data)

            except ValueError:
                break

            while True:
                line = file.readline()
                if '$$$$' in line:
                    break

    allinone = dict()
    allinone['atom_counts'] = atom_counts
    allinone['bond_counts'] = bond_counts
    allinone['atoms_all'] = atoms_all
    allinone['coordinates'] = coordinates_all
    allinone['bond_data'] = bond_data
    
    return allinone


def max_distance(coordinates):
    for n, molecule in enumerate(coordinates, start=1):
        maxdist = 0
        atom1 = 1
        atom2 = 1
        
        for i in range(len(molecule)):
            for j in range(i+1, len(molecule)):
                distance = math.sqrt(sum((molecule[i][k] - molecule[j][k])**2 for k in range(3)))
                
                if distance > maxdist:
                    maxdist = distance
                    atom1 = i+1
                    atom2 = j+1
                    
        if args.verbose:
            print('{}: Maxdist between atms {a1} {a2}; {dist:.4f}'.format(n, a1=atom1, a2=atom2, dist=maxdist))


def maxbonds_counters(bond_tuples, atom_counts, atoms_all):
    glob_pairs = list()

    for mol_bonds, count, atoms in zip(bond_tuples, atom_counts, atoms_all):
        maxbonds = [0] * count
        
        # values = tuple: (position of 1st atom, position of 2nd atom, bond type)
        for first_atom, second_atom, bond in mol_bonds:
            if bond > maxbonds[first_atom-1]:   # if bonds > maxbonds[first_atom-1]
                maxbonds[first_atom-1] = bond     # maxbonds[first-1] = bonds
            if bond > maxbonds[second_atom-1]:
                maxbonds[second_atom-1] = bond

        pairs_tuples = list(zip(atoms, maxbonds))
        glob_pairs.append(pairs_tuples)

    # COUNTERS
    final_stat = Counter()   # mol_counter = Counter([])
    # final_stat += mol_counter (to same u elem_Counter)
    for n, molecule_pairs in enumerate(glob_pairs, start=1):
        mol_counter = Counter()
        
        for i in molecule_pairs:
            mol_counter[i] += 1
        for key in mol_counter:
            final_stat[key] += mol_counter[key]        
        
        if args.verbose:
            print('{}: {}'.format(n, mol_counter))
    
    print("\n### FINAL STATISTICS OF ATOM TYPES: ###")
    print_final = pprint.PrettyPrinter(indent=2)
    print_final.pprint(final_stat)


def molecular_mass(elements_lists):
    help_reader = {}
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        
        for row in reader:  # help_reader is created in order to reduce opening of the csv file
            help_reader[row['Symbol']] = float(row['Atomic Weight'])

    for n, molecule in enumerate(elements_lists, start=1):
        """
        example: molecule = ['N', 'O', 'O', 'O', 'H', 'H', 'H', 'H', 'H']
        Repetition of elements. The summarizing element_sum Counter() is created
        in order to reduce looping over keys in help_reader. """ 
        
        element_sum = Counter()
        for element in molecule:
            element_sum[element] += 1
        
        molecular_mass = 0
            
        for element, count in element_sum.items():
            molecular_mass += help_reader[element] * count
            """for PT_element, atomic_weight in help_reader.items():
                if element == PT_element:
                    molecular_mass += atomic_weight * count"""
                    
        fin_molecular_mass = molecular_mass * 1.66053904e-27
                        
        if args.verbose:
            print('{}: Weight = {:.3e} kg'.format(n, fin_molecular_mass))
            
stack = open_read_file(args.filename)
max_distance(stack['coordinates'])
molecular_mass(stack['atoms_all'])
maxbonds_counters(stack['bond_data'], stack['atom_counts'], stack['atoms_all'])

print("\n--- %s seconds ---" % (time.time() - start_time))

# using elements_lists without Counter reduction: 6.1 s
# with Counter reduction: 4.4 s (both in Wing IDE)
