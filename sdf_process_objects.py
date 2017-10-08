import argparse
import csv
import math
import pprint
from collections import Counter
from collections import namedtuple

parser = argparse.ArgumentParser(prog='SDF4ever',
                                 description='SDF format processing.\n',
                                 epilog='End of help block. Now try it yourself. Good luck!')
parser.add_argument('-v', '--verbose', action='store_true', help='Prints computed data in detail.')
parser.add_argument('filename')
parser.add_argument('--maxweight', type=float,
                    help='Prints names of molecules with maximum molecular weight given in the argument')
args = parser.parse_args()
if args.verbose:
    print('Chatty output turned on.')


class Atom:
    def __init__(self, coordinates, element):
        self.coords = coordinates  # renaming arguments inside of the class; coordinates = tuple
        self.element = element

    def __str__(self):
        return f'Atom {self.element}, coordinates {self.coords}'


class Molecule:
    def __init__(self, name, atoms, bonds):
        self.name = name
        self.atoms = atoms
        self.bonds = bonds

    def __str__(self):
        return f'Molecule {self.name}: atoms:{self.atoms}, \nbonds:{self.bonds}'

    def maxdistance(self):
        maxdist = 0
        atom1 = 1    # initialization required by PEP8
        atom2 = 1

        for m, atom_a in enumerate(self.atoms):
            for n, atom_b in enumerate(self.atoms[1:], start=m+1):
                distance = math.sqrt(sum((atom_a.coords[k] - atom_b.coords[k])**2 for k in range(3)))

                if distance > maxdist:
                    maxdist = distance
                    atom1 = m+1
                    atom2 = n+1-m    # m subtraction because of start=m+1

        if args.verbose:
            print(f'Maxdist between atms {atom1} {atom2}; {maxdist:.4f}')

    def molecular_mass(self):
        molecular_mass = sum(relative_atomic_masses[atom.element] for atom in self.atoms)
        #final_molecular_mass = molecular_mass * 1.66053904e-27

        if args.verbose:
            print(f'Relative molecular mass = {molecular_mass:.3f}')  # {final_molecular_mass:.3e} kg

    def maximum_bonds_statistics(self):
        maxbonds = [0] * len(self.atoms)

        for i in self.bonds:
            if i.bond > maxbonds[i.first_atom-1]:
                maxbonds[i.first_atom-1] = i.bond
            if i.bond > maxbonds[i.second_atom-1]:
                maxbonds[i.second_atom-1] = i.bond

        # molecule_maxbonds == list of tuples created for each molecule [('N', 1), ('0', 2)]
        # contains information about elements and their maximum bonds
        molecule_maxbonds = list(zip([atom.element for atom in self.atoms], maxbonds))
        molecule_statistics = Counter(i for i in molecule_maxbonds)

        if args.verbose:
            print(molecule_statistics)

        return molecule_statistics


class MoleculeSet:
    def __init__(self, filename):
        self.filename = filename

        # data loaded in __init__ -> no need to call open_read_file function
        Bonds = namedtuple('Bonds', ['first_atom', 'second_atom', 'bond'])
        molecules = list()

        with open(self.filename) as file:
            while True:
                try:
                    name = file.readline()[:10]
                    for i in range(2):
                        file.readline()
                    info = file.readline()

                    atoms = list()  # list of Atom objects included in a molecule
                    for i in range(int(info[0:3])):
                        line = file.readline()
                        element = line[31:34].strip()
                        coordinates = tuple(float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                        atom = Atom(coordinates, element)
                        atoms.append(atom)

                    bonds = list()  # stores information about atoms in a molecule and their valence
                    for i in range(int(info[3:6])):
                        line = file.readline()
                        first_at, second_at, bond = (int(line[l: r]) for l, r in [(0, 3), (3, 6), (6, 9)])
                        bonds.append(Bonds(first_at, second_at, bond))

                    molecule = Molecule(name, atoms, bonds)
                    molecules.append(molecule)
                    self.molecules = molecules  # new attribute added // another way?

                except ValueError:
                    break

                while True:
                    line = file.readline()
                    if '$$$$' in line:
                        break

    def statistics(self):
        final_statistics = Counter()

        for n, molecule in enumerate(self.molecules, start=1):
            if args.verbose:
                print(str(n) + ':')
            molecule.maxdistance()
            molecule.molecular_mass()
            final_statistics += molecule.maximum_bonds_statistics()

            # if molecule.molecular_mass <= float(args.maxweight):


        print("\n### FINAL STATISTICS OF ATOM TYPES: ###")
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(final_statistics)


def get_atomic_masses():
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            relative_atomic_masses[row['Symbol']] = float(row['Atomic Weight'])

if __name__ == '__main__':
    relative_atomic_masses = dict()
    get_atomic_masses()

    sdf_set = MoleculeSet(args.filename)
    sdf_set.statistics()
    print(args.maxweight)

# @classmethod ??
