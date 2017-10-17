import argparse
import csv
import math
import pprint
from collections import Counter
from collections import namedtuple
from typing import List

parser = argparse.ArgumentParser(prog='SDF4ever',
                                 description='SDF format processing.\n',
                                 epilog='End of help block. Now try it yourself. Good luck!')
parser.add_argument('-v', '--verbose', action='store_true', help='Prints computed data in detail.')
parser.add_argument('filename')
parser.add_argument('--weight', type=str,
                    help='Prints names and atom statistics of molecules with molecular mass given in the argument.\nUSAGE: type e.g. "60:70 blablabla"')
parser.add_argument('--onlyatomtypes', type=str, help='Returns molecules containing given elements only.')
"""
parser.add_argument('--maxweight', type=float,
                    help='Prints names of molecules with maximum molecular weight given in the argument')
parser.add_argument('--minweight', type=float,
                    help='Prints names of molecules with minimum molecular weight given in the argument')
"""
args = parser.parse_args()
if args.verbose:
    print('Chatty output turned on.')

Bonds = namedtuple('Bonds', ['first_atom', 'second_atom', 'bond'])


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

        for m, atom_a in enumerate(self.atoms):
            for n, atom_b in enumerate(self.atoms[1:], start=m+1):
                distance = math.sqrt(sum((atom_a.coords[k] - atom_b.coords[k])**2 for k in range(3)))

                if distance > maxdist:
                    maxdist = distance

        return maxdist

    @property
    def molecular_mass(self):
        molecular_mass = sum(relative_atomic_masses[atom.element] for atom in self.atoms)
        # final_molecular_mass = molecular_mass * 1.66053904e-27

        return molecular_mass

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

        return molecule_statistics


class MoleculeSet:
    def __init__(self, molecules: List[Molecule]):
        self.molecules = molecules

    @classmethod
    def load_from_file(cls, filename):
        molecules = []

        with open(filename) as file:
            while True:
                try:
                    name = file.readline()[:10].strip()
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

                except ValueError:
                    break

                while True:
                    line = file.readline()
                    if '$$$$' in line:
                        break

        return MoleculeSet(molecules)

    def get_statistics(self):
        final_statistics = Counter()
        for n, molecule in enumerate(self.molecules, start=1):
            final_statistics += molecule.maximum_bonds_statistics()

        return final_statistics

    def print_statistics(self, final_statistics):
        if args.verbose:
            if len(self.molecules) == 500:  # sets with different number of molecules? --> def print_filtered_mols() ?
                pass
            else:
                for n, molecule in enumerate(self.molecules, start=1):
                    print(str("%3d" % n)+':', molecule.name, end=" ")
                    print(f'Relative molecular mass = {molecule.molecular_mass:.3f}', end=' ')
                    print(set([Atom.element for Atom in molecule.atoms]))
                print()

        print("### FINAL STATISTICS OF ATOM TYPES: ###")
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(final_statistics)
        print()

    def filter_molecules_by_minweight(self):
        molecules = []
        names = []

        for molecule in self.molecules:
            if molecule.molecular_mass >= args.minweight:
                names.append(molecule.name)
                molecules.append(molecule)

        if len(names) == 0:
            print('No molecule with those conditions')
        else:
            print('Molecules with required minimum weight:', str(len(names)))
            print_names = pprint.PrettyPrinter(indent=3, compact=True, width=75)
            print_names.pprint(names)
        print()

        return MoleculeSet(molecules)

    def print_filtered_molecules(self):
        # names and relative molecular masses; print_statistics() - only final_statistics of atom types ?
        # implement in order to print molecular mass of filtered molecules only
        # currently: molecules printed out in prints_statistics() under verbose
        # problem: if args.verbose: program prints mol mass of all molecules in SDF and I don't want it to do it
        pass

    def filter_molecules_by_weight(self):
        molecules = []

        try:
            # (float(i) for i in args.weight.split(':')) - not possible because of ''
            minimum, maximum = args.weight.split(':')
            if minimum == '': minimum = '0'
            if maximum == '': maximum = '100000'  ###

            for molecule in self.molecules:
                molecular_mass = molecule.molecular_mass
                if molecular_mass >= float(minimum) and molecular_mass <= float(maximum):
                    molecules.append(molecule)

        except ValueError:
            print('Wrong separator. Usage of ":" required.')

        return MoleculeSet(molecules)

    def filter_molecules_by_atom_types(self):
        requested_elements = [element.strip().upper() for element in args.onlyatomtypes.split(',')]
        molecules = []

        for molecule in self.molecules:
            molecule_elements = set(Atom.element for Atom in molecule.atoms)

            for n, element in enumerate(molecule_elements, start=1):
                if element in requested_elements:
                    if n == len(molecule_elements):
                        molecules.append(molecule)
                    else:
                        continue
                else:
                    break

        return MoleculeSet(molecules)


def get_atomic_masses():
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            relative_atomic_masses[row['Symbol']] = float(row['Atomic Weight'])

if __name__ == '__main__':
    relative_atomic_masses = dict()
    get_atomic_masses()

    sdf_set = MoleculeSet.load_from_file(args.filename)
    sdf_set.get_statistics()
    sdf_set.print_statistics(sdf_set.get_statistics())   ###

    set2 = sdf_set.filter_molecules_by_weight()
    set2.print_statistics(set2.get_statistics())

    set3 = sdf_set.filter_molecules_by_atom_types()
    set3.print_statistics(set3.get_statistics())

    set4 = sdf_set.filter_molecules_by_weight().filter_molecules_by_atom_types()
    set4.print_statistics(set4.get_statistics())


"""
DUTY TO DO

"""