import argparse
import csv
import math
import pprint
from collections import Counter
from collections import namedtuple
from typing import List
import sys

parser = argparse.ArgumentParser(prog='SDF4ever',
                                 description='SDF format processing.\n',
                                 epilog='End of help block. Now try it yourself. Good luck!')
parser.add_argument('filename')
parser.add_argument('--weight', type=str,
                    help='Returns molecule set of molecules with molecular mass given in the argument.\nUSAGE: type e.g."60:70"')
parser.add_argument('--onlyatomtypes', type=str, help='Returns molecule set containing molecules of given elements or their subset.')
args = parser.parse_args()

Bonds = namedtuple('Bonds', ['first_atom', 'second_atom', 'bond'])


class Atom:
    def __init__(self, coordinates, element):
        self.coords = coordinates
        self.element = element

    def __str__(self):
        return f'Atom {self.element}, coordinates {self.coords}'


class Molecule:
    def __init__(self, name, atoms, bonds, moleculestring):
        self.name = name
        self.atoms = atoms
        self.bonds = bonds
        self.moleculestring = moleculestring

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
                    moleculestring = ""
                    name = file.readline()[:10].strip()
                    moleculestring += (name + '\n')

                    for i in range(2):
                        line = file.readline()
                        moleculestring += line
                    info = file.readline()
                    moleculestring += info

                    atoms = list()  # list of Atom objects included in a molecule
                    for i in range(int(info[0:3])):
                        line = file.readline()
                        moleculestring += line
                        element = line[31:34].strip()
                        coordinates = tuple(float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                        atom = Atom(coordinates, element)
                        atoms.append(atom)

                    bonds = list()  # stores information about atoms in a molecule and their valence
                    for i in range(int(info[3:6])):
                        line = file.readline()
                        moleculestring += line
                        first_at, second_at, bond = (int(line[l: r]) for l, r in [(0, 3), (3, 6), (6, 9)])
                        bonds.append(Bonds(first_at, second_at, bond))

                    moleculestring += '$$$$\n'
                    molecule = Molecule(name, atoms, bonds, moleculestring)
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

    def print_statistics(self, print_detailed=False):
        if print_detailed:
            for n, molecule in enumerate(self.molecules, start=1):
                print(str("%3d" % n) + ':', molecule.name, end=" ")
                print(f'Relative molecular mass = {molecule.molecular_mass:.3f}', end=' ')
                print(set([Atom.element for Atom in molecule.atoms]))
            print()

        print("### FINAL STATISTICS OF ATOM TYPES: ###")
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(self.get_statistics())
        print()

        return MoleculeSet(self.molecules)

    def filter_molecules_by_weight(self, minweight, maxweight):
        molecules = []

        for molecule in self.molecules:
            molecular_mass = molecule.molecular_mass
            if molecular_mass >= float(minweight) and molecular_mass <= float(maxweight):
                molecules.append(molecule)

        return MoleculeSet(molecules)

    def filter_molecules_by_atom_types(self, requested_elements):
        molecules = []

        for molecule in self.molecules:
            molecule_elements = set(Atom.element for Atom in molecule.atoms)

            if molecule_elements <= requested_elements:
                molecules.append(molecule)

        return MoleculeSet(molecules)

    def create_new_SDFfile(self, filename):
        with open(filename + '.sdf', 'w') as file:
            for Molecule in self.molecules:
                file.write(Molecule.moleculestring)


def get_atomic_masses():
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            relative_atomic_masses[row['Symbol']] = float(row['Atomic Weight'])


def check_arguments():
    minimum = None
    maximum = None
    requested_elements = None

    if args.weight:
        if ':' in args.weight:
            try:
                minimum, maximum = args.weight.split(':')
                if minimum == '':
                    minimum = 0
                elif maximum == '':
                    maximum = float('inf')
                else:
                    minimum = float(minimum)
                    maximum = float(maximum)

            except ValueError:
                print("--weight argument: wrong input. Check it and try again.\nEnd of program.")
                sys.exit()
        else:
            print('Usage of ":" required. Check it and try again.\nEnd of program.')
            sys.exit()

    if args.onlyatomtypes:
        for i in [element.strip() for element in args.onlyatomtypes.split(',')]:
            if i not in relative_atomic_masses.keys():
                print('--onlyatomtypes argument: wrong input. Check it and try again.\nEnd of program.')
                sys.exit()
        else:
            requested_elements = set([element.strip() for element in args.onlyatomtypes.split(',')])

    return minimum, maximum, requested_elements

if __name__ == '__main__':
    relative_atomic_masses = dict()
    get_atomic_masses()

    minimum, maximum, requested_elements = check_arguments()

    sdf_set = MoleculeSet.load_from_file(args.filename)
    sdf_set.get_statistics()
    sdf_set.print_statistics()

    set2 = sdf_set.filter_molecules_by_weight(minimum, maximum)
    set2.print_statistics(set2.get_statistics())

    set3 = sdf_set.filter_molecules_by_atom_types(requested_elements)
    set3.print_statistics(print_detailed=True)
    #set3.create_new_SDFfile('atom_types_'+''.join(sorted(list(requested_elements)))+'_'+args.filename[:-4])

    # set4 = sdf_set.filter_molecules_by_weight().filter_molecules_by_atom_types()
    # set4.print_statistics(print_detailed=True)


"""
DUTY TO DO
1) atribut Molecule - moleculestring - nacit radky puv souboru jako string
2) nacist do slovniku jmeno molekuly: zapamatovat cisla radku od-do dict[NSC1] = (60, 70)
vymyslet si vlastnosti molekul, co by se daly z SDFka vyvodit
filtr molekul - vola jednotlive podfukce, ktere kontroluji jednotlive vlastosti def check_weight), vraci bool
"""