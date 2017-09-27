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
args = parser.parse_args()
if args.verbose:
    print('Chatty output turned on.')


class Atom:
    def __init__(self, coordinates, element):
        self.coords = coordinates  # renaming arguments inside of the class; coordinates = tuple
        self.element = element

    def __str__(self):
        return 'Atom {}, coordinates {x:.2f}, {y}, {z}'.format(self.element, x=self.x, y=self.y, z=self.z)


class Molecule:
    def __init__(self, atom_count, bond_count, atoms, bonds):
        self.atom_count = atom_count
        self.bond_count = bond_count
        self.atoms = atoms
        self.bonds = bonds

    def __str__(self):
        return 'Molecule: atom & bond count:{}, {}, \natms:{}, bnds:{}'.format(self.atom_count, self.bond_count,
                                                                               self.atoms, self.bonds)

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
            print('Maxdist between atms {a1} {a2}; {dist:.4f}'.format(a1=atom1, a2=atom2, dist=maxdist))

    def molecular_mass(self):
        element_sum = Counter()

        for atom in self.atoms:
            """
            example: elements in a molecule = ['N', 'O', 'O', 'O', 'H', 'H', 'H', 'H', 'H']
            Repetition of elements. The summarizing element_sum Counter() is created
            in order to reduce looping over keys in relative_atomic_masses. 
            """
            element_sum[atom.element] += 1

        molecular_mass = 0

        for element, count in element_sum.items():
            molecular_mass += relative_atomic_masses[element] * count  # getting value through dictionary
        fin_molecular_mass = molecular_mass * 1.66053904e-27

        if args.verbose:
            print('Weight = {:.3e} kg'.format(fin_molecular_mass))

    def maxbonds_counters(self):
        maxbonds = [0] * self.atom_count

        for i in self.bonds:
            if i.bond > maxbonds[i.first_atom-1]:
                maxbonds[i.first_atom-1] = i.bond
            if i.bond > maxbonds[i.second_atom-1]:
                maxbonds[i.second_atom-1] = i.bond

        # molecule_maxbonds == list of tuples created for each molecule [('N', 1), ('0', 2)]
        # contains information about elements and their maximum bonds
        molecule_maxbonds = list(zip([atom.element for atom in self.atoms], maxbonds))

        # COUNTERS
        molecule_statistics = Counter()
        for i in molecule_maxbonds:
            molecule_statistics[i] += 1

        for key in molecule_statistics:
            final_statistics[key] += molecule_statistics[key]

        if args.verbose:
            print(molecule_statistics)


def open_read_file(filename):
    # Coordinates = namedtuple('Coordinates', ['x', 'y', 'z'])
    Bonds = namedtuple('Bonds', ['first_atom', 'second_atom', 'bond'])
    molecules = list()

    with open(filename) as file:
        while True:
            try:
                for i in range(3):
                    file.readline()
                info = file.readline()

                atom_count = int(info[0:3])
                bond_count = int(info[3:6])

                atoms = list()  # list of Atom objects included in a molecule
                for i in range(atom_count):
                    line = file.readline()
                    element = line[31:34].strip()
                    coordinates = tuple(float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                    atom = Atom(coordinates, element)
                    atoms.append(atom)

                bonds = list()  # stores information about atoms in a molecule and their valence
                for i in range(bond_count):
                    line = file.readline()
                    first_at, second_at, bond = (int(line[l: r]) for l, r in [(0, 3), (3, 6), (6, 9)])
                    bonds.append(Bonds(first_at, second_at, bond))

                molecule = Molecule(atom_count, bond_count, atoms, bonds)
                molecules.append(molecule)

            except ValueError:
                break

            while True:
                line = file.readline()
                if '$$$$' in line:
                    break

    return molecules


def get_atomic_masses():
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            relative_atomic_masses[row['Symbol']] = float(row['Atomic Weight'])

    # periodic table storing - other possibility:
    # relative_atomic_masses as local in get_atomic_masses, return it

if __name__ == '__main__':
    final_statistics = Counter()
    relative_atomic_masses = dict()

    stack = open_read_file(args.filename)
    get_atomic_masses()

    for n, molecule in enumerate(stack, start=1):
        if args.verbose:
            print(str(n)+':')
        molecule.maxdistance()
        molecule.molecular_mass()
        molecule.maxbonds_counters()

    print("\n### FINAL STATISTICS OF ATOM TYPES: ###")
    print_final = pprint.PrettyPrinter(indent=2)
    print_final.pprint(final_statistics)
