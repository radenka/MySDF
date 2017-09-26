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
    def __init__(self, atom_count, bond_count, atom_list, bond_list):
        self.atom_count = atom_count
        self.bond_count = bond_count
        self.atom_list = atom_list
        self.bond_list = bond_list

    def __str__(self):
        return 'Molecule: atom & bond count:{}, {}, \natms:{}, bnds:{}'.format(self.atom_count, self.bond_count,
                                                                               self.atom_list, self.bond_list)

    def maxdistance(self):

        maxdist = 0
        atom1 = 1    # initialization required by PEP8
        atom2 = 1
        # print(len(self.atom_list)) # output: always 24
        for i in range(len(self.atom_list)):
            # print(self.atom_list[i].coords[0]) # output: always -4.4647
            for j in range(i+1, len(self.atom_list)):  # namedtuple makes it ugly
                distance = math.sqrt(sum((self.atom_list[i].coords[k] - self.atom_list[j].coords[k])**2 for k in range(3)))

                if distance > maxdist:
                    maxdist = distance
                    atom1 = i+1
                    atom2 = j+1

        if args.verbose:
            print('Maxdist between atms {a1} {a2}; {dist:.4f}'.format(a1=atom1, a2=atom2, dist=maxdist))


def open_read_file(filename):
    # Coordinates = namedtuple('Coordinates', ['x', 'y', 'z'])
    Bonds = namedtuple('Bonds', ['first_at', 'second_at', 'valence'])
    molecules = list()

    with open(filename) as file:
        while True:
            try:
                for i in range(3):
                    file.readline()
                info = file.readline()

                Molecule.atom_count = int(info[0:3])
                Molecule.bond_count = int(info[3:6])
                # print(Molecule.atom_count) # OK

                Molecule.atom_list = list()  # list of Atom objects included in a molecule
                for i in range(Molecule.atom_count):
                    line = file.readline()
                    Atom.element = line[31:34].strip()
                    # Atom.coords = Coordinates((float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)]))
                    """x, y, z = (float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                    Atom.coords = Coordinates(x, y, z)"""
                    x, y, z = (float(line[l: r]) for l, r in [(3, 10), (13, 20), (23, 30)])
                    Atom.coords = (x, y, z)
                    # print(Atom.element)  # OK
                    # print(Atom.coords)   # OK
                    Molecule.atom_list.append(Atom)


                Molecule.bond_list = list()  # stores information about atoms in a molecule and their valence
                for i in range(Molecule.bond_count):
                    line = file.readline()
                    # originally: first, second, bond = (int(line[l: r]) for l, r in [(0, 3), (3, 6), (6, 9)])
                    first_at, second_at, valence = (int(line[l: r]) for l, r in [(0, 3), (3, 6), (6, 9)])
                    Molecule.bond_list.append(Bonds(first_at, second_at, valence))  # visually more complicated than before
                                                                                    # previous version - better?
                # print(Molecule.bond_list) # OK
                # print(Molecule.atom_list[0].element) # output: H H H
                molecules.append(Molecule)

            except ValueError:
                break

            while True:
                line = file.readline()
                if '$$$$' in line:
                    break

    return molecules

if __name__ == '__main__':
    stack = open_read_file(args.filename)
    #for molecule in stack:
        #print(molecule.atom_list[0].coords)

"""help_reader = dict()
with open('Periodic_Table_Of_Elements.csv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:
        help_reader[row['Symbol']] = float(row['Atomic Weight'])"""

"""help_reader = dict()
def get_molecular_mass_using_objects(external_file):

    with open(external_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            help_reader[row['Symbol']] = float(row['Atomic Weight'])

get_molecular_mass_using_objects('Periodic_Table_Of_Elements.csv')"""
