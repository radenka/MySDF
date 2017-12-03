import argparse
import csv
import math
import pprint
from collections import Counter
from collections import namedtuple
from collections import defaultdict
from typing import List
import sys
import statistics

parser = argparse.ArgumentParser(prog='SDF4ever',
                                 description='SDF format processing.\n'
                                             'EXIT STATUS:'
                                             '\n\t1: Inconsistency of external CSV and SDF file.'
                                             ' Files do not contain identical molecules.'
                                             '\n\t2: Weight argument error: wrong input. '
                                             'Check presence of ":" and numeric range correctness.'
                                             '\n\t3: Only-atom-types argument error: wrong input. '
                                             'Check symbols of given elements and their delimiter.'
                                             '\n\t4: One of the "--csvfilter" or "--externalcsv" arguments is missing.'
                                             'Program requires to use both of them at once.'
                                             '\n\t5: CSVfilter argument error: wrong input. '
                                             'Check numeric range correctness and presence of delimiters.',
                                 epilog='End of help block. Now try it yourself. Good luck!')
parser.add_argument('filename')
parser.add_argument('--weight', type=str,
                    help='Returns molecule set of molecules with molecular mass given in the argument.\nUSAGE: type e.g.--weight "60:70"')
parser.add_argument('--onlyatomtypes', type=str, help='Returns molecule set containing molecules of given elements or their subset.')
parser.add_argument('--externalcsv', help='Allows user to filter molecules in sdf file in dependence on their characteristics given in external csv.'
                                          'Adding optional argument --csvfilter required.')
parser.add_argument('--csvfilter', nargs='+')
args = parser.parse_args()

Bonds = namedtuple('Bonds', ['first_atom', 'second_atom', 'bond'])
Arguments = namedtuple('Arguments', ['minweight', 'maxweight', 'requested_elements', 'csv_arguments', 'csv_contradict_args'])


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

    def __getitem__(self, i):
        return self.atoms[i]

    def __len__(self):
        return len(self.atoms)

    def __str__(self):
        return f'Molecule {self.name}: atoms:{self.atoms}, \nbonds:{self.bonds}'

    def maxdistance(self):
        maxdist = 0

        for m, atom_a in enumerate(self): # add .atoms everywhere next to 'self'
            for n, atom_b in enumerate(self[1:], start=m+1):
                distance = math.sqrt(sum((atom_a.coords[k] - atom_b.coords[k])**2 for k in range(3)))

                maxdist = max(distance, maxdist)

        return maxdist

    @property
    def molecular_mass(self):
        molecular_mass = sum(relative_atomic_masses[atom.element] for atom in self)
        # final_molecular_mass = molecular_mass * 1.66053904e-27

        return molecular_mass

    def maximum_bonds_statistics(self):
        maxbonds = [0] * len(self)

        for i in self.bonds:
            if i.bond > maxbonds[i.first_atom-1]:
                maxbonds[i.first_atom-1] = i.bond
            if i.bond > maxbonds[i.second_atom-1]:
                maxbonds[i.second_atom-1] = i.bond

        # molecule_maxbonds == list of tuples created for each molecule [('N', 1), ('0', 2)]
        # contains information about elements and their maximum bonds
        molecule_maxbonds = list(zip([atom.element for atom in self], maxbonds))
        molecule_statistics = Counter(i for i in molecule_maxbonds)

        return molecule_statistics

    def check_external_properties(self, csv_arguments, csv_contradict_args):
        checkups = []

        if len(csv_contradict_args) > 0:
            contradicts_checkups = []
            for m_property in csv_contradict_args:
                contradicts_checkups.append(m_property not in self.external_properties.values())
            checkups.append(all(contradicts_checkups))

        for attr, value in csv_arguments.items():
            # tuple indicates numeric range
            if type(value) == tuple:
                checkups.append(
                    self.external_properties[attr] >= float(value[0]) and self.external_properties[attr] <= float(value[1])
                )
            else:
                if len(csv_contradict_args) > 0:
                    print("CSV argument error: Can't use --csvfilter '^property' and 'type=property' at once."
                          "\nEnd of program.")
                    sys.exit(5)
                checkups.append(self.external_properties[attr] == value)

        return all(checkups)

    def check_weight(self, minweight, maxweight):
        return self.molecular_mass >= float(minweight) and self.molecular_mass <= float(maxweight)

    def check_elements(self, requested_elements):
        molecule_elements = set(Atom.element for Atom in self)
        return molecule_elements <= requested_elements


class MoleculeSet:
    def __init__(self, molecules: List[Molecule]):
        self.molecules = molecules

    def __getitem__(self, i):
        return self.molecules[i]

    def __len__(self):
        return len(self.molecules)

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

                except ValueError:
                    break

                while True:
                    line = file.readline()
                    moleculestring += line
                    if 'END' in line:
                        molecule = Molecule(name, atoms, bonds, moleculestring)
                        molecules.append(molecule)
                    if '$$$$' in line:
                        break

        return MoleculeSet(molecules)

    def get_statistics(self):
        numeric_properties = defaultdict(list)
        categorical_properties = defaultdict(list)
        atom_types_statistics = Counter()

        for molecule in self:
            atom_types_statistics += molecule.maximum_bonds_statistics()
            numeric_properties['molecular_mass'].append(molecule.molecular_mass)
            numeric_properties['maximum_distance'].append(molecule.maxdistance())  # maximum_distance or molecule_size?

            if args.externalcsv:
                for m_property, value in molecule.external_properties.items():
                    if type(value) == float:
                        numeric_properties[m_property].append(value)
                    else:
                        categorical_properties[m_property].append(value)

        np_to_print = defaultdict(list)   # np == numeric properties
        cp_to_print = {}   # cp == categorical properties
        for key, values in numeric_properties.items():
            np_to_print[key].append(statistics.mean(values))
            try:
                np_to_print[key].append(statistics.stdev(values))
            except statistics.StatisticsError:
                np_to_print[key].append(0)

        for key, values in categorical_properties.items():
            cp_to_print[key] = Counter(values)

        return atom_types_statistics, np_to_print, cp_to_print

    def print_statistics(self, print_detailed=False):
        print('MoleculeSet:')
        if print_detailed:
            for n, molecule in enumerate(self, start=1):
                print(str("%3d" % n) + ':', molecule.name, end=" ")
                print(f'Relative molecular mass = {molecule.molecular_mass:.3f}', end=' ')
                print(set([Atom.element for Atom in molecule]))
            print()

        atom_types_statistics, numeric_properties, categorical_properties = self.get_statistics()
        print('### NUMERIC AND CATEGORICAL PROPERTIES')
        print('{:>35} {:>10}'.format('MEAN', 'SD'))

        for key, values in numeric_properties.items():
            print('{:25s}{:10.3f} {:10.3f}'.format(key.capitalize(), values[0], values[1]))

        if len(categorical_properties) > 0:
            print('{:>46}'.format('ABSOLUTE (RELATIVE) FREQUENCY'))
            for m_property, counter in categorical_properties.items():
                amount = sum(count for count in counter.values())
                print(m_property.capitalize())
                for item, frequency in counter.items():
                    print('\t{:17s}{:>3d}{:10.3f}'.format(item, frequency, frequency/amount))

        print("\n### ATOM TYPES:")
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(atom_types_statistics)
        print(2*'\n')


    def add_external_properties(self):
        molecule_storage = {molecule.name: molecule for molecule in self}

        with open(args.externalcsv) as csvfile:
            line = csvfile.readline()
            properties = [m_property.strip() for m_property in line.split(';')][1:]
            # CSV CHECKUPS
            if len(self) != len(csvfile.readlines()):
                print(
                    f"External CSV Error: File {args.filename} and {args.externalcsv} do not contain same number of molecules."
                    f"\nEnd of program.")
                sys.exit(1)

            csvfile.seek(0)
            reader = csv.DictReader(csvfile, delimiter=';')
            for row in reader:
                if row['name'] not in molecule_storage.keys():
                    print(f"External CSV Error: File {args.filename} doesn't contain molecule {row['name']}."
                          f"Check it and try again.\nEnd of program.")
                    sys.exit(1)
                else:
                    # ADDING NEW MOLECULE ATTRIBUTE
                    external_properties = {}
                    for m_property in properties:
                        try:
                            external_properties[m_property] = float(row[m_property])
                        except ValueError:
                            external_properties[m_property] = row[m_property]
                    molecule_storage[row['name']].external_properties = external_properties

    def filter_molecules_by_properties(self, arguments):
        molecules = []

        for molecule in self:
            if args.externalcsv and not molecule.check_external_properties(arguments.csv_arguments, arguments.csv_contradict_args):
                continue
            if args.weight and not molecule.check_weight(arguments.minweight, arguments.maxweight):
                continue
            if args.onlyatomtypes and not molecule.check_elements(arguments.requested_elements):
                continue
            molecules.append(molecule)

        return MoleculeSet(molecules)

    def create_new_SDFfile(self, filename):
        with open(filename + '.sdf', 'w') as file:
            for molecule in self:
                file.write(molecule.moleculestring)
                file.write('$$$$\n')


def get_atomic_masses():
    with open('Periodic_Table_Of_Elements.csv') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            relative_atomic_masses[row['Symbol']] = float(row['Atomic Weight'])


def check_arguments():
    minimum = None
    maximum = None
    requested_elements = None
    csv_arguments = {}
    csv_contradict_args = None

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
                sys.exit(2)
        else:
            print('Usage of ":" required. Check it and try again.\nEnd of program.')
            sys.exit(2)

    if args.onlyatomtypes:
        for i in [element.strip() for element in args.onlyatomtypes.split(',')]:
            if i not in relative_atomic_masses.keys():
                print('--onlyatomtypes argument: wrong input. Check it and try again.\nEnd of program.')
                sys.exit(3)
        else:
            requested_elements = set([element.strip() for element in args.onlyatomtypes.split(',')])

    if (args.externalcsv and not args.csvfilter) or (args.csvfilter and not args.externalcsv):
        print("External CSV Error: Missing '--csvfilter' or '--externalcsv' argument. Check them and try again."
              "\nEnd of program.")
        sys.exit(4)

    if args.csvfilter:
        for arg in args.csvfilter:
            if '=' in arg:
                argument, value = arg.split('=')

                try:
                    mini, maxi = value.split(':')
                    if mini == '':
                        mini = 0
                    elif maxi == '':
                        maxi = float('inf')
                    else: pass  ###

                    try:
                        mini = float(mini)
                        maxi = float(maxi)
                    except ValueError:
                        # exception for typing error in numeric range input (e.g. 50:f70) (string in digits)
                        print(
                            "--csvfilter argument Error: wrong input. Problem with numeric range. Check it and try again.\nEnd of program."
                        )
                        sys.exit(5)
                    value = mini, maxi

                except ValueError:
                    # value.split(':') returns number different to 2
                    if len(value.split(':')) == 1:
                        pass
                    else:
                        # numeric range contains 2 colons or more
                        print("--csvfilter argument Error: wrong input. Problem with numeric range. Check it and try again.\nEnd of program.")
                        sys.exit(5)

                csv_arguments[argument] = value

            elif arg[0] == '^':
                csv_contradict_args = [item.strip() for item in arg[1:].split(',')]

            else:
                print("--csvfilter argument Error: wrong arguments. Check them and try again.\nEnd of program.")
                sys.exit(5)

    return Arguments(minimum, maximum, requested_elements, csv_arguments, csv_contradict_args)

if __name__ == '__main__':
    relative_atomic_masses = dict()
    get_atomic_masses()

    arguments = check_arguments()

    sdf_set = MoleculeSet.load_from_file(args.filename)
    if args.externalcsv:
        sdf_set.add_external_properties()
    sdf_set.print_statistics()   # calling MoleculeSet.get_statistics() inside the print_stat() function

    set2 = sdf_set.filter_molecules_by_properties(arguments)
    set2.print_statistics(print_detailed=False)
    #set2.create_new_SDFfile(args.filename[:-4]+'_filter_by_props_external_csv')

# rozdelit atom, Molecule, MoleculeSet do  souboru -> nacitat jako moduly
# python data model
# from collectons import choice
