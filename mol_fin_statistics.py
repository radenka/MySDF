# zpracovani sdf souboru, extrahovani informaci o atomech a jejich vazbach
# na zaklade struktury sfd souboru vytvorit dvojice atom: max_vazba
# ulozit informace o vaznosti / napr. C 1 (jednovazny uhlik): 2, H 1: 12 apod.

with open('ChI_BI_project.sdf') as file:

    final_stat = dict()
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
            for i in range(atom_count):
                atoms.append(file.readline()[31:34].strip())  # vytvoreni seznamu prvku ['N', 'O', ...]

            maxbonds = [0] * atom_count   # seznam pro zapisovani max vazeb prvku

            for i in range(bond_count):
                line = file.readline()

                bonds = int(line[6:9])
                first = int(line[0:3])    # ulozeni pozice prvniho atomu vazby
                second = int(line[3:6])   # ulozeni pozice druheho atomu vazby

                # hodnota vazby se prepise pouze kdyz je aktualne ctena hodnota vazby je vyssi nez ulozena
                # je treba od pozice 'first' odecist 1, aby pozice atomu odpovidala pozici v 'maxbonds'
                if bonds > maxbonds[first-1]:
                    maxbonds[first-1] = bonds

                # v maxbonds jsou aktualne ulozene pozice atomu v prvnim sloupci, doplnuji druhy sloupec
                # princip prepisovani je stejny (viz vyse)
                if bonds > maxbonds[second-1]:
                    maxbonds[second-1] = bonds

            pairs_tuples = list(zip(atoms, maxbonds))
            mol_stat = dict()     # vytvari set(pairs_tuples); key = kolikrat se v molekule dany typ atomu vyskytuje

            
            for i in pairs_tuples:
                # pokud dany typ atomu neni klicem ve slovniku final, tak se tento klic vytvori, prirazena hodnota 1
                if i not in mol_stat.keys():
                    mol_stat[i] = 1
                else:
                    mol_stat[i] += 1

            # print('sample ' + str(n) + ': ' + str(mol_stat))

            for key in mol_stat:
                ## pokud atom neni klicem ve final_stat, tak se tento klic vytvori, prirazeni prislusne hodnoty
                if key not in final_stat:  
                    final_stat[key] = mol_stat[key]  # dict[key] = dict.get(key)   # puvodne final_stat[key] = mol_stat.get(key)
                ## v opacnem pripade se k aktualni hodnote pricte pocet atomu, ktery se v molekule vyskytuje
                else:
                    final_stat[key] += mol_stat[key]  

        except ValueError:
            break

        while True:
            line = file.readline()
            if '$$$$' in line:
                n += 1
                break

    print(final_stat)

