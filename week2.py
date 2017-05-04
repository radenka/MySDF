# zpracovani sdf souboru, extrahovani informaci o atomech a jejich vazbach
# ukol: na zaklade struktury sfd souboru vytvorit dvojice atom: vazba pro vsechny atomu v molekule, 
     # pote ulozit informace o vaznosti / napr. C 1 (jednovazny uhlik): 2, H 1: 12 apod.

with open('Projekt ChI-BI.sdf') as file:
     
     n = 1
     while True:
          try:
               for i in range(3):         # zahozeni prvnich tri radku souboru 
                    file.readline()
          
          
               info = file.readline()     
               atom_count = int(info[0:2])   # prvni pozice ctvrteho radku: pocet atomu
               
               ###
               bond_count = int(info[3:5])   # druha pozice: pocet vazeb
          
               atoms = list() 
               for i in range(atom_count):                   ## for i in range(atom_count): atom = file.readline().split()[3] ; atoms.append(atom)
                    atoms.append(file.readline()[])    # vytvoreni seznamu prvku ['N', 'O', ...]
          
               maxbonds = [0 for i in range(atom_count)]     # seznam pro zapisovani vazeb prvku (jak to nascitavat?)
          
               bonds = list()  ###
               for i in range(bond_count):
                    line = file.readline().split()
                    bonds.append(int(line[2]))   # list vazeb jednotlivych dvojic atomu (pozice atomu dane dvojice jsou na prvni a druhe pozici v radku)
          
                    first = int(line[0])    # ulozeni pozice prvniho atomu vazby
                    second = int(line[1])   # ulozeni pozice druheho atomu vazby
          
                    if bonds[i] > maxbonds[first-1]:  # hodnota vazby se prepise pouze v pripade, ze aktualne ctena hodnota vazby je vyssi nez vazba ulozena
                         maxbonds[first-1] = bonds[i]   # kvuli indexovani je potreba od pozice 'first' odecist 1, aby pozice atomu odpovidala pozici v 'maxbonds'
                         
                    if bonds[i] > maxbonds[second-1]:   # v maxbonds jsou aktualne ulozene pozice atomu v prvnim sloupci, doplnuji druhy sloupec
                         maxbonds[second-1] = bonds[i]    # princip prepisovani je stejny (viz vyse)              
               
               print('sample ' + str(n) + ': ' + str(list(zip(atoms, maxbonds))) )
               # print(len(maxbonds))
               
          except IndexError:
               if bond_count == 0:
                    print('File error. Invalid data for sample ' + str(n) + '.')
               else:
                    print("No more samples to read. End of programme.")
                    break
          
          while True:
               line = file.readline()
               if '$$$$' in line:
                    n += 1
                    break
               
               # slovniky a tuples
               # vytvorit statistiku pro kazdu molekulu, + pro cely soubor