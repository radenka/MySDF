with open("sdf.txt") as file:
    atom_block = False
    bond_block = False
    atoms = [];    nums1 = [];    nums2 = [];    bonds = []
    
    for line in file: 
        data = line.strip().split()
   
        if len(data) == 11:
            atom_count = int(data[0])
            bond_count = int(data[1])
            atom_block = True
            
        elif atom_block:
            if len(data) != 16:
                atom_block = False
                bond_block = True
            atoms.append(data[3])
        
        elif bond_block:
            if len(data) != 7:
                bond_block = False
                break
            nums1.append(int(data[0]))
            nums2.append(int(data[1]))
            bonds.append(int(data[2]))             
    