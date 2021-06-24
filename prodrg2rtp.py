# define bond, angle and dihedral parameters

def format_atoms(line1):
    if line1[0] == ";":
        formated_string = "; Not required"
    else:
        formated_string = '{:>5} {:>6} {:>12} {:>6}'.format(line1[4], line1[1], line1[6], line1[5])
    return formated_string

def format_bonds(line1):
    if line1[0] == ";":
        formated_string = "; Not required"
    else:
        formated_string = '{:>5} {:>6} {:>12} {:>6}'.format(line1[-2], line1[-1], line1[3], line1[4])
        formated_string = formated_string +" "
        # find similarity and append
        for ab in bonds:
            if 0.99*float(line1[3]) < float(ab.avg_dis) < 1.01*float(line1[3]) and 0.99*float(line1[4]) < float(ab.forc_constant) < 1.01*float(line1[4]):
                formated_string += ab.bond_type +"; "                
    return formated_string

def format_angles(line1):
    if line1[0] == ";":
        formated_string = ";  ai    aj    ak   gromos type"
    else:
        formated_string = '{:>5} {:>6} {:>6} {:>12} {:>6}'.format(line1[-3], line1[-2], line1[-1], line1[4], line1[5])
        formated_string = formated_string +" "
        # find similarity and append
        for ab in angles:
            if 0.99*float(line1[4]) < float(ab.avg_ang) < 1.01*float(line1[4]) and 0.99*float(line1[5]) < float(ab.forc_constant) < 1.01*float(line1[5]):
                formated_string += ab.angle_type +"; " 
    return formated_string       

def format_pairs(line1):
    if line1[0] == ";":
        formated_string = ";  ai    aj"
    else:
        formated_string = '{:>5} {:>6}'.format(line1[-2], line1[-1])
    return formated_string
       
def format_dihedrals(line1):
    if line1[0] == ";":
        formated_string = ";  ai    aj    ak    al   gromos type"
    else:
        formated_string = '{:>5} {:>6} {:>6} {:>6} {:>12} {:>6}'.format(line1[-4], line1[-3], line1[-2], line1[-1], line1[5], line1[6])
        # find similarity and append
        for ab in angles:
            if 0.9*float(line1[5]) < float(ab.avg_ang) < 1.1*float(line1[5]) and 0.9*float(line1[6]) < float(ab.forc_constant) < 1.1*float(line1[6]):
                formated_string += ab.dih_type +"; " 
    return formated_string
 

class Bond:

    def __init__(self, bond_type = None, avg_dis = '0', forc_constant = '0'):
        self.bond_type = bond_type
        self.avg_dis = avg_dis
        self.forc_constant = forc_constant
    def __str__(self):
        return "Type:" + self.bond_type + " Average distance:" + self.avg_dis + "Force constant:" + self.forc_constant


class Angle:
    def __init__(self, angle_type, avg_ang, forc_constant):
        self.angle_type = angle_type
        self.avg_ang = avg_ang
        self.forc_constant = forc_constant
    def __str__(self):
        return "Type:" + self.angle_type + " Average angle:" + self.avg_ang + " Force constant:" + self.forc_constant



class Improper:
    def __init__(self, dih_type, avg_ang, forc_constant):
        self.dih_type = dih_type
        self.avg_ang = avg_ang
        self.forc_constant = forc_constant


class Dihedral:
    def __init__(self, dih_type, avg_ang, forc_constant):
        self.dih_type = dih_type
        self.avg_ang = avg_ang
        self.forc_constant = forc_constant
    def __str__(self):
        return "Type:" + self.dih_type + " Average angle:" + self.avg_ang + " Force constant:" + self.forc_constant


class Pair:
    pass

# Read atom types
r0 = open("atomtypes.atp", "r")
atoms = {}
for y in r0:
    line = y.split()
    atoms.__setitem__(line[0], line[1])
r0.close()    

# Read atom bonded parameters
types = []
angles = []
bonds = []
dihedrals = []
previous = False
pre_angle = False
pre_bond = False
pre_dihedral = False
r1 = open("ffbonded.itp", "r")
for x in r1:
    line = x.split()
    if len(line) >= 1:
        if line[0] == "#define":
            previous = True
            types.append(line[1])
            if "ga" in line[1]:
                pre_angle = True
                a = Angle(line[1], line[2], line[3])
                angles.append(a)
            elif "gb" in line[1]:
                pre_bond = True
                b = Bond(line[1], line[2], line[3])
                bonds.append(b)
            elif "gd" in line[1]:
                pre_dihedral = True
                d = Dihedral(line[1], line[2], line[3])
                dihedrals.append(d)
        elif previous:
            previous = False
            x.replace(',', '')
            line = x.split()
            for y in line:
                if y not in atoms:
                    line.remove(y)
for i in angles:
    print (i)
for i in bonds:
    print (i)
for i in dihedrals:
    print (i)
    
# Open itp file.
in_file = input("input file name:")
out_file = input("Output file name:")
block = None
print_to_output = """ ; This programe is written by Dr Tushar Ranjan Moharana. 
Kindly report your query and trouble to tusharranjanmoharana@gmail.com"""  
r2 = open(out_file, "w") 
r2.write(print_to_output)
r3 = open(in_file, "r")
for a in r3:     
    line1 = a.split()
    if len(line1) < 2:
        continue
    if line1[0] == "[":
        if line1[1] == "atoms":
            block = "atoms"
            print_to_output = "[ atoms ] \n"            
            r2.write(print_to_output)
        elif line1[1] == "bonds":
            block = "bonds"
            print_to_output = "[ bonds ] \n"
            r2.write(print_to_output)
        elif line1[1] == "pairs":
            block = "pairs"
            print_to_output = "[ exclusions ] \n"
            r2.write(print_to_output)
        elif line1[1] == "angles":
            block = "angles"
            print_to_output = "[ angles ] \n"
            r2.write(print_to_output)
        elif line1[1] == "dihedrals":            
            block = "dihedrals"
            print_to_output = "[ impropers ] \n"
            r2.write(print_to_output)
    else:
        if block == "atoms":
            print_to_output = format_atoms(line1)
            r2.write(print_to_output +"\n")
        elif block == "bonds":
            print_to_output = format_bonds(line1)
            r2.write(print_to_output +"\n")
        elif block == "angles":
            print_to_output = format_angles(line1)
            r2.write(print_to_output +"\n")
        elif block == "pairs":
            print_to_output = format_pairs(line1)
            r2.write(print_to_output +"\n")
        elif block == "dihedrals":
            print_to_output = format_dihedrals(line1)
            r2.write(print_to_output +"\n")
        else:
            r2.write("-----The End-----")
      
r3.close()
r2.close()