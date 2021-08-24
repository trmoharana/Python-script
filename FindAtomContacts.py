# This script find all the atom within a given distance from a group of atoms
# This script ignores hydrogen and ligand atoms. If you want hydrogen or ligand, remove line 38 and 39
from Bio.PDB import *
r_cutoff = input("Cutoff distance in Angstrom:")
in_file = input("input file name:")
out_file = input("output file name:")
parser = PDBParser(PERMISSIVE=1)
r1 = open(out_file, "w")
output = str("Neighbour distance:" + r_cutoff + "\n")
r1.write(output)
r0 = open(in_file, "r")
for x in r0:
    line = x.split()
    pdb_file = line[0]
    protein = parser.get_structure("", pdb_file)
    allatoms = []
    for atom in protein.get_atoms():
        allatoms.append(atom)
    all_atoms = NeighborSearch(allatoms)
    models = list(protein.get_models())
    chains = list(models[0].get_chains())
    residues = list(chains[0].get_residues())
    for y in residues:
        residue_id = y.get_full_id()
        resseq = residue_id[3][1]
        line1 = set(line)
        if (y.get_resname() != "ALA") and (y.get_resname() != "CYS") and (y.get_resname() != "ASP") and (y.get_resname() != "GLU") and (y.get_resname() != "PHE") and (y.get_resname() != "GLY") and (y.get_resname() != "HIS") and (y.get_resname() != "ILE") and (y.get_resname() != "LYS") and (y.get_resname() != "LEU") and (y.get_resname() != "MET") and (y.get_resname() != "ASN") and (y.get_resname() != "PRO") and (y.get_resname() != "GLN") and (y.get_resname() != "ARG") and (y.get_resname() != "SER") and (y.get_resname() != "THR") and (y.get_resname() != "SEC") and (y.get_resname() != "VAL") and (y.get_resname() != "TRP") and (y.get_resname() != "TYR"):
            continue
        if str(resseq) in line1:
            atoms = list(y.get_atoms())
            for z in atoms:
                if (z.get_name() != "CA") and (z.get_name() != "CB") and (z.get_name() != "CG") and (z.get_name() != "CD") and (z.get_name() != "NE2") and (z.get_name() != "OE1") and (z.get_name() != "SG") and (z.get_name() != "CG1") and (z.get_name() != "CG2") and (z.get_name() != "ND2") and (z.get_name() != "OD1") and (z.get_name() != "CD1") and (z.get_name() != "CD2") and (z.get_name() != "OG1") and (z.get_name() != "NE") and (z.get_name() != "CZ") and (z.get_name() != "NH1") and (z.get_name() != "NH2") and (z.get_name() != "ND1") and (z.get_name() != "CE1") and (z.get_name() != "OG") and (z.get_name() != "OD2") and (z.get_name() != "OE2") and (z.get_name() != "CZ") and (z.get_name() != "CE") and (z.get_name() != "") and (z.get_name() != "NZ") and (z.get_name() != "CE2") and (z.get_name() != "CE3") and (z.get_name() != "NE1") and (z.get_name() != "CZ2") and (z.get_name() != "CZ3") and (z.get_name() != "CH2") and (z.get_name() != "OH"):
                    continue
                neighbors = list(all_atoms.search(z.get_coord(), int(r_cutoff), "A"))
                output = "\t".join([str(line[0]), str(resseq), y.get_resname(), z.get_name()])
                r1.write(output)
                for w in neighbors:
                    if (w.get_name() != "N") and (w.get_name() != "C") and (w.get_name() != "O") and (w.get_name() != "CA") and (w.get_name() != "CB") and (w.get_name() != "CG") and (w.get_name() != "CD") and (w.get_name() != "NE2") and (w.get_name() != "OE1") and (w.get_name() != "SG") and (w.get_name() != "CG1") and (w.get_name() != "CG2") and (w.get_name() != "ND2") and (w.get_name() != "OD1") and (w.get_name() != "CD1") and (w.get_name() != "CD2") and (w.get_name() != "OG1") and (w.get_name() != "NE") and (w.get_name() != "CZ") and (w.get_name() != "NH1") and (w.get_name() != "NH2") and (w.get_name() != "ND1") and (w.get_name() != "CE1") and (w.get_name() != "OG") and (w.get_name() != "OD2") and (w.get_name() != "OE2") and (w.get_name() != "CZ") and (w.get_name() != "CE") and (w.get_name() != "") and (w.get_name() != "NZ") and (w.get_name() != "CE2") and (w.get_name() != "CE3") and (w.get_name() != "NE1") and (w.get_name() != "CZ2") and (w.get_name() != "CZ3") and (w.get_name() != "CH2") and (w.get_name() != "OH"):
                        continue
                    parent_res = w.get_parent()
                    parent_res_id = parent_res.get_full_id()
                    if (resseq != parent_res_id[3][1]) and (parent_res.get_resname() != "HOH"): 
                        output = "".join(["\t", parent_res.get_resname(), w.get_name()])
                        r1.write(output)
                r1.write("\n")
r1.close()
r0.close()    
