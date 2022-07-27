"""
Simply removes solvent from a PDB file
"""

import argparse

# I/O
parser = argparse.ArgumentParser(description="Remove solvent from a pdb file")
parser.add_argument("-i", "--input", dest="input", required=True, help="Input pdb file")
parser.add_argument("-o", "--output", dest="output", required=True, help="Output file")
args = parser.parse_args()

# Arguments
inputfile = args.input
outputfile = args.output

# Initialize variables
ATOM_list = []


def read_pdb(pdbfile):
    """
    Read through a pdb file and save the essential protein model lines
    """
    with open(pdbfile, "r") as pdb:
        for line in pdb:
            if line.startswith(("ATOM", "TER", "MODEL", "ENDMDL", "END")):
                ATOM_list.append(line)
    return ATOM_list


def write_pdb(outputfile):
    """
    Write out a pdb file with only the essential protein model lines
    """
    with open(outputfile, "w") as out:
        for line in ATOM_list:
            out.write(line)


ATOM_list = read_pdb(inputfile)
write_pdb(outputfile)
