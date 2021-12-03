#!/usr/bin/python

# Hailey Wallace
# October 25, 2021
# Parses out all functional groups that COMBS recognizes from ligand of interest

import sys, argparse, os, re, subprocess
import pandas as pd
from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, ChemicalFeatures, rdchem, rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
from IPython.display import SVG

import time
from datetime import timedelta
start_time = time.time()

####################
#os.remove("ligand_CG_coords.txt")

# Parser and commandline options
parser = argparse.ArgumentParser()
parser.add_argument('--ligand', '-l', type=str, required=True)
parser.add_argument('--coords', '-c', type=str, required=False)
parser.add_argument('--output', '-o', type=str, required=False)
parser.add_argument('--image', '-i', type=str, required=False, default=False)

if len(sys.argv)==1:
    help = """
*********************************************************
*       Ligand Matcher parses COMBS-compatible          *
*        functional groups from input ligand            *
*      __________________________________________       *
*                     options:                          *
*     --ligand, -l: input mol2 ligand (required)        *
*     --ouput, -o: output name for ligand.txt file      *
*                   "ligand.txt" if not specified       *
*     --coords, -c: output coordinate file              *
*                   "none" for no coords/ligand file;   *
*                   will only output number of matches  *
*     --image, -i: ouputs png images with               *
*                   highlighted functional groups       *
*                                                       *
*********************************************************
    """
    print(help)
    parser.exit()
args = parser.parse_args()

# If options are specified:
if args.coords == "none":
    print('Input ligand:', args.ligand, '\n'"No coordinates file.")
elif args.coords == None:
    print('Input ligand:', args.ligand, '\n'"Output coordinate file: ligand_CG_coords.txt")
    PATH = './ligand_CG_coords.txt'
    if os.path.isfile(PATH):
        print("File already exists. Moving previous ligand_CG_coords.txt to ./old_files")
        PATH2 = './old_files/'
        isExist = os.path.exists(PATH2)
        # Directory does exist
        if isExist:
            os.rename(PATH, './old_files/ligand_CG_coords.txt')
        # Create a new directory because it does not exist
        if not isExist:
            os.mkdir(PATH2)
            os.rename(PATH, './old_files/ligand_CG_coords.txt')
else:
    print('Input ligand:', args.ligand, '\n'"Output coordinate file:", args.coords)

########################################
# Functional groups that COMBS recognizes
########################################
conh2 = Chem.MolFromSmarts("[C,c]C(=O)N([H])")
bb_cco = Chem.MolFromSmarts("[N,n,C,c][C,c;X3](=O)[!O]")
ph = Chem.MolFromSmarts("c1ccccc1")
bb_cnh = Chem.MolFromSmarts("[*][N,n;X3][H]")
ccn = Chem.MolFromSmarts("[C,c][C,c]N([H])([H])[H]")
ccoh = Chem.MolFromSmarts("[C,c][CX4]O[H]")
coh = Chem.MolFromSmarts("[CX4]O[H]")
coo = Chem.MolFromSmarts("CC(=O)O")
csc = Chem.MolFromSmarts("C([H])([H])SC([H])([H])[H]")
csh = Chem.MolFromSmarts("CS[H]")
gn = Chem.MolFromSmarts("[C,c][N,n]([H])[C,c](~N([H])[H])N([H])[H]")
hid = Chem.MolFromSmarts("c1c([H])[n;X2]c([H])n1([H])")
hie = Chem.MolFromSmarts("c1c([H])n([H])c([H])[n;X2]1")
hip = Chem.MolFromSmarts("c1c([H])n([H])c([H])n1([H])")
indole = Chem.MolFromSmarts("c21c(c(c(c(c1n(c(c2)[H])[H])[H])[H])[H])[H]")
phenol = Chem.MolFromSmarts("c1(c(c(c(c(c1)[H])[H])O[H])[H])[H]")
isopropyl = Chem.MolFromSmarts("C([H])([H])([H])C([H])C([H])([H])[H]")
pro = Chem.MolFromSmarts("C1([H])([H])C([H])([H])**C1([H])([H])")
ch3 = Chem.MolFromSmarts("[C,c][C,c]C([H])([H])[H]")

func_groups = [conh2, bb_cco, ph, bb_cnh, ccn, ccoh, coh, coo, csc, csh, gn, hid, hie, hip, indole, phenol, isopropyl, pro, ch3]

########################################
# Setting up file and SMARTS strings
########################################
input = args.ligand

# Import mol2 ligand file
file = Chem.MolFromMol2File(input)
#file = Chem.AddHs(start, addCoords=True, addResidueInfo=True)

# SMARTS pattern of your input mol2
ligand = Chem.MolToSmarts(file)
try:
  print("Ligand SMARTS pattern:", ligand)
except NameError:
  print("Mol2 is not defined")

########################################
# Find matching substructures, output to new outfile
########################################
# To return the variable name
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

# Iterate functional group SMARTS patterns through your ligand
for combs_groups in func_groups:
    substruct = file.GetSubstructMatches(combs_groups)

    # assign new variable with old variable name, aka the functional group name
    pre_var= str(namestr(combs_groups, globals())).split("'")
    var = pre_var[1]

    # if -c "none", do not make text file
    if args.coords == "none":
        pass

    # if no argument, make default: "ligand_CG_coords.txt"
    elif args.coords == None:
        f = open("ligand_CG_coords.txt", "a")
        for types in substruct:
            coords = Chem.MolFragmentToCXSmiles(file,types)
            print(var, types, coords, file=f)
        f.close()

    # if specified, make custom-named text file
    else:
        output = args.coords
        f = open(output, "a")
        for types in substruct:
            coords = Chem.MolFragmentToCXSmiles(file,types)
            print(var, types, coords, file=f)
        f.close()

    # Count COMBS functional groups in ligand
    if len(substruct) != 0:
        print(var, "found", len(substruct), "matches.")

########################################
# Optional -- Image creation --
# print substructure match to a JPG of your ligand structure only if -i is specified
########################################
# If -i is specified
if args.image:
    # Image options
    image = args.image
    IPythonConsole.drawOptions.addAtomIndices = True
    IPythonConsole.molSize = 300,300

    # Iterate over substructure matches
    for combs_groups_img in func_groups:
        substruct2 = file.GetSubstructMatches(combs_groups_img)
        if len(substruct2) != 0:
            for group in substruct2:
                i = 0
                while os.path.exists(image+"%s.png" % i):
                    i += 1

                # Make a label for the images using functional group variable names
                name = str(namestr(combs_groups_img, globals())).split("'")
                label = str(name[1])

                # Atom indices for highlight
                highlights = list(group)
                print("Atom indices",highlights,"highlighted in",image+"%s.png" % i,"for",label)

                # Drawing the image -- still has Hydrogens because of the H atom indices used in highlight
                d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(file, highlightAtoms=highlights, legend=label)
                d.FinishDrawing()
                d.WriteDrawingText(image+"%s.png" % i)

########################################
# Parse the previously output file for matching atom names
# Precise atom names are crucial for COMBS
########################################
# Extract coordinates and match with group name
#def func_coord(coords):
#    for line in coords:
        # group name
#        a = line.split(" ")
#        name = str(a[0])
        # x, y, z coords
#        a0 = str(line).split("|")
#        a1 = a0[1]+','+name
#        a2 = re.sub('[()]','', a1)
#        xyz = re.sub(';', ','+name+'\n', a2)
#        print(xyz+'\n')

# Open previously output file, assign it to variable
#if args.coords == "none":
#    pass
#elif args.coords == None:
#    f = open("ligand_match_coords.txt","r")
#else:
#    f = open(args.coords, "r")


# Proceed if variable is assigned (if file was created)
#try:
#    f
#except NameError:
#    print("...Finished!")
#else:
#    print(f)
#    process = subprocess.run(['./ligand_support.sh','-l',input,'-i',f])


if args.coords == "none":
    elapsed_time_secs = time.time() - start_time
    msg = "...Finished in %s seconds." % timedelta(seconds=round(elapsed_time_secs))
    print(msg)

elif args.coords == None:
    if args.output == None:
        process = subprocess.run(['./ligand_support.sh','-l',input,'-i','ligand_CG_coords.txt'])
        elapsed_time_secs = time.time() - start_time
        msg = "Finished making ligand.txt in %s seconds" % timedelta(seconds=round(elapsed_time_secs))
        print(msg)
    else:
        process = subprocess.run(['./ligand_support.sh','-l',input,'-i','ligand_CG_coords.txt','-o',args.output])
        elapsed_time_secs = time.time() - start_time
        msg = "%s seconds" % timedelta(seconds=round(elapsed_time_secs))
        print('Finished making', args.output, msg)

else:
    if args.output == None:
        process = subprocess.run(['./ligand_support.sh','-l',input,'-i',args.coords])
        elapsed_time_secs = time.time() - start_time
        msg = "Finished making ___ in %s seconds" % timedelta(seconds=round(elapsed_time_secs))
        print(msg)
    else:
        process = subprocess.run(['./ligand_support.sh','-l',input,'-i',args.coords,'-o',args.output])
        elapsed_time_secs = time.time() - start_time
        msg = "%s seconds" % timedelta(seconds=round(elapsed_time_secs))
        print('Finished making', args.output, msg)
