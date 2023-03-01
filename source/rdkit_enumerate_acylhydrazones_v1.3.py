
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions, Draw
import numpy as np
import pandas as pd
import sys


def smi2cansmi(smi):
    '''convert SMILES to RDKit-canonicalized form'''
    mol = Chem.MolFromSmiles(smi)
    if mol:
        for atm in mol.GetAtoms():
            atm.SetIsotope(0)
        return Chem.MolToSmiles(mol)
    else:
        return None

def smi2cansmi2(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        return Chem.MolToSmiles(mol)
    else:
        return None

# load aldehyde reactants
incsv = sys.argv[1] # "./data/aldehyde_reactants_cln.csv"
outcsv = open( sys.argv[2], "a" ) # output CSV to write for reactants and products
df = pd.read_csv( incsv, low_memory=False)
df = df.loc[ df['match_aldehyde'] == True ]
aldehyde_smiles = df['neutral_smiles'].tolist()
aldehyde_smiles = [ smi2cansmi2(smi) for smi in aldehyde_smiles ]
aldehyde_IDs = df['aldehyde_ID'].tolist()

# store aldehyde names in dictionary
aldehyde_dict = dict( zip(aldehyde_smiles, aldehyde_IDs) )

# 5 core hydrazide reactants
raw_cores = [ "C1=C(C=CC(=C1)N2CCC(NC2=O)=O)C(=O)NN",
              "C1(=CC=CC(=C1)N2CCC(NC2=O)=O)C([N][N])=O",
              "C1=C(C=CC(=C1)N2CCC(NC2=O)=O)OCC(=O)[N][N]",
              "C1=CC(=CC(=C1)N2CCC(NC2=O)=O)OCC([N][N])=O",
              "C1(CCN(C(N1)=O)N)=O" ]

# RDKit-canonicalize SMILES for hydrazide core molecules (also de-Kekulizes them)
cores = [ smi2cansmi2( smi ) for smi in raw_cores ]

# run chemical transformation
# hydrazide core + aldehyde --> acylhydrazone
# SMARTS for hydrazide substructure = "[C:1]([N:3][N:4])=[O:2]"
# SMARTS for aldehyde functional group = "[CX3H1:5](=O)[#6:6]"

# specify rxn for cis (rxn1) and trans (rxn2) products--couldn't get this to dump both with single reaction
rxn1 = rdChemReactions.ReactionFromSmarts("[C:1]([N:3][N:4])=[O:2].[CX3H1:5](=O)[#6:6]>>[N:3](\[N:4]=[C:5](/[#6:6])[H])[C:1]=[O:2]")
rxn2 = rdChemReactions.ReactionFromSmarts("[C:1]([N:3][N:4])=[O:2].[CX3H1:5](=O)[#6:6]>>[N:3](\[N:4]=[C:5](\[#6:6])[H])[C:1]=[O:2]")

# create lists to store reactants and products
mol_list = []

print( "hydrazide_name,hydrazide_smiles,aldehyde_name,aldehyde_smiles,product_name,product_smiles", file=outcsv )

# loop through each core hydrazide as a reactant (r1)
for core_idx, core in enumerate(cores):

    # build hydrazide core reactant molecules from SMILES
    r1 = Chem.MolFromSmiles( core )
    r1.SetProp( "_Name", "core-"+str(core_idx+1) )
    # build aldehyde reactant molecules from SMILES
    for a in aldehyde_smiles:

        # look up aldehyde reactant name, e.g., "aldehyde-36"
        aldehyde_name = aldehyde_dict[a]
        r2 = Chem.MolFromSmiles( a )
        r2.SetProp( "_Name", aldehyde_dict[a] )
        # run reaction/transformation
        products1 = rxn1.RunReactants( (r1, r2) )
        products2 = rxn2.RunReactants( (r1, r2) )

        # loop through products in case multiple products are generated (multiple aldehyde sites on r2)
        for p in range(len(products1)):

            # set names for cis and trans products
            products1[p][0].SetProp( "_Name", r1.GetProp("_Name")+"_"+r2.GetProp("_Name")+"_cis-"+str(p+1) )
            products2[p][0].SetProp( "_Name", r1.GetProp("_Name")+"_"+r2.GetProp("_Name")+"_trans-"+str(p+1) )

            # https://github.com/rdkit/rdkit/issues/1596
            products1[p][0].UpdatePropertyCache()
            products2[p][0].UpdatePropertyCache()

            # print reactants and products (names and SMILES)
            print(  r1.GetProp("_Name")+","+core+","+r2.GetProp("_Name")+","+a+","+products1[p][0].GetProp("_Name")+","+Chem.MolToSmiles( products1[p][0] ), file=outcsv )
            print(  r1.GetProp("_Name")+","+core+","+r2.GetProp("_Name")+","+a+","+products2[p][0].GetProp("_Name")+","+Chem.MolToSmiles( products2[p][0] ), file=outcsv )

            mol_list.extend( [r1, r2, products1[p][0]] )
            mol_list.extend( [r1, r2, products2[p][0]] )


outcsv.close()

# draw molecules to grid 

# align on phenyldihydrouracil core structure (hydrouracil or phenyldihydrouracil)
#template = Chem.MolFromSmiles("C1=CC=CC(=C1)N2CCC(NC2=O)=O")
#template = Chem.MolFromSmiles("C1=C(C=CC(=C1)N2CCC(N(C2=O)[H])=O)C(N)=O")
template = Chem.MolFromSmiles("C1(CCNC(N1)=O)=O")

# build 2D structure for template
AllChem.Compute2DCoords(template)

# try to align molecules on template when drawing structures
trunc_mol_list = mol_list[0:9] + mol_list[-9:]
for m in trunc_mol_list:
    try:
        AllChem.GenerateDepictionMatching2DStructure( m, template )
    except:
        #print( "error aligning: {}".format( m.GetProp("_Name") ) )
        pass

img = Draw.MolsToGridImage( trunc_mol_list, molsPerRow=3, subImgSize=(400,400), legends=[ m.GetProp("_Name") for m in trunc_mol_list ]  )
img.save( 'acylhydrazones_example_6_rxns.png' )

