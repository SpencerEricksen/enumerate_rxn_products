
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import SaltRemover
import pandas as pd
import sys


def alt_frag_remover( m ):
    '''alternative FRAGMENT REMOVER to remove some non-standard counterions
       that are not in RDKit's default salt fragment list: $RDBASE/Data/Salts.txt'''
    # split m into mol fragments, keep fragment with highest num atoms
    mols = list(Chem.GetMolFrags( m, asMols=True ))
    if (mols):
        mols.sort(reverse=True, key=lambda x: x.GetNumAtoms() )
        mol = mols[0]
    else:
        mol = None
    return mol


def neutralize_atoms(mol):
    '''Noel O'Boyle (Vincent Scalfani adapted code from RDKit)
       https://www.rdkit.org/docs/Cookbook.html'''
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


# load Weiping's aldehyde reactant table, dumped as CSV from Excel file
incsv = sys.argv[1]  # "./data/aldehyde_reactants.csv"

df = pd.read_csv( incsv, low_memory=False )
df = df.dropna(axis=0)

# add rdkit mol objects for each molecule
PandasTools.AddMoleculeColumnToFrame( df, 'Smiles', 'rdkit_mol', includeFingerprints=False )

# remove salts
_saltRemover = SaltRemover.SaltRemover()
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: _saltRemover.StripMol(x) if x is not None else None)
# apply alternative salt remover for optimal cleanliness
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: alt_frag_remover(x) if x is not None else None )

# generate clean smiles
df['smiles'] = df['rdkit_mol'].apply( lambda x: Chem.MolToSmiles(x) if x is not None else None )

# neutralize to help with protonation ambiguity (amines/carboxylates) for substructure
# matching and rxns later
df['rdkit_mol'] = df['rdkit_mol'].apply( lambda x: neutralize_atoms(x) if x is not None else None )

# replace smiles with clean smiles
df['neutral_smiles'] = df['rdkit_mol'].apply( lambda x: Chem.MolToSmiles(x) if x is not None else None )

# add flag to indicate whether an aldehyde functional group is recognized in each molecule from Weiping's table
frag_aldehyde = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
df['match_aldehyde'] = df['rdkit_mol'].apply( lambda x: x.HasSubstructMatch(frag_aldehyde) if x is not None else None )

# rename old SMILES column
df.rename( columns={'Smiles':'raw_smiles'}, inplace=True)
# drop the rdkit mol object column
df.drop( columns='rdkit_mol', inplace=True )

# write new CSV with clean SMILES
df.to_csv('aldehyde_reactants_cln.csv', index=False )

