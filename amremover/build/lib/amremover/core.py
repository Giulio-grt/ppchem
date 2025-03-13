import re
from rdkit import Chem

def remove_atom_mapping(smiles: str) -> str:
    """
    Remove atom mapping from a reaction SMILES string.
    """
    return re.sub(r"(?<=[^\\*])(:\\d+)]", "]", smiles)

def canonicalize_smiles(smiles: str) -> str:
    """
    Canonicalize a SMILES string using RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    return ""

def remove_atom_mapping_and_canonicalize_rxn_smiles(rxn_smiles: str) -> str:
    """
    Remove atom mapping and canonicalize a reaction SMILES string.
    """
    reactants, reagents, products = rxn_smiles.split(">")
    reactants = ".".join([canonicalize_smiles(remove_atom_mapping(s)) for s in reactants.split(".")])
    products = ".".join([canonicalize_smiles(remove_atom_mapping(s)) for s in products.split(".")])
    return f"{reactants}>>{products}"
