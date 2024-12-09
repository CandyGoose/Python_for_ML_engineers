from rdkit import Chem


def substructure_search(molecules, substructure):
    sub_mol = Chem.MolFromSmiles(substructure)
    if not sub_mol:
        raise ValueError("Invalid substructure provided.")

    result = []
    for smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(sub_mol):
            result.append(smiles)
    return result
