from os import getenv

from fastapi import FastAPI, HTTPException
from rdkit import Chem

from models import Molecule, SubstructureSearch

app = FastAPI()

molecule_store = {}

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

@app.post("/add")
def add_molecule(molecule: Molecule):
    if molecule.identifier in molecule_store:
        raise HTTPException(status_code=400, detail="Identifier already exists.")
    mol = Chem.MolFromSmiles(molecule.smiles)
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecule_store[molecule.identifier] = molecule.smiles
    return {"message": "Molecule added successfully.", "identifier": molecule.identifier}


@app.get("/get/{identifier}")
def get_molecule(identifier: str):
    if identifier not in molecule_store:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    return {"identifier": identifier, "smiles": molecule_store[identifier]}


@app.put("/update/{identifier}")
def update_molecule(identifier: str, molecule: Molecule):
    if identifier not in molecule_store:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    mol = Chem.MolFromSmiles(molecule.smiles)
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecule_store[identifier] = molecule.smiles
    return {"message": "Molecule updated successfully.", "identifier": identifier}


@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str):
    if identifier not in molecule_store:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    del molecule_store[identifier]
    return {"message": "Molecule deleted successfully.", "identifier": identifier}


@app.get("/list")
def list_molecules():
    return {"total_molecules": len(molecule_store), "molecules": molecule_store}


@app.post("/search")
def search_substructure(search: SubstructureSearch):
    sub_mol = Chem.MolFromSmiles(search.substructure)
    if not sub_mol:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string.")
    matching_molecules = [
        {"identifier": id, "smiles": smiles}
        for id, smiles in molecule_store.items()
        if Chem.MolFromSmiles(smiles).HasSubstructMatch(sub_mol)
    ]
    return {"total_found": len(matching_molecules), "matching_molecules": matching_molecules}
