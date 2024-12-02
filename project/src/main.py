from fastapi import FastAPI, HTTPException, Depends
from sqlalchemy.orm import Session
from rdkit import Chem
from models import MoleculeDB, Molecule, SubstructureSearch, SessionLocal, engine
from os import getenv

app = FastAPI()

MoleculeDB.metadata.create_all(bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add")
def add_molecule(molecule: Molecule, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == molecule.identifier).first()
    if db_molecule:
        raise HTTPException(status_code=400, detail="Identifier already exists.")

    mol = Chem.MolFromSmiles(molecule.smiles)
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")

    db_molecule = MoleculeDB(identifier=molecule.identifier, smiles=molecule.smiles)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)

    return {"message": "Molecule added successfully.", "identifier": molecule.identifier}


@app.get("/get/{identifier}")
def get_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    return {"identifier": db_molecule.identifier, "smiles": db_molecule.smiles}


@app.put("/update/{identifier}")
def update_molecule(identifier: str, molecule: Molecule, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")

    mol = Chem.MolFromSmiles(molecule.smiles)
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")

    db_molecule.smiles = molecule.smiles
    db.commit()
    db.refresh(db_molecule)

    return {"message": "Molecule updated successfully.", "identifier": identifier}


@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")

    db.delete(db_molecule)
    db.commit()

    return {"message": "Molecule deleted successfully.", "identifier": identifier}


@app.get("/list")
def list_molecules(db: Session = Depends(get_db)):
    molecules = db.query(MoleculeDB).all()
    return {
        "total_molecules": len(molecules),
        "molecules": [{"identifier": mol.identifier, "smiles": mol.smiles} for mol in molecules]
    }


@app.post("/search")
def search_substructure(search: SubstructureSearch, db: Session = Depends(get_db)):
    sub_mol = Chem.MolFromSmiles(search.substructure)
    if not sub_mol:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string.")

    matching_molecules = db.query(MoleculeDB).filter(MoleculeDB.smiles.like(f"%{search.substructure}%")).all()

    return {"total_found": len(matching_molecules), "matching_molecules": [
        {"identifier": mol.identifier, "smiles": mol.smiles} for mol in matching_molecules
    ]}
