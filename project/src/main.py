import json
import logging
from os import getenv

import redis
from fastapi import FastAPI, HTTPException, Depends
from rdkit import Chem
from sqlalchemy.orm import Session

from models import MoleculeDB, Molecule, SubstructureSearch, SessionLocal, engine

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

redis_client = redis.Redis(host='redis', port=6379, db=0)

MoleculeDB.metadata.create_all(bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 300):
    redis_client.setex(key, expiration, json.dumps(value))


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

    logger.info(f"Added new molecule with identifier: {molecule.identifier}")

    return {"message": "Molecule added successfully.", "identifier": molecule.identifier}


@app.get("/get/{identifier}")
def get_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")

    logger.info(f"Retrieved molecule with identifier: {identifier}")

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

    logger.info(f"Updated molecule with identifier: {identifier}")

    return {"message": "Molecule updated successfully.", "identifier": identifier}


@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")

    db.delete(db_molecule)
    db.commit()

    logger.info(f"Deleted molecule with identifier: {identifier}")

    return {"message": "Molecule deleted successfully.", "identifier": identifier}


@app.get("/list")
def list_molecules(limit: int = 100, offset: int = 0, db: Session = Depends(get_db)):
    logger.info(f"Fetching {limit} molecules starting from offset {offset}.")

    cache_key = f"list_molecules:{limit}:{offset}"

    cached_result = get_cached_result(cache_key)
    if cached_result:
        logger.info(f"Cache hit for listing molecules: {cache_key}")
        return cached_result

    def molecule_iterator(limit, offset):
        query = db.query(MoleculeDB).offset(offset).limit(limit)
        for molecule in query.all():
            yield {"identifier": molecule.identifier, "smiles": molecule.smiles}

    total_molecules = db.query(MoleculeDB).count()
    molecules = list(molecule_iterator(limit, offset))

    response_data = {
        "total_molecules": total_molecules,
        "molecules": molecules
    }

    set_cache(cache_key, response_data)

    return response_data


@app.post("/search")
def search_substructure(search: SubstructureSearch, db: Session = Depends(get_db)):
    sub_mol = Chem.MolFromSmiles(search.substructure)
    if not sub_mol:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string.")

    cache_key = f"search:{search.substructure}"

    cached_result = get_cached_result(cache_key)
    if cached_result:
        logger.info(f"Cache hit for substructure search: {cache_key}")
        return {"source": "cache", "data": cached_result}

    matching_molecules = db.query(MoleculeDB).filter(MoleculeDB.smiles.like(f"%{search.substructure}%")).all()

    search_result = {
        "total_found": len(matching_molecules),
        "matching_molecules": [
            {"identifier": mol.identifier, "smiles": mol.smiles} for mol in matching_molecules
        ]
    }

    set_cache(cache_key, search_result)

    logger.info(f"Substructure search for {search.substructure} returned {len(matching_molecules)} results.")

    return search_result
