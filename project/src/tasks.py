import os

from rdkit import Chem
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from celery_worker import celery
from models import MoleculeDB

SQLALCHEMY_DATABASE_URL = os.getenv("DATABASE_URL", "postgresql://user:password@db/db")


@celery.task
def perform_substructure_search(substructure_smiles: str):
    engine = create_engine(SQLALCHEMY_DATABASE_URL, echo=True)
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    db = SessionLocal()

    sub_mol = Chem.MolFromSmiles(substructure_smiles)
    if not sub_mol:
        return {"error": "Invalid substructure SMILES string."}

    matching_molecules = db.query(MoleculeDB).filter(MoleculeDB.smiles.like(f"%{substructure_smiles}%")).all()

    result = {
        "total_found": len(matching_molecules),
        "matching_molecules": [
            {"identifier": mol.identifier, "smiles": mol.smiles} for mol in matching_molecules
        ]
    }

    db.close()
    return result
