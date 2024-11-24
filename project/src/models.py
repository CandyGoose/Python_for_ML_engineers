from pydantic import BaseModel


class Molecule(BaseModel):
    identifier: str
    smiles: str


class SubstructureSearch(BaseModel):
    substructure: str
