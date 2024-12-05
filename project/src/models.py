from sqlalchemy import create_engine, Column, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from pydantic import BaseModel
import os

SQLALCHEMY_DATABASE_URL = os.getenv("DATABASE_URL", "postgresql+psycopg2://user:password@db:5432/db")

engine = create_engine(SQLALCHEMY_DATABASE_URL, echo=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()


class MoleculeDB(Base):
    __tablename__ = "molecules"

    identifier = Column(String, primary_key=True, index=True)
    smiles = Column(String, index=True)


class Molecule(BaseModel):
    identifier: str
    smiles: str

    class Config:
        orm_mode = True


class SubstructureSearch(BaseModel):
    substructure: str
