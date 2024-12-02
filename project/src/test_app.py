import pytest
from fastapi.testclient import TestClient
from src.main import app

client = TestClient(app)

def test_add_and_search_complex_molecules():
    molecule_data = [
        {"identifier": "mol1", "smiles": "O"},
        {"identifier": "mol2", "smiles": "C"},
        {"identifier": "mol3", "smiles": "CCO"},
        {"identifier": "mol4", "smiles": "c1ccccc1"},
        {"identifier": "mol5", "smiles": "CC(=O)O"},
        {"identifier": "mol6", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
    
    for molecule in molecule_data:
        response = client.post("/add", json=molecule)
        assert response.status_code == 200, f"Failed to add {molecule['identifier']}"
    
    search_data = {"substructure": "C(=O)O"}
    response = client.post("/search", json=search_data)
    assert response.status_code == 200
    assert response.json()['total_found'] == 2

def test_add_empty_smiles():
    molecule_data = {"identifier": "mol4", "smiles": "C"}
    response = client.post("/add", json=molecule_data)
    assert response.status_code == 400
    assert response.json()['detail'] == "Identifier already exists."

def test_substructure_search():
    search_data = {"substructure": "CC"}
    response = client.post("/search", json=search_data)
    assert response.status_code == 200
    assert response.json()['total_found'] == 3

def test_aspirin_substructure_search():
    search_data = {"substructure": "CC(=O)Oc1ccccc1C(=O)O"}
    response = client.post("/search", json=search_data)
    assert response.status_code == 200
    assert response.json()['total_found'] == 1

def test_empty_database_search():
    search_data = {"substructure": "CCOO"}
    response = client.post("/search", json=search_data)
    assert response.status_code == 200
    assert response.json()['total_found'] == 0
