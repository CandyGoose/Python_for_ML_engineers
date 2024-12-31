import pytest
from src.utils import substructure_search


def test_substructure_search_valid():
    molecules = [
        "CCO",
        "CCCO",
        "CCCCO",
    ]
    substructure = "O"

    result = substructure_search(molecules, substructure)
    assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed"


def test_substructure_search_invalid_substructure():
    molecules = [
        "CCO",
        "CCCO",
    ]
    substructure = "invalid_smiles"

    with pytest.raises(ValueError, match="Invalid substructure provided."):
        substructure_search(molecules, substructure)


def test_substructure_search_empty_molecules():
    molecules = []
    substructure = "O"

    result = substructure_search(molecules, substructure)
    assert result == [], "Substructure search should return empty list for empty molecules list"


def test_substructure_search_some_match():
    molecules = [
        "CCO",
        "CCCO",
        "CCCCO",
        "CCN",
    ]
    substructure = "O"

    result = substructure_search(molecules, substructure)
    assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed to find matching molecules"
