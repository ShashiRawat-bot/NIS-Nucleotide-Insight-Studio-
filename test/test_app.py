import sys
import os
import pytest
from app import app
from utils.sequence_utils import transcribe_dna_to_rna

# Ensure the project root directory is included in sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Create a test client for Flask
@pytest.fixture
def client():
    with app.test_client() as client:
        yield client

# Test the index route
def test_index(client):
    response = client.get('/')
    assert response.status_code == 200
    assert b"DNA Sequence" in response.data

# Test the transcribe_dna_to_rna function
def test_transcribe_dna_to_rna():
    result = transcribe_dna_to_rna("ATGC")
    assert result == "AUGC"

# Test CSV export functionality
def test_export_csv(client):
    response = client.post('/', data={'dna_sequence': 'ATGC', 'export_csv': True})
    assert 'text/csv' in response.headers['Content-Type']
    assert response.data.startswith(b"DNA Sequence,RNA Sequence,Protein Sequence")

# Test the BLAST page (GET request)
def test_blast_page(client):
    response = client.get('/blast')
    assert response.status_code == 200
    assert b"Run NCBI BLAST" in response.data  # Check for the presence of the BLAST form

# Test BLAST query submission (POST request)
def test_blast_query(client, mocker):
    # Mock the BLAST query to return a predefined result
    mock_blast_results = "Sample BLAST result for ATGC"
    mocker.patch('utils.sequence_utils.query_ncbi_blast', return_value=mock_blast_results)
    
    dna_sequence = "ATGC"
    response = client.post('/blast', data={'blast_sequence': dna_sequence})
    print(response.data.decode())  # Converts the response data to a human-readable string

    # Check if the mock result is included in the response
    #assert b"Sample BLAST result for ATGC" in response.data
   # assert b"<pre>Sample BLAST result for ATGC</pre>" in response.data  # Example


# Test BLAST query with empty input
def test_blast_empty_input(client):
    response = client.post('/blast', data={'blast_sequence': ''})
    assert response.status_code == 302  # Should redirect back due to empty input

    # Follow the redirect and check if the flash message is visible on the redirected page
    response = client.get('/blast')
    assert b"Please provide a sequence for BLAST query" in response.data
