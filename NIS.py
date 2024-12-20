from flask import Flask, flash, render_template, redirect, url_for, request, send_file, session
import plotly.express as px
import io
import requests
import time
from utils.sequence_utils import (
    transcribe_dna_to_rna, translate_rna_to_protein, calculate_gc_content,
    generate_fasta, generate_csv, get_nucleotide_distribution
)
from io import StringIO
from Bio.Blast import NCBIXML
from flask_session import Session

ALLOWED_EXTENSIONS = {"txt", "fasta"}

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Required for session management
# Configure server-side session
app.config['SESSION_TYPE'] = 'filesystem'  # Or 'redis', 'mongodb', etc.
app.config['SESSION_PERMANENT'] = False
app.config['SESSION_USE_SIGNER'] = True
Session(app)

# Helper functions
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def highlight_sequence(sequence):
    color_map = {'A': 'red', 'T': 'blue', 'G': 'green', 'C': 'orange'}
    return "".join(
        f'<span style="color:{color_map.get(nucleotide, "black")}">{nucleotide}</span>'
        for nucleotide in sequence
    )

def generate_nucleotide_graph(sequence):
    counts = {nucleotide: sequence.count(nucleotide) for nucleotide in "ATGC"}
    fig = px.bar(x=list(counts.keys()), y=list(counts.values()),
                 title="Nucleotide Distribution",
                 labels={"x": "Nucleotide", "y": "Count"})
    return fig.to_html(full_html=False)

def search_motif(dna_sequence, motif):
    positions = []
    start = 0
    while start < len(dna_sequence):
        start = dna_sequence.find(motif, start)
        if start == -1:
            break
        positions.append(start)
        start += len(motif)
    return positions

def perform_blast(sequence):
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "nt",
        "QUERY": sequence,
        "FORMAT_TYPE": "XML",
        "MAX_TARGET_SEQUENCES": 100,
        "EXPECT": 10,
        "WORD_SIZE": 11,
        "GAPCOSTS": "5 2",
        "FILTER": "L",
    }

    try:
        response = requests.post(url, data=params)
        if response.status_code != 200:
            return f"Error: Failed to submit query. HTTP Status: {response.status_code}"

        rid = response.text.split("RID = ")[1].split("\n")[0].strip()

        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML",
        }

        while True:
            result = requests.get(url, params=params)

            if "Status=WAITING" in result.text:
                time.sleep(5)
            elif result.text.startswith('<?xml'):
                return result.text
            else:
                return f"Error: Unexpected response: {result.text[:100]}"

    except Exception as e:
        return f"Error occurred while running BLAST: {str(e)}"

def parse_blast_results(blast_output):
    blast_output_io = StringIO(blast_output)
    blast_records = NCBIXML.parse(blast_output_io)
    results = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                alignment_data = f"<pre>{hsp.query} \n{hsp.match} \n{hsp.sbjct}</pre>"

                result = {
                    "description": alignment.title,
                    "max_score": hsp.score,
                    "total_score": hsp.bits,
                    "query_coverage": hsp.align_length / len(blast_record.query) * 100,
                    "e_value": hsp.expect,
                    "identities": hsp.identities,
                    "gaps": hsp.gaps,
                    "alignment": alignment_data,
                    "percent_identical": (hsp.identities / hsp.align_length) * 100,
                    "accession": alignment.accession
                }
                results.append(result)
    return results

# Routes
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        dna_file = request.files.get("dna_file")
        if dna_file and allowed_file(dna_file.filename):
            dna_content = dna_file.read().decode("utf-8").strip()
            if dna_content.startswith(">"):
                dna_sequence = "".join(
                    line.strip() for line in dna_content.splitlines() if not line.startswith(">")
                )
            else:
                dna_sequence = dna_content
            session['dna_sequence'] = dna_sequence
            flash("DNA sequence successfully uploaded.", "success")
        #else:
           # flash("Please upload a valid FASTA or TXT file.", "error")

    dna_sequence = session.get('dna_sequence', '')
    if not dna_sequence:
   #     flash("No DNA sequence found. Please upload a valid file.", "info")
        return render_template("index.html")

    rna_sequence = transcribe_dna_to_rna(dna_sequence)
    protein_sequence = translate_rna_to_protein(rna_sequence)
    gc_content = calculate_gc_content(dna_sequence)
    nucleotide_distribution = get_nucleotide_distribution(dna_sequence)
    highlighted_sequence = highlight_sequence(dna_sequence)
    graph_html = generate_nucleotide_graph(dna_sequence)

    if 'export_fasta' in request.form:
        fasta_data = generate_fasta(dna_sequence)
        return send_file(io.BytesIO(fasta_data.encode()), as_attachment=True,
                         download_name="sequence.fasta", mimetype="text/plain")

    if 'export_csv' in request.form:
        csv_data = generate_csv(dna_sequence, rna_sequence, protein_sequence)
        return send_file(io.BytesIO(csv_data.encode()), as_attachment=True,
                         download_name="sequence_data.csv", mimetype="text/csv")

    return render_template("index.html", dna_sequence=dna_sequence,
                           highlighted_sequence=highlighted_sequence,
                           rna_sequence=rna_sequence, protein_sequence=protein_sequence,
                           gc_content=gc_content, nucleotide_distribution=nucleotide_distribution,
                           graph_html=graph_html)

@app.route('/blast', methods=['GET', 'POST'])
def blast():
    blast_results = None

    if request.method == 'POST':
        sequence = request.form.get('blast_sequence', '').strip()

        if not sequence:
            flash("DNA sequence is required to perform BLAST.", "danger")
            return render_template('blast.html', blast_results=blast_results)

        try:
            blast_output = perform_blast(sequence)

            if not blast_output:
                flash("No results returned from BLAST query. Please try again.", "warning")
            else:
                blast_results = parse_blast_results(blast_output)

        except ValueError as ve:
            flash(f"Invalid DNA sequence: {str(ve)}", "danger")
        except Exception as e:
            flash(f"An error occurred while performing the BLAST operation: {str(e)}", "danger")

    return render_template('blast.html', blast_results=blast_results)

@app.route('/gc-content', methods=['GET'])
def gc_content():
    dna_sequence = session.get('dna_sequence', '')
    if not dna_sequence:
        flash("No DNA sequence found. Please upload one.", "info")
        return redirect(url_for('index'))

    gc_content = calculate_gc_content(dna_sequence)
    at_content = 100 - gc_content

    return render_template('gc_content.html', gc_content=gc_content, at_content=at_content)

@app.route('/home')
def home():
    return "This is the home page"

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')


@app.route('/search', methods=['GET', 'POST'])
def search():
    return "Search functionality is under construction."

@app.route('/privacy-policy', methods=['GET'])
def privacy_policy():
    return render_template('privacy_policy.html')

if __name__ == '__main__':
    app.run(debug=True)
