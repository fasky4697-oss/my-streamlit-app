import streamlit as st
from Bio import Entrez, SeqIO
import re

# ------------------------------
# ตั้งค่า Entrez
# ------------------------------
Entrez.email = "fasky4697@gmail.com"  # <— แก้เป็นอีเมลของคุณ

# ------------------------------
# ตารางแปลรหัส / โค้ดอน
# ------------------------------
codon_table = {
    "ATA":"I","ATC":"I","ATT":"I","ATG":"M",
    "ACA":"T","ACC":"T","ACG":"T","ACT":"T",
    "AAC":"N","AAT":"N","AAA":"K","AAG":"K",
    "AGC":"S","AGT":"S","AGA":"R","AGG":"R",
    "CTA":"L","CTC":"L","CTG":"L","CTT":"L",
    "CCA":"P","CCC":"P","CCG":"P","CCT":"P",
    "CAC":"H","CAT":"H","CAA":"Q","CAG":"Q",
    "CGA":"R","CGC":"R","CGG":"R","CGT":"R",
    "GTA":"V","GTC":"V","GTG":"V","GTT":"V",
    "GCA":"A","GCC":"A","GCG":"A","GCT":"A",
    "GAC":"D","GAT":"D","GAA":"E","GAG":"E",
    "GGA":"G","GGC":"G","GGG":"G","GGT":"G",
    "TCA":"S","TCC":"S","TCG":"S","TCT":"S",
    "TTC":"F","TTT":"F","TTA":"L","TTG":"L",
    "TAC":"Y","TAT":"Y","TAA":"*","TAG":"*",
    "TGC":"C","TGT":"C","TGA":"*","TGG":"W"
}
start_codons = {"ATG"}
stop_codons = {"TAA", "TAG", "TGA"}

# ------------------------------
# state (เก็บข้อมูลใช้งานร่วมกันระหว่างปุ่ม)
# ------------------------------
state = {
    "genome_seq": None,
    "genbank_record": None,
    "genome_id": None,
    "current_frame": None,
    "last_protein_header": None, # เพิ่มสำหรับเก็บ header โปรตีนล่าสุด
    "last_protein_seq": None,    # เพิ่มสำหรับเก็บลำดับโปรตีนล่าสุด
}

# ------------------------------
# Utilities
# ------------------------------
def translate_dna(seq):
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        protein.append(codon_table.get(codon, "?"))
    return "".join(protein)

def highlight_codons(seq):
    """ไฮไลท์ ATG (เขียว) และ stop codons (แดง) โดยแบ่งเป็น triplet"""
    html = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            html.append(codon)
            break
        if codon in start_codons:
            html.append(f"<span style='color:green;font-weight:bold'>{codon}</span>")
        elif codon in stop_codons:
            html.append(f"<span style='color:red;font-weight:bold'>{codon}</span>")
        else:
            html.append(codon)
    return "".join(html)

def find_codons_positions(seq, codon_set, start=0):
    """คืนรายการตำแหน่ง absolute และ frame สำหรับ codon ในลำดับ seq"""
    positions = []
    n = len(seq)
    for i in range(0, n - 2):
        codon = seq[i:i+3]
        if codon in codon_set:
            pos = start + i           # absolute (0-based)
            frame = pos % 3
            positions.append((pos, frame, codon))
    return positions

def fetch_genome_and_gb(bacteria_name):
    term = f"{bacteria_name}[orgn] AND complete genome[title]"
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        return None, None, None

    genome_id = ids[0]

    # FASTA (ซีเควนซ์)
    fh = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
    fasta_rec = SeqIO.read(fh, "fasta")
    fh.close()
    genome_seq = str(fasta_rec.seq).upper()

    # GenBank (annotation)
    try:
        gbh = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
        gb_rec = SeqIO.read(gbh, "genbank")
        gbh.close()
    except Exception:
        gb_rec = None

    return genome_id, genome_seq, gb_rec

def get_overlapping_cds_annotations(gb_record, start, end):
    results = []
    if gb_record is None:
        return results
    for feat in gb_record.features:
        if feat.type != "CDS":
            continue
        f_start = int(feat.location.start)
        f_end = int(feat.location.end)
        if not (end <= f_start or start >= f_end):  # overlap
            q = feat.qualifiers
            product = "; ".join(q.get("product", [])) if "product" in q else ""
            gene = "; ".join(q.get("gene", [])) if "gene" in q else ""
            pid = "; ".join(q.get("protein_id", [])) if "protein_id" in q else ""
            note = "; ".join(q.get("note", [])) if "note" in q else ""
            results.append({
                "start": f_start, "end": f_end, "strand": feat.location.strand,
                "gene": gene, "product": product, "protein_id": pid, "note": note
            })
    return results

# ------------------------------
# Streamlit UI Components
# ------------------------------
st.title("Genome Sequence Translation")

species_name = st.text_input("Enter the bacteria species name (e.g., Escherichia coli):", "Escherichia coli")

if st.button("Fetch Genome"):
    st.text("Fetching genome...")
    genome_id, genome_seq, genbank_record = fetch_genome_and_gb(species_name)
    
    if not genome_seq:
        st.error("No genome found for this species.")
    else:
        state["genome_id"] = genome_id
        state["genome_seq"] = genome_seq
        state["genbank_record"] = genbank_record
        st.success(f"Genome fetched successfully! Genome length: {len(genome_seq)} bases.")
        
        # Display first 120 bases of genome with codon highlights
        st.markdown(f"**Sample Sequence (First 120 bases):**")
        st.text(highlight_codons(genome_seq[:120]))
        
        start_pos, end_pos = st.slider("Select a range of bases to analyze", min_value=1, max_value=len(genome_seq), value=(1, 3000))
        
        selected_range = genome_seq[start_pos-1:end_pos]
        st.text(f"Selected Range: {selected_range[:120]}... (Length: {end_pos - start_pos + 1} bases)")
        
        # Find start and stop codons in the selected range
        starts = find_codons_positions(selected_range, start_codons, start=start_pos-1)
        stops = find_codons_positions(selected_range, stop_codons, start=start_pos-1)
        
        st.markdown(f"**Start Codons:** {starts}")
        st.markdown(f"**Stop Codons:** {stops}")
        
        # Translate the selected range to protein sequence
        translated_seq = translate_dna(selected_range)
        st.markdown(f"**Translated Protein Sequence:** {translated_seq[:400]}... (Length: {len(translated_seq)} aa)")
        
        # Annotate protein if available
        annotations = get_overlapping_cds_annotations(genbank_record, start_pos-1, end_pos)
        if annotations:
            st.markdown("**Annotations from GenBank (Overlapping CDS):**")
            for ann in annotations[:3]:
                st.markdown(f"- Gene: {ann['gene']} | Product: {ann['product']} | Protein ID: {ann['protein_id']} | Position: {ann['start']} - {ann['end']}")
        else:
            st.markdown("No overlapping annotations found.")
        
        st.download_button("Download Protein Sequence", data=translated_seq, file_name="protein_sequence.faa", mime="text/fasta")

        st.download_button("Download Genome Sequence", data=genome_seq, file_name=f"{genome_id}.fna", mime="text/fasta")

