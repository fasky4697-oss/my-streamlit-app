import streamlit as st
from Bio import Entrez, SeqIO
import re
import os

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
    "last_protein_header": None,  # เพิ่มสำหรับเก็บ header โปรตีนล่าสุด
    "last_protein_seq": None,  # เพิ่มสำหรับเก็บลำดับโปรตีนล่าสุด
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
            pos = start + i  # absolute (0-based)
            frame = pos % 3
            positions.append((pos, frame, codon))
    return positions

def nearest_positions(target, positions, limit=10):
    return sorted(positions, key=lambda x: abs(x[0] - target))[:limit]

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

def heuristic_function_prediction(aa):
    """ทายฟังก์ชันอย่างหยาบจาก motif/คุณสมบัติ"""
    hints = []
    if not aa:
        return ["(ไม่มีกรดอะมิโน)"]

    # P-loop NTPase (Walker A): GxxxxGKS/T
    if re.search(r"G....GK[ST]", aa):
        hints.append("อาจเป็น P-loop NTPase (Walker A) เช่น ATPase/GTPase")

    # Metalloprotease HExH
    if re.search(r"H.EH", aa):
        hints.append("อาจเป็น metalloprotease (HExH motif)")

    # Serine hydrolase (GxSxG)
    if re.search(r"G.S.G", aa):
        hints.append("อาจเป็นเซรีนไฮดรอลาเซ (GxSxG motif)")

    # Transmembrane helix (ช่วง hydrophobic ยาว)
    hydrophobic = set("AILMFWVY")
    run = 0; max_run = 0
    for ch in aa:
        if ch in hydrophobic:
            run += 1; max_run = max(max_run, run)
        else:
            run = 0
    if max_run >= 18:
        hints.append("อาจเป็นโปรตีนเยื่อหุ้ม (มี transmembrane helix)")

    # Signal peptide (หยาบๆ)
    if len(aa) >= 25 and sum(1 for c in aa[:25] if c in hydrophobic) >= 10:
        hints.append("อาจมีสัญญาณหลั่ง/นำส่ง (signal peptide)")

    if len(aa) < 50:
        hints.append("สั้นมาก (อาจเป็น ORF สั้น/เปปไทด์เล็ก)")

    if not hints:
        hints.append("ฟังก์ชันไม่ชัด ควรทำ BLASTp/HMMER เพื่อเทียบฐานข้อมูล")
    return hints

def to_fasta(header, seq, width=70):
    lines = [f">{header}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i+width])
    return "\n".join(lines)

# ------------------------------
# Streamlit UI Components
# ------------------------------
st.title('Genome & Protein Translation Tool')

species_name = st.text_input('กรอกชื่อแบคทีเรีย (เช่น Escherichia coli)', value="
