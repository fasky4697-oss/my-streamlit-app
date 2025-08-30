# app.py - This file will contain the Streamlit application code

from Bio import Entrez, SeqIO
import streamlit as st
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
# Utilities - Adapted for Streamlit
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
            pos = start + i      # absolute (0-based)
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
    try:
        fh = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        fasta_rec = SeqIO.read(fh, "fasta")
        fh.close()
        genome_seq = str(fasta_rec.seq).upper()
    except Exception:
        genome_seq = None

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

# Keep file saving functions separate to be called by Streamlit UI
def save_protein_to_file(filename, header, protein_seq):
    if header and protein_seq:
        fasta_content = to_fasta(header, protein_seq)
        return fasta_content
    return None

def save_genome_to_file(filename, header, genome_seq):
    if header and genome_seq:
        fasta_content = to_fasta(header, genome_seq)
        return fasta_content
    return None

# Adapt analyze_current_range
def analyze_current_range(gseq, start_pos, end_pos):
    s0, e0 = start_pos - 1, end_pos  # Convert to 0-based [s0, e0)
    selected = gseq[s0:e0]

    # Find start/stop in the range
    starts = find_codons_positions(selected, start_codons, start=s0)
    stops  = find_codons_positions(selected, stop_codons, start=s0)

    return selected, starts, stops

# Adapt annotate_and_show to return content for Streamlit
def annotate_and_show_streamlit(s0, e0, aa, gb_record, label_prefix="ช่วงที่ใช้"):
    header = f"range_{s0+1}_{e0}_frame{(s0%3)}"

    ann = get_overlapping_cds_annotations(gb_record, s0, e0)
    hints = heuristic_function_prediction(aa)

    output_parts = []

    output_parts.append(
        f"<div><b>{label_prefix}:</b> {s0+1:,}–{e0:,} | frame {(s0%3)} | "
        f"ความยาวโปรตีน: {len(aa)} aa</div>"
    )
    output_parts.append(f"<pre style='font-size:13px'>{aa[:400]}{'...' if len(aa)>400 else ''}</pre>")


    if ann:
        ann_html = "<b>Annotation (GenBank CDS ที่ทับซ้อน):</b><ul>"
        for a in ann[:3]:
            ann_html += (
                f"<li><b>{a.get('gene','').strip() or '(unknown gene)'}:</b> "
                f"{a.get('product','').strip() or '(unknown product)'} "
                f"[protein_id: {a.get('protein_id','')}] "
                f"(ตำแหน่ง {a['start']+1}–{a['end']}, strand: {a['strand']})</li>"
            )
        ann_html += "</ul>"
        output_parts.append(ann_html)
    else:
        hints_html = "<b>คาดเดาฟังก์ชัน (heuristic):</b><ul>" + "".join(f"<li>{h}</li>" for h in hints) + "</ul>"
        output_parts.append(hints_html)

    return "".join(output_parts), header, aa


# ------------------------------
# Streamlit UI
# ------------------------------

st.title("Genome Sequence Translation")

# Initialize session state
if 'genome_seq' not in st.session_state:
    st.session_state.genome_seq = None
if 'genbank_record' not in st.session_state:
    st.session_state.genbank_record = None
if 'genome_id' not in st.session_state:
    st.session_state.genome_id = None
if 'last_protein_header' not in st.session_state:
    st.session_state.last_protein_header = None
if 'last_protein_seq' not in st.session_state:
    st.session_state.last_protein_seq = None
if 'current_range_start' not in st.session_state:
    st.session_state.current_range_start = 1
if 'current_range_end' not in st.session_state:
    st.session_state.current_range_end = 3000
if 'start_codon_options' not in st.session_state:
    st.session_state.start_codon_options = [("— ไม่เลือก —", None)]
if 'stop_codon_options' not in st.session_state:
    st.session_state.stop_codon_options = [("— ไม่เลือก —", None)]


# Input: Bacteria name
col_input, col_dropdown = st.columns([2, 1])
with col_input:
    custom_species_name = st.text_input(
        'กรอกชื่อแบคทีเรียเอง:',
        value="Escherichia coli",
        key='custom_species_name_input'
    )
with col_dropdown:
    st.markdown("<br>", unsafe_allow_html=True)  # Add some vertical space for alignment
    species_name_from_dropdown = st.selectbox(
        'หรือเลือกจากรายการ:',
        options=[
            "",
            "Escherichia coli", "Staphylococcus aureus", "Salmonella enterica",
            "Bacillus subtilis", "Pseudomonas aeruginosa", "Mycobacterium tuberculosis",
            "Vibrio cholerae", "Klebsiella pneumoniae", "Acinetobacter baumannii"
        ],
        index=0,
        key='species_name_dropdown'
    )

# Use the value from the text input unless a dropdown value is selected
species_name = species_name_from_dropdown if species_name_from_dropdown else custom_species_name

col1, col2 = st.columns(2)

with col1:
    fetch_button = st.button("ดึงจีโนม", key='fetch_button')
with col2:
    reset_button = st.button("ล้างผล", key='reset_button')

if reset_button:
    st.session_state.genome_seq = None
    st.session_state.genbank_record = None
    st.session_state.genome_id = None
    st.session_state.last_protein_header = None
    st.session_state.last_protein_seq = None
    st.session_state.current_range_start = 1
    st.session_state.current_range_end = 3000
    st.session_state.start_codon_options = [("— ไม่เลือก —", None)]
    st.session_state.stop_codon_options = [("— ไม่เลือก —", None)]
    st.rerun() # Rerun to clear the UI

if fetch_button and species_name:
    with st.spinner(f"กำลังดึงจีโนมของ {species_name} จาก NCBI..."):
        genome_id, genome_seq, genbank_record = fetch_genome_and_gb(species_name)
        if genome_seq:
            st.session_state.genome_id = genome_id
            st.session_state.genome_seq = genome_seq
            st.session_state.genbank_record = genbank_record
            st.session_state.current_range_start = 1
            st.session_state.current_range_end = min(3000, len(genome_seq))
            st.success(f"✅ พบจีโนมความยาว {len(genome_seq):,} เบส | NCBI ID: {genome_id}")
            preview = genome_seq[:120]
            st.markdown("<b>ตัวอย่าง 120 nt แรก:</b>", unsafe_allow_html=True)
            st.markdown(f"<pre style='font-size:13px'>{highlight_codons(preview)}</pre>", unsafe_allow_html=True)
        else:
            st.error("❌ ไม่พบจีโนมสำหรับคำค้นนี้")
            st.session_state.genome_seq = None
            st.session_state.genbank_record = None
            st.session_state.genome_id = None

# Range slider and analysis
if st.session_state.genome_seq:
    genome_length = len(st.session_state.genome_seq)
    st.session_state.current_range_start, st.session_state.current_range_end = st.slider(
        'ช่วงเบส:',
        min_value=1,
        max_value=genome_length,
        value=(st.session_state.current_range_start, st.session_state.current_range_end),
        step=1,
        key='range_slider'
    )

    selected_seq, starts_in_range, stops_in_range = analyze_current_range(
        st.session_state.genome_seq,
        st.session_state.current_range_start,
        st.session_state.current_range_end
    )

    st.markdown(f"<b>ช่วงที่เลือก:</b> {st.session_state.current_range_start:,}–{st.session_state.current_range_end:,} (ยาว {len(selected_seq):,} nt)", unsafe_allow_html=True)
    st.markdown(f"<pre style='font-size:13px'>{highlight_codons(selected_seq)}</pre>", unsafe_allow_html=True)


    # Update start/stop dropdown options
    s_opts = [("— ไม่เลือก —", None)] + [
        (f"ATG @ {p+1} (frame {f})", (p, f)) for (p, f, c) in starts_in_range[:50]
    ]
    t_opts = [("— ไม่เลือก —", None)] + [
        (f"{c} @ {p+1} (frame {f})", (p, f, c)) for (p, f, c) in stops_in_range[:50]
    ]

    # Add nearby codons if the range is sparse
    msg = []
    if not starts_in_range:
        all_starts = find_codons_positions(st.session_state.genome_seq, start_codons, start=0)
        near_s = nearest_positions(st.session_state.current_range_start - 1, all_starts, 10)
        s_opts.extend([(f"แนะนำ: ATG @ {p+1} (frame {f})", (p, f)) for (p, f, c) in near_s])
        if near_s: # Only add message if suggestions were added
            msg.append("⚠️ ไม่พบ start codon ในช่วงที่เลือก — มีตำแหน่งใกล้เคียงให้เลือกในเมนู")


    if not stops_in_range:
        all_stops  = find_codons_positions(st.session_state.genome_seq, stop_codons, start=0)
        near_t = nearest_positions(st.session_state.current_range_end, all_stops, 10)
        t_opts.extend([(f"แนะนำ: {c} @ {p+1} (frame {f})", (p, f, c)) for (p, f, c) in near_t])
        if near_t: # Only add message if suggestions were added
            msg.append("⚠️ ไม่พบ stop codon ในช่วงที่เลือก — มีตำแหน่งใกล้เคียงให้เลือกในเมนู")


    if msg:
        st.warning("\n".join(msg))


    st.session_state.start_codon_options = s_opts
    st.session_state.stop_codon_options = t_opts

    col3, col4, col5 = st.columns(3)

    with col3:
        start_selection = st.selectbox(
            "เลือก Start:",
            options=st.session_state.start_codon_options,
            format_func=lambda x: x[0], # Display the text part of the tuple
            key='start_dropdown'
        )
    with col4:
        stop_selection = st.selectbox(
            "เลือก Stop:",
            options=st.session_state.stop_codon_options,
            format_func=lambda x: x[0], # Display the text part of the tuple
            key='stop_dropdown'
        )
    with col5:
        # Add a small space or use markdown for alignment
        st.markdown("<br>", unsafe_allow_html=True) # Add some vertical space
        use_picks_button = st.button("แปลจาก Start/Stop ที่เลือก", key='use_picks_button')


    col6, col7 = st.columns(2)
    with col6:
        translate_selected_button = st.button("แปลจากช่วงที่เลือก", key='translate_selected_button')
    with col7:
        if st.session_state.last_protein_seq:
            fasta_content = save_protein_to_file(st.session_state.last_protein_header + ".faa", st.session_state.last_protein_header, st.session_state.last_protein_seq)
            st.download_button(
                label="บันทึกโปรตีน (FASTA)",
                data=fasta_content,
                file_name=st.session_state.last_protein_header + ".faa",
                mime="text/plain",
                key='save_protein_button'
            )
        else:
            st.button("บันทึกโปรตีน (FASTA)", key='save_protein_button_disabled', disabled=True)


    if translate_selected_button:
        s0 = st.session_state.current_range_start - 1
        e0 = st.session_state.current_range_end
        # Adjust to be divisible by 3
        adj_len = ((e0 - s0) // 3) * 3
        e0_adj = s0 + adj_len
        if adj_len <= 0:
            st.warning("ช่วงที่เลือกสั้นเกินไปสำหรับการแปลรหัส")
        else:
            aa = translate_dna(st.session_state.genome_seq[s0:e0_adj])
            output_html, header, protein_seq = annotate_and_show_streamlit(s0, e0_adj, aa, st.session_state.genbank_record, label_prefix="แปลตามช่วงที่เลือก")
            st.session_state.last_protein_header = header
            st.session_state.last_protein_seq = protein_seq
            st.markdown(output_html, unsafe_allow_html=True)

    if use_picks_button:
        s_pick = start_selection # (pos, frame) or None
        t_pick = stop_selection  # (pos, frame, codon) or None

        if not s_pick and not t_pick:
            st.warning("โปรดเลือกอย่างน้อย Start หรือ Stop")
        else:
            s0, e0 = None, None
            label = "ช่วงที่เลือก"

            # Both selected
            if s_pick and t_pick:
                s_pos, s_frame = s_pick
                t_pos, t_frame, t_cod = t_pick
                if s_frame != t_frame or t_pos <= s_pos:
                    st.error("Start/Stop ไม่สอดคล้องกัน (frame/ลำดับ)")
                else:
                    s0, e0 = s_pos, t_pos + 3
                    label = "ช่วง Start–Stop ที่เลือก"

            # Only Start selected
            elif s_pick and not t_pick:
                s_pos, s_frame = s_pick
                t_pos = None
                for i in range(s_pos, len(st.session_state.genome_seq)-2, 3):
                    if (i % 3) == s_frame and st.session_state.genome_seq[i:i+3] in stop_codons:
                        t_pos = i
                        break
                if t_pos is None:
                    # No stop found - use end of range slider, adjusted for frame
                    s0 = s_pos
                    e0 = st.session_state.current_range_end
                    e0 = s0 + ((e0 - s0)//3)*3 # Adjust length
                    st.warning("ไม่พบ Stop codon หลัง Start ที่เลือก ใช้ขอบช่วงที่เลือกแทน")
                else:
                    s0, e0 = s_pos, t_pos + 3
                label = "ช่วงจาก Start → Stop/ขอบช่วง"

            # Only Stop selected
            elif (not s_pick) and t_pick:
                t_pos, t_frame, t_cod = t_pick
                s_pos = None
                for i in range(t_pos, -1, -3):
                    if (i % 3) == t_frame and st.session_state.genome_seq[i:i+3] in start_codons:
                        s_pos = i
                        break
                if s_pos is None:
                    # No start found - use start of range, adjusted for frame
                    s0 = st.session_state.current_range_start - 1
                    while (s0 % 3) != t_frame and s0 < t_pos:
                        s0 += 1
                    if s0 >= t_pos:
                        st.error("ไม่พบ Start codon ก่อน Stop ที่เลือกใน Frame เดียวกัน")
                        s0 = None # Indicate failure
                    else:
                        st.warning("ไม่พบ Start codon ก่อน Stop ที่เลือก ใช้ขอบช่วงที่เลือกแทน")
                else:
                    s0 = s_pos
                e0 = t_pos + 3
                label = "ช่วงจากขอบช่วง/Start → Stop"

            if s0 is not None and e0 is not None and e0 > s0:
                aa = translate_dna(st.session_state.genome_seq[s0:e0])
                if len(aa) == 0:
                    st.warning("ช่วงที่เลือกสั้นเกินไปสำหรับการแปลรหัส")
                else:
                    output_html, header, protein_seq = annotate_and_show_streamlit(s0, e0, aa, st.session_state.genbank_record, label_prefix=label)
                    st.session_state.last_protein_header = header
                    st.session_state.last_protein_seq = protein_seq
                    st.markdown(output_html, unsafe_allow_html=True)
            elif s0 is not None and e0 is not None and e0 <= s0:
                st.error("ตำแหน่งสิ้นสุดต้องมากกว่าตำแหน่งเริ่มต้น")

    if st.session_state.genome_seq and st.session_state.genome_id:
        fasta_genome_content = save_genome_to_file(st.session_state.genome_id + ".fna", st.session_state.genome_id, st.session_state.genome_seq)
        st.download_button(
            label="บันทึกจีโนม (FASTA)",
            data=fasta_genome_content,
            file_name=st.session_state.genome_id + ".fna",
            mime="text/plain",
            key='save_genome_button'
        )
    else:
        st.button("บันทึกจีโนม (FASTA)", key='save_genome_button_disabled', disabled=True)
