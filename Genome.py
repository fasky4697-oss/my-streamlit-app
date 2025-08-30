

# @title #ดึง genome ของ Bacteria จากNCBI แปลงเป็นโปรตีนได้ เลือกช่วงที่จะแปลงได้

from Bio import Entrez, SeqIO
from IPython.display import display, HTML
import ipywidgets as widgets
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

def save_protein_to_file(filename, header, protein_seq):
    if header and protein_seq:
        with open(filename, "w") as f:
            f.write(to_fasta(header, protein_seq))
        print(f"✅ บันทึกผลเป็นไฟล์ {filename}")
    else:
        print("❌ ไม่มีข้อมูลโปรตีนที่จะบันทึก")

def save_genome_to_file(filename, header, genome_seq):
     if header and genome_seq:
        with open(filename, "w") as f:
            f.write(to_fasta(header, genome_seq))
        print(f"✅ บันทึกจีโนมเป็นไฟล์ {filename}")
     else:
        print("❌ ไม่มีข้อมูลจีโนมที่จะบันทึก")


# ------------------------------
# Widgets (กล่องตอบโต้ + ปุ่ม)
# ------------------------------
species_box = widgets.Combobox(
    placeholder='กรอกชื่อแบคทีเรีย เช่น Escherichia coli',
    options=[
        "Escherichia coli", "Staphylococcus aureus", "Salmonella enterica",
        "Bacillus subtilis", "Pseudomonas aeruginosa", "Mycobacterium tuberculosis",
        "Vibrio cholerae", "Klebsiella pneumoniae", "Acinetobacter baumannii"
    ],
    description='ชื่อแบค:',
    ensure_option=False,
    value='Escherichia coli',
    layout=widgets.Layout(width='450px')
)

fetch_button = widgets.Button(description="ดึงจีโนม", button_style="info")
analyze_button = widgets.Button(description="วิเคราะห์ช่วง", button_style="")
translate_selected_button = widgets.Button(description="แปลจากช่วงเลือก", button_style="warning")
use_picks_button = widgets.Button(description="แปลจาก Start/Stop ที่เลือก", button_style="success")
save_protein_button = widgets.Button(description="บันทึกโปรตีน", button_style="primary", disabled=True) # เพิ่มปุ่มบันทึก
save_genome_button = widgets.Button(description="บันทึกลำดับเบส (FASTA)", button_style="primary", disabled=True) # เปลี่ยนชื่อปุ่มและ disabled เริ่มต้น
reset_button = widgets.Button(description="ล้างผล", button_style="danger")

info_html = widgets.HTML("")
output_html = widgets.HTML("")

range_slider = widgets.IntRangeSlider(
    value=[1, 3000],
    min=1, max=3000, step=1,
    description='ช่วงเบส:',
    continuous_update=False, readout=True, layout=widgets.Layout(width='95%')
)

start_dropdown = widgets.Dropdown(options=[("— ไม่เลือก —", None)], description="เลือก Start:")
stop_dropdown = widgets.Dropdown(options=[("— ไม่เลือก —", None)], description="เลือก Stop:")

# ------------------------------
# Event handlers
# ------------------------------
def on_fetch_clicked(b):
    name = (species_box.value or "").strip()
    if not name:
        info_html.value = "<span style='color:red'>❌ โปรดกรอกชื่อแบคทีเรีย</span>"
        return
    info_html.value = "<em>⏳ กำลังดึงจีโนมจาก NCBI...</em>"
    output_html.value = ""
    gid, gseq, grecord = fetch_genome_and_gb(name)
    if not gseq:
        info_html.value = "<span style='color:red'>❌ ไม่พบจีโนมสำหรับคำค้นนี้</span>"
        return
    state["genome_id"] = gid
    state["genome_seq"] = gseq
    state["genbank_record"] = grecord
    state["last_protein_header"] = None # ล้างข้อมูลโปรตีนเมื่อดึงจีโนมใหม่
    state["last_protein_seq"] = None
    save_protein_button.disabled = True
    save_genome_button.disabled = False # เปิดใช้งานปุ่มบันทึกจีโนม

    n = len(gseq)
    range_slider.min = 1
    range_slider.max = n
    range_slider.value = (1, min(3000, n))

    preview = gseq[:120]
    info_html.value = f"✅ พบจีโนมความยาว {n:,} เบส | NCBI ID: {gid}"
    output_html.value = (
        "<b>ตัวอย่าง 120 nt แรก:</b>"
        f"<pre style='font-size:13px'>{highlight_codons(preview)}</pre>"
    )

    # รีเซ็ตเมนูเลือก
    start_dropdown.options = [("— ไม่เลือก —", None)]
    stop_dropdown.options = [("— ไม่เลือก —", None)]

def analyze_current_range():
    gseq = state["genome_seq"]
    if not gseq:
        output_html.value = "<span style='color:red'>❌ ยังไม่ได้ดึงจีโนม</span>"
        return None

    s, e = range_slider.value
    s0, e0 = s-1, e              # แปลงเป็น 0-based [s0, e0)
    selected = gseq[s0:e0]

    # หา start/stop ในช่วง
    starts = find_codons_positions(selected, start_codons, start=s0)
    stops  = find_codons_positions(selected, stop_codons, start=s0)

    # เติม dropdown (แสดงสูงสุด 50 รายการเพื่อความลื่น)
    s_opts = [("— ไม่เลือก —", None)] + [
        (f"ATG @ {p+1} (frame {f})", (p, f)) for (p, f, c) in starts[:50]
    ]
    t_opts = [("— ไม่เลือก —", None)] + [
        (f"{c} @ {p+1} (frame {f})", (p, f, c)) for (p, f, c) in stops[:50]
    ]
    start_dropdown.options = s_opts
    stop_dropdown.options = t_opts

    # หากไม่พบในช่วง ให้แนะนำตำแหน่งใกล้เคียง
    msg = []
    if not starts or not stops:
        all_starts = find_codons_positions(gseq, start_codons, start=0)
        all_stops  = find_codons_positions(gseq, stop_codons, start=0)
        near_s = nearest_positions(s0, all_starts, 10) if not starts else []
        near_t = nearest_positions(e0, all_stops, 10)  if not stops else []
        if near_s:
            start_dropdown.options = list(start_dropdown.options) + [
                (f"แนะนำ: ATG @ {p+1} (frame {f})", (p, f)) for (p, f, c) in near_s
            ]
            msg.append("⚠️ ไม่พบ start codon ในช่วงที่เลือก — มีตำแหน่งใกล้เคียงให้เลือกในเมนู")
        if near_t:
            stop_dropdown.options = list(stop_dropdown.options) + [
                (f"แนะนำ: {c} @ {p+1} (frame {f})", (p, f, c)) for (p, f, c) in near_t
            ]
            msg.append("⚠️ ไม่พบ stop codon ในช่วงที่เลือก — มีตำแหน่งใกล้เคียงให้เลือกในเมนู")

    html = (
        f"<b>ช่วงที่เลือก:</b> {s:,}–{e:,} (ยาว {e-s+1:,} nt)"
        f"<pre style='font-size:13px'>{highlight_codons(selected)}</pre>"
    )
    if msg:
        html += "<div>" + "<br>".join(msg) + "</div>"
    output_html.value = html
    # ไม่ต้อง reset โปรตีนล่าสุดเมื่อวิเคราะห์ช่วงเฉยๆ

    return (s0, e0)

def on_analyze_clicked(b):
    analyze_current_range()

def annotate_and_show(s0, e0, aa, label_prefix="ช่วงที่ใช้"):
    header = f"range_{s0+1}_{e0}_frame{(s0%3)}"
    # ไม่บันทึกไฟล์ตรงนี้แล้ว
    # filename = f"protein_{s0+1}_{e0}_frame{s0%3}.faa"
    # save_protein_to_file(filename, header, aa)

    # เก็บโปรตีนล่าสุดไว้ใน state
    state["last_protein_header"] = header
    state["last_protein_seq"] = aa
    save_protein_button.disabled = False # เปิดใช้งานปุ่มบันทึก

    ann = get_overlapping_cds_annotations(state["genbank_record"], s0, e0)
    hints = heuristic_function_prediction(aa)
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
    else:
        ann_html = "<b>คาดเดาฟังก์ชัน (heuristic):</b><ul>" + "".join(f"<li>{h}</li>" for h in hints) + "</ul>"

    output_html.value = (
        f"<div><b>{label_prefix}:</b> {s0+1:,}–{e0:,} | frame {(s0%3)} | "
        f"ความยาวโปรตีน: {len(aa)} aa</div>"
        f"<pre style='font-size:13px'>{aa[:400]}{'...' if len(aa)>400 else ''}</pre>"
        + ann_html
    )

def on_translate_selected_clicked(b):
    gseq = state["genome_seq"]
    if not gseq:
        output_html.value = "<span style='color:red'>❌ ยังไม่ได้ดึงจีโนม</span>"
        return
    rng = analyze_current_range()
    if rng is None:
        # analyze_current_range จะแสดง error message ถ้าไม่มี genome_seq
        return
    s0, e0 = rng
    # ปรับให้หาร 3 ลงตัว
    adj_len = ((e0 - s0) // 3) * 3
    e0_adj = s0 + adj_len
    aa = translate_dna(gseq[s0:e0_adj])
    annotate_and_show(s0, e0_adj, aa, label_prefix="แปลตามช่วงที่เลือก")

def on_use_picks_clicked(b):
    gseq = state["genome_seq"]
    if not gseq:
        output_html.value = "<span style='color:red'>❌ ยังไม่ได้ดึงจีโนม</span>"
        return

    s_pick = start_dropdown.value  # (pos, frame) or None
    t_pick = stop_dropdown.value   # (pos, frame, codon) or None

    if not s_pick and not t_pick:
        output_html.value = "<span style='color:red'>โปรดเลือกอย่างน้อย Start หรือ Stop</span>"
        return

    # เลือกทั้งคู่: ต้อง frame เดียว และ stop > start
    if s_pick and t_pick:
        s_pos, s_frame = s_pick
        t_pos, t_frame, t_cod = t_pick
        if s_frame != t_frame or t_pos <= s_pos:
            output_html.value = "<span style='color:red'>Start/Stop ไม่สอดคล้องกัน (frame/ลำดับ)</span>"
            return
        s0, e0 = s_pos, t_pos + 3
        aa = translate_dna(gseq[s0:e0])
        annotate_and_show(s0, e0, aa, label_prefix="ช่วง Start–Stop ที่เลือก")
        return

    # เลือก start อย่างเดียว: หา stop ถัดไปใน frame เดียวกัน
    if s_pick and not t_pick:
        s_pos, s_frame = s_pick
        t_pos = None
        for i in range(s_pos, len(gseq)-2, 3):
            if (i % 3) == s_frame and gseq[i:i+3] in stop_codons:
                t_pos = i
                break
        if t_pos is None:
            # ไม่พบ stop — ใช้ขอบช่วงจากสไลเดอร์
            s, e = range_slider.value
            s0 = s_pos
            e0 = e
            # ปรับความยาวหาร 3 ลงตัว
            e0 = s0 + ((e0 - s0)//3)*3
        else:
            s0, e0 = s_pos, t_pos + 3
        aa = translate_dna(gseq[s0:e0])
        annotate_and_show(s0, e0, aa, label_prefix="ช่วงจาก Start → Stop/ขอบช่วง")
        return

    # เลือก stop อย่างเดียว: ย้อนหา start ก่อนหน้าใน frame เดียวกัน
    if (not s_pick) and t_pick:
        t_pos, t_frame, t_cod = t_pick
        s_pos = None
        for i in range(t_pos, -1, -3):
            if (i % 3) == t_frame and gseq[i:i+3] in start_codons:
                s_pos = i
                break
        if s_pos is None:
            # ไม่พบ start — ใช้ขอบช่วง (ปรับให้เข้ากับ frame)
            s, e = range_slider.value
            s0 = s-1
            while (s0 % 3) != t_frame and s0 < t_pos:
                s0 += 1
        else:
            s0 = s_pos
        e0 = t_pos + 3
        aa = translate_dna(gseq[s0:e0])
        annotate_and_show(s0, e0, aa, label_prefix="ช่วงจากขอบช่วง/Start → Stop")
        return

def on_save_protein_clicked(b):
    header = state.get("last_protein_header")
    seq = state.get("last_protein_seq")
    if header and seq:
         filename = f"{header}.faa" # ตั้งชื่อไฟล์ตาม header
         save_protein_to_file(filename, header, seq)
    else:
        output_html.value = "<span style='color:red'>❌ ไม่มีข้อมูลโปรตีนที่แปลล่าสุดให้บันทึก</span>"

def on_save_genome_clicked(b): # เปลี่ยนชื่อฟังก์ชันกลับ
    gid = state.get("genome_id")
    gseq = state.get("genome_seq")
    if gid and gseq:
        filename = f"{gid}.fna" # ตั้งชื่อไฟล์ตาม genome ID
        save_genome_to_file(filename, gid, gseq) # เรียกใช้ฟังก์ชันบันทึกไฟล์
    else:
        output_html.value = "<span style='color:red'>❌ ยังไม่ได้ดึงข้อมูลจีโนม</span>"


def on_range_change(change):
    if change["name"] == "value" and state["genome_seq"]:
        analyze_current_range()

def on_reset_clicked(b):
    info_html.value = ""
    output_html.value = ""
    start_dropdown.options = [("— ไม่เลือก —", None)]
    stop_dropdown.options = [("— ไม่เลือก —", None)]
    state["last_protein_header"] = None # ล้างข้อมูลโปรตีนล่าสุด
    state["last_protein_seq"] = None
    save_protein_button.disabled = True
    save_genome_button.disabled = True # ปิดใช้งานปุ่มบันทึกจีโนมเมื่อ reset


# ------------------------------
# Bind events & display UI
# ------------------------------
fetch_button.on_click(on_fetch_clicked)
analyze_button.on_click(on_analyze_clicked)
translate_selected_button.on_click(on_translate_selected_clicked)
use_picks_button.on_click(on_use_picks_clicked)
save_protein_button.on_click(on_save_protein_clicked)
save_genome_button.on_click(on_save_genome_clicked) # ผูก event ใหม่กับฟังก์ชันใหม่
reset_button.on_click(on_reset_clicked)
range_slider.observe(on_range_change)

ui = widgets.VBox([
    widgets.HBox([species_box, fetch_button, reset_button]),
    info_html,
    range_slider,
    widgets.HBox([analyze_button, translate_selected_button, save_protein_button, save_genome_button]), # เปลี่ยนปุ่มใน UI
    widgets.HBox([start_dropdown, stop_dropdown, use_picks_button]),
    output_html
])

display(ui)
