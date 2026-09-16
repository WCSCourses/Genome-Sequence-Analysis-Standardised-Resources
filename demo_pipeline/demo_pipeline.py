import os
import sys

# Define the reference genome sequence (500 bp)
REF_SEQ = (
    "ATGCGTACGTATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA"
    "CGGCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG"
    "AATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATTTTAATT"
    "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCC"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    "TCGATCGATCG"
)

# Define mock reads (paired-end)
# format: (name, R1_seq, R1_qual, R2_seq, R2_qual)
MOCK_READS = [
    # 1. Perfect match, properly paired
    # R1: bases 1-20 (ATGCGTACGTATCGATCGAT)
    # R2: bases 51-70 RC (GCTAGCTAGCCGATCGATCG -> RC of CGATCGATCGGCTAGCTAGC)
    ("read1", "ATGCGTACGTATCGATCGAT", "IIIIIIIIIIIIIIIIIIII", "GCTAGCTAGCCGATCGATCG", "IIIIIIIIIIIIIIIIIIII"),
    
    # 2. Match with mismatch/SNP in R1, properly paired
    # R1: bases 1-20, last base T -> C (ATGCGTACGTATCGATCGAC)
    # R2: bases 51-70 RC (GCTAGCTAGCCGATCGATCG)
    ("read2", "ATGCGTACGTATCGATCGAC", "IIIIIIIIIIIIIIIIIIII", "GCTAGCTAGCCGATCGATCG", "IIIIIIIIIIIIIIIIIIII"),
    
    # 3. Match with low-quality tail in R1, properly paired
    # R1: bases 141-160 (TATATATATATATATATATA)
    # R2: bases 211-230 RC (CCGGCCGGCCGGCCGGCCGG)
    ("read3", "TATATATATATATATATATA", "IIIIIIIIIIIIII!!!!!!", "CCGGCCGGCCGGCCGGCCGG", "IIIIIIIIIIIIIIIIIIII"),
    
    # 4. Unmapped reads (all Ns or random junk)
    ("read4", "NNNNNNNNNNNNNNNNNNNN", "!!!!!!!!!!!!!!!!!!!!", "NNNNNNNNNNNNNNNNNNNN", "!!!!!!!!!!!!!!!!!!!!"),
    
    # 5. R1 mapped, R2 unmapped (single mapped, mate unmapped)
    # R1: bases 351-370 (GGGGCCCCGGGGCCCCGGGG)
    # R2: unmapped
    ("read5", "GGGGCCCCGGGGCCCCGGGG", "IIIIIIIIIIIIIIIIIIII", "NNNNNNNNNNNNNNNNNNNN", "!!!!!!!!!!!!!!!!!!!!")
]

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def calculate_qscore(qual_str):
    if not qual_str:
        return 0
    return sum(ord(c) - 33 for c in qual_str) / len(qual_str)

def generate_datasets():
    print("[+] Generating mock reference genome and reads...")
    os.makedirs("data", exist_ok=True)
    
    # 1. Reference FASTA
    with open("data/ref.fa", "w") as f:
        f.write(">chr1\n")
        # wrap sequence to 70 bases per line
        for i in range(0, len(REF_SEQ), 70):
            f.write(REF_SEQ[i:i+70] + "\n")
            
    # 2. FASTQ R1 and R2
    with open("data/reads_R1.fastq", "w") as f1, open("data/reads_R2.fastq", "w") as f2:
        for name, r1, q1, r2, q2 in MOCK_READS:
            f1.write(f"@{name}/1\n{r1}\n+\n{q1}\n")
            f2.write(f"@{name}/2\n{r2}\n+\n{q2}\n")
            
    print("    - Created data/ref.fa")
    print("    - Created data/reads_R1.fastq")
    print("    - Created data/reads_R2.fastq")

def run_qc():
    print("[+] Running Quality Control analysis (Simulating FastQC)...")
    
    def parse_fastq(filepath):
        stats = []
        with open(filepath, "r") as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline() # plus line
                qual = f.readline().strip()
                
                gc_content = sum(1 for c in seq if c in "GCgc") / len(seq) * 100 if seq else 0
                avg_q = calculate_qscore(qual)
                stats.append({
                    "name": header[1:],
                    "length": len(seq),
                    "gc": gc_content,
                    "qscore": avg_q,
                    "seq": seq,
                    "qual": qual
                })
        return stats

    r1_stats = parse_fastq("data/reads_R1.fastq")
    r2_stats = parse_fastq("data/reads_R2.fastq")
    
    # Save text report
    with open("qc_report.txt", "w") as f:
        f.write("=== QUALITY CONTROL REPORT ===\n\n")
        f.write("Reads File 1: data/reads_R1.fastq\n")
        f.write(f"Total Reads: {len(r1_stats)}\n")
        f.write(f"Avg Read Length: {sum(r['length'] for r in r1_stats)/len(r1_stats):.1f} bp\n")
        f.write(f"Avg GC Content: {sum(r['gc'] for r in r1_stats)/len(r1_stats):.1f}%\n")
        f.write(f"Avg Phred Score: {sum(r['qscore'] for r in r1_stats)/len(r1_stats):.1f}\n\n")
        
        f.write("Reads File 2: data/reads_R2.fastq\n")
        f.write(f"Total Reads: {len(r2_stats)}\n")
        f.write(f"Avg Read Length: {sum(r['length'] for r in r2_stats)/len(r2_stats):.1f} bp\n")
        f.write(f"Avg GC Content: {sum(r['gc'] for r in r2_stats)/len(r2_stats):.1f}%\n")
        f.write(f"Avg Phred Score: {sum(r['qscore'] for r in r2_stats)/len(r2_stats):.1f}\n")
        
    print("    - QC Analysis complete. Saved qc_report.txt")
    return r1_stats, r2_stats

def align_read(ref_seq, read_seq, max_mismatches=3):
    if "NNNNNNNNNNNNNNNNNNNN" in read_seq:
        return 0, "*", -1
        
    best_dist = len(read_seq) + 1
    best_pos = -1
    best_strand = "+"
    
    # Scan forward strand
    for i in range(len(ref_seq) - len(read_seq) + 1):
        ref_sub = ref_seq[i:i+len(read_seq)]
        dist = sum(1 for a, b in zip(ref_sub, read_seq) if a != b and a != 'N' and b != 'N')
        if dist < best_dist:
            best_dist = dist
            best_pos = i + 1  # 1-based index
            best_strand = "+"
            
    # Scan reverse complement strand
    rc_read = reverse_complement(read_seq)
    for i in range(len(ref_seq) - len(rc_read) + 1):
        ref_sub = ref_seq[i:i+len(rc_read)]
        dist = sum(1 for a, b in zip(ref_sub, rc_read) if a != b and a != 'N' and b != 'N')
        if dist < best_dist:
            best_dist = dist
            best_pos = i + 1
            best_strand = "-"
            
    if best_dist <= max_mismatches:
        return best_pos, best_strand, best_dist
    else:
        return 0, "*", -1

def run_alignment(r1_stats, r2_stats):
    print("[+] Aligning reads to reference (Simulating BWA mem)...")
    
    alignments = []
    
    for idx, (r1, r2) in enumerate(zip(r1_stats, r2_stats)):
        read_name = r1["name"].split("/")[0]
        
        # Align R1
        r1_pos, r1_strand, r1_dist = align_read(REF_SEQ, r1["seq"])
        # Align R2
        r2_pos, r2_strand, r2_dist = align_read(REF_SEQ, r2["seq"])
        
        # Calculate bitwise flag
        # Bit 0 (1): paired
        # Bit 1 (2): mapped in proper pair (both mapped, R1 forward, R2 reverse, and close together)
        # Bit 2 (4): R1 unmapped
        # Bit 3 (8): R2 unmapped
        # Bit 4 (16): R1 on reverse strand
        # Bit 5 (32): R2 on reverse strand
        # Bit 6 (64): R1 first in pair
        # Bit 7 (128): R2 second in pair
        
        r1_mapped = (r1_pos > 0)
        r2_mapped = (r2_pos > 0)
        
        r1_proper = r1_mapped and r2_mapped and (r1_strand == "+") and (r2_strand == "-") and (r2_pos > r1_pos)
        r2_proper = r1_proper
        
        # Build Flags
        r1_flag = 1
        r2_flag = 1
        
        if r1_proper:
            r1_flag |= 2
            r2_flag |= 2
            
        if not r1_mapped:
            r1_flag |= 4
        if not r2_mapped:
            r1_flag |= 8
            
        if not r2_mapped:
            r2_flag |= 4
        if not r1_mapped:
            r2_flag |= 8
            
        if r1_mapped and r1_strand == "-":
            r1_flag |= 16
        if r2_mapped and r2_strand == "-":
            r2_flag |= 16
            
        if r2_mapped and r2_strand == "-":
            r1_flag |= 32
        if r1_mapped and r1_strand == "-":
            r2_flag |= 32
            
        r1_flag |= 64
        r2_flag |= 128
        
        alignments.append({
            "name": read_name,
            "read_num": 1,
            "flag": r1_flag,
            "pos": r1_pos if r1_mapped else 0,
            "strand": r1_strand,
            "dist": r1_dist,
            "seq": r1["seq"],
            "qual": r1["qual"],
            "mapped": r1_mapped,
            "proper": r1_proper
        })
        
        alignments.append({
            "name": read_name,
            "read_num": 2,
            "flag": r2_flag,
            "pos": r2_pos if r2_mapped else 0,
            "strand": r2_strand,
            "dist": r2_dist,
            "seq": r2["seq"],
            "qual": r2["qual"],
            "mapped": r2_mapped,
            "proper": r2_proper
        })
        
    # Write SAM file
    with open("aligned.sam", "w") as f:
        # Header
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        f.write(f"@SQ\tSN:chr1\tLN:{len(REF_SEQ)}\n")
        f.write("@PG\tID:miniBWA\tPN:mini-bwa-mem-sim\tVN:1.0\n")
        
        for a in alignments:
            pos_str = str(a["pos"])
            cigar = f"{len(a['seq'])}M" if a["mapped"] else "*"
            mqual = "60" if a["mapped"] else "0"
            rname = "chr1" if a["mapped"] else "*"
            
            f.write(f"{a['name']}\t{a['flag']}\t{rname}\t{pos_str}\t{mqual}\t{cigar}\t*\t0\t0\t{a['seq']}\t{a['qual']}\n")
            
    print("    - Alignment complete. Saved aligned.sam")
    return alignments

def run_sorting_flagstat(alignments):
    print("[+] Sorting SAM alignments and running metrics (Simulating samtools)...")
    
    # Sort alignments: unmapped go to end, mapped sorted by position
    mapped_alignments = [a for a in alignments if a["mapped"]]
    unmapped_alignments = [a for a in alignments if not a["mapped"]]
    
    sorted_mapped = sorted(mapped_alignments, key=lambda x: x["pos"])
    sorted_all = sorted_mapped + unmapped_alignments
    
    # Write Sorted SAM
    with open("aligned_sorted.sam", "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:chr1\tLN:{len(REF_SEQ)}\n")
        f.write("@PG\tID:miniSamtools\tPN:mini-samtools-sim\tVN:1.0\n")
        
        for a in sorted_all:
            pos_str = str(a["pos"])
            cigar = f"{len(a['seq'])}M" if a["mapped"] else "*"
            mqual = "60" if a["mapped"] else "0"
            rname = "chr1" if a["mapped"] else "*"
            f.write(f"{a['name']}\t{a['flag']}\t{rname}\t{pos_str}\t{mqual}\t{cigar}\t*\t0\t0\t{a['seq']}\t{a['qual']}\n")
            
    # Calculate flagstat stats
    total = len(alignments)
    mapped = sum(1 for a in alignments if a["mapped"])
    unmapped = total - mapped
    properly_paired = sum(1 for a in alignments if a["proper"])
    singletons = sum(1 for a in alignments if a["mapped"] and not a["proper"]) # mapped but mate isn't or isn't proper
    
    # Write flagstat.txt
    with open("flagstat.txt", "w") as f:
        f.write(f"{total} + 0 in total (QC-passed reads + QC-failed reads)\n")
        f.write(f"{mapped} + 0 mapped ({mapped/total*100:.2f}% : N/A)\n")
        f.write(f"{properly_paired} + 0 properly paired ({properly_paired/total*100:.2f}% : N/A)\n")
        f.write(f"{singletons} + 0 singletons ({singletons/total*100:.2f}% : N/A)\n")
        
    print("    - Sorting complete. Saved aligned_sorted.sam")
    print("    - Flagstat complete. Saved flagstat.txt")
    
    return {
        "total": total,
        "mapped": mapped,
        "unmapped": unmapped,
        "properly_paired": properly_paired,
        "singletons": singletons
    }, sorted_all

def decode_flag(flag):
    props = []
    if flag & 1: props.append("Paired-end")
    if flag & 2: props.append("Proper pair")
    if flag & 4: props.append("Read unmapped")
    if flag & 8: props.append("Mate unmapped")
    if flag & 16: props.append("Reverse strand")
    if flag & 32: props.append("Mate on reverse")
    if flag & 64: props.append("First in pair")
    if flag & 128: props.append("Second in pair")
    return ", ".join(props)

def generate_visual_dashboard(qc1, qc2, stats, sorted_alignments):
    print("[+] Generating beautiful HTML dashboard...")
    
    # Build HTML table for alignments
    align_rows = ""
    for a in sorted_alignments:
        status_badge = '<span class="badge success">Mapped</span>' if a["mapped"] else '<span class="badge danger">Unmapped</span>'
        flag_desc = decode_flag(a["flag"])
        strand_char = a["strand"] if a["mapped"] else "-"
        pos_val = str(a["pos"]) if a["mapped"] else "-"
        mismatch_val = str(a["dist"]) if a["mapped"] else "-"
        
        align_rows += f"""
        <tr>
            <td><strong>{a['name']}/{a['read_num']}</strong></td>
            <td>{status_badge}</td>
            <td><code>{a['flag']}</code> <small style="color: #666;">({flag_desc})</small></td>
            <td>{pos_val}</td>
            <td><code>{strand_char}</code></td>
            <td>{mismatch_val}</td>
            <td><code class="seq">{a['seq']}</code></td>
        </tr>
        """

    # Create a visual reference mapping visualization
    # We will slice the 500bp reference and show where reads sit
    mapping_viz = ""
    ref_chunks = []
    chunk_size = 50
    for i in range(0, len(REF_SEQ), chunk_size):
        chunk_ref = REF_SEQ[i:i+chunk_size]
        # Check which reads map in this range
        read_tracks = []
        for a in sorted_alignments:
            if a["mapped"] and i <= a["pos"] - 1 < i + chunk_size:
                offset = (a["pos"] - 1) - i
                read_tracks.append((offset, a["name"] + f"/{a['read_num']}", a["seq"], a["strand"]))
        
        # Render tracks
        track_html = ""
        for offset, name, seq, strand in read_tracks:
            spacer = "&nbsp;" * offset
            arrow = "&rarr;" if strand == "+" else "&larr;"
            track_html += f'<div class="read-track">{spacer}<span class="read-span" title="{name}">{arrow} {seq}</span></div>'
            
        ref_chunks.append(f"""
        <div class="ref-block">
            <div class="ref-header">Positions {i+1} - {i+chunk_size}</div>
            <div class="ref-seq"><code>{chunk_ref}</code></div>
            <div class="tracks">
                {track_html if track_html else '<div class="no-reads">No reads mapped in this window</div>'}
            </div>
        </div>
        """)
    mapping_viz = "\n".join(ref_chunks)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome Sequence Analysis Dashboard</title>
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --accent: #38bdf8;
            --text: #f8fafc;
            --text-muted: #94a3b8;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --border: #334155;
        }}
        body {{
            background-color: var(--bg-color);
            color: var(--text);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            margin: 0;
            padding: 24px;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            border-bottom: 2px solid var(--border);
            padding-bottom: 16px;
            margin-bottom: 24px;
        }}
        h1 {{
            margin: 0;
            font-size: 2.2rem;
            color: var(--accent);
            letter-spacing: -0.025em;
        }}
        h2 {{
            color: var(--accent);
            font-size: 1.4rem;
            margin-top: 32px;
            margin-bottom: 16px;
            border-bottom: 1px solid var(--border);
            padding-bottom: 8px;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 20px;
            margin-bottom: 24px;
        }}
        .card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            transition: transform 0.2s;
        }}
        .card:hover {{
            transform: translateY(-2px);
        }}
        .card-title {{
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            margin-bottom: 8px;
        }}
        .card-value {{
            font-size: 2rem;
            font-weight: 700;
            color: var(--text);
        }}
        .card-value.success {{ color: var(--success); }}
        .card-value.warning {{ color: var(--warning); }}
        .card-value.danger {{ color: var(--danger); }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 12px;
            background-color: var(--card-bg);
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid var(--border);
        }}
        th, td {{
            padding: 12px 16px;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background-color: rgba(51, 65, 85, 0.5);
            font-weight: 600;
            color: var(--accent);
        }}
        tr:last-child td {{
            border-bottom: none;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            font-size: 0.75rem;
            font-weight: 600;
            border-radius: 9999px;
        }}
        .badge.success {{ background-color: rgba(16, 185, 129, 0.2); color: var(--success); }}
        .badge.danger {{ background-color: rgba(239, 68, 68, 0.2); color: var(--danger); }}
        
        code.seq {{
            font-family: monospace;
            background-color: rgba(15, 23, 42, 0.6);
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.9rem;
            word-break: break-all;
        }}
        
        .viz-container {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
            margin-top: 24px;
        }}
        .ref-block {{
            margin-bottom: 20px;
            border-bottom: 1px dashed var(--border);
            padding-bottom: 12px;
        }}
        .ref-block:last-child {{
            border-bottom: none;
            margin-bottom: 0;
            padding-bottom: 0;
        }}
        .ref-header {{
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-bottom: 4px;
        }}
        .ref-seq {{
            font-family: monospace;
            color: #a7f3d0;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
        }}
        .tracks {{
            margin-top: 8px;
            padding-left: 8px;
            border-left: 2px solid var(--accent);
        }}
        .read-track {{
            font-family: monospace;
            font-size: 0.9rem;
            white-space: pre;
            line-height: 1.2;
            margin: 2px 0;
        }}
        .read-span {{
            background-color: rgba(56, 189, 248, 0.15);
            color: var(--accent);
            padding: 1px 4px;
            border-radius: 3px;
            border: 1px solid rgba(56, 189, 248, 0.3);
            cursor: help;
        }}
        .no-reads {{
            font-size: 0.8rem;
            color: var(--text-muted);
            font-style: italic;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Genome Sequence Analysis Pipeline Results</h1>
            <p style="color: var(--text-muted); margin: 4px 0 0 0;">Interactive Dashboard | Simulated NGS Pipeline</p>
        </header>

        <h2>1. QC Metrics Summary</h2>
        <div class="grid">
            <div class="card">
                <div class="card-title">Total Raw Reads (PE)</div>
                <div class="card-value">{len(qc1) + len(qc2)}</div>
            </div>
            <div class="card">
                <div class="card-title">Average Read Quality (Phred)</div>
                <div class="card-value success">
                    {(sum(r['qscore'] for r in qc1) + sum(r['qscore'] for r in qc2)) / (len(qc1) + len(qc2)):.1f}
                </div>
            </div>
            <div class="card">
                <div class="card-title">Average GC %</div>
                <div class="card-value">
                    {(sum(r['gc'] for r in qc1) + sum(r['gc'] for r in qc2)) / (len(qc1) + len(qc2)):.1f}%
                </div>
            </div>
        </div>

        <h2>2. Alignment Summary (Flagstat)</h2>
        <div class="grid">
            <div class="card">
                <div class="card-title">Mapped Reads</div>
                <div class="card-value success">{stats['mapped']} <span style="font-size: 1rem; color: var(--text-muted);">({stats['mapped']/stats['total']*100:.1f}%)</span></div>
            </div>
            <div class="card">
                <div class="card-title">Properly Paired</div>
                <div class="card-value success">{stats['properly_paired']} <span style="font-size: 1rem; color: var(--text-muted);">({stats['properly_paired']/stats['total']*100:.1f}%)</span></div>
            </div>
            <div class="card">
                <div class="card-title">Unmapped</div>
                <div class="card-value danger">{stats['unmapped']} <span style="font-size: 1rem; color: var(--text-muted);">({stats['unmapped']/stats['total']*100:.1f}%)</span></div>
            </div>
        </div>

        <h2>3. Detailed Read Alignment Table (SAM Output)</h2>
        <table>
            <thead>
                <tr>
                    <th>Read ID</th>
                    <th>Status</th>
                    <th>Flag (Decimal & Bits)</th>
                    <th>Position (chr1)</th>
                    <th>Strand</th>
                    <th>Mismatches</th>
                    <th>Sequence</th>
                </tr>
            </thead>
            <tbody>
                {align_rows}
            </tbody>
        </table>

        <h2>4. Alignment Visualization Map</h2>
        <p style="color: var(--text-muted); margin-top: -8px; font-size: 0.9rem;">Shows how read sequences stack and align to the reference chromosome sections.</p>
        <div class="viz-container">
            {mapping_viz}
        </div>
    </div>
</body>
</html>
"""
    with open("pipeline_dashboard.html", "w") as f:
        f.write(html_content)
    print("    - HTML Dashboard complete. Saved pipeline_dashboard.html")

def main():
    print("=== STARTING BIOINFORMATICS PIPELINE DEMONSTRATION ===")
    generate_datasets()
    qc1, qc2 = run_qc()
    alignments = run_alignment(qc1, qc2)
    stats, sorted_alignments = run_sorting_flagstat(alignments)
    generate_visual_dashboard(qc1, qc2, stats, sorted_alignments)
    print("=== PIPELINE RUN COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
