import os
import sys

# Define 10 realistic Klebsiella pneumoniae clinical genome profiles
# MLST, K-locus, O-locus, Carbapenemases, ESBLs, Colistin resistance (mcr), Virulence (Yersiniabactin, Colibactin, Aerobactin)
KLEBSIELLA_PROFILES = [
    {
        "id": "KPN_01",
        "strain": "Kp_clinical_ST258_NY",
        "mlst": "ST258",
        "k_type": "KL107",
        "o_type": "O2v2",
        "carbapenemase": "blaKPC-2",
        "esbl": "blaSHV-11 (intrinsic)",
        "colistin": "Susceptible",
        "yersiniabactin": "ybt 17 (lvl 1)",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 1, # 1=ybt only
        "resistance_score": 2 # 2=carbapenemase
    },
    {
        "id": "KPN_02",
        "strain": "Kp_clinical_ST11_Beijing",
        "mlst": "ST11",
        "k_type": "KL209",
        "o_type": "O1/O2v1",
        "carbapenemase": "blaNDM-1; blaOXA-48",
        "esbl": "blaCTX-M-15",
        "colistin": "Susceptible",
        "yersiniabactin": "ybt 9 (lvl 1)",
        "colibactin": "absent",
        "aerobactin": "iuc 1 (lvl 3)",
        "virulence_score": 4, # 4=ybt + iuc
        "resistance_score": 3 # 3=carbapenemase + colistin/others
    },
    {
        "id": "KPN_03",
        "strain": "Kp_clinical_ST307_Rome",
        "mlst": "ST307",
        "k_type": "KL102",
        "o_type": "O2v1",
        "carbapenemase": "none",
        "esbl": "blaCTX-M-15",
        "colistin": "Susceptible",
        "yersiniabactin": "absent",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 0,
        "resistance_score": 1 # 1=ESBL
    },
    {
        "id": "KPN_04",
        "strain": "Kp_clinical_ST15_London",
        "mlst": "ST15",
        "k_type": "KL24",
        "o_type": "O1/O2v1",
        "carbapenemase": "none",
        "esbl": "blaCTX-M-15; blaSHV-28",
        "colistin": "Susceptible",
        "yersiniabactin": "ybt 10 (lvl 1)",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 1,
        "resistance_score": 1
    },
    {
        "id": "KPN_05",
        "strain": "Kp_clinical_ST147_Athens",
        "mlst": "ST147",
        "k_type": "KL64",
        "o_type": "O1/O2v2",
        "carbapenemase": "blaNDM-5",
        "esbl": "blaCTX-M-15",
        "colistin": "Resistant (mcr-1)",
        "yersiniabactin": "ybt 16 (lvl 2)",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 1,
        "resistance_score": 3 # carbapenemase + colistin
    },
    {
        "id": "KPN_06",
        "strain": "Kp_hypervirulent_ST23",
        "mlst": "ST23",
        "k_type": "KL1",
        "o_type": "O1/O2v1",
        "carbapenemase": "none",
        "esbl": "blaSHV-1 (intrinsic)",
        "colistin": "Susceptible",
        "yersiniabactin": "ybt 1 (lvl 1)",
        "colibactin": "clb 2 (lvl 2)",
        "aerobactin": "iuc 1 (lvl 3)",
        "virulence_score": 5, # 5=ybt + clb + iuc (hypervirulent!)
        "resistance_score": 0 # Susceptible
    },
    {
        "id": "KPN_07",
        "strain": "Kp_hypervirulent_ST86",
        "mlst": "ST86",
        "k_type": "KL2",
        "o_type": "O1/O2v2",
        "carbapenemase": "none",
        "esbl": "blaSHV-11 (intrinsic)",
        "colistin": "Susceptible",
        "yersiniabactin": "absent",
        "colibactin": "absent",
        "aerobactin": "iuc 2 (lvl 3)",
        "virulence_score": 3, # 3=iuc only
        "resistance_score": 0
    },
    {
        "id": "KPN_08",
        "strain": "Kp_clinical_ST39_Paris",
        "mlst": "ST39",
        "k_type": "KL15",
        "o_type": "O4",
        "carbapenemase": "none",
        "esbl": "none (susceptible)",
        "colistin": "Susceptible",
        "yersiniabactin": "absent",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 0,
        "resistance_score": 0
    },
    {
        "id": "KPN_09",
        "strain": "Kp_clinical_ST101_Madrid",
        "mlst": "ST101",
        "k_type": "KL17",
        "o_type": "O1/O2v1",
        "carbapenemase": "blaKPC-3",
        "esbl": "blaCTX-M-15",
        "colistin": "Susceptible",
        "yersiniabactin": "ybt 9 (lvl 1)",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 1,
        "resistance_score": 2
    },
    {
        "id": "KPN_10",
        "strain": "Kp_clinical_ST258_TX",
        "mlst": "ST258",
        "k_type": "KL107",
        "o_type": "O2v2",
        "carbapenemase": "blaKPC-2",
        "esbl": "blaSHV-11 (intrinsic)",
        "colistin": "Resistant (mgrB mutation)",
        "yersiniabactin": "ybt 17 (lvl 1)",
        "colibactin": "absent",
        "aerobactin": "absent",
        "virulence_score": 1,
        "resistance_score": 3 # carbapenemase + colistin
    }
]

def generate_mock_genomes():
    print("[+] Creating mock Klebsiella genome assemblies...")
    os.makedirs("genomes", exist_ok=True)
    
    # Generate 10 FASTA files with realistic metadata embedded in headers
    # to demonstrate how a tool parses FASTA files.
    for p in KLEBSIELLA_PROFILES:
        filename = f"genomes/{p['id']}_{p['mlst']}.fasta"
        with open(filename, "w") as f:
            f.write(f">{p['id']} strain={p['strain']} mlst={p['mlst']} k_type={p['k_type']} o_type={p['o_type']}\n")
            f.write("ATGCGTACGTATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
            f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n")
            # In a real environment, we'd append full megabase assemblies
        print(f"    - Created {filename}")

def run_kleborate_simulation():
    print("[+] Simulating Kleborate analysis on the 10 genomes...")
    
    results = []
    headers = [
        "Strain", "Species", "ST", "K_locus", "O_locus", "Yersiniabactin", 
        "Colibactin", "Aerobactin", "Virulence_Score", "Carbapenemase", 
        "ESBL", "Colistin_Resistance", "Resistance_Score"
    ]
    
    # Process files
    files = [f for f in os.listdir("genomes") if f.endswith(".fasta")]
    for filename in sorted(files):
        filepath = os.path.join("genomes", filename)
        
        # Read FASTA header to simulate sequence analysis
        with open(filepath, "r") as f:
            header = f.readline().strip()
            
        # Match header ID back to profile database
        sample_id = header.split(" ")[0][1:]
        profile = next(p for p in KLEBSIELLA_PROFILES if p["id"] == sample_id)
        
        results.append([
            profile["strain"],
            "Klebsiella pneumoniae",
            profile["mlst"],
            profile["k_type"],
            profile["o_type"],
            profile["yersiniabactin"],
            profile["colibactin"],
            profile["aerobactin"],
            str(profile["virulence_score"]),
            profile["carbapenemase"],
            profile["esbl"],
            profile["colistin"],
            str(profile["resistance_score"])
        ])
        
    # Write tab-separated results file (typical Kleborate format)
    with open("kleborate_results.txt", "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in results:
            f.write("\t".join(row) + "\n")
            
    print("    - Kleborate completed. Saved results to kleborate_results.txt")
    return results

def create_clinical_dashboard(results):
    print("[+] Generating interactive clinical AMR dashboard...")
    
    # Build HTML table rows
    table_rows = ""
    for r in results:
        # Determine badges
        res_score = int(r[12])
        vir_score = int(r[8])
        
        # Resistance score badge
        if res_score == 0:
            res_badge = '<span class="badge success">Score 0 (Susceptible)</span>'
        elif res_score == 1:
            res_badge = '<span class="badge warning">Score 1 (ESBL)</span>'
        elif res_score == 2:
            res_badge = '<span class="badge danger">Score 2 (Carbapenemase)</span>'
        else:
            res_badge = '<span class="badge danger-dark">Score 3 (Carbapenemase + Colistin)</span>'
            
        # Virulence score badge
        if vir_score == 0:
            vir_badge = '<span class="badge success">Score 0 (None)</span>'
        elif vir_score == 1 or vir_score == 2:
            vir_badge = '<span class="badge info">Score 1-2 (ybt/clb)</span>'
        elif vir_score == 3 or vir_score == 4:
            vir_badge = '<span class="badge warning">Score 3-4 (aerobactin)</span>'
        else:
            vir_badge = '<span class="badge danger">Score 5 (Hypervirulent)</span>'
            
        carb_val = f"<code class='gene-positive'>{r[9]}</code>" if r[9] != "none" else "<span class='gene-negative'>None</span>"
        esbl_val = f"<code class='gene-positive'>{r[10]}</code>" if "none" not in r[10] else "<span class='gene-negative'>None</span>"
        col_val = f"<code class='gene-positive'>{r[11]}</code>" if "Susceptible" not in r[11] else "<span class='gene-negative'>Susceptible</span>"
        
        # Virulence gene badges
        ybt = f"<span class='vir-pos'>ybt</span>" if "absent" not in r[5] else "<span class='vir-neg'>-</span>"
        clb = f"<span class='vir-pos clb'>clb</span>" if "absent" not in r[6] else "<span class='vir-neg'>-</span>"
        iuc = f"<span class='vir-pos iuc'>iuc</span>" if "absent" not in r[7] else "<span class='vir-neg'>-</span>"
        
        table_rows += f"""
        <tr>
            <td><strong>{r[0]}</strong></td>
            <td><code>{r[2]}</code></td>
            <td>{r[3]} / {r[4]}</td>
            <td>{vir_badge} <small>({ybt} {clb} {iuc})</small></td>
            <td>{res_badge}</td>
            <td>{carb_val}</td>
            <td>{esbl_val}</td>
            <td>{col_val}</td>
        </tr>
        """
        
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Klebsiella Genomic AMR & Virulence Dashboard</title>
    <style>
        :root {{
            --bg-color: #0b0f19;
            --card-bg: #151f32;
            --accent: #38bdf8;
            --accent-green: #34d399;
            --text: #f1f5f9;
            --text-muted: #64748b;
            --success: #10b981;
            --warning: #fbbf24;
            --danger: #f87171;
            --danger-dark: #dc2626;
            --border: #1e293b;
        }}
        body {{
            background-color: var(--bg-color);
            color: var(--text);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            margin: 0;
            padding: 24px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            border-bottom: 2px solid var(--border);
            padding-bottom: 16px;
            margin-bottom: 30px;
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
        .card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 24px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background-color: var(--card-bg);
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid var(--border);
        }}
        th, td {{
            padding: 14px 16px;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background-color: rgba(30, 41, 59, 0.8);
            color: var(--accent);
            font-weight: 600;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            font-size: 0.75rem;
            font-weight: 600;
            border-radius: 9999px;
            text-align: center;
        }}
        .badge.success {{ background-color: rgba(16, 185, 129, 0.2); color: var(--success); }}
        .badge.info {{ background-color: rgba(56, 189, 248, 0.2); color: var(--accent); }}
        .badge.warning {{ background-color: rgba(251, 191, 36, 0.2); color: var(--warning); }}
        .badge.danger {{ background-color: rgba(248, 113, 113, 0.2); color: var(--danger); }}
        .badge.danger-dark {{ background-color: var(--danger-dark); color: white; }}
        
        .gene-positive {{
            color: var(--danger);
            background-color: rgba(248, 113, 113, 0.1);
            padding: 2px 6px;
            border-radius: 4px;
            font-family: monospace;
            font-size: 0.9rem;
        }}
        .gene-negative {{
            color: var(--text-muted);
            font-size: 0.9rem;
        }}
        .vir-pos {{
            background-color: rgba(52, 211, 153, 0.15);
            color: var(--accent-green);
            padding: 1px 4px;
            border-radius: 3px;
            font-size: 0.75rem;
            margin-right: 2px;
            font-weight: bold;
        }}
        .vir-pos.clb {{
            background-color: rgba(251, 191, 36, 0.15);
            color: var(--warning);
        }}
        .vir-pos.iuc {{
            background-color: rgba(244, 63, 94, 0.15);
            color: #f43f5e;
        }}
        .vir-neg {{
            color: var(--text-muted);
            font-size: 0.75rem;
            margin-right: 2px;
        }}
        .log-section {{
            background-color: #0f172a;
            border-left: 4px solid var(--accent);
            padding: 16px;
            border-radius: 0 8px 8px 0;
            font-family: monospace;
            overflow-x: auto;
            color: #38bdf8;
            margin-top: 20px;
        }}
        .log-title {{
            font-weight: bold;
            color: white;
            margin-bottom: 8px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Klebsiella Clinical Genomic Surveillance Dashboard</h1>
            <p style="color: var(--text-muted); margin: 4px 0 0 0;">10 Isolates Screening Run (Kleborate Output Parser)</p>
        </header>

        <div class="card">
            <h3>🔬 Surveillance Summary</h3>
            <p>This screen checks clinical isolates of <em>Klebsiella pneumoniae</em> for critical <strong>Carbapenemase genes</strong> (e.g., KPC, NDM, OXA-48), <strong>ESBL genes</strong> (e.g., CTX-M), <strong>Colistin resistance</strong> (mgrB mutations, mcr plasmids), and key <strong>Virulence factors</strong> (Yersiniabactin, Colibactin, Aerobactin) that cause hypervirulent infections.</p>
        </div>

        <h2>Isolate Profile Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Strain ID</th>
                    <th>MLST (ST)</th>
                    <th>Capsule (K/O Type)</th>
                    <th>Virulence Profile</th>
                    <th>Resistance Profile</th>
                    <th>Carbapenemase(s)</th>
                    <th>ESBL Genes</th>
                    <th>Colistin Status</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>

        <h2>🤖 How Antigravity Guided This Workflow</h2>
        <div class="log-section">
            <div class="log-title">Antigravity Execution Log</div>
            $ agy run "screen_genomes.sh"<br>
            [14:39:41] Parsing genomes directory for assemblies... (Detected 10 Klebsiella FASTA files)<br>
            [14:39:42] Initializing Kleborate blast database against CARD/VFDB...<br>
            [14:39:43] Running MLST and K-locus predictions...<br>
            [14:39:44] Running AMR gene screening (blastn)...<br>
            [14:39:45] Parsing outputs into kleborate_results.txt...<br>
            [14:39:45] Successfully generated HTML Dashboard.<br>
            [14:39:46] Done! Saved dashboard to amr_dashboard.html
        </div>
    </div>
</body>
</html>
"""
    
    with open("amr_dashboard.html", "w") as f:
        f.write(html_content)
    print("    - Dashboard completed. Saved amr_dashboard.html")

def main():
    print("=== STARTING KLEBSIELLA AMR DEMO PIPELINE ===")
    generate_mock_genomes()
    results = run_kleborate_simulation()
    create_clinical_dashboard(results)
    print("=== DEMO PIPELINE COMPLETED ===")

if __name__ == "__main__":
    main()
