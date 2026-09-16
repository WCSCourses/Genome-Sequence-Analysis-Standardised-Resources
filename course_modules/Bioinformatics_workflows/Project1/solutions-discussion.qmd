## Why the numbers differ

Since both students ran the *same* task on (nominally) the same data, but got meaningfully different counts, the most likely explanations are:

1. **Region actually used differs.** The brief says `-r 6` (all of chromosome 6, but the bam is pre-subset to ~50 Mb). v2's write-up explicitly narrows this to `6:25000000-50000000` (25 Mb) — if v2 really restricted to that half-range while v1 used the full ~50 Mb bam content, that alone would shift every SNP/indel count.
2. **Tool version differences.** BCFtools, FreeBayes, and GATK all changed default heuristics (e.g. genotype likelihood models, indel-realignment behavior) across versions — even a minor version bump changes calls.
3. **FreeBayes has run-to-run variability.** Its haplotype-based calling isn't fully deterministic under all threading/complex-region settings, which could explain why FreeBayes numbers (esp. ts/tv: 1.49 vs 1.98) diverge more than GATK's (which is much more deterministic).
4. **Counting method differs.** v1 counts with `bcftools query -i 'type="SNP"' | wc -l`, while v2 doesn't show its exact commands — if v2 used `bcftools stats` summary counts instead, multiallelic sites can be tallied differently between the two approaches.
5. **`isec` collapse flag differs** (`-c all` in v1 vs default/`-c both` in v2) — this directly changes what counts as "shared" vs "private," so Task 4/5 numbers are expected to diverge for methodological, not random, reasons.

Point 1 (region) and point 5 (isec flag) are the most likely dominant causes — the GATK numbers, which are nearly identical between the two (8,384 vs 8,393 SNPs), suggest the *data* is basically the same; it's the **downstream commands** that push BCFtools/FreeBayes apart.

## Side-by-side comparison

### Task 2 — Raw callset stats

| Caller | Metric | v1 (docx) | v2 (md) | Assessment |
|---|---|---|---|---|
| BCFtools | SNPs / Indels / ts-tv | 20,862 / 1,811 / 1.88 | 22,125 / 1,073 / 1.84 | Close in ts/tv; indel counts diverge ~40% — worth double-checking counting method |
| FreeBayes | SNPs / Indels / ts-tv | 11,034 / 778 / **1.49** | 9,051 / 794 / **1.98** | v2's ts/tv is far more biologically plausible (WGS human ts/tv ≈ 2.0). v1's 1.49 looks like an error |
| GATK | SNPs / Indels / ts-tv | 8,384 / 948 / 2.14 | 8,393 / 964 / 2.14 | Nearly identical — good agreement, suggests same underlying data |

### Task 3 — Filtered (QUAL>20) stats

| Caller | Metric | v1 (docx) | v2 (md) | Assessment |
|---|---|---|---|---|
| BCFtools | SNPs / Indels / ts-tv | 9,684 / 1,054 / 2.17 | 10,182 / 906 / 2.14 | Both show the expected ts/tv improvement; magnitudes differ but direction agrees |
| FreeBayes | SNPs / Indels / ts-tv | 7,571 / 635 / 2.27 | 7,674 / 643 / 2.26 | Very close — good agreement here despite raw-stage divergence |
| GATK | SNPs / Indels / ts-tv | 8,384 / 948 / 2.14 (unchanged) | 8,393 / 964 / 2.14 (0% reduction) | Both correctly identify the GATK QUAL≥30 explanation — this is the core learning point of Task 3, and both nail it |

### Task 4 — Intersection (isec)

| | v1 (docx) | v2 (md) | Assessment |
|---|---|---|---|
| Method | Single-step: `isec -i 'QUAL>20' -c all -n =3` | Two-step: filter first, then isec with default (exact) **and** `-c both` | v2's approach is more rigorous — `-c all` conflates SNPs and indels at the same position, which isn't standard for this comparison |
| Intersection result | ts/tv "2.25–2.29" (no SNP/indel breakdown given) | Exact: 6,728 SNPs / **0 indels** / ts-tv 2.27. Collapsed: 7,088 SNPs / 504 indels / ts-tv 2.23 | v2 catches and explains the classic "0 indels" representation artifact — this is exactly what Task 4's "explore the output" is testing for, and v1 misses it entirely |
| Verdict | Correct qualitative conclusion, but shallower exploration | More correct methodology, more complete answer | **v2 stronger** |

### Task 5 — Private variants

| Caller | v1 ts/tv (only) | v2 ts/tv, exact / collapsed | Assessment |
|---|---|---|---|
| BCFtools private | 2.12 | 1.96 / 2.03 | Same ballpark, v2 gives fuller breakdown (SNP/indel counts too) |
| FreeBayes private | 2.01 | 1.77 / 1.77 | Direction agrees (private < intersection) in both |
| GATK private | 1.66 | 1.55 / 1.62 | Direction agrees; GATK private consistently lowest in both — correct, since GATK's higher specificity means its unique calls are more likely artifacts |
| Verdict | Correct trend, minimal detail | Correct trend, much richer detail (counts + two methods) | **v2 stronger** |

## Bottom line

Both correctly identify every qualitative trend the project is testing for (GATK's QUAL≥30 floor, ts/tv rising with filtering, intersection having the best ts/tv, private calls being the worst). Where they diverge is **rigor of method** — v2 uses the technically correct `isec` flag, surfaces the indel-representation artifact the task is implicitly probing for, and its raw ts/tv values are more biologically sound. v1's FreeBayes raw ts/tv (1.49) is the one number I'd flag as a likely mistake rather than a legitimate run-to-run difference.