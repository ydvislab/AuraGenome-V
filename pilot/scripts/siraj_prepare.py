"""Second dataset (Siraj et al. 2026, saturation mutagenesis MPRA in K562 and HepG2): map every substitution of
the 128 final elements onto GRCh38 so AlphaGenome can score it in genomic context.

Primary background (declared in docs/internal/preregistration_v2_siraj_zh.md): the reference-allele oligo of each
element (allele R, or RR for two-variant haplotype elements). The data use GRCh37 coordinates: the fine-mapped
variant is mapped to GRCh38 with the Ensembl REST API, and the 200 bp m0 sequence is then located exactly in GRCh38
around it (forward or reverse complement, UCSC REST API); elements whose sequence is not found exactly once are
excluded. Substitutions are then ordinary GRCh38 SNVs. No measured effect is analysed here: effects and
standard errors are only copied into the variant table.

Data: Zenodo 10.5281/zenodo.15297965, data_satmut.zip (CC-BY-NC-4.0). Inputs and outputs stay under the gitignored
pilot/data/; never commit them to the public repository.

usage: python siraj_prepare.py <pilot_dir> [--skip-mapping | --ref-norm]
  --ref-norm  effects relative to the measured reference oligo, written to <pilot_dir>/data/siraj_ref/ (after review)
Inputs:  <pilot_dir>/data/raw/siraj2026/data/satmut/{satmut.txt.gz, satmut_sequences.txt.gz, satmut_final_128elements.txt}
Outputs: <pilot_dir>/data/siraj/variants_grch38.tsv   one row per substitution x cell type
         <pilot_dir>/data/siraj/elements.tsv           one row per element: window, strand, mapping status
         <pilot_dir>/data/siraj/proxy_variants.tsv     the same substitutions in proxy_select.py's format: one
                                                       experiment per element x cell type, with the metadata-matched
                                                       cell line and its related biosamples (declared in the
                                                       pre-registration, chosen from track metadata only)
"""
import os, sys, json, time, urllib.request
import numpy as np
import pandas as pd
from scipy.stats import norm

pilot_dir = sys.argv[1]
skip_mapping = "--skip-mapping" in sys.argv[2:]
src = os.path.join(pilot_dir, "data", "raw", "siraj2026", "data", "satmut")
out = os.path.join(pilot_dir, "data", "siraj")
os.makedirs(out, exist_ok=True)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
rc = lambda s: s.translate(COMP)[::-1]
FLANK = 400

MATCH = {"K562": "K562", "HEPG2": "HepG2"}
RELATED = {"K562": ["common myeloid progenitor, CD34-positive", "hematopoietic multipotent progenitor cell", "HL-60",
                    "myelocyte", "bone marrow", "bone marrow cell"],
           "HEPG2": ["liver", "left lobe of liver", "right lobe of liver", "hepatocyte", "HuH-7", "HuH-7.5"]}
P_SIG, MIN_EFFECT = 1e-5, 0.1      # significance classes as for Kircher (build_variant_table.py); accuracy only


def to_proxy_format(v):
    v = v.assign(locus=v.element.str.replace(":", "_"))
    sig = (v.pvalue < P_SIG) & (v.effect.abs() >= MIN_EFFECT)
    return pd.DataFrame(dict(element=v.locus + "." + v.cell_type, pred_element=v.locus, cell_type=v.cell_type,
                             chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt, effect=v.effect, se=v.se,
                             pvalue=v.pvalue, truth_class=(sig * pd.Series(1, index=v.index).where(v.effect > 0, -1)),
                             model_biosample_keywords=v.cell_type.map(MATCH),
                             related_biosamples=v.cell_type.map(lambda c: ";".join(RELATED[c]))))


if skip_mapping:
    v = pd.read_csv(os.path.join(out, "variants_grch38.tsv"), sep="\t")
    to_proxy_format(v).to_csv(os.path.join(out, "proxy_variants.tsv"), sep="\t", index=False)
    sys.exit(0)

if "--ref-norm" in sys.argv[2:]:
    # review round 2 (docs/internal/siraj_refnorm_zh.md): the published log2Skew subtracts the central mean of every
    # oligo of the background, so a few labels carry information about the whole locus. Here each effect is the
    # substitution's log2FC minus that of the measured unmutated oligo (m0) of the same background and cell type,
    # which an eight-variant experiment would measure too; the standard errors are combined.
    v = pd.read_csv(os.path.join(out, "variants_grch38.tsv"), sep="\t")
    d = pd.read_csv(os.path.join(src, "satmut.txt.gz"), sep="\t", low_memory=False,
                    usecols=["ID", "cell_type", "sat_ref", "oligomut", "log2FC", "log2FC_SE"])
    ref = d[d.oligomut == "m0"].set_index(["sat_ref", "cell_type"])[["log2FC", "log2FC_SE"]]
    m = v.merge(d, left_on=["oligo_id", "cell_type"], right_on=["ID", "cell_type"], how="left", validate="1:1")
    r = ref.loc[list(zip(m.sat_ref, m.cell_type))].to_numpy()
    m["effect"] = m.log2FC - r[:, 0]
    m["se"] = np.sqrt(m.log2FC_SE ** 2 + r[:, 1] ** 2)
    m["pvalue"] = 2 * norm.sf(np.abs(m.effect / m.se))
    assert m.effect.notna().all() and len(m) == len(v)
    out_ref = os.path.join(pilot_dir, "data", "siraj_ref")
    os.makedirs(out_ref, exist_ok=True)
    m = m[v.columns]
    m.to_csv(os.path.join(out_ref, "variants_grch38.tsv"), sep="\t", index=False)
    to_proxy_format(m).to_csv(os.path.join(out_ref, "proxy_variants.tsv"), sep="\t", index=False)
    print(f"{len(m)} rows written to {out_ref}")
    sys.exit(0)

final = [l.strip() for l in open(os.path.join(src, "satmut_final_128elements.txt"))][1:]
d = pd.read_csv(os.path.join(src, "satmut.txt.gz"), sep="\t", low_memory=False)
d = d[d.sat_ref_parent.isin(final) & d.allele.isin(["R", "RR"]) & d.indel.isna()]
seqs = pd.read_csv(os.path.join(src, "satmut_sequences.txt.gz"), sep="\t", usecols=["sequence", "final_id"])
seq_of = dict(zip(seqs.final_id, seqs.sequence))


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 ** attempt)
    raise RuntimeError(url)


def hg38(chrom, start, end):
    """GRCh38 sequence, 0-based half-open, upper case."""
    return get(f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{chrom};start={start};end={end}")["dna"].upper()


def to_grch38(chrom, pos):
    """GRCh37 -> GRCh38 position of one base (Ensembl REST), or None when it does not map uniquely."""
    m = get(f"https://rest.ensembl.org/map/human/GRCh37/{chrom}:{pos}..{pos}:1/GRCh38?content-type=application/json")
    m = [x["mapped"] for x in m["mappings"] if x["mapped"]["seq_region_name"] == chrom]
    return m[0]["start"] if len(m) == 1 else None


elements, rows = [], []
for parent, g in d.groupby("sat_ref_parent", sort=True):
    m0 = g[g.oligomut == "m0"].ID.iloc[0]
    oligo = seq_of[m0].upper()
    chrom, pos37 = str(g.chr.iloc[0]), int(g.pos.iloc[0])    # the fine-mapped variant (var1), GRCh37
    pos = to_grch38(chrom, pos37)
    if pos is None:
        elements.append(dict(element=parent, background=g.sat_ref.iloc[0], chrom=chrom, var1_pos37=pos37,
                             status="unmapped"))
        continue
    ref = hg38(chrom, pos - 1 - FLANK, pos + FLANK)
    hits = [(i, "+") for i in range(len(ref) - len(oligo) + 1) if ref[i:i + len(oligo)] == oligo]
    hits += [(i, "-") for i in range(len(ref) - len(oligo) + 1) if ref[i:i + len(oligo)] == rc(oligo)]
    info = dict(element=parent, background=g.sat_ref.iloc[0], chrom=chrom, var1_pos37=pos37, var1_pos=pos, oligo_len=len(oligo),
                hits=len(hits))
    if len(hits) != 1:
        elements.append(dict(info, status="not_found" if not hits else "ambiguous"))
        continue
    off, strand = hits[0]
    start = pos - 1 - FLANK + off                               # 0-based start of the oligo in GRCh38
    elements.append(dict(info, status="ok", start=start, end=start + len(oligo), strand=strand))
    s = g[g.oligomut != "m0"].copy()
    i = s.mut_pos.astype(int) - 1                               # 0-based position within the oligo
    s["bg"] = [oligo[j] for j in i]
    s["gpos"] = start + 1 + (i if strand == "+" else len(oligo) - 1 - i)   # 1-based GRCh38 position
    s["gref"] = s.bg if strand == "+" else s.bg.str.translate(COMP)
    s["galt"] = s.mut_base.str.upper() if strand == "+" else s.mut_base.str.upper().str.translate(COMP)
    s = s[s.bg != s.mut_base.str.upper()]
    rows.append(pd.DataFrame(dict(element=parent, cell_type=s.cell_type, chrom="chr" + chrom, pos=s.gpos, ref=s.gref,
                                  alt=s.galt, oligo_id=s.ID, effect=s.log2Skew, se=s.log2Skew_SE,
                                  pvalue=s.log2Skew_pval)))
    print(parent, strand, len(s), flush=True)

el = pd.DataFrame(elements)
el.to_csv(os.path.join(out, "elements.tsv"), sep="\t", index=False)
v = pd.concat(rows)
# the oligo matched GRCh38 exactly, so every substitution's reference base is the genome base
v.to_csv(os.path.join(out, "variants_grch38.tsv"), sep="\t", index=False)
to_proxy_format(v).to_csv(os.path.join(out, "proxy_variants.tsv"), sep="\t", index=False)
print(el.status.value_counts().to_dict(), f"{len(v)} rows, {v.element.nunique()} elements, "
      f"{v.drop_duplicates(['element', 'pos', 'alt']).shape[0]} substitutions")
