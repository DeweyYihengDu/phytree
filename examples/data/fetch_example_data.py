"""Download a small public 16S rRNA example set from NCBI (Entrez).

A "microbial tree of life" sampler: well-known model organisms and type
strains spanning the major bacterial phyla, plus four archaea that form a
natural outgroup.  All sequences are public RefSeq 16S rRNA. Writes:

    examples/data/tol_16S.fasta     unaligned 16S sequences
    examples/data/tol_metadata.csv  per-tip metadata (domain, phylum, ...)

Run:  python examples/data/fetch_example_data.py
"""
import csv
import os
import time

from Bio import Entrez, SeqIO

Entrez.email = "yiheng.du@slu.edu"          # required by NCBI
HERE = os.path.dirname(os.path.abspath(__file__))

# organism -> (tip label, domain, phylum)
ORGANISMS = [
    ("Escherichia coli",              "Escherichia_coli",        "Bacteria", "Pseudomonadota"),
    ("Salmonella enterica",           "Salmonella_enterica",     "Bacteria", "Pseudomonadota"),
    ("Pseudomonas aeruginosa",        "Pseudomonas_aeruginosa",  "Bacteria", "Pseudomonadota"),
    ("Vibrio cholerae",               "Vibrio_cholerae",         "Bacteria", "Pseudomonadota"),
    ("Helicobacter pylori",           "Helicobacter_pylori",     "Bacteria", "Campylobacterota"),
    ("Bacillus subtilis",             "Bacillus_subtilis",       "Bacteria", "Bacillota"),
    ("Staphylococcus aureus",         "Staphylococcus_aureus",   "Bacteria", "Bacillota"),
    ("Streptococcus pneumoniae",      "Streptococcus_pneumoniae","Bacteria", "Bacillota"),
    ("Mycobacterium tuberculosis",    "Mycobacterium_tuberculosis","Bacteria","Actinomycetota"),
    ("Bifidobacterium longum",        "Bifidobacterium_longum",  "Bacteria", "Actinomycetota"),
    ("Bacteroides fragilis",          "Bacteroides_fragilis",    "Bacteria", "Bacteroidota"),
    ("Synechocystis sp. PCC 6803",    "Synechocystis_PCC6803",   "Bacteria", "Cyanobacteriota"),
    ("Thermus thermophilus",          "Thermus_thermophilus",    "Bacteria", "Deinococcota"),
    ("Aquifex aeolicus",              "Aquifex_aeolicus",        "Bacteria", "Aquificota"),
    ("Methanocaldococcus jannaschii", "Methanocaldococcus_jannaschii", "Archaea", "Euryarchaeota"),
    ("Halobacterium salinarum",       "Halobacterium_salinarum", "Archaea", "Euryarchaeota"),
    ("Saccharolobus solfataricus",    "Saccharolobus_solfataricus", "Archaea", "Thermoproteota"),
    ("Pyrococcus furiosus",           "Pyrococcus_furiosus",     "Archaea", "Euryarchaeota"),
]


def fetch_16s(organism: str):
    """esearch a near-full-length 16S, return (accession, sequence) or None."""
    term = (f'{organism}[ORGN] AND 16S ribosomal RNA[Title] '
            f'AND 1200:1600[SLEN] AND biomol_rrna[PROP]')
    h = Entrez.esearch(db="nucleotide", term=term, retmax=3)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        # relax: drop the rRNA biomol filter
        term2 = f'{organism}[ORGN] AND 16S ribosomal RNA[Title] AND 1200:1600[SLEN]'
        h = Entrez.esearch(db="nucleotide", term=term2, retmax=3)
        ids = Entrez.read(h)["IdList"]
        h.close()
    if not ids:
        return None
    h = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
    rec = SeqIO.read(h, "fasta")
    h.close()
    return rec.id, str(rec.seq)


def main():
    fasta_path = os.path.join(HERE, "tol_16S.fasta")
    meta_path = os.path.join(HERE, "tol_metadata.csv")
    rows, fasta = [], []
    for organism, label, domain, phylum in ORGANISMS:
        try:
            res = fetch_16s(organism)
        except Exception as e:                      # network / parse hiccup
            print(f"  ! {organism}: {e}")
            res = None
        if res is None:
            print(f"  - {organism}: no hit, skipped")
            continue
        acc, seq = res
        fasta.append((label, seq))
        rows.append({"name": label, "organism": organism, "domain": domain,
                     "phylum": phylum, "accession": acc, "length": len(seq)})
        print(f"  + {label:30s} {acc:14s} {len(seq)} bp  ({domain}/{phylum})")
        time.sleep(0.4)                             # be nice to NCBI

    with open(fasta_path, "w") as f:
        for name, seq in fasta:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")
    with open(meta_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "organism", "domain", "phylum",
                                          "accession", "length"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(fasta)} sequences -> {fasta_path}")
    print(f"wrote metadata           -> {meta_path}")


if __name__ == "__main__":
    main()
