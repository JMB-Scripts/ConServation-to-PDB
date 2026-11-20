# ==============================================================================
# CS-to-PDB : Sequence Conservation to Protein Structure Mapper
# ==============================================================================
#
# OVERVIEW
# --------
# CS-to-PDB is a standalone application that calculates residue conservation
# from a Multiple Sequence Alignment (MSA) and maps evolutionary information
# onto three-dimensional protein structures by storing conservation scores in
# the PDB B-factor field.
#
# The generated structures can be visualized directly in PyMOL or ChimeraX
# using identical color palettes. The program also produces an interactive
# HTML alignment viewer to compare sequence conservation with structural data,
# with the secondary structure (DSSP) drawn above the alignment.
#
#
# INPUT FILES
# -----------
# • Alignment (FASTA / A3M):
#     - Accepts standard aligned FASTA ('-' or '.' for gaps). Upper and lower
#       case are both accepted: the sequences are read case-insensitively.
#     - Accepts A3M format (common in HH-suite/AlphaFold): lowercase letters of
#       the aligned sequences are insertions relative to the reference and are
#       removed. The reference (first) sequence is always kept whole.
#     - NOTE: The FIRST sequence in the alignment must exactly match the
#       sequence of the selected chain of the structure (gaps excluded).
#
# • Structure (PDB / mmCIF):
#     - Only the first model is used (NMR ensembles: model 0).
#     - Choose the chain, or let "Auto" pick the chain whose sequence matches
#       the reference sequence of the alignment.
#
# • Secondary structure:
#     - Computed with DSSP (mkdssp, found automatically in PATH or in CCP4;
#       the path can also be given in the interface).
#     - Without mkdssp, the HELIX/SHEET records of the PDB file (or
#       _struct_conf/_struct_sheet_range of the mmCIF) are used.
#
#
# FEATURES
# --------
# • Calculates residue conservation using multiple scoring methods:
#     - BLOSUM62 similarity
#     - Shannon Entropy
#     - Property Conservation
#     - Jensen-Shannon Divergence (JSD)
#     - Rate4Site-like weighted conservation
#
# • Optional gap penalty for poorly aligned regions.
#
# • Colormap picker borrowed from MolNympheas: every Matplotlib colormap
#   (+ cmcrameri / cmocean if installed), grouped, searchable, with gradient
#   thumbnails, badges, reverse toggle and range slider. The same final map is
#   used in all outputs.
#
#
# OUTPUT FILES  (in the chosen folder; default = folder of the structure)
# ------------
# Files are named after the structure and the metric, e.g. PCPE1_entropy.*
# so runs with different metrics never overwrite each other or the input.
# • [name]_[metric].pdb             Model 0 with scores as B-factors
# • [name]_[metric].txt             Residue-by-residue table (score + DSSP)
# • [name]_[metric].pml             PyMOL visualization script
# • [name]_[metric].cxc             ChimeraX visualization script
# • [name]_[metric]_alignment.html  Alignment viewer with DSSP track
# • [name]_[metric]_settings.json   Every setting of the run; reload it with
#                             "Load settings…" to redo or tweak the run
#
#
# USAGE
# -----
#   python CS-to-pdb-plugin-v10.py          (standalone)
#   PyMOL: run CS-to-pdb-plugin-v10.py
#   Requires: biopython, matplotlib, numpy and one of PySide6 / PyQt6 /
#   PyQt5 / PySide2.  Optional: mkdssp (DSSP 4.x, 3.x or 2.x).
#
#
# VERSIONING
# ----------
# v9  : MolNympheas colormap picker (Tkinter).
# v10 : Qt interface (PySide6 / PyQt6 / PyQt5 / PySide2).
#       DSSP secondary structure above the alignment (HTML) and in the table.
#       Model 0 only, chain choice (Auto = chain matching the alignment).
#       Detailed sequence-mismatch report (first differing residue + context).
#       Output folder choice. JSON settings file written with the results and
#       reloadable. Alignment read case-insensitively. Duplicate sequence IDs
#       removed for every format. Alignment numbering follows the structure.
#
# ==============================================================================
import sys
import os
import re
import io
import json
import glob
import math
import shutil
import tempfile
import datetime
import importlib
import subprocess
import warnings
import html as htmllib
from collections import Counter
from urllib.parse import quote

from Bio import SeqIO, PDB
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
try:
    from Bio import BiopythonDeprecationWarning
    warnings.simplefilter('ignore', BiopythonDeprecationWarning)
except ImportError:
    pass

import matplotlib
import matplotlib.colors as mcolors
import numpy as np

# Optional scientific colormap packages: they register their maps in
# matplotlib on import (cmc.batlow, cmo.thermal, ...).
for _extra in ("cmcrameri.cm", "cmocean"):
    try:
        importlib.import_module(_extra)
    except Exception:
        pass

VERSION = 10

# ============================================================
# Constants for the new scoring metrics
# ============================================================

# Physicochemical grouping used by the "property" metric. Six broad classes;
# feel free to refine (e.g. splitting "special" or separating His out of
# "positive") if your protein family needs finer granularity.
PROPERTY_GROUPS = {
    'A': 'aliphatic', 'V': 'aliphatic', 'L': 'aliphatic', 'I': 'aliphatic', 'M': 'aliphatic',
    'F': 'aromatic', 'W': 'aromatic', 'Y': 'aromatic',
    'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar', 'C': 'polar',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'D': 'negative', 'E': 'negative',
    'G': 'special', 'P': 'special',
}

# Standard amino acid background frequencies (Robinson & Robinson, 1991),
# used as the reference distribution for the Jensen-Shannon Divergence metric.
BACKGROUND_FREQ = {
    'A': 0.078, 'R': 0.051, 'N': 0.041, 'D': 0.052, 'C': 0.024,
    'Q': 0.034, 'E': 0.059, 'G': 0.083, 'H': 0.025, 'I': 0.062,
    'L': 0.092, 'K': 0.056, 'M': 0.024, 'F': 0.044, 'P': 0.043,
    'S': 0.059, 'T': 0.055, 'W': 0.014, 'Y': 0.034, 'V': 0.072,
}
_bg_total = sum(BACKGROUND_FREQ.values())
BACKGROUND_FREQ = {aa: f / _bg_total for aa, f in BACKGROUND_FREQ.items()}

METRIC_LABELS = {
    "blosum": "BLOSUM62",
    "entropy": "Shannon Entropy",
    "property": "Property Conservation",
    "jsd": "Jensen-Shannon Divergence",
    "rate4site": "Rate4Site-like (proxy)",
}

# ============================================================
# Scoring helper functions
# ============================================================

def shannon_entropy_score(valid_residues):
    """Classic Shannon-entropy conservation score, 0-100 (100 = fully conserved)."""
    if not valid_residues:
        return 0
    counts = Counter(valid_residues)
    total = len(valid_residues)
    H = sum(-(c / total) * math.log2(c / total) for c in counts.values())
    return (1 - H / math.log2(20)) * 100


def blosum_column_score(valid_residues, matrix):
    """Average pairwise BLOSUM62 self/cross similarity for a column, 0-100."""
    valid = [aa for aa in valid_residues if aa in matrix.alphabet]
    if not valid:
        return 0
    raw_score = sum(matrix[a, a] for a in valid)
    for i, aa1 in enumerate(valid):
        for aa2 in valid[i + 1:]:
            raw_score += 2 * matrix[aa1, aa2]
    max_score = len(valid) ** 2 * max(matrix[a, a] for a in valid)
    return (raw_score / max_score) * 100 if max_score > 0 else 0


def property_conservation_score(valid_residues):
    """
    Groups residues into physicochemical classes (see PROPERTY_GROUPS) and
    scores a column by the fraction belonging to the majority class.
    A column that is 100% hydrophobic (even if the exact residue varies,
    e.g. L/I/V) scores near 100.
    """
    groups = [PROPERTY_GROUPS[aa] for aa in valid_residues if aa in PROPERTY_GROUPS]
    if not groups:
        return 0
    counts = Counter(groups)
    majority_fraction = counts.most_common(1)[0][1] / len(groups)
    return majority_fraction * 100


def jsd_conservation_score(valid_residues, background_freq=BACKGROUND_FREQ):
    """
    Jensen-Shannon Divergence between a column's amino-acid distribution and
    a fixed background distribution (Capra & Singh, 2007 style, without the
    sliding-window smoothing they apply). JSD is bounded in [0, 1] bit for
    two distributions, so the result is simply scaled by 100.
    """
    valid = [aa for aa in valid_residues if aa in background_freq]
    if not valid:
        return 0
    counts = Counter(valid)
    total = len(valid)
    col_freq = {aa: counts.get(aa, 0) / total for aa in background_freq}

    mixture = {aa: 0.5 * (col_freq[aa] + background_freq[aa]) for aa in background_freq}

    def kl_divergence(p, q):
        return sum(p[aa] * math.log2(p[aa] / q[aa]) for aa in p if p[aa] > 0 and q[aa] > 0)

    jsd = 0.5 * kl_divergence(col_freq, mixture) + 0.5 * kl_divergence(background_freq, mixture)
    return max(0.0, min(jsd, 1.0)) * 100


def compute_henikoff_weights(alignment):
    """
    Henikoff & Henikoff (1994) position-based sequence weighting.
    Down-weights sequences that are over-represented in the alignment
    (near-duplicates contribute less), which is the key ingredient the
    "rate4site-like" proxy borrows from real evolutionary-rate estimators.
    Returns one weight per sequence, normalized to sum to len(alignment).
    """
    n_seq = len(alignment)
    aln_len = alignment.get_alignment_length()
    weights = [0.0] * n_seq

    for col in range(aln_len):
        column = alignment[:, col]
        counts = Counter(column)
        distinct = len(counts)
        if distinct == 0:
            continue
        for i, aa in enumerate(column):
            weights[i] += 1.0 / (distinct * counts[aa])

    total = sum(weights)
    if total == 0:
        return [1.0] * n_seq
    return [w * n_seq / total for w in weights]


def rate4site_like_score(column, weights):
    """
    SIMPLIFIED PROXY for site-specific evolutionary rate, NOT the actual
    Rate4Site algorithm (which requires ML rate estimation on a phylogenetic
    tree). Here we compute a sequence-weighted Shannon entropy: residues from
    over-represented (near-identical) sequences count less, giving a rate
    estimate that is less biased by uneven taxon sampling than raw entropy.
    Returned as a conservation score, 0-100 (100 = most conserved / slowest).
    """
    weighted_counts = {}
    total_weight = 0.0
    for aa, w in zip(column, weights):
        if aa == "-":
            continue
        weighted_counts[aa] = weighted_counts.get(aa, 0.0) + w
        total_weight += w

    if total_weight == 0:
        return 0

    H = 0.0
    for c in weighted_counts.values():
        p = c / total_weight
        if p > 0:
            H += -p * math.log2(p)
    return (1 - H / math.log2(20)) * 100


def compute_scores(alignment, metric="blosum", gap_penalty=False):
    """
    Dispatches to the appropriate per-column scoring function.
    metric one of: "blosum", "entropy", "property", "jsd", "rate4site".
    All metrics are normalized to a common 0-100 scale so downstream
    B-factor mapping / color spectra work unchanged regardless of choice.
    """
    scores = []
    aln_len = alignment.get_alignment_length()

    matrix = substitution_matrices.load("BLOSUM62") if metric == "blosum" else None
    # Only needed for the rate4site-like proxy; cheap to skip otherwise.
    henikoff_weights = compute_henikoff_weights(alignment) if metric == "rate4site" else None

    for col in range(aln_len):
        column = alignment[:, col]
        valid = [aa for aa in column if aa != "-"]

        if metric == "entropy":
            score = shannon_entropy_score(valid)
        elif metric == "property":
            score = property_conservation_score(valid)
        elif metric == "jsd":
            score = jsd_conservation_score(valid)
        elif metric == "rate4site":
            score = rate4site_like_score(column, henikoff_weights)
        else:  # "blosum" (default/fallback)
            score = blosum_column_score(valid, matrix)

        if gap_penalty:
            score *= (1 - column.count("-") / len(column))
        scores.append(score)
    return scores



# ============================================================
# Colormap catalogue (same families and badges as MolNympheas)
# ============================================================

MOLNYMPHEAS_COLORS = ['#1d3557', '#457b9d', '#a8d5ba', '#f1faee', '#fbb4b9', '#b56576']

MPL_GROUPS = [
    ("Perceptually uniform", ["viridis", "plasma", "inferno", "magma", "cividis"]),
    ("Sequential", ["Greys", "Purples", "Blues", "Greens", "Oranges", "Reds", "YlOrBr", "YlOrRd", "OrRd", "PuRd",
                    "RdPu", "BuPu", "GnBu", "PuBu", "YlGnBu", "PuBuGn", "BuGn", "YlGn"]),
    ("Sequential (2)", ["binary", "gist_yarg", "gist_gray", "gray", "bone", "pink", "spring", "summer", "autumn",
                        "winter", "cool", "Wistia", "hot", "afmhot", "gist_heat", "copper"]),
    ("Diverging", ["PiYG", "PRGn", "BrBG", "PuOr", "RdGy", "RdBu", "RdYlBu", "RdYlGn", "Spectral", "coolwarm",
                   "bwr", "seismic"]),
    ("Diverging · dark centre", ["berlin", "managua", "vanimo"]),
    ("Cyclic", ["twilight", "twilight_shifted", "hsv"]),
    ("Qualitative", ["Pastel1", "Pastel2", "Paired", "Accent", "Dark2", "Set1", "Set2", "Set3", "tab10", "tab20",
                     "tab20b", "tab20c"]),
    ("Miscellaneous", ["turbo", "cubehelix", "terrain", "ocean", "gist_earth", "gist_stern", "gnuplot", "gnuplot2",
                       "CMRmap", "brg", "gist_rainbow", "rainbow", "jet", "nipy_spectral", "gist_ncar", "flag", "prism"]),
]
CMAP_ALIASES = {"Grays", "grey", "gist_grey", "gist_yerg"}
CVD_SAFE = {"viridis", "plasma", "inferno", "magma", "cividis"}
NEW_310 = {"berlin", "managua", "vanimo"}


class CmapEntry:
    def __init__(self, name, group, qual=False):
        self.name, self.group, self.qual = name, group, qual
        self.badges = []
        if name in CVD_SAFE or name.startswith("cmc."):
            self.badges.append("CVD-safe")
        if name in NEW_310:
            self.badges.append("New · mpl 3.10")
        if name.startswith("cmc."):
            self.badges.append("cmcrameri")
        if name.startswith("cmo."):
            self.badges.append("cmocean")


def registered_cmaps():
    try:
        return list(matplotlib.colormaps)
    except AttributeError:  # matplotlib < 3.5
        import matplotlib.pyplot as plt
        return list(plt.colormaps())


def base_cmap(name):
    """The colormap as registered, before reverse/range are applied."""
    if name.lower() == "molnympheas":
        return mcolors.LinearSegmentedColormap.from_list("molnympheas", MOLNYMPHEAS_COLORS)
    try:
        return matplotlib.colormaps[name]
    except (AttributeError, KeyError):
        try:
            import matplotlib.cm as mcm
            return mcm.get_cmap(name)
        except (ValueError, AttributeError):
            return matplotlib.colormaps["viridis"]


def is_qualitative(cm):
    return isinstance(cm, mcolors.ListedColormap) and cm.N <= 24


def build_catalogue():
    cat = [CmapEntry("molnympheas", "Custom")]
    reg = {n for n in registered_cmaps() if not n.endswith("_r")}
    used = set()
    for group, names in MPL_GROUPS:
        for n in names:
            if n in reg:
                cat.append(CmapEntry(n, group, is_qualitative(base_cmap(n))))
                used.add(n)
    extra = {"cmc.": "Scientific · cmcrameri", "cmo.": "Scientific · cmocean"}
    for n in sorted(reg - used - CMAP_ALIASES, key=str.lower):
        group = next((g for p, g in extra.items() if n.startswith(p)), "Other")
        cat.append(CmapEntry(n, group, is_qualitative(base_cmap(n))))
    order = ["Custom"] + [g for g, _ in MPL_GROUPS] + list(extra.values()) + ["Other"]
    cat.sort(key=lambda e: order.index(e.group) if e.group in order else 99)
    return cat


def make_cmap(name, reverse=False, lo=0.0, hi=1.0):
    """
    Final colormap used by every output. The range [lo, hi] is taken on the
    (possibly reversed) map, exactly as MolNympheas does, then stretched back
    over 0-100 conservation.
    """
    base = base_cmap(name)
    if is_qualitative(base):
        colors = list(base.colors)
        return mcolors.ListedColormap(colors[::-1] if reverse else colors, name=name)
    if reverse:
        base = base.reversed()
    if lo <= 0.0 and hi >= 1.0:
        return base
    return mcolors.ListedColormap(base(np.linspace(lo, hi, 256)), name=f"{name}_{lo:.2f}_{hi:.2f}")


def palette_label(name, reverse=False, lo=0.0, hi=1.0):
    extras = []
    if reverse:
        extras.append("reversed")
    if lo > 0.0 or hi < 1.0:
        extras.append(f"range {lo:.2f}–{hi:.2f}")
    return name + (f" ({', '.join(extras)})" if extras else "")


def palette_tag(name, reverse=False):
    """Safe identifier for PyMOL colour names (cmc.batlow -> cmc_batlow)."""
    return re.sub(r"\W", "_", name) + ("_r" if reverse else "")


def cmap_hex(cmap, t):
    r, g, b, _ = cmap(float(min(max(t, 0.0), 1.0)))
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def cmap_samples(cmap, k=48, hard=False):
    """List of hex colours used to draw previews."""
    if hard and hasattr(cmap, "colors"):
        return [mcolors.to_hex(c) for c in cmap.colors]
    return [cmap_hex(cmap, i / (k - 1)) for i in range(k)]


# ============================================================
# Alignment reading
# ============================================================

class MappingError(ValueError):
    """Error whose message is meant to be shown as-is (monospace) to the user."""


def read_alignment(path):
    """
    Reads an aligned FASTA or an A3M/A2M file.
    • Aligned FASTA (all sequences the same length): read case-insensitively,
      '.' and '~' are gaps.
    • A3M (lengths differ, or .a3m/.a2m extension): lowercase letters and '.'
      of the aligned sequences are insertions and are removed. The reference
      (first) sequence is kept whole and upper-cased.
    Duplicate IDs are removed (first occurrence kept) in both cases.
    Returns (alignment, notes) where notes is a list of human-readable remarks.
    """
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise MappingError("Alignment file is empty or not in FASTA / A3M format.")
    notes = []
    raw = [str(r.seq) for r in records]
    lengths = {len(s) for s in raw}
    is_a3m = path.lower().endswith((".a3m", ".a2m")) or len(lengths) > 1

    if is_a3m:
        ref = raw[0].replace(".", "-").replace("~", "-").upper()
        seqs = [ref]
        fixed = 0
        for s in raw[1:]:
            s = "".join(c for c in s if not (c.islower() or c == "."))
            s = s.replace("~", "-").upper()
            if len(s) != len(ref):
                fixed += 1
                s = s.ljust(len(ref), "-")[:len(ref)]
            seqs.append(s)
        if any(c.islower() for c in raw[0]):
            notes.append("The reference sequence contained lowercase letters: they were kept as residues.")
        if fixed:
            notes.append(f"{fixed} A3M sequence(s) did not match the reference length after removing "
                         "insertions and were padded/truncated.")
        fmt = "A3M"
    else:
        seqs = [s.replace(".", "-").replace("~", "-").upper() for s in raw]
        fmt = "aligned FASTA"

    out, seen, dupes = [], set(), 0
    for rec, s in zip(records, seqs):
        if rec.id in seen:
            dupes += 1
            continue
        seen.add(rec.id)
        out.append(SeqRecord(Seq(s), id=rec.id, description=rec.description))
    if dupes:
        notes.append(f"{dupes} sequence(s) with a duplicate ID were removed (first occurrence kept).")
    return MultipleSeqAlignment(out), fmt, notes


def ungapped(seq):
    return "".join(c for c in seq if c != "-")


# ============================================================
# Structure: model 0, chains, sequence check
# ============================================================

def load_model0(path):
    low = path.lower()
    parser = PDB.MMCIFParser(QUIET=True) if low.endswith((".cif", ".mmcif")) else PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("struct", path)
    try:
        model = next(iter(structure))
    except StopIteration:
        raise MappingError("The structure file contains no model.")
    return structure, model


def chain_residues(chain):
    return [r for r in chain if PDB.is_aa(r, standard=True)]


def residues_sequence(residues):
    return "".join(seq1(r.get_resname()) for r in residues)


def res_label(r):
    icode = r.id[2].strip()
    return f"{r.get_resname().capitalize()}{r.id[1]}{icode}"


def list_chains(path):
    """[(chain_id, n_residues, sequence)] for the protein chains of model 0."""
    _, model = load_model0(path)
    out = []
    for ch in model:
        res = chain_residues(ch)
        if res:
            out.append((ch.id, len(res), residues_sequence(res)))
    return out


def describe_mismatch(chain_id, residues, struct_seq, msa_seq):
    """Plain-text report: first differing residue with context."""
    n = min(len(struct_seq), len(msa_seq))
    i = next((k for k in range(n) if struct_seq[k] != msa_seq[k]), n)
    lines = [f"Chain {chain_id}: structure {len(struct_seq)} residues, "
             f"alignment reference {len(msa_seq)} residues (gaps removed)."]
    if i < n:
        lines.append(f"First difference at position {i + 1} (structure residue {chain_id}:{res_label(residues[i])}): "
                     f"structure '{struct_seq[i]}'  vs  alignment '{msa_seq[i]}'.")
    elif len(struct_seq) > len(msa_seq):
        lines.append(f"The alignment reference stops after {n} residues; the structure continues at "
                     f"{chain_id}:{res_label(residues[n])}.")
    else:
        lines.append(f"The structure stops after {n} residues (last: {chain_id}:{res_label(residues[n - 1]) if n else '-'}); "
                     f"the alignment reference has {len(msa_seq) - n} more.")
    if 0 < i < len(residues) and residues[i].id[1] - residues[i - 1].id[1] > 1:
        lines.append(f"Residue numbering jumps from {residues[i - 1].id[1]} to {residues[i].id[1]} here: "
                     "the structure probably has unmodelled residues that the alignment reference still contains.")
    a = max(0, i - 15)
    b = i + 16
    lines.append("")
    lines.append(f"  from position {a + 1}:")
    lines.append("  PDB  " + struct_seq[a:b])
    lines.append("  MSA  " + msa_seq[a:b])
    lines.append("       " + " " * (i - a) + "^")
    return "\n".join(lines), i


def pick_chain(model, msa_seq, choice):
    """
    Returns (chain_id, residues). choice is a chain ID or 'auto'.
    Raises MappingError with a detailed report when nothing matches.
    """
    chains = [(ch.id, chain_residues(ch)) for ch in model]
    chains = [(cid, res) for cid, res in chains if res]
    if not chains:
        raise MappingError("No protein chain found in model 0 of the structure.")
    if choice and choice != "auto":
        sel = [(cid, res) for cid, res in chains if cid == choice]
        if not sel:
            raise MappingError(f"Chain '{choice}' not found in model 0 "
                               f"(available: {', '.join(c for c, _ in chains)}).")
        cid, res = sel[0]
        seq = residues_sequence(res)
        if seq != msa_seq:
            report, _ = describe_mismatch(cid, res, seq, msa_seq)
            raise MappingError("Sequence mismatch between the structure and the alignment.\n\n" + report +
                               "\n\nThe first sequence of the alignment (gaps removed) must be identical to the "
                               "sequence of the selected chain.")
        return cid, res
    reports = []
    for cid, res in chains:
        seq = residues_sequence(res)
        if seq == msa_seq:
            return cid, res
        reports.append(describe_mismatch(cid, res, seq, msa_seq))
    reports.sort(key=lambda t: -t[1])  # chain matching longest first
    raise MappingError("No chain of model 0 matches the first sequence of the alignment.\n\n" +
                       "\n\n".join(r for r, _ in reports))


# ============================================================
# Secondary structure (DSSP, or records of the input file)
# ============================================================

SS_NAMES = {"H": "α-helix", "G": "3₁₀-helix", "I": "π-helix", "P": "PPII helix", "E": "β-strand",
            "B": "β-bridge", "T": "turn", "S": "bend", "-": "coil"}

_MKDSSP_GLOBS = [
    "~/programs/CCP4/*/bin/mkdssp", "~/programs/ccp4*/bin/mkdssp",
    "/Applications/ccp4-*/bin/mkdssp", "/Applications/CCP4*/bin/mkdssp",
    "/opt/xtal/ccp4-*/bin/mkdssp", "/opt/ccp4*/bin/mkdssp", "/usr/local/ccp4*/bin/mkdssp",
    "/opt/homebrew/bin/mkdssp", "/usr/local/bin/mkdssp", "/usr/bin/mkdssp",
    "~/miniconda3/bin/mkdssp", "~/miniforge3/bin/mkdssp", "~/anaconda3/bin/mkdssp", "~/mambaforge/bin/mkdssp",
]


def find_mkdssp(user_path=""):
    """User-given path first, then PATH, then the usual CCP4 / conda / brew locations."""
    if user_path:
        p = os.path.expanduser(user_path)
        return p if os.path.isfile(p) and os.access(p, os.X_OK) else None
    for name in ("mkdssp", "dssp"):
        p = shutil.which(name)
        if p:
            return p
    for pattern in _MKDSSP_GLOBS:
        hits = sorted(glob.glob(os.path.expanduser(pattern)), reverse=True)  # newest version first
        for p in hits:
            if os.access(p, os.X_OK):
                return p
    return None


def parse_dssp_classic(text):
    """{(chain, resseq, icode): ss} from a classic DSSP file."""
    ss, started = {}, False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 17 or line[13] == "!":
            continue
        try:
            resseq = int(line[5:10])
        except ValueError:
            continue
        key = (line[11], resseq, line[10].strip())
        ss[key] = line[16] if line[16] != " " else "-"
    return ss


def run_dssp(model, exe):
    """
    Runs mkdssp on model 0. mkdssp 4.x refuses PDB files that do not start
    with a HEADER line (AlphaFold, PDBIO output...), so one is added together
    with a dummy CRYST1. Falls back to the 2.x/3.x '-i/-o' syntax.
    """
    tmp = tempfile.mkdtemp(prefix="cs2pdb_")
    try:
        buf = io.StringIO()
        pio = PDB.PDBIO()
        pio.set_structure(model)
        pio.save(buf)
        inp = os.path.join(tmp, "model.pdb")
        out = os.path.join(tmp, "model.dssp")
        with open(inp, "w") as f:
            f.write("HEADER    CS-TO-PDB                               01-JAN-00   XXXX              \n")
            f.write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n")
            f.write(buf.getvalue())
        attempts = [[exe, "--output-format", "dssp", inp, out], [exe, "-i", inp, "-o", out]]
        last_err = ""
        for cmdline in attempts:
            if os.path.exists(out):
                os.remove(out)
            try:
                proc = subprocess.run(cmdline, capture_output=True, text=True, timeout=300)
            except (OSError, subprocess.TimeoutExpired) as e:
                last_err = str(e)
                continue
            if os.path.exists(out) and os.path.getsize(out) > 0:
                with open(out) as f:
                    ss = parse_dssp_classic(f.read())
                if ss:
                    return ss
            last_err = (proc.stderr or proc.stdout or "").strip()[:300]
        raise RuntimeError(last_err or "mkdssp produced no output")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def ss_from_file_records(path, chain_id, residues):
    """HELIX/SHEET records (PDB) or _struct_conf/_struct_sheet_range (mmCIF)."""
    ranges = []  # (ss, beg_resseq, end_resseq)
    if path.lower().endswith((".cif", ".mmcif")):
        d = MMCIF2Dict(path)

        def col(key):
            v = d.get(key, [])
            return v if isinstance(v, list) else [v]
        conf_map = {"HELX_RH_3T_P": "G", "HELX_RH_PI_P": "I", "HELX_LH_PP_P": "P", "STRN": "E", "TURN_TY1_P": "T",
                    "TURN_P": "T", "BEND": "S"}
        for t, c, b, e in zip(col("_struct_conf.conf_type_id"), col("_struct_conf.beg_auth_asym_id"),
                              col("_struct_conf.beg_auth_seq_id"), col("_struct_conf.end_auth_seq_id")):
            if c == chain_id:
                ranges.append((conf_map.get(t, "H" if t.startswith("HELX") else "-"), b, e))
        for c, b, e in zip(col("_struct_sheet_range.beg_auth_asym_id"), col("_struct_sheet_range.beg_auth_seq_id"),
                           col("_struct_sheet_range.end_auth_seq_id")):
            if c == chain_id:
                ranges.append(("E", b, e))
    else:
        helix_class = {"1": "H", "3": "I", "5": "G"}
        with open(path, errors="replace") as f:
            for line in f:
                if line.startswith("HELIX ") and len(line) >= 38 and line[19] == chain_id:
                    ranges.append((helix_class.get(line[38:40].strip(), "H"), line[21:25], line[33:37]))
                elif line.startswith("SHEET ") and len(line) >= 38 and line[21] == chain_id:
                    ranges.append(("E", line[22:26], line[33:37]))
    ss = ["-"] * len(residues)
    for s, b, e in ranges:
        try:
            b, e = int(str(b).strip()), int(str(e).strip())
        except ValueError:
            continue
        for k, r in enumerate(residues):
            if b <= r.id[1] <= e:
                ss[k] = s
    return ss if ranges else None


def secondary_structure(model, chain_id, residues, structure_path, mkdssp_path=""):
    """Returns (list of one SS letter per residue or None, source description, warning or '')."""
    exe = find_mkdssp(mkdssp_path)
    warning = ""
    if exe:
        try:
            table = run_dssp(model, exe)
            ss = [table.get((chain_id, r.id[1], r.id[2].strip()), "-") for r in residues]
            return ss, f"DSSP ({os.path.basename(exe)})", ""
        except Exception as e:
            warning = f"mkdssp failed ({e}); secondary structure taken from the file records instead."
    elif mkdssp_path:
        warning = f"mkdssp not found at '{mkdssp_path}'."
    try:
        ss = ss_from_file_records(structure_path, chain_id, residues)
    except Exception:
        ss = None
    if ss:
        return ss, "HELIX/SHEET records of the structure file", warning
    return None, "none", warning or "No mkdssp and no secondary structure records: no SS track."


# ============================================================
# Output writers
# ============================================================

def generate_pml(out_pdb, out_pml, cmap, tag, step_value, chain_id):
    num_steps = int(100 / step_value) + 1
    intervals = np.linspace(0, 100, num_steps)

    color_defs, color_names = "", []
    for i, val in enumerate(intervals):
        rgba = cmap(val / 100)
        cname = f"c_{tag}_{i}"
        color_defs += f"set_color {cname}, [{rgba[0]:.3f}, {rgba[1]:.3f}, {rgba[2]:.3f}]\n"
        color_names.append(cname)
    intervals_list = ", ".join(f"{v:.0f}" for v in intervals)
    colors_list = ", ".join(color_names)
    sel = f"conservation_struct and chain {chain_id}" if chain_id.strip() else "conservation_struct"

    with open(out_pml, "w") as f:
        f.write(f"""load {os.path.basename(out_pdb)}, conservation_struct
bg_color white
hide all
show surface, conservation_struct

{color_defs}
# chains that were not scored stay grey
color grey80, conservation_struct
# spectrum command uses spaces
spectrum b, {" ".join(color_names)}, {sel}, minimum=0, maximum=100

# ramp_new command MUST use commas
ramp_new conservation_ramp, conservation_struct, [{intervals_list}], [{colors_list}]
zoom conservation_struct
""")


def generate_cxc(out_pdb, out_cxc, cmap, step_value, chain_id):
    """ChimeraX script using exactly the same colours."""
    num_steps = int(100 / step_value) + 1
    intervals = np.linspace(0, 100, num_steps)
    palette_string = ":".join(cmap_hex(cmap, v / 100.0) for v in intervals)
    spec = f" /{chain_id}" if chain_id.strip() else ""
    with open(out_cxc, "w") as f:
        f.write(f"""# ChimeraX visualization script
open {os.path.basename(out_pdb)}
hide all
show surface
color lightgray
color byattribute bfactor{spec} palette {palette_string} range 0,100 key true
""")


def _svg_url(svg):
    return 'url("data:image/svg+xml,' + quote(svg) + '")'


SS_COLORS = {"H": "#c8414b", "G": "#9b59b6", "I": "#e67e22", "P": "#7f8c8d", "E": "#2f74c0", "B": "#2f74c0"}


def _ss_css():
    rules = []
    for k in ("H", "G", "I", "P"):
        wave = (f"<svg xmlns='http://www.w3.org/2000/svg' width='16' height='18'>"
                f"<path d='M0 9 Q4 2 8 9 T16 9' stroke='{SS_COLORS[k]}' stroke-width='2.6' fill='none'/></svg>")
        rules.append(f".ss{k}{{background:{_svg_url(wave)} repeat-x left center;}}")
    body = (f"<svg xmlns='http://www.w3.org/2000/svg' width='16' height='18'>"
            f"<rect x='0' y='6' width='16' height='6' fill='{SS_COLORS['E']}'/></svg>")
    head = (f"<svg xmlns='http://www.w3.org/2000/svg' width='16' height='18'>"
            f"<polygon points='0,2 16,9 0,16' fill='{SS_COLORS['E']}'/></svg>")
    bridge = (f"<svg xmlns='http://www.w3.org/2000/svg' width='16' height='18'>"
              f"<rect x='3' y='6' width='10' height='6' fill='{SS_COLORS['B']}'/></svg>")
    coil = ("<svg xmlns='http://www.w3.org/2000/svg' width='16' height='18'>"
            "<rect x='0' y='8.5' width='16' height='1' fill='#b8b8b8'/></svg>")
    rules += [f".ssE{{background:{_svg_url(body)} repeat-x left center;}}",
              f".ssEend{{background:{_svg_url(head)} no-repeat left center;}}",
              f".ssB{{background:{_svg_url(bridge)} no-repeat center;}}",
              f".ssC{{background:{_svg_url(coil)} repeat-x left center;}}",
              ".ssT{color:#8a8597;font-size:11px;text-align:center;font-family:Arial,sans-serif;}"]
    legend = {"wave": {k: _svg_url((
        f"<svg xmlns='http://www.w3.org/2000/svg' width='32' height='18'>"
        f"<path d='M0 9 Q4 2 8 9 T16 9 T24 9 T32 9' stroke='{SS_COLORS[k]}' stroke-width='2.6' fill='none'/></svg>"))
        for k in ("H", "G", "I")},
        "arrow": _svg_url(f"<svg xmlns='http://www.w3.org/2000/svg' width='32' height='18'>"
                          f"<rect x='0' y='6' width='18' height='6' fill='{SS_COLORS['E']}'/>"
                          f"<polygon points='16,2 32,9 16,16' fill='{SS_COLORS['E']}'/></svg>")}
    for k in ("H", "G", "I"):
        rules.append(f".lg{k}{{background:{legend['wave'][k]} no-repeat center;}}")
    rules.append(f".lgE{{background:{legend['arrow']} no-repeat center;}}")
    return "\n".join(rules), legend


def _ss_columns(alignment, ss):
    """
    One SS letter per alignment column (None outside the reference). Gap
    columns of the reference inside an element are bridged so the drawing
    stays continuous.
    """
    ref = str(alignment[0].seq)
    cols, k = [], 0
    for c in ref:
        if c == "-":
            cols.append(None)
        else:
            cols.append(ss[k] if ss is not None else None)
            k += 1
    filled = list(cols)
    for i, v in enumerate(cols):
        if v is None:
            prev = next((cols[j] for j in range(i - 1, -1, -1) if cols[j] is not None), None)
            nxt = next((cols[j] for j in range(i + 1, len(cols)) if cols[j] is not None), None)
            if prev is not None and prev == nxt and prev not in ("T", "S", "B"):
                filled[i] = prev
    return filled


def generate_alignment_html(alignment, scores, output_html, cmap, palette, metric,
                            residues, chain_id, ss, ss_source):
    """
    ESPript-like HTML alignment viewer: sticky headers and names, infinite
    horizontal scrolling, DSSP track and structure numbering above the alignment.
    """
    esc = htmllib.escape
    ss_css, _ = _ss_css()

    def score_to_hex(score):
        return cmap_hex(cmap, score / 100.0)

    def text_on(score):
        r, g, b, _ = cmap(float(min(max(score / 100.0, 0.0), 1.0)))
        lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
        return "#ffffff" if lum < 0.45 else "#000000"

    grad = ",".join(f"{score_to_hex(i)} {i}%" for i in range(101))
    ref = str(alignment[0].seq)
    ss_cols = _ss_columns(alignment, ss) if ss is not None else None

    # structure residue for every reference column
    col_res, k = [], 0
    for c in ref:
        if c == "-":
            col_res.append(None)
        else:
            col_res.append(residues[k])
            k += 1

    out = [f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Alignment viewer</title>
<style>
body{{font-family:Arial,sans-serif;margin:30px;background-color:#fdfdfd;}}
h1{{margin-bottom:5px;}}
.info{{margin-bottom:20px;line-height:1.5;}}
.legend{{width:700px;height:20px;border:1px solid #999;background:linear-gradient(to right,{grad});}}
.scale{{width:700px;display:flex;justify-content:space-between;margin-bottom:18px;font-size:12px;}}
.sslegend{{display:flex;flex-wrap:wrap;gap:18px;font-size:12px;color:#444;margin-bottom:22px;align-items:center;}}
.sslegend span.k{{display:inline-block;width:32px;height:18px;vertical-align:middle;margin-right:6px;}}
.table-container{{max-width:100%;max-height:70vh;overflow:auto;border:1px solid #ccc;background:#fff;
    margin-bottom:30px;box-shadow:0 4px 6px rgba(0,0,0,0.05);}}
table{{border-collapse:collapse;font-family:Consolas,monospace;font-size:15px;}}
td{{padding:0;}}
.name{{font-weight:bold;padding-right:10px;padding-left:10px;max-width:20ch;overflow:hidden;text-overflow:ellipsis;
    white-space:nowrap;min-width:180px;position:sticky;left:0;background-color:#fdfdfd;z-index:1;
    border-right:2px solid #ddd;}}
tr.ssrow td, tr.numrow td{{position:sticky;background-color:#fdfdfd;z-index:2;}}
tr.ssrow td{{top:0;height:20px;}}
tr.numrow td{{top:20px;border-bottom:2px solid #ddd;}}
tr.ssrow td.name, tr.numrow td.name{{z-index:3;font-weight:normal;color:#666;font-size:12px;font-family:Arial,sans-serif;}}
.num{{height:22px;font-size:11px;color:#555;min-width:16px;max-width:16px;width:16px;}}
.num span{{position:absolute;left:50%;bottom:3px;transform:translateX(-50%);white-space:nowrap;}}
.res{{width:16px;min-width:16px;height:20px;text-align:center;}}
.ss{{width:16px;min-width:16px;}}
.ref{{font-weight:bold;border-radius:2px;}}
{ss_css}
</style>
</head>
<body>
<h1>Alignment Conservation</h1>
<div class="info">
<b>Metric:</b> {esc(metric)}<br>
<b>Palette:</b> {esc(palette)}<br>
<b>Structure chain:</b> {esc(chain_id)} &nbsp;·&nbsp; <b>Secondary structure:</b> {esc(ss_source)}
</div>
<div class="legend"></div>
<div class="scale"><span>0</span><span>50</span><span>100</span></div>
"""]
    if ss_cols is not None:
        out.append(f"""<div class="sslegend">
<span><span class="k lgH"></span>α-helix (H)</span>
<span><span class="k lgG"></span>3₁₀-helix (G)</span>
<span><span class="k lgI"></span>π-helix (I)</span>
<span><span class="k lgE"></span>β-strand (E)</span>
<span><span class="k ssB"></span>β-bridge (B)</span>
<span><span class="k ssT" style="line-height:18px">T</span>turn</span>
</div>
""")

    out.append("<div class='table-container'>\n<table>")

    # Secondary-structure row
    if ss_cols is not None:
        out.append("<tr class='ssrow'><td class='name'>" + esc(ss_source.split(" (")[0]) + "</td>")
        for i, s in enumerate(ss_cols):
            r = col_res[i]
            where = f"{chain_id}:{res_label(r)} · " if r is not None else ""
            if s is None:
                out.append("<td class='ss'></td>")
                continue
            tip = esc(f"{where}{SS_NAMES.get(s, s)} ({s})", quote=True)
            if s in ("H", "G", "I", "P"):
                cls = f"ss{s}"
            elif s == "E":
                nxt = ss_cols[i + 1] if i + 1 < len(ss_cols) else None
                cls = "ssE" if nxt == "E" else "ssEend"
            elif s == "B":
                cls = "ssB"
            elif s == "T":
                out.append(f"<td class='ss ssT' title='{tip}'>T</td>")
                continue
            else:
                cls = "ssC"
            out.append(f"<td class='ss {cls}' title='{tip}'></td>")
        out.append("</tr>")

    # Numbering row (structure numbering, every 10)
    out.append("<tr class='numrow'><td class='name'>" + esc(f"Chain {chain_id}") + "</td>")
    for r in col_res:
        if r is not None and r.id[1] % 10 == 0 and not r.id[2].strip():
            out.append(f"<td class='num'><span>{r.id[1]}</span></td>")
        else:
            out.append("<td class='num'></td>")
    out.append("</tr>")

    # Alignment
    for seq_index, record in enumerate(alignment):
        out.append("<tr>")
        out.append(f"<td class='name' title='{esc(record.id, quote=True)}'>{esc(record.id[:20])}</td>")
        seq = str(record.seq)
        for i, aa in enumerate(seq):
            if aa == "-":
                out.append("<td class='res'>-</td>")
                continue
            score = scores[i]
            r = col_res[i]
            where = f"{chain_id}:{res_label(r)}" if r is not None else "insertion vs reference"
            tip = f"{where}&#10;Score {score:.1f}"
            if ss_cols is not None and ss_cols[i] is not None and r is not None:
                tip += f"&#10;{SS_NAMES.get(ss_cols[i], ss_cols[i])}"
            cls = "res ref" if seq_index == 0 else "res"
            out.append(f"<td class='{cls}' style='background:{score_to_hex(score)};color:{text_on(score)};' title='{tip}'>{esc(aa)}</td>")
        out.append("</tr>")

    out.append("</table>\n</div>\n</body>\n</html>\n")
    with open(output_html, "w", encoding="utf-8") as f:
        f.write("".join(out))


# ============================================================
# Settings (JSON) and main processing
# ============================================================

DEFAULT_SETTINGS = {
    "alignment": "", "structure": "", "chain": "auto", "output_folder": "",
    "metric": "blosum", "gap_penalty": False, "step": 10, "mkdssp": "",
    "palette": {"name": "molnympheas", "reverse": False, "lo": 0.0, "hi": 1.0},
}


def load_settings(path):
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    s = json.loads(json.dumps(DEFAULT_SETTINGS))
    for k in s:
        if k in data:
            if k == "palette" and isinstance(data[k], dict):
                s[k].update({kk: data[k][kk] for kk in s[k] if kk in data[k]})
            else:
                s[k] = data[k]
    return s


def output_base(settings):
    """(folder, base name) of the output files: <structure name>_<metric>."""
    struct = settings.get("structure", "")
    folder = settings.get("output_folder", "")
    if not folder and struct:
        folder = os.path.dirname(os.path.abspath(struct))
    stem = os.path.splitext(os.path.basename(struct))[0] if struct else "structure"
    return folder, f"{stem}_{settings.get('metric', 'blosum')}"


def process_structure(settings):
    """
    Runs the full mapping from a settings dict (see DEFAULT_SETTINGS).
    Returns a dict describing the run (output dir, chain used, SS source, notes).
    """
    s = settings
    aln_path, struct_path = s["alignment"], s["structure"]
    if not aln_path or not os.path.isfile(aln_path):
        raise MappingError("Please choose an alignment file.")
    if not struct_path or not os.path.isfile(struct_path):
        raise MappingError("Please choose a structure file.")
    output_dir, name = output_base(s)
    os.makedirs(output_dir, exist_ok=True)

    alignment, aln_format, notes = read_alignment(aln_path)
    ref_seq = str(alignment[0].seq)
    msa_seq = ungapped(ref_seq)

    structure, model = load_model0(struct_path)
    if len(structure) > 1:
        notes.append(f"The structure has {len(structure)} models: only model 0 was used.")
    chain_id, residues = pick_chain(model, msa_seq, s.get("chain", "auto"))

    scores = compute_scores(alignment, s["metric"], bool(s["gap_penalty"]))
    mapped = [scores[i] for i, aa in enumerate(ref_seq) if aa != "-"]
    for residue, score in zip(residues, mapped):
        for atom in residue:
            atom.set_bfactor(score)

    ss, ss_source, ss_warning = secondary_structure(model, chain_id, residues, struct_path, s.get("mkdssp", ""))
    if ss_warning:
        notes.append(ss_warning)

    pal = s["palette"]
    cmap = make_cmap(pal["name"], pal["reverse"], pal["lo"], pal["hi"])
    step = int(s["step"])

    base = os.path.join(output_dir, name)
    pdb_path = base + ".pdb"
    pio = PDB.PDBIO()
    pio.set_structure(model)  # model 0 only
    pio.save(pdb_path)

    with open(base + ".txt", "w") as f:
        f.write(f"# Metric: {METRIC_LABELS.get(s['metric'], s['metric'])}\n")
        f.write(f"# Gap/Insertion penalty applied: {bool(s['gap_penalty'])}\n")
        f.write(f"# Secondary structure: {ss_source}\n")
        f.write("Chain\tResID\tResName\tScore\tSS\n")
        for k, (residue, score) in enumerate(zip(residues, mapped)):
            resid = f"{residue.id[1]}{residue.id[2].strip()}"
            f.write(f"{chain_id}\t{resid}\t{residue.resname}\t{score:.3f}\t{ss[k] if ss else '-'}\n")

    generate_pml(pdb_path, base + ".pml", cmap, palette_tag(pal["name"], pal["reverse"]), step, chain_id)
    generate_cxc(pdb_path, base + ".cxc", cmap, step, chain_id)
    generate_alignment_html(alignment, scores, base + "_alignment.html", cmap,
                            palette_label(pal["name"], pal["reverse"], pal["lo"], pal["hi"]),
                            METRIC_LABELS.get(s["metric"], s["metric"]),
                            residues, chain_id, ss, ss_source)

    record = json.loads(json.dumps(s))
    record["alignment"] = os.path.abspath(aln_path)
    record["structure"] = os.path.abspath(struct_path)
    record["output_folder"] = os.path.abspath(output_dir)
    record = {"tool": "CS-to-PDB", "version": VERSION,
              "created": datetime.datetime.now().isoformat(timespec="seconds"), **record,
              "results": {"chain_used": chain_id, "alignment_format": aln_format,
                          "sequences": len(alignment), "columns": alignment.get_alignment_length(),
                          "residues_mapped": len(residues), "secondary_structure": ss_source, "notes": notes}}
    with open(base + "_settings.json", "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2, ensure_ascii=False)

    return {"dir": output_dir, "chain": chain_id, "ss_source": ss_source, "notes": notes,
            "format": aln_format, "n_seq": len(alignment)}


# ============================================================
# Qt bindings (PyMOL's Qt, PySide6, PyQt6, PyQt5, PySide2)
# ============================================================
QtCore = QtGui = QtWidgets = None
try:
    from pymol.Qt import QtCore, QtGui, QtWidgets  # inside PyMOL
except Exception:
    for _pkg in ("PySide6", "PyQt6", "PyQt5", "PySide2"):
        try:
            QtCore = importlib.import_module(_pkg + ".QtCore")
            QtGui = importlib.import_module(_pkg + ".QtGui")
            QtWidgets = importlib.import_module(_pkg + ".QtWidgets")
            break
        except ImportError:
            QtWidgets = None
if QtWidgets is None:
    raise ImportError("CS-to-PDB needs one of PySide6, PyQt6, PyQt5 or PySide2 (pip install PySide6).")

Qt = QtCore.Qt
Signal = getattr(QtCore, "pyqtSignal", None) or getattr(QtCore, "Signal")


def _e(owner, scope, name):
    """Enum lookup that works with both scoped (Qt6) and flat (Qt5) enums."""
    sc = getattr(owner, scope, None)
    if sc is not None and hasattr(sc, name):
        return getattr(sc, name)
    return getattr(owner, name)


def _exec(obj, *args):
    fn = getattr(obj, "exec", None) or getattr(obj, "exec_")
    return fn(*args)


def _pos(ev):
    return ev.position().toPoint() if hasattr(ev, "position") else ev.pos()


AlignLeft = _e(Qt, "AlignmentFlag", "AlignLeft")
AlignVCenter = _e(Qt, "AlignmentFlag", "AlignVCenter")
AlignBottom = _e(Qt, "AlignmentFlag", "AlignBottom")
AlignCenter = _e(Qt, "AlignmentFlag", "AlignCenter")
AlignHCenter = _e(Qt, "AlignmentFlag", "AlignHCenter")
AlignTop = _e(Qt, "AlignmentFlag", "AlignTop")
PointingHand = _e(Qt, "CursorShape", "PointingHandCursor")
WaitCursor = _e(Qt, "CursorShape", "WaitCursor")
NoPen = _e(Qt, "PenStyle", "NoPen")
NoBrush = _e(Qt, "BrushStyle", "NoBrush")
ElideRight = _e(Qt, "TextElideMode", "ElideRight")
Transparent = _e(Qt, "GlobalColor", "transparent")
UserRole = _e(Qt, "ItemDataRole", "UserRole")
NoItemFlags = _e(Qt, "ItemFlag", "NoItemFlags")
PopupFlag = _e(Qt, "WindowType", "Popup")
ScrollOff = _e(Qt, "ScrollBarPolicy", "ScrollBarAlwaysOff")
RichText = _e(Qt, "TextFormat", "RichText")
K = lambda n: _e(Qt, "Key", n)
Antialias = _e(QtGui.QPainter, "RenderHint", "Antialiasing")
FixedFont = _e(QtGui.QFontDatabase, "SystemFont", "FixedFont")
StateSelected = _e(QtWidgets.QStyle, "StateFlag", "State_Selected")
StateHover = _e(QtWidgets.QStyle, "StateFlag", "State_MouseOver")
SPExpanding = _e(QtWidgets.QSizePolicy, "Policy", "Expanding")
SPFixed = _e(QtWidgets.QSizePolicy, "Policy", "Fixed")
ScrollPerPixel = _e(QtWidgets.QAbstractItemView, "ScrollMode", "ScrollPerPixel")
AbsSpacing = _e(QtGui.QFont, "SpacingType", "AbsoluteSpacing")
WDemi = _e(QtGui.QFont, "Weight", "DemiBold")
WMedium = _e(QtGui.QFont, "Weight", "Medium")
MsgWarning = _e(QtWidgets.QMessageBox, "Icon", "Warning")
MsgInfo = _e(QtWidgets.QMessageBox, "Icon", "Information")
QRectF, QPointF, QSize, QPoint = QtCore.QRectF, QtCore.QPointF, QtCore.QSize, QtCore.QPoint

# ============================================================
# Design tokens (MolNympheas "Water Lily" look)
# ============================================================
WIN = "#FFFFFF"
SURFACE = "#F2EFF7"
LINE = "#E6E1EE"
INK = "#1F1B2D"
MUTED = "#6E6883"
ACCENT = "#D63F79"
ACCENT_SOFT = "#FBE6EE"
OK = "#178A5B"
WARN = "#B07800"

QSS = f"""
QMainWindow#csRoot, QWidget#csBody, QScrollArea#csScroll {{ background:{WIN}; }}
QWidget#csBody QLabel {{ color:{INK}; background:transparent; }}
QLabel[cls="lbl"] {{ color:{MUTED}; }}
QLabel[cls="hint"] {{ color:{MUTED}; font-size:11px; }}
QLabel[cls="ok"] {{ color:{OK}; font-size:11px; }}
QLabel[cls="warn"] {{ color:{WARN}; font-size:11px; }}
QLabel[cls="section"] {{ color:{INK}; font-size:13px; font-weight:600; padding-top:4px; }}
QFrame[cls="hline"] {{ background:{LINE}; max-height:1px; min-height:1px; }}
QFrame[cls="help"] {{ background:{ACCENT_SOFT}; border-radius:10px; }}
QFrame[cls="footer"] {{ border-top:1px solid {LINE}; background:{WIN}; }}
QLineEdit[cls="fld"], QComboBox[cls="fld"] {{ background:{SURFACE}; border:1px solid transparent; border-radius:8px;
    padding:4px 8px; color:{INK}; selection-background-color:{ACCENT}; min-height:20px; }}
QLineEdit[cls="fld"]:focus, QComboBox[cls="fld"]:focus {{ border-color:{ACCENT}; background:{WIN}; }}
QComboBox[cls="fld"]::drop-down {{ border:0; width:22px; }}
QComboBox QAbstractItemView {{ background:{WIN}; border:1px solid {LINE}; selection-background-color:{SURFACE};
    selection-color:{INK}; color:{INK}; outline:0; }}
QPushButton[cls="btn"] {{ border:0; border-radius:8px; background:{SURFACE}; padding:6px 12px; color:{INK}; }}
QPushButton[cls="btn"]:hover {{ background:{LINE}; }}
QPushButton[cls="btn"]:checked {{ background:{ACCENT_SOFT}; color:{ACCENT}; }}
QPushButton[cls="primary"] {{ border:0; border-radius:8px; background:{ACCENT}; padding:8px 18px; color:white;
    font-weight:600; }}
QPushButton[cls="primary"]:hover {{ background:#C2346B; }}
QPushButton[cls="link"] {{ border:0; background:transparent; color:{MUTED}; text-decoration:underline; padding:2px 4px; }}
QPushButton[cls="link"]:hover {{ color:{INK}; }}
QCheckBox {{ color:{INK}; spacing:6px; }}
QFrame[cls="popup"] {{ background:{WIN}; border:1px solid {LINE}; border-radius:10px; }}
QLineEdit[cls="search"] {{ border:0; border-bottom:1px solid {LINE}; background:{WIN}; padding:6px 8px; color:{INK}; }}
QLineEdit[cls="search"]:focus {{ border-bottom-color:{ACCENT}; }}
QListWidget {{ border:0; background:transparent; outline:0; }}
QToolTip {{ background:{INK}; color:{WIN}; border:0; padding:4px 6px; }}
QScrollBar:vertical {{ background:transparent; width:8px; margin:2px; }}
QScrollBar::handle:vertical {{ background:{LINE}; border-radius:3px; min-height:24px; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height:0; }}
"""


def _dpr():
    app = QtWidgets.QApplication.instance()
    try:
        return max(1.0, float(app.devicePixelRatio()))
    except Exception:
        return 2.0


def _new_pixmap(w, h):
    d = _dpr()
    pm = QtGui.QPixmap(int(math.ceil(w * d)), int(math.ceil(h * d)))
    pm.setDevicePixelRatio(d)
    pm.fill(Transparent)
    return pm


def gradient_pixmap(colors, w, h, hard=False, radius=3):
    pm = _new_pixmap(w, h)
    p = QtGui.QPainter(pm)
    p.setRenderHint(Antialias)
    path = QtGui.QPainterPath()
    path.addRoundedRect(QRectF(0, 0, w, h), radius, radius)
    p.setClipPath(path)
    if hard or len(colors) == 1:
        seg = w / float(len(colors))
        for i, c in enumerate(colors):
            p.fillRect(QRectF(i * seg, 0, seg + 0.6, h), QtGui.QColor(c))
    else:
        g = QtGui.QLinearGradient(0, 0, w, 0)
        for i, c in enumerate(colors):
            g.setColorAt(i / (len(colors) - 1), QtGui.QColor(c))
        p.fillRect(QRectF(0, 0, w, h), QtGui.QBrush(g))
    p.setClipping(False)
    p.setPen(QtGui.QPen(QtGui.QColor(31, 27, 45, 22), 1))
    p.setBrush(NoBrush)
    p.drawRoundedRect(QRectF(0.5, 0.5, w - 1, h - 1), radius, radius)
    p.end()
    return pm


def ui_font(px, weight=None):
    f = QtGui.QFont(QtWidgets.QApplication.font())
    f.setPixelSize(px)
    if weight is not None:
        f.setWeight(weight)
    return f


def mono_font(px):
    f = QtGui.QFontDatabase.systemFont(FixedFont)
    f.setPixelSize(px)
    return f


def mk_label(text, cls=None, wrap=False):
    lab = QtWidgets.QLabel(text)
    if cls:
        lab.setProperty("cls", cls)
    lab.setWordWrap(wrap)
    return lab


def mk_btn(text, cls="btn"):
    b = QtWidgets.QPushButton(text)
    b.setProperty("cls", cls)
    b.setCursor(PointingHand)
    return b


def mk_field(text="", placeholder=""):
    f = QtWidgets.QLineEdit(text)
    f.setProperty("cls", "fld")
    f.setPlaceholderText(placeholder)
    return f


def mk_combo(items=()):
    c = QtWidgets.QComboBox()
    c.setProperty("cls", "fld")
    for it in items:
        if isinstance(it, tuple):
            c.addItem(it[0], it[1])
        else:
            c.addItem(it, it)
    return c


def hline():
    f = QtWidgets.QFrame()
    f.setProperty("cls", "hline")
    return f


def hbox(*items, spacing=8):
    lay = QtWidgets.QHBoxLayout()
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(spacing)
    for it in items:
        if it is None:
            lay.addStretch(1)
        elif isinstance(it, int):
            lay.addSpacing(it)
        elif isinstance(it, QtWidgets.QLayout):
            lay.addLayout(it)
        else:
            lay.addWidget(it)
    return lay


def _draw_badge(p, right, cy, text):
    f = ui_font(9, WDemi)
    p.setFont(f)
    bw = QtGui.QFontMetrics(f).horizontalAdvance(text.upper()) + 10
    br = QRectF(right - bw, cy - 8, bw, 16)
    p.setPen(NoPen)
    p.setBrush(QtGui.QColor(ACCENT_SOFT))
    p.drawRoundedRect(br, 4, 4)
    p.setPen(QtGui.QColor(ACCENT))
    p.drawText(br, AlignCenter, text.upper())
    return bw


# ============================================================
# Colormap picker (ported from MolNympheas v6)
# ============================================================

class PickerDelegate(QtWidgets.QStyledItemDelegate):
    def sizeHint(self, opt, idx):
        d = idx.data(UserRole) or {}
        return QSize(100, 26 if d.get("header") else 28)

    def paint(self, p, opt, idx):
        d = idx.data(UserRole) or {}
        r = QRectF(opt.rect)
        p.save()
        p.setRenderHint(Antialias)
        if d.get("header"):
            f = ui_font(10, WDemi)
            f.setLetterSpacing(AbsSpacing, 0.8)
            p.setFont(f)
            p.setPen(QtGui.QColor(MUTED))
            p.drawText(r.adjusted(8, 0, -8, -4), AlignLeft | AlignBottom, d["name"].upper())
            p.restore()
            return
        if opt.state & StateSelected or opt.state & StateHover:
            p.setPen(NoPen)
            p.setBrush(QtGui.QColor(SURFACE))
            p.drawRoundedRect(r.adjusted(2, 1, -2, -1), 7, 7)
        x, right = r.left() + 10, r.right() - 8
        pm = d.get("pix")
        if pm is not None:
            ph = pm.height() / pm.devicePixelRatio()
            p.drawPixmap(QPointF(x, r.center().y() - ph / 2), pm)
            x += pm.width() / pm.devicePixelRatio() + 10
        for badge in reversed(d.get("badges", [])):
            right -= _draw_badge(p, right, r.center().y(), badge) + 6
        cur = d.get("current")
        f = ui_font(13, WDemi if cur else WMedium)
        p.setFont(f)
        p.setPen(QtGui.QColor(ACCENT if cur else INK))
        tw = int(right - x)
        p.drawText(QRectF(x, r.top(), tw, r.height()), AlignLeft | AlignVCenter,
                   QtGui.QFontMetrics(f).elidedText(d["name"], ElideRight, tw))
        p.restore()


class SearchEdit(QtWidgets.QLineEdit):
    nav = Signal(str)

    def keyPressEvent(self, e):
        k = e.key()
        if k in (K("Key_Down"), K("Key_Up")):
            self.nav.emit("down" if k == K("Key_Down") else "up")
        elif k in (K("Key_Return"), K("Key_Enter")):
            self.nav.emit("enter")
        elif k == K("Key_Escape"):
            self.nav.emit("esc")
        else:
            super().keyPressEvent(e)


class PickerPopup(QtWidgets.QFrame):
    chosen = Signal(str)

    def __init__(self, parent, entries_fn, search, width=440):
        super().__init__(parent, PopupFlag)
        self.setProperty("cls", "popup")
        self.setStyleSheet(QSS)
        self.entries_fn = entries_fn
        self.setFixedSize(width, 420)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)
        self.search = SearchEdit()
        self.search.setProperty("cls", "search")
        self.search.setPlaceholderText(search)
        self.search.textChanged.connect(self.rebuild)
        self.search.nav.connect(self._nav)
        lay.addWidget(self.search)
        self.list = QtWidgets.QListWidget()
        self.list.setItemDelegate(PickerDelegate(self.list))
        self.list.setMouseTracking(True)
        self.list.setVerticalScrollMode(ScrollPerPixel)
        self.list.setHorizontalScrollBarPolicy(ScrollOff)
        self.list.itemClicked.connect(self._clicked)
        self.list.itemActivated.connect(self._clicked)
        lay.addWidget(self.list)
        self.empty = mk_label("", "hint", wrap=True)
        self.empty.hide()
        lay.addWidget(self.empty)

    def open_below(self, widget):
        self.search.blockSignals(True)
        self.search.clear()
        self.search.blockSignals(False)
        self.rebuild()
        g = widget.mapToGlobal(QPoint(0, widget.height() + 4))
        scr = QtWidgets.QApplication.screenAt(g) if hasattr(QtWidgets.QApplication, "screenAt") else None
        if scr is not None:
            av = scr.availableGeometry()
            x = min(g.x(), av.right() - self.width() - 8)
            y = g.y() if g.y() + self.height() < av.bottom() else widget.mapToGlobal(QPoint(0, 0)).y() - self.height() - 4
            g = QPoint(max(av.left() + 8, x), max(av.top() + 8, y))
        self.move(g)
        self.show()
        self.search.setFocus()

    def rebuild(self):
        q = self.search.text().strip().lower()
        self.list.clear()
        cur_row, group = -1, None
        for d in self.entries_fn():
            if q and q not in (d["name"] + " " + d["group"] + " " + " ".join(d.get("badges", []))).lower():
                continue
            if d["group"] != group:
                group = d["group"]
                h = QtWidgets.QListWidgetItem()
                h.setData(UserRole, {"header": True, "name": group})
                h.setFlags(NoItemFlags)
                self.list.addItem(h)
            it = QtWidgets.QListWidgetItem()
            it.setData(UserRole, d)
            self.list.addItem(it)
            if d.get("current"):
                cur_row = self.list.count() - 1
        has = self.list.count() > 0
        self.list.setVisible(has)
        self.empty.setVisible(not has)
        self.empty.setText(f"Nothing matches “{htmllib.escape(q)}”. Try a family such as “diverging”.")
        if cur_row < 0:
            cur_row = self._step(-1, 1)
        if cur_row >= 0:
            self.list.setCurrentRow(cur_row)
            self.list.scrollToItem(self.list.item(cur_row))

    def _selectable(self, row):
        it = self.list.item(row)
        return it is not None and not (it.data(UserRole) or {}).get("header")

    def _step(self, row, d):
        r = row + d
        while 0 <= r < self.list.count():
            if self._selectable(r):
                return r
            r += d
        return row if self._selectable(row) else -1

    def _nav(self, what):
        if what == "esc":
            self.close()
        elif what == "enter":
            it = self.list.currentItem()
            if it is not None:
                self._clicked(it)
        else:
            r = self._step(self.list.currentRow(), 1 if what == "down" else -1)
            if r >= 0:
                self.list.setCurrentRow(r)

    def _clicked(self, it):
        d = it.data(UserRole) or {}
        if d.get("header"):
            return
        self.close()
        self.chosen.emit(d["id"])


class PickerButton(QtWidgets.QPushButton):
    """Flat button that shows the current colormap with its gradient and badges."""

    def __init__(self):
        super().__init__()
        self.setCursor(PointingHand)
        self.setFixedHeight(34)
        self.setMinimumWidth(260)
        self.setSizePolicy(SPExpanding, SPFixed)
        self.pix, self.name, self.suffix, self.badges = None, "", "", []

    def set_content(self, name, pix=None, suffix="", badges=()):
        self.name, self.pix, self.suffix, self.badges = name, pix, suffix, list(badges)
        self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(Antialias)
        r = QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        p.setBrush(QtGui.QColor(SURFACE))
        p.setPen(QtGui.QPen(QtGui.QColor(ACCENT if self.hasFocus() else (LINE if self.underMouse() else SURFACE)), 1))
        p.drawRoundedRect(r, 8, 8)
        x = 10.0
        if self.pix is not None:
            ph = self.pix.height() / self.pix.devicePixelRatio()
            p.drawPixmap(QPointF(x, (self.height() - ph) / 2), self.pix)
            x += self.pix.width() / self.pix.devicePixelRatio() + 10
        right = self.width() - 30.0
        for b in reversed(self.badges):
            right -= _draw_badge(p, right, self.height() / 2, b) + 6
        f = ui_font(13, WMedium)
        p.setFont(f)
        p.setPen(QtGui.QColor(INK))
        text = self.name + (("  ·  " + self.suffix) if self.suffix else "")
        tw = int(right - x)
        p.drawText(QRectF(x, 0, tw, self.height()), AlignLeft | AlignVCenter,
                   QtGui.QFontMetrics(f).elidedText(text, ElideRight, max(20, tw)))
        # chevron
        p.setPen(QtGui.QPen(QtGui.QColor(MUTED), 1.6))
        cx, cy = self.width() - 17, self.height() / 2
        p.drawPolyline(QtGui.QPolygonF([QPointF(cx - 4, cy - 2), QPointF(cx, cy + 2), QPointF(cx + 4, cy - 2)]))
        p.end()

    def enterEvent(self, e):
        self.update()
        super().enterEvent(e)

    def leaveEvent(self, e):
        self.update()
        super().leaveEvent(e)


class RangeSlider(QtWidgets.QWidget):
    """Two-handle slider drawn over the colormap gradient. Double-click resets."""
    changed = Signal(float, float)

    def __init__(self):
        super().__init__()
        self.setFixedHeight(22)
        self.setMinimumWidth(200)
        self.a, self.b = 0.0, 1.0
        self.colors, self.hard = ["#000000", "#FFFFFF"], False
        self._drag = None
        self.setCursor(PointingHand)
        self.setToolTip("Drag the handles to keep only part of the map · double-click to reset")

    def set_values(self, a, b):
        self.a, self.b = a, b
        self.update()

    def set_colors(self, colors, hard):
        self.colors, self.hard = colors, hard
        self.update()

    def _x(self, t):
        return 7 + t * (self.width() - 14)

    def _t(self, x):
        return min(max((x - 7) / float(max(1, self.width() - 14)), 0.0), 1.0)

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(Antialias)
        if not self.isEnabled():
            p.setOpacity(0.35)
        tw = self.width() - 14
        p.drawPixmap(QPoint(7, 7), gradient_pixmap(self.colors, tw, 8, self.hard, 4))
        shade = QtGui.QColor(WIN)
        shade.setAlphaF(0.72)
        p.setPen(NoPen)
        p.setBrush(shade)
        p.drawRect(QRectF(7, 7, self._x(self.a) - 7, 8))
        p.drawRect(QRectF(self._x(self.b), 7, self.width() - 7 - self._x(self.b), 8))
        if self.isEnabled():
            for t in (self.a, self.b):
                p.setBrush(QtGui.QColor(WIN))
                p.setPen(QtGui.QPen(QtGui.QColor(INK), 2))
                p.drawEllipse(QPointF(self._x(t), 11), 6, 6)
        p.end()

    def mousePressEvent(self, e):
        x = _pos(e).x()
        self._drag = "a" if abs(x - self._x(self.a)) <= abs(x - self._x(self.b)) else "b"
        self.mouseMoveEvent(e)

    def mouseMoveEvent(self, e):
        if not self._drag:
            return
        t = round(self._t(_pos(e).x()), 2)
        if self._drag == "a":
            self.a = max(0.0, min(t, self.b - 0.05))
        else:
            self.b = min(1.0, max(t, self.a + 0.05))
        self.update()
        self.changed.emit(self.a, self.b)

    def mouseReleaseEvent(self, e):
        self._drag = None

    def mouseDoubleClickEvent(self, e):
        self.a, self.b = 0.0, 1.0
        self.update()
        self.changed.emit(self.a, self.b)


class ResultBar(QtWidgets.QWidget):
    """Final map exactly as exported, with a tick at every step value."""

    def __init__(self):
        super().__init__()
        self.setFixedHeight(34)
        self.setMinimumWidth(200)
        self.colors, self.hard, self.step = ["#000000", "#FFFFFF"], False, 10

    def set_data(self, colors, hard, step):
        self.colors, self.hard, self.step = colors, hard, step
        self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(Antialias)
        w = self.width() - 14
        p.drawPixmap(QPoint(7, 2), gradient_pixmap(self.colors, w, 14, self.hard, 4))
        p.setPen(QtGui.QPen(QtGui.QColor(MUTED), 1))
        p.setFont(ui_font(10))
        for v in range(0, 101, max(1, self.step)):
            x = 7 + v / 100 * w
            p.drawLine(QPointF(x, 18), QPointF(x, 22))
            if v % 50 == 0:
                if v == 0:
                    p.drawText(QRectF(x - 3, 21, 40, 13), AlignLeft | AlignTop, str(v))
                elif v == 100:
                    p.drawText(QRectF(x - 37, 21, 40, 13), _e(Qt, "AlignmentFlag", "AlignRight") | AlignTop, str(v))
                else:
                    p.drawText(QRectF(x - 20, 21, 40, 13), AlignHCenter | AlignTop, str(v))
        p.end()


# ============================================================
# Main window
# ============================================================

HELP_TEXT = (
    "<b>BLOSUM</b>: evolutionary substitution rates (BLOSUM62). High score = frequently conserved amino acids.<br>"
    "<b>Entropy</b>: Shannon entropy of the column. High score = low diversity (highly conserved).<br>"
    "<b>Property</b>: amino acids grouped by physicochemical class (aliphatic, aromatic, polar, positive, "
    "negative, special). High score = column dominated by one class, even if the exact residue varies.<br>"
    "<b>JSD</b>: Jensen-Shannon divergence between the column and a background distribution. "
    "High score = strongly biased / conserved column.<br>"
    "<b>Rate4Site-like</b>: SIMPLIFIED PROXY only — sequence-weighted entropy (Henikoff weights), "
    "not the true ML/phylogenetic Rate4Site.<br>"
    "<b>Penalize gaps</b>: lowers each column score in proportion to its fraction of gaps.<br><br>"
    "<b>Chain · Auto</b>: uses the chain of model 0 whose sequence is identical to the first sequence of the "
    "alignment.<br>"
    "<b>Secondary structure</b>: computed with mkdssp when found (PATH, CCP4, conda, Homebrew, or the path "
    "given here); otherwise read from the HELIX/SHEET records of the file.<br>"
    "<b>Settings</b>: every run writes <i>name_metric</i>_settings.json next to the results; "
    "<i>Load settings…</i> restores it."
)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setObjectName("csRoot")
        self.setWindowTitle(f"CS-to-PDB v{VERSION}")
        self.setStyleSheet(QSS)
        self.catalogue = build_catalogue()
        self.cat_by_name = {e.name: e for e in self.catalogue}
        self.cm = dict(DEFAULT_SETTINGS["palette"])
        self._pix_cache = {}
        self._build()
        self.sync_cmap_ui()
        self.update_dssp_status()
        self.resize(780, 900)

    # ---------------------------------------------------------------- layout
    def _build(self):
        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("csScroll")
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(ScrollOff)
        scroll.setFrameShape(_e(QtWidgets.QFrame, "Shape", "NoFrame"))
        body = QtWidgets.QWidget()
        body.setObjectName("csBody")
        scroll.setWidget(body)
        grid = QtWidgets.QGridLayout(body)
        grid.setContentsMargins(20, 16, 20, 16)
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(1, 1)
        grid.setColumnMinimumWidth(0, 110)
        self.grid, self._row = grid, 0

        # Inputs
        self._section("Inputs")
        self.aln = mk_field(placeholder="Aligned FASTA or A3M — the first sequence is the reference")
        self._file_row("Alignment", self.aln, self.browse_alignment)
        self.pdb = mk_field(placeholder="PDB or mmCIF — model 0 is used")
        self._file_row("Structure", self.pdb, self.browse_structure)
        self.chain = mk_combo([("Auto (chain matching the alignment)", "auto")])
        self.chain_info = mk_label("", "hint", wrap=True)
        self._row_widgets("Chain", hbox(self.chain, self.chain_info, None))
        self.aln.editingFinished.connect(self.refresh_chains)
        self.pdb.editingFinished.connect(self.refresh_chains)

        # Output
        self._section("Output")
        self.folder = mk_field(placeholder="Default: folder of the structure")
        self._file_row("Folder", self.folder, self.browse_folder)
        self.out_hint = mk_label("", "hint", wrap=True)
        self._row_widgets("", hbox(self.out_hint))
        for w in (self.folder, self.pdb):
            w.textChanged.connect(self.update_out_hint)

        # Scoring
        self._section("Scoring")
        self.metric_cb = mk_combo([(METRIC_LABELS[k], k) for k in ("blosum", "entropy", "property", "jsd", "rate4site")])
        self.metric_cb.currentIndexChanged.connect(lambda _: self.update_out_hint())
        self.gap = QtWidgets.QCheckBox("Penalize gaps")
        self.help_btn = mk_btn("What do they mean?", "link")
        self.help_btn.setCheckable(True)
        self._row_widgets("Metric", hbox(self.metric_cb, 12, self.gap, None, self.help_btn))
        self.help = QtWidgets.QFrame()
        self.help.setProperty("cls", "help")
        hl = QtWidgets.QVBoxLayout(self.help)
        hl.setContentsMargins(12, 10, 12, 10)
        ht = mk_label(HELP_TEXT, wrap=True)
        ht.setTextFormat(RichText)
        hl.addWidget(ht)
        self.help.hide()
        self.help_btn.toggled.connect(self.help.setVisible)
        self.grid.addWidget(self.help, self._row, 0, 1, 3)
        self._row += 1

        # Palette
        self._section("Colour palette")
        self.cmap_btn = PickerButton()
        self.cmap_popup = PickerPopup(self, self.cmap_entries, search=f"Search {len(self.catalogue)} colormaps")
        self.cmap_popup.chosen.connect(self.choose_cmap)
        self.cmap_btn.clicked.connect(lambda: self.cmap_popup.open_below(self.cmap_btn))
        self.rev_btn = mk_btn("⇄  Reverse")
        self.rev_btn.setCheckable(True)
        self.rev_btn.setToolTip("Reverse the colormap")
        self.rev_btn.toggled.connect(self.toggle_reverse)
        self._row_widgets("Colormap", hbox(self.cmap_btn, self.rev_btn))
        self.range = RangeSlider()
        self.range.changed.connect(self.on_range)
        self.range_lbl = QtWidgets.QLabel("0.00 – 1.00")
        self.range_lbl.setFont(mono_font(11))
        self.range_lbl.setMinimumWidth(90)
        rl = hbox(self.range, self.range_lbl, spacing=10)
        rl.setStretch(0, 1)
        self._row_widgets("Range", rl)
        self.result = ResultBar()
        self.step = mk_combo([(f"{v} %", v) for v in (5, 10, 20, 25, 50)])
        self.step.setCurrentIndex(1)
        self.step.setFixedWidth(84)
        self.step.setToolTip("Colour stops written in the PyMOL / ChimeraX scripts")
        self.step.currentIndexChanged.connect(lambda _: self.update_result())
        rr = hbox(self.result, 6, mk_label("Step", "lbl"), self.step, spacing=10)
        rr.setStretch(0, 1)
        self._row_widgets("Result", rr)

        # Secondary structure
        self._section("Secondary structure")
        self.dssp = mk_field(placeholder="mkdssp path — leave empty to find it automatically")
        self._file_row("mkdssp", self.dssp, self.browse_dssp)
        self.dssp_status = mk_label("", "hint", wrap=True)
        self._row_widgets("", hbox(self.dssp_status))
        self.dssp.editingFinished.connect(self.update_dssp_status)

        self.grid.setRowStretch(self._row, 1)

        # Footer
        footer = QtWidgets.QFrame()
        footer.setProperty("cls", "footer")
        fl = QtWidgets.QHBoxLayout(footer)
        fl.setContentsMargins(20, 10, 20, 10)
        load_btn = mk_btn("Load settings…")
        load_btn.setToolTip("Restore a previous run from its _settings.json file")
        load_btn.clicked.connect(self.load_settings_dialog)
        self.status = mk_label("", "hint")
        run_btn = mk_btn("Run mapping", "primary")
        run_btn.clicked.connect(self.run)
        fl.addWidget(load_btn)
        fl.addSpacing(10)
        fl.addWidget(self.status, 1)
        fl.addWidget(run_btn)

        central = QtWidgets.QWidget()
        central.setObjectName("csBody")
        v = QtWidgets.QVBoxLayout(central)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)
        v.addWidget(scroll, 1)
        v.addWidget(footer)
        self.setCentralWidget(central)
        self.update_out_hint()

    def _section(self, title):
        if self._row:
            self.grid.addWidget(hline(), self._row, 0, 1, 3)
            self._row += 1
        self.grid.addWidget(mk_label(title, "section"), self._row, 0, 1, 3)
        self._row += 1

    def _row_widgets(self, label, layout):
        lab = mk_label(label, "lbl")
        self.grid.addWidget(lab, self._row, 0, AlignVCenter)
        self.grid.addLayout(layout, self._row, 1, 1, 2)
        self._row += 1

    def _file_row(self, label, field, slot):
        b = mk_btn("Browse…")
        b.clicked.connect(slot)
        self.grid.addWidget(mk_label(label, "lbl"), self._row, 0)
        self.grid.addWidget(field, self._row, 1)
        self.grid.addWidget(b, self._row, 2)
        self._row += 1

    # --------------------------------------------------------------- files
    def _start_dir(self):
        for f in (self.aln.text(), self.pdb.text()):
            if f and os.path.exists(os.path.dirname(f)):
                return os.path.dirname(f)
        return os.getcwd()

    @staticmethod
    def _first(ret):
        return ret[0] if isinstance(ret, tuple) else ret

    def browse_alignment(self):
        f = self._first(QtWidgets.QFileDialog.getOpenFileName(
            self, "Alignment", self._start_dir(), "Alignments (*.fasta *.fa *.fas *.faa *.aln *.a3m *.a2m);;All files (*)"))
        if f:
            self.aln.setText(f)
            self.refresh_chains()

    def browse_structure(self):
        f = self._first(QtWidgets.QFileDialog.getOpenFileName(
            self, "Structure", self._start_dir(), "Structures (*.pdb *.ent *.cif *.mmcif);;All files (*)"))
        if f:
            self.pdb.setText(f)
            self.refresh_chains()

    def browse_folder(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Output folder", self.folder.text() or self._start_dir())
        if d:
            self.folder.setText(d)

    def browse_dssp(self):
        f = self._first(QtWidgets.QFileDialog.getOpenFileName(self, "mkdssp executable", "/"))
        if f:
            self.dssp.setText(f)
            self.update_dssp_status()

    def update_out_hint(self):
        s = {"structure": self.pdb.text().strip(), "output_folder": self.folder.text().strip(),
             "metric": self.metric_cb.currentData() if hasattr(self, "metric_cb") else "blosum"}
        if not s["structure"] and not s["output_folder"]:
            self.out_hint.setText("Files are named after the structure and the metric (e.g. PCPE1_entropy.pdb).")
            return
        folder, name = output_base(s)
        where = "" if s["output_folder"] else f"  (in {folder})"
        self.out_hint.setText(f"→ {name}.pdb · .pml · .cxc · .txt · _alignment.html · _settings.json{where}")

    def update_dssp_status(self):
        user = self.dssp.text().strip()
        exe = find_mkdssp(user)
        if exe:
            self.dssp_status.setProperty("cls", "ok")
            self.dssp_status.setText(f"✓ DSSP: {exe}")
        else:
            self.dssp_status.setProperty("cls", "warn")
            self.dssp_status.setText(("✗ Not found at this path. " if user else "mkdssp not found. ") +
                                     "The HELIX/SHEET records of the structure file will be used instead.")
        self.dssp_status.style().unpolish(self.dssp_status)
        self.dssp_status.style().polish(self.dssp_status)

    def refresh_chains(self, wanted=None):
        """Fills the chain list from model 0 and marks the chain matching the alignment."""
        current = wanted or self.chain.currentData() or "auto"
        path = self.pdb.text().strip()
        self.chain.blockSignals(True)
        self.chain.clear()
        self.chain.addItem("Auto (chain matching the alignment)", "auto")
        self.chain_info.setText("")
        if path and os.path.isfile(path):
            try:
                chains = list_chains(path)
            except Exception as e:
                chains = []
                self.chain_info.setText(f"Cannot read structure: {e}")
            ref = None
            aln = self.aln.text().strip()
            if aln and os.path.isfile(aln):
                try:
                    ref = ungapped(str(read_alignment(aln)[0][0].seq))
                except Exception:
                    ref = None
            matches = []
            for cid, n, seq in chains:
                ok = ref is not None and seq == ref
                if ok:
                    matches.append(cid)
                self.chain.addItem(f"{cid}  ·  {n} aa" + ("  ✓ matches alignment" if ok else ""), cid)
            if ref is not None and chains:
                self.chain_info.setText(f"Auto → chain {matches[0]}" if matches else
                                        "No chain matches the alignment exactly")
        idx = self.chain.findData(current)
        self.chain.setCurrentIndex(idx if idx >= 0 else 0)
        self.chain.blockSignals(False)

    # ------------------------------------------------------------ colormaps
    def _cm_pix(self, e, rev, w, h):
        key = (e.name, rev, w, h)
        if key not in self._pix_cache:
            self._pix_cache[key] = gradient_pixmap(cmap_samples(make_cmap(e.name, rev), 40, e.qual), w, h, e.qual)
        return self._pix_cache[key]

    def cmap_entries(self):
        return [{"id": e.name, "group": e.group, "name": e.name, "badges": e.badges,
                 "pix": self._cm_pix(e, self.cm["reverse"], 88, 12), "current": e.name == self.cm["name"]}
                for e in self.catalogue]

    def sync_cmap_ui(self):
        e = self.cat_by_name.get(self.cm["name"])
        if e is None:
            e = self.catalogue[0]
            self.cm["name"] = e.name
        rev = self.cm["reverse"]
        self.cmap_btn.set_content(e.name, self._cm_pix(e, rev, 130, 14),
                                  suffix="reversed" if rev else "", badges=e.badges)
        self.range.set_colors(cmap_samples(make_cmap(e.name, rev), 40, e.qual), e.qual)
        if e.qual:
            self.cm["lo"], self.cm["hi"] = 0.0, 1.0
        self.range.set_values(self.cm["lo"], self.cm["hi"])
        self.range.setEnabled(not e.qual)
        self.range.setToolTip("A qualitative set has no range: its colours are taken in order." if e.qual else
                              "Drag the handles to keep only part of the map · double-click to reset")
        self.range_lbl.setText("whole set" if e.qual else f"{self.cm['lo']:.2f} – {self.cm['hi']:.2f}")
        self.rev_btn.blockSignals(True)
        self.rev_btn.setChecked(rev)
        self.rev_btn.blockSignals(False)
        self.update_result()

    def update_result(self):
        e = self.cat_by_name.get(self.cm["name"], self.catalogue[0])
        cm = make_cmap(self.cm["name"], self.cm["reverse"], self.cm["lo"], self.cm["hi"])
        self.result.set_data(cmap_samples(cm, 64, e.qual), e.qual, int(self.step.currentData() or 10))

    def choose_cmap(self, name):
        self.cm["name"] = name
        self.sync_cmap_ui()

    def toggle_reverse(self, on):
        self.cm["reverse"] = bool(on)
        self.sync_cmap_ui()

    def on_range(self, a, b):
        self.cm["lo"], self.cm["hi"] = round(a, 2), round(b, 2)
        self.range_lbl.setText(f"{a:.2f} – {b:.2f}")
        self.update_result()

    # ------------------------------------------------------------ settings
    def current_settings(self):
        return {
            "alignment": self.aln.text().strip(), "structure": self.pdb.text().strip(),
            "chain": self.chain.currentData() or "auto", "output_folder": self.folder.text().strip(),
            "metric": self.metric_cb.currentData(), "gap_penalty": self.gap.isChecked(),
            "step": int(self.step.currentData() or 10), "mkdssp": self.dssp.text().strip(),
            "palette": dict(self.cm),
        }

    def apply_settings(self, s):
        self.aln.setText(s["alignment"])
        self.pdb.setText(s["structure"])
        self.folder.setText(s["output_folder"])
        i = self.metric_cb.findData(s["metric"])
        self.metric_cb.setCurrentIndex(max(0, i))
        self.gap.setChecked(bool(s["gap_penalty"]))
        i = self.step.findData(int(s["step"]))
        self.step.setCurrentIndex(i if i >= 0 else 1)
        self.dssp.setText(s["mkdssp"])
        self.cm = {"name": s["palette"]["name"], "reverse": bool(s["palette"]["reverse"]),
                   "lo": float(s["palette"]["lo"]), "hi": float(s["palette"]["hi"])}
        self.refresh_chains(wanted=s["chain"])
        self.sync_cmap_ui()
        self.update_dssp_status()
        self.update_out_hint()

    def load_settings_dialog(self):
        f = self._first(QtWidgets.QFileDialog.getOpenFileName(
            self, "Load settings", self._start_dir(), "CS-to-PDB settings (*.json);;All files (*)"))
        if not f:
            return
        try:
            s = load_settings(f)
        except Exception as e:
            self._message(MsgWarning, "Cannot load settings", str(e))
            return
        self.apply_settings(s)
        missing = [p for p in (s["alignment"], s["structure"]) if p and not os.path.isfile(p)]
        if missing:
            self._message(MsgWarning, "Settings loaded",
                          "Settings restored, but these files no longer exist:\n\n" + "\n".join(missing))
        else:
            self.status.setText(f"Settings loaded from {os.path.basename(f)}")

    # ------------------------------------------------------------------ run
    def _message(self, icon, title, text, mono=False):
        box = QtWidgets.QMessageBox(self)
        box.setIcon(icon)
        box.setWindowTitle(title)
        if mono:
            box.setTextFormat(RichText)
            box.setText("<pre style='font-family:Menlo,Consolas,monospace;font-size:12px'>" +
                        htmllib.escape(text) + "</pre>")
        else:
            box.setText(text)
        _exec(box)

    def run(self):
        s = self.current_settings()
        self.status.setText("Running…")
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(WaitCursor))
        QtWidgets.QApplication.processEvents()
        try:
            info = process_structure(s)
        except MappingError as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            self.status.setText("Failed.")
            self._message(MsgWarning, "CS-to-PDB", str(e), mono=True)
            return
        except Exception as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            self.status.setText("Failed.")
            self._message(MsgWarning, "CS-to-PDB", f"{type(e).__name__}: {e}")
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        self.status.setText(f"Done · chain {info['chain']} · {info['dir']}")
        text = (f"Files written to:\n{info['dir']}\n\n"
                f"Chain: {info['chain']}   ·   {info['n_seq']} sequences ({info['format']})\n"
                f"Secondary structure: {info['ss_source']}")
        if info["notes"]:
            text += "\n\nNotes:\n• " + "\n• ".join(info["notes"])
        self._message(MsgInfo, "CS-to-PDB", text)


_WINDOWS = []


def main():
    app = QtWidgets.QApplication.instance()
    owns_app = app is None
    if owns_app:
        app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    _WINDOWS.append(win)  # keep a reference when running inside PyMOL
    win.show()
    win.raise_()
    if owns_app:
        sys.exit(_exec(app))


if __name__ in ("__main__", "pymol"):
    main()
