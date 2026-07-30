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
# HTML alignment viewer to compare sequence conservation with structural data.
#
#
# INPUT FILES
# -----------
# • Alignment (FASTA / A3M):
#     - Accepts standard aligned FASTA (using '-' for gaps).
#     - Accepts A3M format (common in HH-suite/AlphaFold). The script
#       automatically filters out lowercase unaligned insertions to ensure 
#       strict mapping to the reference consensus sequence.
#     - NOTE: The FIRST sequence in the alignment must exactly match the 
#       sequence of the provided 3D structure (excluding gaps/insertions).
#
# • Structure (PDB / mmCIF):
#     - Accepts standard .pdb or .cif coordinate files.
#     - The tool reads the atomic coordinates and writes the computed 0-100 
#       conservation scores directly into the B-factor column.
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
# • Maps conservation scores onto PDB or mmCIF structures.
#
# • Supports all Matplotlib colormaps, including custom palettes
#   (e.g. MolNympheas), ensuring identical coloring across all outputs.
#
# • Generates:
#     - PDB structure with conservation values stored as B-factors
#     - Residue conservation table (.txt)
#     - PyMOL visualization script (.pml)
#     - ChimeraX visualization script (.cxc)
#     - Interactive HTML alignment viewer
#
#
# OUTPUT FILES
# ------------
# • [prefix].pdb              Protein structure with mapped conservation scores
# • [prefix].txt              Residue-by-residue conservation table
# • [prefix].pml              PyMOL visualization script
# • [prefix].cxc              ChimeraX visualization script
# • [prefix]_alignment.html   Interactive conservation-colored alignment viewer
#
#
# APPLICATIONS
# ------------
# • Identification of conserved functional regions
# • Structural interpretation of sequence conservation
# • Comparative analysis of protein families
# • Preparation of publication-quality figures
# • Interactive exploration of sequence–structure relationships
#
# ==============================================================================
