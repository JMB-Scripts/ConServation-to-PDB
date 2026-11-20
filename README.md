# CS-to-PDB

**CS-to-PDB** is a standalone Python application that calculates residue conservation from a Multiple Sequence Alignment (MSA) and maps evolutionary information onto three-dimensional protein structures by storing conservation scores in the PDB B-factor field. The resulting structures can be visualized directly in **PyMOL** and **ChimeraX**, while an interactive **HTML alignment viewer** provides a convenient way to compare sequence conservation with structural data.

---

30-07-2026 v8
Add a3m alignement file input  

20-11-2026 v10
Add load previous session
Cosmetic in the gui 
Secondary Strucutures from the model are diplayed on the top the alignement sequence
---

## Features

- Multiple conservation metrics:
  - BLOSUM62 similarity
  - Shannon Entropy
  - Property Conservation
  - Jensen-Shannon Divergence (JSD)
  - Rate4Site-like weighted conservation
- Optional gap penalty for poorly aligned regions.
- Supports **PDB** and **mmCIF** structures.
- Maps conservation scores directly into the **B-factor** field.
- Generates visualization scripts for **PyMOL** and **ChimeraX**.
- Generates an interactive **HTML alignment viewer**.
- Supports all **Matplotlib** colormaps, including custom palettes (e.g. MolNympheas).
  
<img width="892" height="1040" alt="image" src="https://github.com/user-attachments/assets/5874cc75-3743-41e7-8053-d2906b7cc965" />


- Uses identical color palettes for:
  - Protein structures
  - HTML alignment viewer

<img width="474" height="229" alt="image" src="https://github.com/user-attachments/assets/e1dca1e9-7fb6-4a46-8ef3-d5f17f114e39" />


  - PyMOL
    
<img width="308" height="265" alt="image" src="https://github.com/user-attachments/assets/fd4b635a-de3f-40bb-9fd5-38194b063c30" />


  - ChimeraX
  - 
<img width="327" height="247" alt="image" src="https://github.com/user-attachments/assets/cec59c56-9a1c-48a4-8758-2f395290128f" />


- Simple graphical user interface.
- Cross-platform (Linux, macOS and Windows).

---

## Output files

```
Project/
│
├── protein.pdb
├── conservation.txt
├── visualization.pml
├── visualization.cxc
└── alignment.html
└── settings.json
```

| File | Description |
|------|-------------|
| **.pdb** | Protein structure with conservation scores stored in the B-factor field |
| **.txt** | Residue-by-residue conservation table |
| **.pml** | PyMOL visualization script |
| **.cxc** | ChimeraX visualization script |
| **.html** | Interactive conservation-colored alignment viewer |
| **.json** | To reload a session |


## Requirements

- Python 3.9+
- Biopython
- NumPy
- Matplotlib
- pSyde

Install the required packages with:

```bash
pip install biopython numpy matplotlib
```

## License
Do whatever you want with the scripts and keep in mind that it has been code by a stupid biochemist (for that matter me).
