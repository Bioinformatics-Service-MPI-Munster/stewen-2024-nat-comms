# Eph-ephrin signaling couples endothelial cell sorting and arterial specification

This repository contains the code for the single cell analysis of 

Eph-ephrin signaling couples endothelial cell sorting and arterial specification

*Stewen J, Kruse K, Godoi-Filip AT, Jeong H, Adams S, Berkenfeld F, Stehling M, Red-Horse K, Adams RH, Pitulescu ME*

Abstract:

Cell segregation allows the compartmentalization of cells with similar fates during morphogenesis, which can be enhanced by cell fate plasticity in response to local molecular and biomechanical cues. Endothelial tip cells in the growing retina, which lead vessel sprouts, give rise to arterial endothelial cells and thereby mediate arterial growth. Here, we have combined cell type-specific and inducible mouse genetics, flow experiments in vitro, single cell RNA sequencing and biochemistry to show that the balance between ephrin-B2 and its receptor EphB4 is critical for arterial specification, cell sorting and arteriovenous patterning. At the molecular level, elevated ephrin-B2 function after loss of EphB4 enhances signaling responses by the Notch pathway, VEGF and the transcription factor Dach1, which is influenced by endothelial shear stress. Our findings reveal how Eph-ephrin interactions integrate cell segregation and arteriovenous specification in the vasculature, which has potential relevance for human vascular malformations caused by EPHB4 mutations.


## Dependencies

This repository depends on the `anndataview` and `scrna-tools` helper packages from this organisation. Please install them first, before attempting to run any scripts:

### Anndataview

```
git clone https://github.com/Bioinformatics-Service-MPI-Munster/anndataview.git
python -m pip install --user ./anndataview
```

### scrna-tools

```
git clone https://github.com/Bioinformatics-Service-MPI-Munster/scrna-tools.git
python -m pip install --user ./scrna-tools
```

## File content

- `mapping.sh` contains commands to extract BD Rhapsody barcodes from FASTQ files, and map them using STARsolo. You will need to change folder paths accordingly
- `preprocessing.py` is used to import, filter, annotate, and nomalise count data produced by STARsolo and create an AnnData object for downstream analysis
- `embedding_and_clustering.py` contains code to run dimensionality reduction and clustering
- `de.py` can be used to repeat pseudobulk DE analyses from the study. R is needed for the code to work
- `cell_cycle_details.py` contains code to reproduce the cell cycle detail plots from the manuscript
