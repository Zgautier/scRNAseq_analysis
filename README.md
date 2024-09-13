### scRNAseq analysis workflow for the creation of cell type specific signatures ###

The scRNAseq analysis workflow is divided into 3 parts:
- Data pre-processing, including the quality control (QC) and the normalisation.
- Cell type markers identification.
- Curation of cell type gene signatures using TF and CCC analysis.
The workflow does not perform any clustering and annotation of cells, therefore the input datasets need to be already cell type-annotated.
The workflow takes as an input either a count matrix (with a features and a barcodes object) and a metadata file, or a cell-type annotated RData file.
The relevant input files to conduct the analysis are avaialble in the repository.
