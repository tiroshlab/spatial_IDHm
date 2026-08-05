This repository contains R analysis code accompanying:

Hoefflin, Greenwald, Galili Darnell, Mount, et al.
Spatial analysis reveals the evolving organization of IDH-mutant glioma (2026) Cancer Cell

The scripts in this repository document the main spatial transcriptomics and spatial proteomics analyses performed for the study. They are provided as an analysis record and as a starting point for researchers who wish to adapt the approaches to related datasets.

This repository is not intended to be a standalone R package or a fully automated end-to-end pipeline. The scripts reflect the working analysis environment used for the study and may require modification of file paths and other project-specific settings before use.

The core computational framework including metaprogram generation, spatial coherence, spatial associations, and consensus interactions was developed and more extensively documented for our earlier glioblastoma study:
•	Greenwald, Galili Darnell, Hoefflin, et al. Integrative spatial analysis reveals a multi-layered organization of glioblastoma. Cell 187, 2485–2501.e26 (2024). https://doi.org/10.1016/j.cell.2024.03.029
•	Code and documentation: https://github.com/tiroshlab/Spatial_Glioma
•	Associated data resource: https://doi.org/10.5281/zenodo.12624860
The current repository applies and extends that framework to IDH-mutant gliomas.

Data availability
The data used by these scripts are available from the following repositories:
•	10x Visium data: https://zenodo.org/records/18380571
Includes IDH-mutant glioma and GBM Visium datasets, a sample list, and a sample-ID conversion table for the external GBM cohort.
•	Processed CODEX data: https://zenodo.org/records/21335411
Includes cell-level metadata, expression and Nimbus-score matrices, colocalization results, sample metadata, and the associated QuPath project.
•	CODEX imaging data: https://doi.org/10.6019/S-BIAD2840
The Visium and CODEX workflows can be used independently.

Visium spatial transcriptomics
1_Vis_PerSamp_QC_LeidenClustering.R	Per-sample quality control, dimensionality reduction, Leiden clustering, and cluster gene programs
2_Vis_PerSampleNMF.R	Per-sample non-negative matrix factorization
3_Vis_Metaprograms.R	Generation and annotation of recurrent metaprograms
4_Vis_SpotAnn_SampleComp.R	Spot annotation and sample composition analyses
5_Vis_StateCoherence.R	Spatial coherence analyses
6_Vis_SpatialRelationships.R	Colocalization, adjacency, and proximity analyses
7_Vis_ConsensusInteractions.R	Definition of recurrent consensus interactions
8_Vis_InteractionTypes.R	Classification and comparison of interaction types

CODEX
9_CODEX_Import_functions_objects.R	Data import, shared functions, annotations, and analysis objects
10_CODEX_Expression_heatmaps.R	Marker expression and cell state heatmaps
11_CODEX_Composition_analysis.R	Cell composition analyses
12_CODEX_TME_composition_per_grade.R	Tumour-microenvironment composition across grades
13_CODEX_Spatial_maps_by_grade.R	Spatial maps organized by tumour grade
14_CODEX_Colocalization.R	Cell type colocalization and neighborhood analyses
15_CODEX_Junction_analysis.R	Analyses of anatomical and tumour-state junctions

Scripts 10–15 use functions and objects initialized in script 9.

Use
1.	Clone or download this repository.
2.	Download the relevant Visium and/or CODEX data from the repositories listed above.
3.	Update project-specific file paths in the scripts to match your local directory structure.
4.	Install the R packages loaded by the scripts.
5.	Run the scripts relevant to the analysis of interest.
The scripts do not need to be run as a single uninterrupted workflow. Researchers starting from processed matrices or cell tables may begin with the corresponding downstream analysis.
Reproducibility notes
•	The repository contains the analysis code used in the study but does not include every temporary or intermediate object generated during analysis.
•	Some scripts contain hard-coded paths or settings inherited from the original computing environment and must be adapted before use.
•	Package versions and computing environments may affect results, visual appearance, or clustering outcomes.
•	Large input datasets are hosted externally.
Citation
When using this code or the associated datasets, please cite the IDH-mutant glioma study. When reusing the general spatial analysis framework, please also cite the 2024 Cell study listed above.

Support
For questions about the shared methodological framework (metaprogram generation, spatial statistics, CNA inference), please refer to https://github.com/tiroshlab/Spatial_Glioma, which contains the fuller documentation. This repository is provided primarily as an analysis record. Maintenance and user support may be limited. 
