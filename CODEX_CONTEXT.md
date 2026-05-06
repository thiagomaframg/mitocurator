# MitoCurator context

MitoCurator is a diagnostic and curation-oriented pipeline for mitochondrial genomes assembled from short reads, long reads, or hybrid data.

## Main design decisions

- The pipeline must accept FASTA or GenBank as input.
- The user must provide full/absolute paths in config.yaml.
- The initial version focuses on diagnostic mode and gene-by-gene reports.
- The pipeline should support modular commands and an all-in-one `run` mode.
- The start gene for rotation is user-defined.
- The pipeline should be general for metazoan mitogenomes.
- The genetic code must be user-defined:
  - 5 for invertebrate mitochondrial code
  - 2 for vertebrate mitochondrial code
- The pipeline should eventually compare annotations from MitoFinder, MiTFi, ARWEN, tRNAscan-SE, and BLAST/TBLASTN-based reference transfer.
- MitoFinder may require Python 2.7, so it must be treated as an external tool.
- MitoFinder execution modes:
  - python_interpreter
  - conda_env
  - wrapper
- The tool must not assume the user has a conda environment called mitofinder_py2.
- If the user has `/usr/bin/python2.7`, they should be able to configure that executable.
- The number of threads must be configurable per step.
- The report must include commands, tool paths, thread usage, and reproducibility information.

## Current v0.1-dev status

Implemented:
- config.example.yaml
- environment.yml
- CLI skeleton
- check-tools
- diagnose
- rotate
- run
- basic gene_qc.tsv
- problematic_features.tsv
- intergenic_regions.tsv
- diagnostic_summary.md

Known issue fixed manually:
- `arwen` is not available in bioconda/conda-forge in this environment; treat ARWEN as an optional external tool, not an environment dependency.

Current limitation:
- Rotation currently requires an annotated GenBank input and only searches existing features.
- FASTA input currently produces only sequence-level diagnosis unless annotation is implemented.
- Annotation comparison is not implemented yet.
- Read mapping is not implemented yet.
- MitoFinder calling is not implemented yet.

## Immediate next tasks

1. Remove ARWEN from environment.yml and document it as optional external tool.
2. Make FASTA input trigger preliminary annotation or at least reference-based BLAST/TBLASTN scan.
3. Improve gene name normalization:
   - ND6, nad6, NAD6, ND-6, NADH6 should be recognized as equivalent.
   - COI/COXI/CO1/cox1 should normalize to COX1.
   - cob/cytb/CYTB should normalize to CYTB.
4. Implement MitoFinder runner:
   - supports python_interpreter, conda_env, wrapper
   - supports `-t mitfi` and `-t arwen`
   - logs exact command
5. Implement read mapping module:
   - HiFi: minimap2 -x map-hifi
   - ONT: minimap2 -x map-ont
   - PacBio CLR: minimap2 -x map-pb
   - Illumina: minimap2 -ax sr
6. Implement coverage_by_gene.tsv from BAM.
7. Implement annotation_comparison.tsv.
8. Implement final HTML/Markdown report.
