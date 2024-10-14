# Genomic analysis of adaptive immune loci in the primate T2T project

The repository contains data and scripts for analysis of immunoglobulin (IG) and T-cell receptor (TR) loci in genomes of four great ape species (bonobo, gorilla, Sumatran orangutan, and Bornean orangutan) assembled by the [Primate T2T project](https://github.com/marbl/Primates?tab=readme-ov-file#code-availability).

## Data

- `data_primate_igtrloci_locus_config.csv`: positions of IG/TR loci in genome assemblies of four great species. 
- `data_primate_igtrloci_fasta/`: FASTA files with haplotype-resolved sequences of IG/TR loci in assemblies of four great species.
- `data_primate_gene_positions/`: positions of germline IG/TR genes across four great species.
- `data_human_t2t_gene_positions/`: positions of germline IG/TR genes in the human T2T assembly.
- `data_locus_alignments/`: positions of long non-overlapping alignment blocks computed across pairs of IG/TR loci.  
- `data_SV_block_positions/`: positions of structural variation blocks computed within IG/TR loci.
- `configs/`: config files containing paths to all data files for each of the loci as well as information about locus lengths and orientations in the assembly.

## Code
### Dependencies:
- BioPython.
- [pyGenomeViz](https://github.com/moshi4/pyGenomeViz), v.0.4.4.

### Visualization of genomic diagrams
`python visualize_genome_diagram.py locus_config.txt output_fname.png`

E.g.:
`python visualize_genome_diagram.py configs/config_IGH.txt IGH_diagram.png`

### Visualization of gene count plots
`python visualize_gene_counts.py locus_config.txt output_fname.png`

E.g.:
`python visualize_gene_counts.py configs/config_IGH.txt IGH_gene_counts.png`

## Citation
Yoo D, et al. [Complete sequencing of ape genomes.](https://www.biorxiv.org/content/10.1101/2024.07.31.605654v1) BioRxiv, 2024
