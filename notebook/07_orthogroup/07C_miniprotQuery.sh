## Author: Aimee Schulz
## Last updated: 04/21/2025
## Generic MiniProt script parameters that was used at scale on SciNet. 
## This script requires an assembly fasta file and the reference/orthogroup consensus sequence peptide file that you want to use for extracting sequence for each orthogroup. 
## For SLURM scripts, contact Aimee Schulz at ajs692@cornell.edu

miniprot --gff-only -t 30 -L 0 {path/to/assembly.fa} {path/to/orthogroupSequencePeptides.fa} > assemblyOutput.gff
