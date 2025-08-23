Here're the notebooks to construct and filter the orthogroups (OG).

With 34 representative high quality long read assemblies (listed in Supplementary Table 2), we constructed orthogroups using orthofinder with 00_buildHelixerOG.

Constructed OGs were filtered based on the (1) # of taxa presence (>8) and (2) max. # of copy per taxa (<100) with 01_OGFilter.

Ancestral protein sequences of the filtered OGs were reconstructed with 00_buildHelixerOG and queried against all assemblies using miniProt with 02_miniprotQuery.

MiniProt alignment results were evaluated with 03_miniProtResult_eval.

Suggested citations for key tools/resources:

- Overall orthogroup pipeline: Hsu SK, Schulz AJ, Hale CO, et al. The Genomic basis of environmental adaptation in Poaceae. Unpublished data.

- Helixer: Stiehler F, Steinborn M, Scholz S, Dey D, Weber APM, Denton AK. 2021. Helixer: cross-species gene annotation of large eukaryotic genomes using deep learning. Bioinformatics 36:5291–5298.

- OrthoFinder: Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol. 20:238.

- MiniProt: Li H. 2023. Protein-to-genome alignment with miniprot. Bioinformatics [Internet] 39. Available from: http://dx.doi.org/10.1093/bioinformatics/btad014

- phangorn: Schliep KP. 2011. phangorn: phylogenetic analysis in R. Bioinformatics 27:592–593.
