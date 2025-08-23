These notebooks run cross-species association models with asreml, associating motif features with environmental PCs.
00: Calculates mono and di-nucleotide frequencies in short introns across species to use as a covariate in the global association models.
01: Runs global association models. We recommend running models from the command line using the asremlR conda environment.
02: Runs orthogroup-specific association models. We recommend running models from the command line using the asremlR conda environment.
03: Prepares data for genome browser visualization of HSF/GARP motif variation between rice and Zizania.
    Use the promoterAnnotation.yaml environment to run Anchorwave and the crossMap.yaml environment to run CrossMap.

Recommended citations for key tools:
- asreml-R: Butler D, Cullis B, Gilmour A, Gogel B. 2009. mixed models for S language environments ASReml-R reference manual ASReml estimates variance components under a general linear mixed model by residual maximum likelihood (REML). Available from: https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-R-Reference-Manual-4.2.pdf
- Anchorwave: Song B, Marco-Sola S, Moreto M, Johnson L, Buckler ES, Stitzer MC. 2022. AnchorWave: Sensitive alignment of genomes with high sequence diversity, extensive structural polymorphism, and whole-genome duplication. Proc. Natl. Acad. Sci. U. S. A. [Internet] 119. Available from: http://dx.doi.org/10.1073/pnas.2113075119
- CrossMap: Zhao H, Sun Z, Wang J, Huang H, Kocher J-P, Wang L. 2014. CrossMap: a versatile tool for coordinate conversion between genome assemblies. Bioinformatics 30:1006–1007.