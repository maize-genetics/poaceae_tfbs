00_assembleShortReads.py uses SRA raw reads as input for a short read assembly pipeline described in Schulz et al 2023.
01-03 are QC steps to flag contaminated samples (Kraken) and any obviously mislabeled samples (by building a matK phylogeny and comparing to known matK sequences)
04_PlotTABASCO.Rmd creates BUSCO-like plots used to assess assembly completeness
05_AssemblyStats.sh: Calculates assembly contiguity statistics.

Recommended citations for key tools/resources:
- Short read assembly pipeline: Schulz AJ, Zhai J, AuBuchon-Elder T, El-Walid M, Ferebee TH, Gilmore EH, Hufford MB, Johnson LC, Kellogg EA, La T, et al. 2023. Fishing for a reelGene: evaluating gene models with evolution and machine learning. bioRxiv [Internet]. Available from: https://www.biorxiv.org/content/10.1101/2023.09.19.558246.abstract
- BOLD database: Ratnasingham S, Wei C, Chan D, Agda J, Agda J, Ballesteros-Mejia L, Boutou HA, El Bastami ZM, Ma E, Manjunath R, et al. 2024. BOLD v4: A centralized bioinformatics platform for DNA-based biodiversity data. Methods Mol. Biol. 2744:403–441.
- Kraken: Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M. 2022. Metagenome analysis using the Kraken software suite. Nat. Protoc. 17:2815–2839.
- TABASCO: Schulz AJ, Zhai J, AuBuchon-Elder T, El-Walid M, Ferebee TH, Gilmore EH, Hufford MB, Johnson LC, Kellogg EA, La T, et al. 2023. Fishing for a reelGene: evaluating gene models with evolution and machine learning. bioRxiv [Internet]. Available from: https://www.biorxiv.org/content/10.1101/2023.09.19.558246.abstract 