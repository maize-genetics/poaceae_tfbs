00_assembleShortReads.py uses SRA raw reads as input for a short read assembly pipeline described in Schulz et al 2023.
01-03 are QC steps to flag contaminated samples (Kraken) and any obviously mislabeled samples (by building a matK phylogeny and comparing to known matK sequences)
04_PlotTABASCO.Rmd creates BUSCO-like plots used to assess assembly completeness