The ConstructJASPAR2024Key.sh script downloads the set of JASPAR 2024 CORE plant motifs and prepare a key grouping them into clusters as defined by JASPAR.

01_SRA_download.sh downloads raw reads from SRA given a list of accession numbers (found in supplementary data). 
We use these raw reads to generate dozens of new short read assemblies.
02_FetchDataForMotifEnrichment.sh downloads data for motif enrichment analyses.

In addition, we downloaded existing assemblies from NCBI, Phytozome, CoGE, and Figshare (from Schulz 2023, Schulz 2025 (in prep)). 
See supplementary data for a full list of assembly sources.

Recommended citations for key tools/resources:
- Raw reads from SRA: Sayers EW, Beck J, Bolton EE, Brister JR, Chan J, Comeau DC, Connor R, DiCuccio M, Farrell CM, Feldgarden M, et al. 2024. Database resources of the National Center for Biotechnology Information. Nucleic Acids Res. 52:D33–D43.
- JASPAR database: Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, Castro-Mondragon JA, Ferenc K, Kumar V, Lemma RB, Lucas J, Chèneby J, Baranasic D, et al. 2023. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res. [Internet]. Available from: http://dx.doi.org/10.1093/nar/gkad1059
- Unmethylated regions: Crisp PA, Marand AP, Noshay JM, Zhou P, Lu Z, Schmitz RJ, Springer NM. 2020. Stable unmethylated DNA demarcates expressed genes and their cis-regulatory space in plant genomes. Proc. Natl. Acad. Sci. U. S. A. 117:23991–24000.
- Accessible chromatin regions: Lu Z, Marand AP, Ricci WA, Ethridge CL, Zhang X, Schmitz RJ. 2019. The prevalence, evolution and chromatin signatures of plant regulatory elements. Nat Plants 5:1250–1259.
- Maize scATAC data: Marand AP, Chen Z, Gallavotti A, Schmitz RJ. 2021. A cis-regulatory atlas in maize at single-cell resolution. Cell 184:3041–3055.e21.
- Maize MOA-seq data: Engelhorn J, Snodgrass SJ, Kok A, Seetharam AS, Schneider M, Kiwit T, Singh A, Banf M, Doan DTH, Khaipho-Burch M, et al. 2025. Genetic variation at transcription factor binding sites largely explains phenotypic heritability in maize. Nat. Genet.:1–10.
- Andropogoneae conserved noncoding sequences: Stitzer MC, Seetharam AS, Scheben A, Hsu S-K, Schulz AJ, AuBuchon-Elder T, El-Walid M, Ferebee TH, Hale CO, La T, et al. 2025. Extensive genome evolution distinguishes maize within a stable tribe of grasses. bioRxiv [Internet]:2025.01.22.633974. Available from: https://www.biorxiv.org/content/10.1101/2025.01.22.633974v1.abstract 

Genome assemblies:
    - NCBI: Sayers EW, Beck J, Bolton EE, Brister JR, Chan J, Comeau DC, Connor R, DiCuccio M, Farrell CM, Feldgarden M, et al. 2024. Database resources of the National Center for Biotechnology Information. Nucleic Acids Res. 52:D33–D43.
    - Phytozome: Goodstein DM, Shu S, Howson R, Neupane R, Hayes RD, Fazo J, Mitros T, Dirks W, Hellsten U, Putnam N, et al. 2012. Phytozome: a comparative platform for green plant genomics. Nucleic Acids Res. 40:D1178–D1186.
    - CoGE: Lyons E, Freeling M. 2008. How to usefully compare homologous plant genes and chromosomes as DNA sequences: How to usefully compare plant genomes. Plant J. 53:661–673.
    - PanAnd long read assemblies: Stitzer MC, Seetharam AS, Scheben A, Hsu S-K, Schulz AJ, AuBuchon-Elder T, El-Walid M, Ferebee TH, Hale CO, La T, et al. 2025. Extensive genome evolution distinguishes maize within a stable tribe of grasses. bioRxiv [Internet]:2025.01.22.633974. Available from: https://www.biorxiv.org/content/10.1101/2025.01.22.633974v1.abstract
    - PanAnd short read assemblies: 
        a. Schulz AJ, Zhai J, AuBuchon-Elder T, El-Walid M, Ferebee TH, Gilmore EH, Hufford MB, Johnson LC, Kellogg EA, La T, et al. 2023. Fishing for a reelGene: evaluating gene models with evolution and machine learning. bioRxiv [Internet]. Available from: https://www.biorxiv.org/content/10.1101/2023.09.19.558246.abstract
        b. Schulz AJ, Hsu SK, Hale CO, et al. The molecular evolution of perenniality across the grasses. Unpublished data.
    - SRA short read assemblies: this study
