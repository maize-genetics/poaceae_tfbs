Here're the notebooks to construct and filter the orthogroups (OG).

With 32 representative high quality long read assemblies, we constructed orthogroups using orthofinder with 00_buildHelixerOG.

Constructed OGs were filtered based on the (1) # of taxa presence (>8) and (2) max. # of copy per taxa (<100) with 01_OGFilter.

Ancestral protein sequences of the filtered OGs were reconstructed with 00_buildHelixerOG and queried against all assemblies using miniProt with 02_miniprotQuery.

MiniProt alignment results were evaluated with 03_miniProtResult_eval.
