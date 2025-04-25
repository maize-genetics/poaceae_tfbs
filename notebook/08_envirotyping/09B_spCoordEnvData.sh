#!/bin/bash
# env. data
# run 12_metadataProcessing.ipynb
# species name -> coordinate
mkdir output/metadataFormalOut
Rscript src/09_envirotying/01_pulling_geo_data.R data/spNameMetadata_20240819.txt output/metadataFormalOut

Rscript src/09_envirotying/02_pulling_envData.R
