# Prepares and aggregates nucleotide frequencies for usage in linear motif model,
# given raw output from fasta-get-markov. 

# Usage: Rscript NucFreqsToCSV.R <inputDir> <outCSVPath>

library(dplyr)
library(stringr)

getwd()
args <- commandArgs(trailingOnly = TRUE)

# Access the arguments
inputDir <- args[1]
outputFile <- args[2]

# Print the arguments to verify they are being passed correctly
cat("Input Directory:", args[1], "\n")
cat("Output File:", args[2], "\n")

#### NUCLEOTIDE FREQS FOR INTRONS < 150 bp
intronfiles.names<-list.files(inputDir,pattern="*.txt", full.names = T, include.dirs = T)
# read files and return them as items in a list()
intronfreqsList <- lapply(intronfiles.names,function(x){
  theData <- read.csv(x,header = F) %>%
    dplyr::slice(-c(1,2,3,8)) %>%
    dplyr::slice(c(1:20)) %>%
    t() 
  
  # Pull out assemblyIDs to facilitate data cleaning
  assemblyID <- basename(tools::file_path_sans_ext(x))
  
  theData <- gsub(pattern = "[^0-9.e-]", "", theData) %>%
    as.numeric() %>%
    append(assemblyID)})

# rbind data frames  into a single data frame  
intronNucFreqs <- do.call(rbind,intronfreqsList) %>%
  as.data.frame() %>%
  dplyr::rename("A" = V1, "C" = V2, "G" = V3, "T" = V4, "AA" = V5, "AC" = V6, "AG" = V7, "AT" = V8, 
         "CA" = V9, "CC" = V10,"CG" = V11, "CT" = V12, "GA" = V13, "GC" = V14, "GG" = V15, 
         "GT" = V16, "TA" = V17, "TC" = V18, "TG" = V19, "TT" = V20, "assemblyID" = V21)
# Write to csv file
cat("Writing to output file:", outputFile, "\n")
write.csv(intronNucFreqs, file = outputFile, row.names = FALSE, quote = F)
