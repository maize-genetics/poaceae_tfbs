options(scipen = 20)
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(fitdistrplus))
args=commandArgs(T)

inputFile <- args[1]
inputBGFile <- args[2]
outputFile <- gsub(pattern = '.txt', 
                   replacement = '_stat.txt', 
                   x = inputFile)

predictedTFBS <- data.table::fread(input = inputFile,
                                   sep = '\t', header = T, quote = '', stringsAsFactors = F) %>% as.data.frame()

background <- data.table::fread(input = inputBGFile,
                                sep = '\t', header = T, quote = '', stringsAsFactors = F) %>% as.data.frame()
extracted_numbers <- lapply(background$sequence_name, function(sequence) {
  # Use regmatches and gregexpr correctly
  regmatches(sequence, gregexpr("(?<=shuf_)\\d+", sequence, perl = TRUE))[[1]]
}) %>% unlist()

fit <- fitdist(as.numeric(table(extracted_numbers)), "norm")



observed_value <- nrow(predictedTFBS)

# two-tailed p-value
p_value <- 2 * pnorm(-abs(observed_value - fit$estimate["mean"]), 
                     mean = 0, 
                     sd = fit$estimate["sd"])

res <- c(inputFile, nrow(predictedTFBS),
         observed_value,
         fit$estimate["mean"], fit$estimate["sd"], p_value)
write.table(x = t(res), file = outputFile, 
            sep = '\t', quote = F, 
            row.names = F, col.names = F)
