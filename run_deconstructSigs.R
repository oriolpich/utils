#!/usr/bin/env Rscript

library(deconstructSigs)
set.seed(666)

args=commandArgs(trailingOnly=TRUE)

# read the cancer type name
input_file=args[1]
signature_file=args[2]

sigs.input <- as.data.frame(t(read.table(input_file, header=T, sep="\t", 
                         stringsAsFactors = FALSE, row.names = 1)))

sigs <- read.csv(signature_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

new_sigs <- sigs[c("SBS1", "SBS2", "SBS4", "SBS5", "SBS13", "SBS17a",
                   "SBS17b", "SBS18", "SBS31", "SBS35", "SBS40",
                   # add artefactual signatures, but should be manually inspected
                  'SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48',
                  'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55',
                  'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60'
                   ) ]

new_sigs = as.data.frame(t(as.matrix(new_sigs)))
flag = 0

for (sample in row.names(sigs.input))
{
  test = whichSignatures(tumor.ref = sigs.input,
                         signatures.ref = new_sigs,
                         sample.id = sample,
                         contexts.needed = TRUE,
                         
                         # if whole genome 
                         tri.counts.method = 'default',
                         signature.cutoff = 0.06,
                         signatures.limit = 6,
  )
  
  a = test$weights  # save the weights for each signature.
  a['SSE']  = round(sqrt(sum(test$diff * test$diff)), digits = 3) 
  # append the results of each sample in to dataframe
  if (flag == 0){total = a; flag=1}
  else{total <- rbind(total, a)}
}

# prepare CSV file
myDF <- cbind(sample_id = rownames(total), total)  # assign row names
rownames(myDF) <- NULL
write.table(myDF, file= paste(basename(input_file), '.decons.out.tsv', sep =''), 
            sep="\t", col.names = TRUE, row.names=FALSE)

drops <- c("SSE")
subset = myDF[ , !(names(myDF) %in% drops)]
write.table(subset, file= paste(basename(input_file), '.exposures.tsv', sep =''), 
            sep="\t", col.names = TRUE, row.names=FALSE)
