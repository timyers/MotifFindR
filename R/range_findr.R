# 12 April 2024
# Description: This script reads in the results generated in step 3
# of `findr` and subsequently saved as a CSV file. The resulting data
# data frame (df) is converted to a `GRanges` (gr) object so that
# additional analyses and/or visualiztion can be done without the
# need to re-run the script from the beginning.

library(GenomicRanges)
library(motifbreakR)

# Read in the data generated in step 3 of `findr` script and
# saved to an RDS file
gr_results <- readRDS("data/data-out/ld/rcc_tfbs_encodemotif_high-ld_r2-0.9_win-10K_granges_object_pval_20240412_113759.rds")


##### Examples of additional analyses #####

# Subset results by rsID
rs7794538 <- gr_results[gr_results$SNP_id == "rs7794538"]
rs7794538

# Visualize results
plotMB(results = gr_results, rsid = "rs7794538", effect = "strong", altAllele = "T")
plotMB(results = rs7794538, rsid = "rs7794538", effect = "strong", altAllele = "T")
