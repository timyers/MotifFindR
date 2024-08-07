# Followed example found in vignette
# https://bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html

library(BiocManager)

BiocManager::install(c("BiocParallel", "motifStack", "BSgenome", "BiocGenerics",
                       "Biostrings", "GenomeInfoDb", "GenomicRanges", "Gviz", "S4Vectors",
                       "rtracklayer", "IRanges", "MotifDb", "BSgenome.Hsapiens.UCSC.hg38",
                       "SNPlocs.Hsapiens.dbSNP.20120608", "SNPlocs.Hsapiens.dbSNP155.GRCh38",
                       "VariantAnnotation", "matrixStats", "BiocStyle", "motifBreakR"))


library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(MotifDb)
library(motifbreakR)
library(BiocParallel)
library(motifStack)
library(Gviz)
library(VariantAnnotation)
library(BiocStyle)
library(data.table)
library(Cairo)
library(gt)


# 0) Load function below

############################################################
##### Start Function to Filter MotifList by organisms  #####

# Function to accept a vector of organisms to modify a 
# MotifList (e.g., c("Hsapiens", "Mmusculus"))
# This function can take a MotifList (e.g. MotifDb) and one   
# or more desired organisms as input and return a subset of the 
# MotifList where the organism(s) match the specified value.

filter_motif_by_organisms <- function(motif_list, organisms) {
  # Ensure that motif_list is a MotifList
  if (!inherits(motif_list, "MotifList")) {
    stop("Input must be a MotifList")
  }
  
  # Ensure organisms is a vector
  if (!is.vector(organisms)) {
    stop("Organisms input must be a vector")
  }
  
  # Create a logical vector indicating whether the organism matches any in the vector
  organism_indicator <- mcols(motif_list)$organism %in% organisms & !is.na(mcols(motif_list)$organism)
  
  # Subset the MotifList using this indicator
  filtered_motif_list <- motif_list[organism_indicator]
  
  return(filtered_motif_list)
}
############ End Function ############


## Usage ##
# Assuming your MotifList object is named `MotifDb`:
# hsapiens_mmusculus_motifs <- filter_motif_by_organisms(MotifDb, c("Hsapiens", "Mmusculus"))

############     End Function     ############
##############################################


######### Begin Script #########

# 1) Input rsID's
  # SNPs below are used for testing purpose. The full set of 39 SNPs
  # are automatically read into variable rcc.snps when project is loaded.
rcc.snps <- c('rs6466948', 'rs6466949')

# Or...

# 1.1) Read in new set of rsID's in high LD with full set of 39 SNPs.
# SNPs in high LD found using `LDlinkR::LDproxy_batch` with population CEU,
# genome build grch38, and window size of 10000 bp.

# Read in CSV file with 13 columns
rcc.snps <- read.csv(file = "data/data-out/ld/high_ld_snps_query39_win-10k_20240411_082835.csv",
                     header = TRUE)

# Create a character vector from the SNPs in high LD with original
# list of 39 SNPs that are listed in column `RS_Number`.  This will
# be the input for step 2 below.

# Create character vector
rcc.snps <- as.character(rcc.snps$RS_Number)

# Or...
# 1.2) Read in CSV file of rsIDs, one column named V1, of six
# RCC SNPs
rcc.snps <- read.csv(file = "/Users/myersta/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Rprojects_OneDrive/winter/MotifFindR/data-raw/rcc_snps_6.txt",
                     header = FALSE
                    )
rcc.snps <- as.character(rcc.snps$V1)

# Orr...
# 1.3) Use full list of 39 SNPs automatically read
# into 'rcc.snps' when project is started.

# Orrr....
# 1.3 Use subset 32 from full list of SNPs
rcc.snps <- read.csv(file = "data-raw/TW_32-SNPs_List.txt",
                     header = FALSE
)
rcc.snps <- as.character(rcc.snps$V1)

# Clean up - remove elements that are not rsID numbers,
# if necessary.
rcc.snps <- rcc.snps[grep("^rs[0-9]+$", rcc.snps)]


# 2) Retrieve info about rsID's
snps.mb <- snps.from.rsid(rsid = rcc.snps,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)


# 3) Find broken Motifs
# 3.1) Determine which Motif databases to use
# 3.1.2) Combine all databases found in 'MotifDb' R package, an annotated
# collection of DNA-binding sequence motifs.
combined_mb <- c("HOCOMOCOv10", "HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C", 
                 "HOCOMOCOv11-secondary-A", "HOCOMOCOv11-secondary-B", "HOCOMOCOv11-secondary-C", 
                 "HOCOMOCOv11-secondary-D", "JASPAR_CORE", "jaspar2022", "SwissRegulon", 
                 "cisbp_1.02", "hPDI", "stamlab", "jolma2013", "UniPROBE"
                 )

# 3.1.2.1) Filter 'MotifDb' by organism, 'Hsapiens'
hsapiens_MotifDb <- filter_motif_by_organisms(MotifDb, "Hsapiens")

# Find motifs predicted to be disrupted by rsID's stored in 'snp.mb',
# using MotifList object selected databases found in 'MotifDb' and
# combined into 'combined_mb'.
results <- motifbreakR(snpList = snps.mb, 
                       filterp = TRUE,
                       pwmList = subset(hsapiens_MotifDb, 
                                        dataSource %in% combined_mb
                       ),
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam()
                      )
# Or....
# 3.1.3) Use 'hocomoco' database for CRISPRa/i project
data("hocomoco")
results <- motifbreakR(snpList = snps.mb, 
                       filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam()
)

# Subset results by rsID
rs6466949 <- results[results$SNP_id == "rs6466949"]
rs6466949

# Calculate p-values
p_results <- calculatePvalue(results, granularity = 1e-6)
p_results


# 3.1.4) 'motifbreakR_motif' is a MotifDb object containing motif information 
# from the motif databases of hocomoco, homer, factorbook, and encodemotif.
# (These can also be loaded individually.)

data("motifbreakR_motif")

# Find motifs predicted to be disrupted by rsID's in 'snps.mb'
# using MotifList object 'motifbreakR_motif'
results <- motifbreakR(snpList = snps.mb, 
                       filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam()
                      )

# 3.1.5) 'encodemotif' is a MotifDb object containing motif information 
# from the known and discovered motifs for ENCODE TF ChIP-seq datasets.

data("encodemotif")

# Find motifs predicted to be disrupted by rsID's in 'snps.mb'
# using MotifList object 'encodemotif'
results <- motifbreakR(snpList = snps.mb, 
                       filterp = TRUE,
                       pwmList = encodemotif,
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam()
                      )


# Subset results by rsID
rs6466949 <- results[results$SNP_id == "rs6466949"]
rs6466949

# Calculate p-values
p_results <- calculatePvalue(results, granularity = 1e-6)
p_results

# 4) Visualize results
plotMB(results = results, rsid = "rs6466949", effect = "strong", altAllele = "C")
plotMB(results = results, rsid = "rs6466949", effect = "strong", altAllele = "A")

plotMB(results = results, rsid = "rs6466948", effect = "strong", altAllele = "A")
plotMB(results = results, rsid = "rs6466948", effect = "strong", altAllele = "G")
plotMB(results = results, rsid = "rs6466948", effect = "weak", altAllele = "A")

##########
# Save 'results' or 'p_results' to an RDS native R data
# format that will preserve complex data structure of
# `GRanges` objects for reloading later.
saveRDS(p_results, "data/data-out/ld/rcc_tfbs_encodemotif_high-ld_r2-0.9_win-10K_granges_object_pval_20240412_113759.rds")
saveRDS(results, "data/data-out/ld/rcc_tfbs_encodemotif_high-ld_r2-0.9_win-10K_granges_object_20240412_113821.rds")
saveRDS(p_results, "data/data-out/rcc_6-snps_hocomoco_no-ld_granges_object_pval_20240711_140121.rds")
saveRDS(p_results, "data/data-out/rcc_39-snps_hocomoco_no-ld_granges_object_pval_20240724_081716.rds")
saveRDS(p_results, "data/data-out/rcc_32-snps_hocomoco_no-ld_granges_object_pval_20240801_111804.rds")

######### End Script #########

##############################################
############ Create Results Table ############

# Turn GRanges 'results' with p-values object into data.table
# 1)
dt_p_results <- as.data.table(p_results)
# 2)
dt_rs6466949 <- as.data.table(rs6466949)


# Save 'dt_p_results' to CSV file using the `fwrite` function
# from the `data.table` package

fwrite(dt_p_results, "data/data-out/rcc_tfbs_encodemotif.csv")
# or
fwrite(dt_p_results, 
       "data/data-out/ld/rcc_tfbs_encodemotif_high-ld_r2-0.9_win-10K_20240411_155434.csv")
# or
fwrite(dt_p_results, 
       "data/data-out/rcc_39-snps_hocomoco_no-ld_granges_object_pval_20240724_081716.csv")
# or
fwrite(dt_p_results, 
       "data/data-out/rcc_32-snps_hocomoco_no-ld_granges_object_pval_20240801_111804.csv")


## Create publication ready table of 'dt_p_results' using 'gt' package

# Create publication-ready table
# 1)
gt_table <- dt_p_results |>
  gt()

# 2)
gt_rs6466949_table <- dt_rs6466949 |>
  gt()

# Add Titles & Labels
gt_table <- gt_table |>
  tab_header(
    title = md("**Renal Cell Carcinoma**"),
    # subtitle = "TF Motifs from multiple motif databases"
    # subtitle = "Known TF Motifs from Encode"
    subtitle = "TF Motifs from HOCOMOCO"
  )

# Format Columns
gt_table <- gt_table |>
  fmt_number(
    columns = c(pctRef, pctAlt, scoreRef, scoreAlt, alleleDiff, alleleEffectSize),
    decimals = 4
  )

## Styling Options

# Formatting Headers
gt_table <- gt_table |>
  tab_style(
    style = list(
      cell_fill(color = "gray"),
      cell_text(weight = "bold", color = "white")
    ),
    locations = cells_column_labels(columns = everything())
  )

# Save to HTML file
# 1)
gtsave(gt_table, "data/data-out/rcc_tfbs_encode.html")
# 2)
gtsave(gt_rs6466949_table, "rcc_rs6466949.html")
# 3)
gtsave(gt_table, "data/data-out/rcc_39-snps_hocomoco_no-ld_granges_object_pval_20240724_081716.html")


##############################################
######### Write Session Info to file #########

# Get current date and time
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create a unique file name using the current date and time
file_name <- paste0("data/SessionInfo_", current_time, ".txt")

# Capture the output of sessionInfo() and write it to a file
writeLines(capture.output(sessionInfo()), file_name)
