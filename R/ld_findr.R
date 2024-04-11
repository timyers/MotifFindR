# Purpose:  Identify SNPs in high linkage disequilibrium (LD)
# with affect target genes.

# Check if the `LDlinkR` is installed, if not, install and load
if (!require("LDlinkR", character.only = TRUE)) {
  # Install the package if it's not installed
  install.packages(LDlinkR)
  # Load the package after installing
  library(LDlinkR)
} else {
  message(sprintf("Package '%s' is already installed.", "LDlinkR"))
}

# 1) Find SNPs in high LD with affect target genes
# Note: "rcc.snps" variable contains full list of 39 SNPs
LDproxy_batch(snp = rcc.snps,
              pop = "CEU",
              token = Sys.getenv("LDLINK_TOKEN"),
              append = TRUE,
              genome_build = "grch38"
             )

# Note:  rs116895223 and rs11456122 are not in the 1000G reference
# panel.  No LD results were returned.

# 2) Read in file created in step #1
# Note: file had to be moved to data/data-out/ld/
ld_snps <- read.table("data/data-out/ld/combined_query_snp_list_grch38.txt",
                      header = TRUE,
                      sep = "\t",
                      row.names = NULL)


# 3) Filter SNPs by R^2
high_ld_snps <- subset(ld_snps, R2 >=0.9)


# 4) Write dataframe "high_ld_snps" to CSV file

# Get current date and time to use for unique filename
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create length of "snp_list" to use in "file_name" below
rcc_snps_length <- length(rcc.snps)

# Create a unique file name using the current date and time
file_name <- paste0("data/data-out/ld/high_ld_snps_query", rcc_snps_length, "_", current_time, ".csv")

# Write the dataframe "snp_sequence" to a CSV file called "file_name"
write.csv(high_ld_snps, file = file_name, row.names = FALSE)
