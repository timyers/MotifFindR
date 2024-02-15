source("renv/activate.R")

if (file.exists("data-raw/rcc_snps.txt")) {
  rcc.snps <- readLines("data-raw/rcc_snps.txt", warn = FALSE)
} else {
  message("The file 'data-raw/rcc_snps.txt' does not exist.")
}
