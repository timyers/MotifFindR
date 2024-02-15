# Purpose: read in a text file in the relative directory "data-raw",
# convert it to an R data file (.rda), and save.


# Reading a tab-separated file into a data frame
my_data <- read.delim("data-raw/rcc_snps.txt", header = FALSE)


# Save "my_data_ object to a file
save(my_data, file = "data-raw/rcc_snps.rda")


# Have the R project automatically load the list of rcc_snps for use
# by the `findr` script by adding the lines below to the .Rprofile file.

# if (file.exists("rcc_snps.rda")) {
#   load("rcc_snps.rda")
# } else {
#  message("The file 'rcc_snps.rda' does not exist.")
# }
