#!/usr/bin/env Rscript

#to run Rscript script.R "$P_wald_thrshold" ==>output tiff plots for QQ 
#plots and Manhattan pltos and snps_significant.txt file includes significant SNPs and stat.txt file 
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is supplied; otherwise, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  # Default output file
  args[2] = "out.txt"
}

# Check and install dependencies
if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman")
}
if (!requireNamespace("extrafont", quietly = TRUE)) {
  install.packages("extrafont")
}

# Load necessary libraries
library(qqman)
library(extrafont)

# Set genome-wide threshold
threshold_genome_wide <- as.numeric(args[1])

# Get the list of files in the current directory ending with .assoc.txt
assoc_files <- list.files(pattern = "\\.assoc\\.txt$")

# Generate names and names_r lists
names <- assoc_files
names_r <- gsub("\\.assoc\\.txt$", "", assoc_files)

# QQ_plot
for (i in 1:length(names)) {
  # Read data from the file
  a <- read.table(names[i], head = TRUE)

  # Rename columns for clarity
  colnames(a) <- c("CHR", "SNP", "BP", "n_miss", "beta", "se", "l_remle", "P")

  # Filter data based on CHR values
  list_chr <- seq(1, 18, 1)
  a1 <- subset(a, (CHR %in% list_chr))

  # Ensure CHR column is numeric
  a1$CHR <- as.numeric(a1$CHR)

  # Order data based on CHR and BP columns
  a1 <- a1[with(a1, order(CHR, BP)), ]

  # Create QQ plot and save as TIFF file
  tiff(paste(c(names_r[i], 'QQ.tiff'), collapse = ''), height = 8, width = 16, units = 'cm',
    compression = "lzw", res = 600, pointsize = 10, family = "Arial")
  qq(a$P, main = names_r[i])
  dev.off()
}

# Manhattan_plots
for (i in 1:length(assoc_files)) {
  a <- read.table(assoc_files[i], head = TRUE)

  # Rename columns for clarity
  x <- colnames(a)
  x[x == 'rs'] <- "SNP"
  x[x == 'chr'] <- "CHR"
  x[x == 'ps'] <- "BP"
  x[x == 'p_wald'] <- "P"
  colnames(a) <- x

  # Filter data based on CHR values
  a1 <- a[a$CHR %in% seq(1, 29, 1), ]

  # Ensure CHR column is numeric
  a1$CHR <- as.numeric(as.character(a1$CHR))

  tiff(paste(c(assoc_files[i], '.tiff'), collapse = ''), height = 8.9, width = 17.5, units = 'cm',
    compression = "lzw", res = 600, pointsize = 5.5, family = "Arial")
  manhattan(a1, col = c("black", "gray"), suggestiveline = 5.0E-05,
    genomewideline = -log10(threshold_genome_wide), ylim = c(0, 18),
    main = assoc_files[i], cex = 1.1, cex.axis = 1.1)
  dev.off()
}

# Extract significant results < args[1] p-value threshold
parameters <- list.files(pattern = '.assoc.txt$')
valid_snps <- data.frame()

for (i in 1:length(parameters)) {
  db <- read.table(parameters[i], head = TRUE)
  sb <- subset(db, p_wald < as.numeric(args[1]))
  sb$parameters <- parameters[i]
  valid_snps <- rbind(valid_snps, sb)
}

write.table(valid_snps, 'snps_significant.txt', sep = '\t')

# Compute lambda
for (i in 1:length(parameters)) {
  db <- read.table(parameters[i], head = TRUE)
  chisq <- qchisq(1 - db$p_wald, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)

  print(c(parameters[i], lambda))
}
#####
# Extract information from log files
log_files <- list.files(pattern = "\\.log\\.txt$")
# make data.frame
stat_df <- data.frame(
  file = character(),
  total_individuals = numeric(),
  analyzed_individuals = numeric(),
  covariates = numeric(),
  phenotypes = numeric(),
  total_snps = numeric(),
  analyzed_snps = numeric(),
  pve_estimate = numeric(),
  se_pve = numeric(),
  lambda = numeric()
)
# Loop through log files
for (log_file in log_files) {
  # Read data from the log file
  log_data <- readLines(log_file)

  # Extract relevant information
   total_individuals <- as.numeric(gsub(".*number of total individuals = (\\d+).*", "\\1", log_data[grep("number of total individuals =", log_data)]))
  analyzed_individuals <- as.numeric(gsub(".*number of analyzed individuals = (\\d+).*", "\\1", log_data[grep("number of analyzed individuals =", log_data)]))
  covariates <- as.numeric(gsub(".*number of covariates = (\\d+).*", "\\1", log_data[grep("number of covariates =", log_data)]))
  phenotypes <- as.numeric(gsub(".*number of phenotypes = (\\d+).*", "\\1", log_data[grep("number of phenotypes =", log_data)]))
  total_snps <- as.numeric(gsub(".*number of total SNPs/var = (\\d+).*", "\\1", log_data[grep("number of total SNPs/var =", log_data)]))
  analyzed_snps <- as.numeric(gsub(".*number of analyzed SNPs/var = (\\d+).*", "\\1", log_data[grep("number of analyzed SNPs/var =", log_data)]))
  pve_estimate <- as.numeric(gsub(".*pve estimate in the null model = (.+).*", "\\1", log_data[grep("pve estimate in the null model =", log_data)]))
  se_pve <- as.numeric(gsub(".*se\\(pve\\) in the null model = (.+).*", "\\1", log_data[grep("se\\(pve\\) in the null model =", log_data)]))
  lambda <- lambda  
  stat_df <- rbind(stat_df, c(log_file, total_individuals, analyzed_individuals, covariates, phenotypes, total_snps, analyzed_snps, 
   pve_estimate, se_pve, lambda))
  }

colnames(stat_df)=c("file","total_individuals","analyzed_individuals","covariates","phenotypes","total_snps","analyzed_snps",
"pve_estimate","se_pve","lambda")
write.table(stat_df, "stat.txt", sep = "\t", quote = FALSE, row.names = FALSE)
