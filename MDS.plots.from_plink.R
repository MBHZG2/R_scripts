#!/usr/bin/env Rscript

#!/usr/bin/env Rscript
#USAGE Rscript script.R plink_output.mds
#plink.mds file output from plink after preforming IBD analysis please read from https://zzz.bwh.harvard.edu/plink/strat.shtml
#MDS scatterplots 

# Check if command-line argument is provided
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  stop("Usage: script.R <mds_file_path>")
}

# Install and load necessary packages

required_packages <- c("scatterplot3d")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}


# Function to read plink.mds file
read_mds <- function(file_path) {
  p <- read.table(file_path, header = T)
  p$FID <- factor(p$FID)
  return(p)
}

# Function to create MDS plots
create_mds_plot <- function(data, dim1, dim2, plot_title, output_file) {
  tiff(output_file, height = 12, width = 12, units = 'cm', compression = "lzw", res = 600, pointsize = 5.5, family = "Arial")
  color_palette <- rainbow(length(levels(data$FID)))
  iris_colors <- color_palette[as.numeric(data$FID)]
  
  plot(data[, dim1], data[, dim2], col = iris_colors, pch = 19, xlab = paste(dim1, "Score"), ylab = paste(dim2, "Score"), main = plot_title)
  legend("topleft", legend = levels(data$FID), fill = color_palette, cex = 0.66)
  dev.off()
}

# Function to create 3D scatter plot
create_3d_scatter_plot <- function(data, dims, output_file) {
  tiff(output_file, height = 12, width = 12, units = 'cm', compression = "lzw", res = 600, pointsize = 5.5, family = "Arial")
  
  data$FID <- factor(data$FID)
  shapes <- 1:length(unique(data$FID))
  shapes <- shapes[as.numeric(data$FID)]
  
  color_palette <- rainbow(length(levels(data$FID)))
  iris_colors <- color_palette[as.numeric(data$FID)]
  
  scatterplot3d(data[, dims], xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2), zlim = c(-0.2, 0.2), pch = shapes, color = iris_colors)
  legend("bottomright", legend = levels(data$FID), col = color_palette, pch = shapes, cex = 0.66)
  dev.off()
}

# Read plink.mds file from command-line argument
mds_file <- commandArgs(trailingOnly = TRUE)[1]
mds_data <- read_mds(mds_file)

# Create MDS plots for C1 vs C2, C1 vs C3, C2 vs C3
create_mds_plot(mds_data, "C1", "C2", "MDS Plot: C1 vs C2", "MDS.plot.C1.C2.tiff")
create_mds_plot(mds_data, "C1", "C3", "MDS Plot: C1 vs C3", "MDS.plot.C1.C3.tiff")
create_mds_plot(mds_data, "C2", "C3", "MDS Plot: C2 vs C3", "MDS.plot.C2.C3.tiff")

# Create 3D scatter plot for C1, C2, C3 dimensions
create_3d_scatter_plot(mds_data, c("C1", "C2", "C3"), "3DScatter1.tiff")
