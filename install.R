# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# List of CRAN packages to be installed
cran_packages <- c(
  "terra",
  "dplyr",
  "sf",
  "s3"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Load the packages
library("terra")
library("dplyr")
library("sf")
library("s3")
