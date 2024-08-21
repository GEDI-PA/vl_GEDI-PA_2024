# Set CRAN mirror
#Add if statement to only install packages if not existent
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of CRAN packages to be installed
cran_packages <- c(
  "s3","foreach", "aws.s3","stringr" ,"optmatch"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Load the packages
#library("s3")
