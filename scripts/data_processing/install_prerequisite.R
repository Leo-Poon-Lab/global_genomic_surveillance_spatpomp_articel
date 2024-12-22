# List of packages to check and install
required_packages <- c(
  "doParallel", "countrycode", "here", "lubridate", "readxl", "writexl",
  "geosphere", "doFuture", "iterators", "maps", "rnaturalearth", "rnaturalearthdata",
  "sf", "WDI", "gravity", "ggrepel", "ggpubr", "slider", "MetBrewer", "cartogram", "ggforce", "packcircles", "writexl", "shadowtext", "geofacet", "RColorBrewer", "english", "naturalsort", "scales", "ggpmisc",
  "ggbump", "pacman", "pomp", "spatPomp"
)

# Function to install missing packages
install_missing_packages <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, repos='https://cloud.r-project.org/')
  } else {
    message("All packages are already installed.")
  }
}

# Check and install required packages
install_missing_packages(required_packages)

# Other packages that are not available on CRAN
if(!require(ggflags)) devtools::install_github("rensa/ggflags")
if(!require(scico)) devtools::install_github("thomasp85/scico")
