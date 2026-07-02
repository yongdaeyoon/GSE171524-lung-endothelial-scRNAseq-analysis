# Required packages for this workflow

packages <- c(
  "Seurat",
  "data.table",
  "Matrix",
  "dplyr",
  "ggplot2",
  "patchwork"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(packages, library, character.only = TRUE))
