# GOreparent
R package for clustering of GO terms into related parent groups

## Installation
```
devtools::install_github("signalbash/goreparent")
```
goreparent requires the Bioconductor packages AnnotationDbi, clusterProfiler, GO.db, GOfuncR, and simplifyEnrichment
If these are not already installed, please run the following:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("AnnotationDbi", "clusterProfiler", "GO.db", "GOfuncR", "simplifyEnrichment"))

```

## Usage
Please see the vignette for full usage details. 

## Basic usage
goreparent clusters GO terms into groups of directly related terms with a common parent. The extent of the search for common parents between GO terms can be customised by the variables max_parents (default=2) - which defines how many generations up the search will be performed from initial 'child' GO terms - and max_from_top (default=3) - which defines how many levels down from the top level GO terms ("biological_process", "molecular_function", "cellular_component") is considered *outside* of the search frame and should enable more specific clusters to be identified.

```
library(goreparent)

enrich_file = system.file("extdata","go_enrich_results.txt", package = "goreparent")
go_input = read.table(enrich_file, sep="\t", header=TRUE)
go_parents = add_go_groups(go_input)

```

```
plot_go_parents(go_parents)
```
[PNG LINK]
```
plot_go_parents(go_parents, collapse=FALSE, n_top=2)
```
[PNG LINK]
