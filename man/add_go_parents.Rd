% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_go_parents.R
\name{add_go_parents}
\alias{add_go_parents}
\title{add parent groups to go_ids}
\usage{
add_go_parents(
  go_results,
  max_parents = 2,
  max_from_top = 3,
  ignore_terms = NULL,
  max_children = NULL
)
}
\arguments{
\item{go_results}{data.frame or vector containing GO ids}

\item{max_parents}{how many levels up from a GO term should be used to find common parents?
default = \code{2}}

\item{max_from_top}{how many levels down from top level GO terms should GO terms be excluded when searching for parent terms?
default = \code{3}
Using this option prevents grouping of GO terms into top levels terms e.g. level 1: "biological_process" level 2: "biological regulation", "behavior"}

\item{ignore_terms}{vector of any GO ids to ignore when searching for common parents.
Can be helpful if you need to force one group into two or more smaller, more specific groups.
default = \code{NULL}}

\item{max_children}{ignore all GO parent terms with at least this many 'children' terms.
default = \code{NULL}}
}
\description{
add parent groups to go_ids
}
\examples{
# vector
add_go_parents(c("GO:0070997","GO:1901214","GO:0051402","GO:0043523","GO:1901215"))
# data.frame
add_go_parents(data.frame(GOIDS = c("GO:0070997","GO:1901214","GO:0051402",
"GO:0043523","GO:1901215")))
}
