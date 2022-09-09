#' Annotate GO terms into groups of related terms
#'
#' Takes enrichGO, chipenrich, or vectors of GO terms and groups them into related terms.
#' Either using "parent" grouping, which finds common parents of terms,
#' or using "cluster" grouping relying on simplifyEnrichment semantic similarity clustering.
#' Additionally optionally adds additional group descriptions based on most common related terms
#' (using both direct parent-child relationships and uncle/cousin relationships)
#' overing all GO terms within the group.
#'
#' @param go_input list of or single object (\code{data.frame}, \code{enrichResult},\code{vector}) containing GO ids.
#' @param FDR minimum FDR/adjusted p-value required for grouping of GO terms. default = \code{0.05}
#' @param n_top number of groups to add descriptive parents to.
#' This cuts down on processing time and only annotates the groups with the n-th highest minimum FDR.
#' default = \code{10}
#' @param descriptive_parent add a more descriptive annotation of groups based on common relatives?
#' @param max_parents if using method = \code{"parent"}:
#' how many levels up from a GO term should be used to find common parents?
#' default = \code{2}
#' @param max_from_top if using method = \code{"parent"}:
#' how many levels down from top level GO terms should GO terms be excluded when searching for parent terms?
#' default = \code{2}
#' Using this option prevents grouping of GO terms into top levels terms e.g. level 1: "biological_process" level 2: "biological regulation", "behavior"
#' @param ignore_terms if using method = \code{"parent"}:
#' vector of any GO ids to ignore when searching for common parents.
#' Can be helpful if you need to force one group into two or more smaller, more specific groups.
#' default = \code{NULL}
#' @param method which method to use for grouping.
#' \code{"parent"}: common parents of GO terms, going up two levels and ignoring any top level terms/terms on the ignore list.
#' \code{"cluster"}: simplifyEnrichment semantic similarity clustering.
#' Note that these terms may not have very direct parent/child relations, and non-descriptive group names,
#' therefore running with \code{descriptive_parent = TRUE} is highly recommended.
#' @param verbose run with status messages? set to \code{FALSE} to supress messages.
#' @param ... arguments to be passed to \code{add_clusters_to_go}
#' @returns data.frame with original GO terms, new groups, and parent terms.
#' Contains the same columns as any data.frame (or object able to be conevrted to a data.frame) passed in as \code{go_input}.
#' @export
#' @import clusterProfiler
#' @import simplifyEnrichment
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term
#' @importFrom stats setNames
#' @examples
#'
#' enrich_file = system.file("extdata","go_enrich_results.txt", package = "goreparent")
#' go_input = read.table(enrich_file, sep="\t", header=TRUE)
#' go_parents = add_go_groups(go_input, descriptive_parent=FALSE)
#' go_clusters = add_go_groups(go_input, descriptive_parent=FALSE, method = "cluster")
#'
add_go_groups = function(go_input, FDR=0.05, n_top=10,
                         descriptive_parent = TRUE,
                         max_parents = 2, max_from_top = 2, ignore_terms=NULL,
                         method = c("parent", "cluster"),
                         verbose=TRUE,
                         ...
                         ){

  # to get R CMD CHECK to be quiet
  .id = NULL


  method = tolower(method)
  stopifnot('method must be one of "parent" or "cluster"'=any(method %in% c("parent", "cluster")))

  if(class(go_input)[1] != "list"){
    go_input = list(go_input)
    unlist_at_end = TRUE
  }else{
    unlist_at_end = FALSE
  }

  if(!class(go_input[[1]])[1] %in% c("data.frame", "character", "enrichResult")){
    go_input = lapply(go_input, function(x) as.data.frame(x))
  }

  datatypes = lapply(go_input, check_go_datatype)
  datatype_check = unique(unlist(lapply(datatypes, "[[", 1)))
  if(length(datatype_check) > 1){
    stop("multiple formats detected in go_input:", "\n",
         paste(datatype_check, collapse = ", "), "\n",
         "please check all elements in list are the same format" )
  }
  datatypes = datatypes[[1]]

  if(datatypes[1] == "enrichResult"){
    go_input = lapply(go_input, function(x) as.data.frame(x))
    datatypes[1] = 'data.frame'
  }
  if(datatypes[1] == "data.frame"){
    signif_col = check_signif_col(go_input[[1]])

    go_input = lapply(go_input, function(x){
      x = x[x[,signif_col] < FDR,]
      x = calc_go_ratios(x)
      return(x)
    })

  }else if(datatypes[1] == "vector"){
    go_input = lapply(go_input, function(x){
      data.frame(ID=x)
    })
    datatypes = c("data.frame", "ID")
  }
  all_go_terms = unique(unlist(lapply(go_input, function(x) x[,datatypes[2]])))

  if(method[1] == "parent"){
    if(verbose) message("grouping ",length(all_go_terms), " GO terms by common parents")
    go_results = add_go_parents(all_go_terms,
                                max_parents=max_parents,
                                max_from_top=max_from_top,
                                ignore_terms=ignore_terms)

    enrich_results_list = lapply(go_input, function(x){
      left_join(x, go_results, by = stats::setNames("ID", datatypes[2]))
    })
    if(verbose) message("defined ", length(unique(go_results$parent_description)), " parental groups")

  }else if(method[1] == "cluster"){
    if(verbose) message("clustering ", length(all_go_terms), " GO terms with simplifyEnrichment")
    clusters = add_clusters_to_go(all_go_terms, ...)
    enrich_results_list = lapply(go_input, function(x){
      left_join(x, clusters, by = setNames("id", datatypes[2]))
    })
    if(verbose) message("generated ", length(unique(clusters$cluster)), " clusters")

  }

  if(descriptive_parent){

    ## merge all
    if(verbose) message("adding descriptions to GO groups")

    go_results = data.table::rbindlist(enrich_results_list, idcol=T) %>%
      as.data.frame() %>%
      rename(set = .id)
    keep_cols = c("set", datatypes[2], "parent_description")
    if(!is.na(n_top)){
      keep_cols = c(keep_cols, check_signif_col(go_results))
    }
    if(any(colnames(go_results) == "parent_id")){
      keep_cols = c(keep_cols, "parent_id")
    }
    go_results = go_results[,c(keep_cols)]

    go_results = go_results[!duplicated(go_results[,datatypes[2]]),]

    if(verbose) message("running on ", length(unique(go_results$parent_description)), " groups...")
    #make descriptive parents
    pretest = ifelse(method[1]=="parent", T, F)
    go_results = make_descriptive_parent(go_results, group_col = "parent_description" ,
                                        pretest = pretest, verbose=verbose, n_top = n_top,
                                        max_parents = max_parents, max_from_top = max_from_top, ignore_terms = ignore_terms)

    enrich_results_list = lapply(enrich_results_list, function(x){
      x$parent_description = go_results$parent_description[match(x[,datatypes[2]], go_results[,datatypes[2]])]
      return(x)
    })
    if(verbose) message("done")

  }

  if(unlist_at_end){
    return(do.call("rbind", enrich_results_list))
  }else{
    return(enrich_results_list)
  }
}


#' Add simplifyEnrichment semantic similarity clustering to GO ids
#'
#' @param goids vector of GO ids to cluster
#' @param ... arguments to be passed to simplifyEnrichment::GO_similarity and simplifyEnrichment::simplifyGO
#' @returns data.frame with GO ids and clusters
#' @export
#' @import simplifyEnrichment
#' @import dplyr
#' @importFrom dplyr %>%
#' @examples
#' enrich_file = system.file("extdata","go_enrich_results.txt", package = "goreparent")
#' go_input = read.table(enrich_file, sep="\t", header=TRUE)
#' clusters = add_clusters_to_go(go_input$ID)
add_clusters_to_go = function(goids, ...){

  # guess GO ontology... stops the check from happening later and stops the warning
  ont = simplifyEnrichment::guess_ont(goids[sample(1:length(goids), min(10, length(goids)))])

  clustered_res = simplifyEnrichment::GO_similarity(goids, ont=ont, ...) %>%
    simplifyEnrichment::simplifyGO(plot=F, ..., verbose=F)

  clustered_res = arrange(clustered_res, .data$cluster, .data$term)
  clustered_res$parent_description = paste0("cluster_", clustered_res$cluster)
  return(clustered_res[,-2])

}
