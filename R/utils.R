#' Convert GO ids to GO descriptions
#'
#' @param id character vector with GO ids in the format "GO:XXXXXXX"
#' @param replace replace any non-GO terms with NA \code{TRUE} or return original values \code{FALSE}
#' @returns character vector with GO descriptions
#' @export
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term
#' @examples
#' id2term("GO:0050804")
#' id2term(c("GO:0050804","GO:0099177"))
#' id2term(c("GO:0050804", NA, "GO:0099177"))
id2term = function(id, replace=T){
  na_index = which(!is.na(id))
  id_placeholder = id[na_index]
  id_placeholder = unlist(lapply(id_placeholder, function(x){
    tryCatch({z=AnnotationDbi::Term(GO.db::GOTERM[[x]])},
             error=function(cond){z=NA})
    }))

  if(replace == F){
    id[!is.na(id_placeholder)] = id_placeholder[!is.na(id_placeholder)]
  }else{
    id[na_index] = id_placeholder
  }

  return(id)
}

#' Check format/datatype of input GO ids
#'
#' Checks if input is an enrichResult (clusterProfiler), or other object which can be converted to a data.frame
#' Checks colname of column which contains GO ids
#' @keywords internal
check_go_datatype = function(input){

  if(class(input)[1] == "character"){
    return("vector")
  }else if(class(input)[1] == "data.frame"){

    go_id_cols = colnames(input)[grep("^GO[:]", input[1,])]
    return(c("data.frame", go_id_cols))

  }else if(class(input)[1] == "enrichResult"){
    return(c("enrichResult", "ID"))
  }else{
    stop("unsure what go datatype is being provided as input... please check")
  }

}
#' Check which column name is used for FDR/pvalues from a list of common terms
#'
#' Checks for FDRs/ adjusted pvalues FIRST
#' @keywords internal
check_signif_col = function(go_results, signif_col = c("FDR", "p.adjust", "pval", "p", "p.val")){

  FDR_or_p = NA
  m = which(!is.na(match(tolower(signif_col), tolower(colnames(go_results)))))
  if(length(m) > 0){
    FDR_or_p = signif_col[m][1]
    FDR_or_p = colnames(go_results)[match(tolower(FDR_or_p), tolower(colnames(go_results)))]
  }
  return(FDR_or_p)
}

#' Check which column name is used for GO ids
#'
#' @keywords internal
#' @importFrom utils capture.output
#' @importFrom utils head
check_goid_col = function(go_results, verbose = TRUE){

  go_id_cols = colnames(go_results)[grep("^GO[:]", go_results[1,])]
  if(length(go_id_cols)>1){
    if(verbose) message("multiple columns with GO ids found... ")
    if(verbose) message(paste(utils::capture.output(utils::head(go_results[,go_id_cols], 4)), collapse="\n"))
    if(any(go_id_cols == "ID")){
      if(verbose) message("using ", '"', "ID", '"')
      go_id_cols = "ID"
    }else if(any(go_id_cols == "go_id")){
      if(verbose) message("using ", '"', "go_id", '"')
      go_id_cols = "go_id"
    }else{
      if(verbose) message("using first column ", '"', go_id_cols[1], '"')
      if(verbose) message('if you want to force use go ids from a particular column, please use the column name "ID" or "go_id"')
    }
  }else if(length(go_id_cols) == 0){
    stop("no columns containing GO ids found.")
  }
  return(go_id_cols[1])
}

#' Makes a vector containing "top level" GO terms,
#'
#' Used to exclude when generating parent terms as ALL terms will be under the umbrella of
#' "BP/MF/CC" etc and this creates groups which are too large and non-specific
#'
#' @keywords internal
#' @import GOfuncR
#' @import dplyr

get_top_level_goterms = function(max_parents = 2, max_from_top = 3, ignore_terms = NULL, max_children=NULL){

  ### INTERNAL EXAMPLES
  # get_top_level_goterms()
  #
  # # using ignore_terms to cut out terms above a term of(dis)interest
  # CIA_complex = "GO:0097361"
  # CIA_parents = GOfuncR::get_parent_nodes(CIA_complex)
  # top_terms = get_top_level_goterms(ignore_terms=CIA_complex)
  # # all terms above CIA complex are added to the top_level_terms
  # all(CIA_parents$parent_go_id %in% top_terms)


  children = GOfuncR::get_child_nodes(c("GO:0008150", "GO:0003674", "GO:0005575"))

  top_level_terms =
    children %>%
    filter(.data$distance < max_from_top) %>% dplyr::select(.data$child_go_id)
  top_level_terms = top_level_terms$child_go_id

  if(!is.null(max_children)){
    big_parents =
      GOfuncR::get_child_nodes(children$child_go_id) %>%
      filter(distance>0) %>%
      count(parent_go_id) %>%
      filter(n>=max_children) %>%
      left_join(children, by=c('parent_go_id' = 'child_go_id')) %>% select(c("parent_go_id", "n", "child_name", "distance"))
    top_level_terms = unique(c(top_level_terms, big_parents$parent_go_id))
  }


  if(!is.null(ignore_terms)){
    parents = GOfuncR::get_parent_nodes(ignore_terms)
    top_level_terms = unique(c(top_level_terms, parents$parent_go_id))
  }
  return(top_level_terms)

}

#' Sort a GO_result with groups (in parent_description) by "significance"
#' and return the top n parent terms
#' Used to prevent plotting/further analysis of groups which are less significant.
#' @keywords internal
#' @import dplyr
get_top_parent_descriptions = function(go_results, n_top=10, group_col = "parent_description", ...){

  sort_by = check_signif_col(go_results, ...)

  colnames(go_results)[which(colnames(go_results)==group_col)] = "GROUP_COLUMN"

  if(!is.na(sort_by) & is.numeric(n_top)){
    x = go_results[order(go_results[,sort_by]),]%>%
      filter(!duplicated(.data$GROUP_COLUMN)) %>%
      head(n_top) %>%
      dplyr::select(.data$GROUP_COLUMN)
    return(x$GROUP_COLUMN)

  }else{
    return(unique(go_results$GROUP_COLUMN))
  }

}

#' apply str_split to string vectors
#' @param string_vector vector of strings
#' @param split_by character delimiter to split strings by
#' @param position which position of the split to return
#' @return vector of strings split by the split_by character at the position
#' @export
#' @import stringr
#' @author Brandon Signal
#' @examples
#' str_split_vector(c("ABC_DEF_HIJ_KLM", "NOP_QRS_TUV_WXY"), "_" , 2)
#'
#' str_split_vector(c("ABC_DEF_HIJ_KLM", "NOP_QRS_TUV"), "_" , 2)
#'
#' str_split_vector(c("ABC_DEF_HIJ_KLM", "NOP_QRS_TUV"), "_" , -1)

str_split_vector = function(string_vector, split_by, position){

  if(position > 0){
    out = unlist(lapply(str_split(string_vector, split_by), "[[", position))
  }else{
    position = (position*-1) - 1
    out = unlist(lapply(str_split(string_vector, split_by), function(x) x[length(x) - position]))
  }

  return(out)

}


#' Recalculate GO ratios from enrichResults to a numeric variable
#' @keywords internal
calc_go_ratios = function(go_results){
  if(any(colnames(go_results) == "GeneRatio")){
    go_results$GeneRatio = as.numeric(str_split_vector(go_results$GeneRatio, "/", 1)) /
      as.numeric(str_split_vector(go_results$GeneRatio, "/", 2))
  }
  if(any(colnames(go_results) == "BgRatio")){
    go_results$BgRatio = as.numeric(str_split_vector(go_results$BgRatio, "/", 1)) /
      as.numeric(str_split_vector(go_results$BgRatio, "/", 2))
  }
  return(go_results)
}

test_var_in_data <- function(data, variable, check_null_etc=TRUE) {


  if(is.null(variable) & check_null_etc){
    return(list(TRUE, variable, variable))
  }else if(is.numeric(variable) & check_null_etc){
    return(list(TRUE, variable, variable))
  }else{
    variable_is_parseable = FALSE
    counter=1
    while(variable_is_parseable==FALSE & counter <= length(variable)){
      workingvar = variable[counter]
      variable_is_parseable = tryCatch({
        data %>% mutate(NEW = eval(rlang::parse_expr(workingvar))) %>% head
        TRUE
      }, error = function(e)
        FALSE
      )
      counter=counter+1

    }
    basevar = unlist(lapply(colnames(data), function(x) x[grep(x, workingvar)]))[1]

    return(list(variable_is_parseable, workingvar, basevar))
  }
}

