#' add parent groups to go_ids
#' @param go_results data.frame or vector containing GO ids
#' @param max_parents how many levels up from a GO term should be used to find common parents?
#' default = \code{2}
#' @param max_from_top how many levels down from top level GO terms should GO terms be excluded when searching for parent terms?
#' default = \code{2}
#' Using this option prevents grouping of GO terms into top levels terms e.g. level 1: "biological_process" level 2: "biological regulation", "behavior"
#' @param ignore_terms vector of any GO ids to ignore when searching for common parents.
#' Can be helpful if you need to force one group into two or more smaller, more specific groups.
#' default = \code{NULL}
#' @export
#' @examples
#' # vector
#' add_go_parents(c("GO:0070997","GO:1901214","GO:0051402","GO:0043523","GO:1901215"))
#' # data.frame
#' add_go_parents(data.frame(GOIDS = c("GO:0070997","GO:1901214","GO:0051402",
#' "GO:0043523","GO:1901215")))

add_go_parents = function(go_results, max_parents = 2, max_from_top = 2, ignore_terms=NULL){

  if(class(go_results)== "character"){
    go_results = data.frame(ID=go_results)
    id_col = "ID"
  }
  if(class(go_results)=="data.frame"){
    id_col = check_goid_col(go_results)

    parents = lowest_common_parent(go_term_list = go_results[,id_col],
                                         max_parents=max_parents,
                                         max_from_top=max_from_top,
                                         ignore_terms = ignore_terms)
    parents = get_meta_parents(go_parents = parents)
    m = match(go_results[,id_col], parents$id)
    go_results$parent_id = parents$parent_id[m]
    go_results$parent_description = parents$parent_description[m]
  }

  return(go_results)
}
#' find the lowest level term that is a direct parent of (or is) two GO terms.
#' Runs for all potential combinations of GO terms given as input in \code{go_term_list}
#' @import GOfuncR
#' @param fill_terms fill in parent terms using parent terms that have already been established in the first use of the parent matching?
#' @keywords internal
#' @import GOfuncR
#' @import dplyr
#' @importFrom reshape2 melt
#' @keywords internal
lowest_common_parent = function(go_term_list, max_parents = 2, max_from_top = 2, fill_terms = T, ignore_terms = NULL){

  top_level_terms = get_top_level_goterms(max_parents, max_from_top, ignore_terms)

  go.2 = go_term_list

  parents = GOfuncR::get_parent_nodes(go.2)
  parents = parents[parents$distance < (max_parents+1),]
  parents = parents[!(parents$parent_go_id %in% top_level_terms),]
  parents = arrange(parents, .data$distance, .data$parent_go_id)

  has_nearest =
    full_join(parents, parents, by = c("parent_go_id", "parent_name")) %>%
    filter(.data$child_go_id.x!=.data$child_go_id.y) %>%
    mutate(min_dist = pmin(.data$distance.x, .data$distance.y)) %>%
    dplyr::select(.data$child_go_id.x, .data$child_go_id.y, .data$parent_go_id)
  colnames(has_nearest) = c("go1", "go2", "parent")

  if(nrow(has_nearest) > 0){
    n_offspring = data.frame(id = unique(has_nearest$parent), n_off=  NA)
    offspring = GOfuncR::get_child_nodes(n_offspring$id[!is.na(n_offspring$id)])
    offspring = table(offspring$parent_go_id) %>% as.data.frame() %>% mutate(Freq = .data$Freq -1)
    n_offspring$n_off = offspring$Freq[match(n_offspring$id, offspring$Var1)]
    has_nearest$n_off = n_offspring$n_off[match(has_nearest$parent, n_offspring$id)]


    has_nearest = has_nearest %>%
      mutate(combo = paste(pmin(.data$go1, .data$go2), pmax(.data$go1, .data$go2), sep="_")) %>%
      arrange(.data$n_off) %>%
      filter(!duplicated(.data$combo)) %>%
      dplyr::select(-.data$n_off) %>%
      filter(.data$go1 !=.data$go2) %>% arrange(.data$go1, .data$go2)
  }


  ## check for redundant children
  nearest_parent_matrix = data.frame(matrix(ncol = length(go.2), nrow=length(go.2)))
  colnames(nearest_parent_matrix) = go.2
  rownames(nearest_parent_matrix) = go.2

  nearest_parent =
    cbind(go1 = rownames(nearest_parent_matrix), nearest_parent_matrix) %>%
    reshape2::melt(id.vars = c("go1"), variable.name = "go2", value.name = "parent") %>%
    mutate(go1 = as.character(.data$go1)) %>%
    mutate(go2 = as.character(.data$go2)) %>%
    mutate(combo = paste(pmin(.data$go1, .data$go2), pmax(.data$go1, .data$go2), sep="_")) %>%
    filter(.data$go1 !=.data$go2) %>%
    filter(!duplicated(.data$combo))

  nearest_parent[match(has_nearest$combo, nearest_parent$combo),]  = has_nearest
  nearest_parent$combo = NULL


  if(fill_terms & any(is.na(nearest_parent$parent))){
    parents_top = unique(nearest_parent$parent)
    parents_top = parents_top[!is.na(parents_top)]

    parents = GOfuncR::get_parent_nodes(go.2)
    parents = parents[parents$distance != 0,]
    parents = parents[(parents$parent_go_id %in% parents_top),]


    parents = arrange(parents, .data$distance, .data$parent_go_id)


    has_nearest =
      full_join(parents, parents, by = c("parent_go_id", "parent_name")) %>%
      filter(.data$child_go_id.x!=.data$child_go_id.y) %>%
      mutate(min_dist = pmin(.data$distance.x, .data$distance.y)) %>%
      dplyr::select(.data$child_go_id.x, .data$child_go_id.y, .data$parent_go_id)
    colnames(has_nearest) = c("go1", "go2", "parent")

    if(nrow(has_nearest) > 0){

      n_offspring = data.frame(id = unique(has_nearest$parent), n_off=  NA)
      offspring = GOfuncR::get_child_nodes(n_offspring$id[!is.na(n_offspring$id)])
      offspring = table(offspring$parent_go_id) %>% as.data.frame() %>% mutate(Freq = .data$Freq -1)
      n_offspring$n_off = offspring$Freq[match(n_offspring$id, offspring$Var1)]

      has_nearest$n_off = n_offspring$n_off[match(has_nearest$parent, n_offspring$id)]

      has_nearest = has_nearest %>%
        mutate(combo = paste(pmin(.data$go1, .data$go2), pmax(.data$go1, .data$go2), sep="_")) %>%
        arrange(.data$n_off) %>%
        filter(!duplicated(.data$combo)) %>%
        dplyr::select(-.data$n_off) %>%
        filter(.data$go1 !=.data$go2) %>% arrange(.data$go1, .data$go2)

      nearest_parent = nearest_parent %>% mutate(combo = paste(pmin(.data$go1, .data$go2), pmax(.data$go1, .data$go2), sep="_"))
      blank_parent = which(is.na(nearest_parent$parent))
      m = match(has_nearest$combo, nearest_parent$combo[blank_parent])
      nearest_parent[blank_parent,][m[!is.na(m)],]  = has_nearest[!is.na(m),]
      nearest_parent$combo = NULL
    }


  }

  return(nearest_parent)

}

#' from a list of all two-child->parent relationships, merge parent terms into 'Meta-parents'
#' that cover as many of the parent terms as possible that have direct parent:child relationships with each other
#' @import GOfuncR
#' @keywords internal
#' @import GOfuncR
#' @import dplyr
#' @param go_parents data.frame with two go terms and their corresponding 'parent'
get_meta_parents = function(go_parents){

  parents = unique(go_parents$parent)
  parents = parents[!is.na(parents)]

  all_go_terms = unique(c(go_parents$go1, go_parents$go2))
  parents = c(all_go_terms[!(all_go_terms %in% unique(c(go_parents$go1[go_parents$parent %in% parents],
                                                        go_parents$go2[go_parents$parent %in% parents])))], parents)

  #id2term(parents)

  children = GOfuncR::get_child_nodes(parents)
  children = children[children$distance > 0,]

  meta_parents = parents[!(parents %in% children$child_go_id)]

  meta_children = GOfuncR::get_child_nodes(meta_parents)

  total_children = table(meta_children$parent_go_id) %>%
    as.data.frame() %>%
    rename(parent_go_id = .data$Var1, n_offspring = .data$Freq)

  gene_to_meta = meta_children[meta_children$child_go_id %in% all_go_terms,] %>%
    left_join(total_children, by='parent_go_id') %>%
    mutate(parent_description = id2term(.data$parent_go_id)) %>%
    arrange(.data$n_offspring, .data$parent_go_id) %>%
    filter(!duplicated(.data$child_go_id))


  ## if parent only has one child, reset to parent to child ID
  single_go = as.data.frame(table(gene_to_meta$parent_go_id)) %>%
    filter(.data$Freq==1) %>%
    mutate(Var1=as.character(.data$Var1)) %>%
    rename(parent_go_id = .data$Var1)
  if(nrow(single_go) > 0){
    m = which(gene_to_meta$parent_go_id %in% single_go$parent_go_id)
    gene_to_meta$parent_go_id[m] = gene_to_meta$child_go_id[m]
    gene_to_meta$parent_description[m] = gene_to_meta$child_name[m]
  }

  gene_to_meta$n_offspring = NULL
  gene_to_meta = gene_to_meta[,c(2,3,1,5)]
  colnames(gene_to_meta) = c("id", "description", "parent_id", "parent_description")
  return(gene_to_meta)
}

#' Make a more descriptive name for groups of go_ids
#'
#' Uses stepparent/children (i.e. once removed related offspring terms) to
#' redefine parent groups to include semi-related terms and therefore increase
#' the specificity of the descriptive parent.
#' @import GOfuncR
#' @import dplyr
#' @param go_results data.frame with two go terms and their corresponding group
#' @param group_col which column contains the grouping variable for go terms. default \code{"parent_description"}
#' @param group_go_id_col which column contains the grouping variable formatted as a parent GO ID for go terms. default \code{"parent_id"}
#' Not required to match a column, but allows pretest filtering if used correctly.
#' @param n_top number of groups to add descriptive parents to.
#' This cuts down on processing time and only annotates the groups with the n-th highest minimum FDR.
#' default = \code{NA}. If set to NA, will run on all groups, unless pre-testing is applicable.
#' @param use_stepchild_ratio use the true stepchild ratio, not the stepchild:total children ratio.
#' default \code{TRUE}
#' Please do not turn off unless debugging.
#' @param pretest pre-test groups of GO ids for high relatedness. default \code{TRUE}
#' This can cut down on processing time, and for currently tested data, should not accidentally exclude any groups,
#' but allows some groups to be tested for which there is not a better descriptive parent than the current term.
#' @param replace_group_col replace the data in the \code{group_col} with the original id and the new description in brackets?
#' if FALSE, will add an additional column to the output, named after the \code{group_col} with "new_" as a prefix. i.e. "new_parent_description"
#' default \code{TRUE}
#' @param verbose display progress messages? default \code{FALSE}
#' @param ... arguments to be passed to \code{get_top_level_goterms} if pretest is \code{TRUE}
#' @export
#' @examples
#'
#' go_results = add_go_parents(data.frame(GOIDS = c("GO:0048167", "GO:0061001",
#' "GO:0060291" ,"GO:0050803" ,"GO:0051963", "GO:0050807", "GO:0051965")))
#' make_descriptive_parent(go_results)

#' go_results = add_go_parents(data.frame(FDR=c(rep(0.001,5), rep(0.005, 7)),
#' ID = c("GO:0070997", "GO:1901214", "GO:0051402", "GO:0043523", "GO:1901215",
#' "GO:0048167", "GO:0061001", "GO:0060291" ,"GO:0050803" ,"GO:0051963",
#' "GO:0050807", "GO:0051965")))
#' make_descriptive_parent(go_results, n_top=1, verbose=TRUE)
#' make_descriptive_parent(go_results, verbose=TRUE, replace_group_col=FALSE)

make_descriptive_parent = function(go_results,
                                   group_col="parent_description",
                                   group_go_id_col = "parent_id",
                                   n_top=NA,
                                   use_stepchild_ratio = TRUE,
                                   pretest=TRUE,
                                   replace_group_col = TRUE,
                                   verbose=FALSE, ...){

  # get top level terms required for pretesting
  if(pretest){
    top_level_terms = get_top_level_goterms(...)
  }

  # check columns are specified
  # and if pretesting is viable
  go_id_cols = colnames(go_results)[grep("^GO[:]", go_results[1,])]
  has_parent_go_id = FALSE
  if(group_go_id_col %in% go_id_cols) has_parent_go_id = TRUE
  if(pretest & has_parent_go_id == FALSE & verbose){
    message("skipping pretesting, not able to compute without go ids for groups")
    message("please set group_go_id_col to a column with group go ids if you think this step should be able to be run")
  }
  if(has_parent_go_id){
    id_col = check_goid_col(go_results[,!(colnames(go_results) == group_go_id_col)], verbose = verbose)
  }else{
    id_col = check_goid_col(go_results, verbose = verbose)
  }
  # split into groups
  all_go_terms_split = split(go_results[,id_col], go_results[,group_col])

  # filter/order by lowest FDRs
  group_order = get_top_parent_descriptions(go_results, n_top=n_top, group_col=group_col)
  skip_desc = names(all_go_terms_split)[!names(all_go_terms_split) %in% group_order]

  all_go_terms_split = all_go_terms_split[group_order]

  if(length(skip_desc) > 0){
    parent_to_new = data.frame(parent = skip_desc, new_parent=skip_desc)
    if(verbose) {
      message("skipping following groups with low significance (set n_top higher to annotate more)...")
      skip_desc = paste0('"', skip_desc, '"')
      message(paste(skip_desc, collapse = "; "))
    }
  }else{
    parent_to_new = NULL
  }

  for(i in seq_along(all_go_terms_split)){
    has_parent_go_id.internal = has_parent_go_id
    all_go_terms = all_go_terms_split[[i]]

    # only run descriptive parenting if there is at least 2 GO terms.
    if(length(all_go_terms) <= 1){
      if(verbose) message("skipping ",'"', names(all_go_terms_split)[i],'"', ": group size <=1")
      parent_to_new = rbind(parent_to_new, data.frame(parent = names(all_go_terms_split)[i], new_parent=all_go_terms))

    # RUN DESCRIPTIVE PARENTING
    }else{
      ## check for parent go id

      if(has_parent_go_id){
        parent_id = go_results[,group_go_id_col][match(all_go_terms[1], go_results[,id_col])]
      }else{
        parent_id = names(all_go_terms_split)[i]
      }
      # get lowest common parents (for later step) and..
      # if no parent go ids, determine if there is ONE single parent term or not.
      go_ids = all_go_terms
      #go_parents = lowest_common_parent(go_ids, ...)
      go_parents = lowest_common_parent(go_ids)
      if(has_parent_go_id == F){
        meta_parents = get_meta_parents(go_parents)
        meta_parents = unique(meta_parents$parent_id)
        if(length(meta_parents) == 1){
          parent_id = meta_parents
          has_parent_go_id.internal = TRUE
        }
      }

      # now try running....
      if(verbose) message("trying group ", '"',names(all_go_terms_split)[i],'"', " with ", length(all_go_terms), " GO terms")

      # run pretesting to see if skipping is possible
      if(pretest & has_parent_go_id.internal){
        try({
          go_paths = do.call('rbind', lapply(all_go_terms, function(x) get_direct_go_path(x, parent_id, top_level_terms)))
          #message("trying..", parent_id)
          if(nrow(go_paths) > 0){
            path_table = table(go_paths$child_name) %>% as.data.frame() %>% arrange(desc(.data$Freq)) %>%
              #mutate(is_cp = ifelse(.data$Var1 %in% id2term(c(go_term_list, parent_id)), T, F)) %>%
              mutate(is_cp = ifelse(.data$Var1 %in% id2term(c(all_go_terms, parent_id)), T, F)) %>%

              filter(.data$is_cp == F)
            if(nrow(path_table) == 0){
              if(verbose) message("skipped")
              all_go_terms = ""
            }
          }else{
            if(verbose) message("skipped")
            all_go_terms = ""
          }
        }, silent = T)
      }

      # THE MEAT OF IT: run descriptive parenting with stepchildren
      if(length(all_go_terms) > 1){
        ## commented out because it's been moved up...
        #go_ids = all_go_terms
        #go_parents = lowest_common_parent(go_ids)
        #meta_parents = get_meta_parents(go_parents)
        parents = unique(go_parents$parent)
        parents = parents[!is.na(parents)]
        all_go_terms = unique(c(go_parents$go1, go_parents$go2))
        parents = c(all_go_terms[!(all_go_terms %in% unique(c(go_parents$go1[go_parents$parent %in% parents],
                                                              go_parents$go2[go_parents$parent %in% parents])))], parents)


        children_new = GOfuncR::get_child_nodes(parents)
        stepchild_new = get_stepchildren(parents)
        ## FUCKED UP HERE
        new_parent = get_top_parents_df(all_go_terms, parents, children_new, stepchild_new, stepchild_ratio = use_stepchild_ratio)

        parent_to_new = rbind(parent_to_new, data.frame(parent = names(all_go_terms_split)[i], new_parent))

        if(verbose){
          if(any(names(all_go_terms_split)[i] != new_parent)){
           message("reannotated as ",'"', gsub("[\n]", "", gsub(";", "+", new_parent)), '"')
          }
        }
      }else{
        if(has_parent_go_id == F & has_parent_go_id.internal == T){
          new_parent = parent_id
        }else{
          new_parent = names(all_go_terms_split)[i]
        }
        parent_to_new = rbind(parent_to_new, data.frame(parent = names(all_go_terms_split)[i], new_parent))
      }
    }
  }#end of i-loop
  parent_to_new$new_parent = id2term(parent_to_new$new_parent, replace = FALSE)

  m = match(go_results[,group_col], parent_to_new$parent)
  if(replace_group_col){
    check_diff = go_results[,group_col] != parent_to_new$new_parent[m]
    go_results[,group_col][check_diff] = paste0(go_results[,group_col][check_diff], "\n(", parent_to_new$new_parent[m][check_diff], ")")
  }else{
    go_results$new_group_description = parent_to_new$new_parent[m]
    colnames(go_results)[which(colnames(go_results) =="new_group_description")]=paste0("new_", group_col)
  }
  return(go_results)

}


