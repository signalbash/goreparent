# add_go_parents_utils.R
# INTERNAL UTILITY FUNCTIONS FOR GO PARENTING

# Old very slow stepchild finding function
# get_stepchildren_old = function(go_terms){
#
#   #inclusion body assembly (parent) -> regulation of inclusion body assembly (child) ->
#   #regulation of cellular component biogenesis (partner) -> regulation of extent of heterochromatin assembly (stepchild)
#
#   children = lapply(go_terms, function(x) unique(unlist(AnnotationDbi::mget(x, GOBPCHILDREN))))
#   partners = lapply(children, function(x) unique(unlist(lapply(x, function(y) unique(unlist(AnnotationDbi::mget(y, GOBPPARENTS)))))))
#   stepchild = lapply(partners, function(x) unique(unlist(lapply(x, function(y) unique(unlist(AnnotationDbi::mget(y, GOBPOFFSPRING)))))))
#   stepchild = mapply(function(x, y) c(y, x), x=stepchild, y=go_terms)
#   if(!is.list(stepchild)){
#     stepchild = split(stepchild, 1)
#     names(stepchild) = NULL
#   }
#   return(stepchild)
# }


#' Get "step-child" terms for GO terms
#'
#' stepchild are found by finding partners (now birth partners) (GO -> all child terms (1 level down) -> all parent terms (1 level up))
#' then all offspring of partner terms (all offspring of birth partners -> step children) (from above: all parent terms (1 level up)/"partners" -> all child terms)
#'
#' @param go_terms vector of GO ids
#' @keywords internal
#' @import GOfuncR
#' @import dplyr
get_stepchildren = function(go_terms, return_partner=FALSE, return_descriptions=FALSE, verbose =FALSE){

  ## INTERNAL EXAMPLES
  # go_ids = c("GO:0061564", "GO:0048813")
  # id2term(go_ids)
  # get_stepchildren(go_ids[1], return_partner=T, return_description=T)[1:10,]
  # get_stepchildren(go_ids) %>% head
  #
  # # will fail if all GO ids have no children
  # get_stepchildren("GO:0048929")
  # # will skip GO ids with no children
  # # and give warning of the skip if verbose == TRUE
  # get_stepchildren(c(go_ids[1], "GO:0048929"), verbose=TRUE)[1:5,]
  #

  children = GOfuncR::get_child_nodes(go_terms)
  children = children[children$distance == 1,]
  colnames(children)[1] = "go_term"
  youngest_child = unique(go_terms)[!go_terms %in% unique(children$go_term)]
  if(verbose) message("failed to find child terms for: ", paste(youngest_child, collapse = ";"))
  if(nrow(children) == 0){
    stop("Can't find any child terms any input go_terms: ", paste(go_terms, collapse = ";"))
  }

  #####
  # #inclusion body assembly (parent) -> regulation of inclusion body assembly (child) ->
  # #regulation of cellular component biogenesis (partner) -> regulation of extent of heterochromatin assembly (stepchild)

  partners = GOfuncR::get_parent_nodes(children$child_go_id)
  partners = partners[partners$distance == 1,]

  stepchild = GOfuncR::get_child_nodes(partners$parent_go_id)

  partners = full_join(partners, children, by = "child_go_id")
  stepchild =
    full_join(stepchild, partners, by="parent_go_id") %>% arrange(.data$go_term) %>%
    mutate(to_stepchild = paste(.data$go_term, "_", .data$child_go_id.x)) %>%
    filter(!duplicated(.data$to_stepchild))%>%
    #dplyr::select(go_term, child_go_id.x, distance) %>%
    rename(child_go_id = .data$child_go_id.x, child_description = .data$child_name.x) %>%
    rename(partner_go_id = .data$parent_go_id, partner_description = .data$parent_name) %>%
      select(-c(.data$distance.x, .data$to_stepchild, ends_with(".y"))) %>%
      relocate(.data$go_term, .data$partner_go_id, .data$partner_description,
               .data$child_go_id, .data$child_description, .data$distance)


  if(!return_partner){
    stepchild = stepchild %>% select(-starts_with("partner"))
  }
  if(!return_descriptions){
    stepchild = stepchild %>% select(-ends_with("description"))
  }


  n_stepchild = table(stepchild$go_term) %>% as.data.frame() %>%
    rename(go_term = .data$Var1, n_stepchild=.data$Freq)

  stepchild = left_join(stepchild, n_stepchild, by='go_term')
  #

  return(stepchild)
}

#' get the direct path between two GO terms with any kind of offspring/ancestor relationship

#' @param child child GO id
#' @param parent parent GO id
#' @param top_level_terms any 'top level terms' to remove from analysis
#' @import GOfuncR
#' @import dplyr
#' @export
#' @examples
#'
#' get_direct_go_path(child = "GO:0021955", parent = "GO:0061564")
#' # returns ALL potential direct pathways
#' get_direct_go_path(child = "GO:0021955", parent = "GO:0031175")
#' # returns empty if parent:child relationship does not exist
#' get_direct_go_path(child = "GO:0021955", parent = "GO:0050890")
get_direct_go_path = function(child, parent, top_level_terms=NULL){

  tmp = GOfuncR::get_parent_nodes(child) %>% filter(!(.data$parent_go_id %in% top_level_terms))
  all_parents = tmp[(tmp$distance < tmp$distance[tmp$parent_go_id ==  parent] )|( tmp$parent_go_id == parent),]

  all_children = GOfuncR::get_child_nodes(parent) %>% filter(.data$child_go_id %in% all_parents$parent_go_id)
  return(all_children)

}

make_parent_df = function(all_go_terms, parents, children_new, stepchild_new){

  children_new.filtered = children_new[children_new$parent_go_id %in% parents,]
  stepchild_new.filtered = stepchild_new[stepchild_new$go_term %in% parents,]

  #children_new = GOfuncR::get_child_nodes(parents)
  #stepchild_new = get_stepchildren_new(parents)
  n_children_new = as.data.frame(table(children_new.filtered$parent_go_id))

  children_new.filtered$contains_go = F
  children_new.filtered$contains_go[children_new.filtered$child_go_id %in% all_go_terms] = T

  stepchild_new.filtered$contains_go = F
  stepchild_new.filtered$contains_go[stepchild_new.filtered$child_go_id %in% all_go_terms] = T

  ###
  contains_go = table(children_new.filtered$parent_go_id, children_new.filtered$contains_go) %>%
    as.data.frame() %>%
    filter(.data$Var2==TRUE)
  contains_go_stepchild =
    table(stepchild_new.filtered$go_term, stepchild_new.filtered$contains_go) %>%
    as.data.frame() %>%
    filter(.data$Var2==TRUE) %>%
    select(-.data$Var2) %>%
    ## code for adding in n_stepchild_new
    #left_join(stepchild_new[!duplicated(stepchild_new$go_term),c("go_term", "n_stepchild_new")], by=c("Var1"="go_term")) %>%
    #rename(n_stepchild = n_stepchild_new) %>%
    rename(contains_go_stepchild=.data$Freq) %>% rename(parents=.data$Var1)

  if("n_stepchild" %in% colnames(stepchild_new.filtered)){
    contains_go_stepchild$n_stepchild = (stepchild_new.filtered$n_stepchild[match(contains_go_stepchild$parents, stepchild_new.filtered$go_term)])
  }


  parent_term_df =
    data.frame(parents, parent_description = id2term(parents)) %>%
    left_join(n_children_new, by= c('parents'='Var1')) %>% rename(n_children=.data$Freq) %>%
    left_join(contains_go[,-2], by= c('parents'='Var1')) %>% rename(contains_go=.data$Freq)%>%
    left_join(contains_go_stepchild, by= 'parents') %>% filter(.data$contains_go > 0)



  return(parent_term_df)

}

#' Using defined lists of parents for go terms, find parent terms that cover all go terms
#'
#' using the stepchild relationships
#' AND order by specificity of parent/stepparent terms on the defined go terms
#' @param stepchild_ratio use total STEPCHILDREN (\code{TRUE}) or total CHILDREN (\code{FALSE}) of parent terms to calculate ratios/specificity of stepparent terms
#' Using total CHILDREN (\code{FALSE}) may result in ratios higher than 1.
#' @keywords internal
get_top_parents_df = function(all_go_terms, parents, children_new, stepchild_new, stepchild_ratio = T){
  #all_go_terms = all_go_terms_split[[i]]
  go_terms = all_go_terms
  #parent_desc_orig = names(all_go_terms_split)[i]
  parent_terms_cat = NULL
  while(length(all_go_terms) > 0){
    #message(paste(all_go_terms, collapse = "; "))
    #message(check %in% all_go_terms)

    parent_term_df_new = make_parent_df(all_go_terms, parents, children_new, stepchild_new)


    if(stepchild_ratio == F){
      parent_term_df_new$n_stepchild=parent_term_df_new$n_children
    }

    parent_term_df_new.ordered =
      parent_term_df_new %>%
      mutate(go_spec_stepchild = case_when(stepchild_ratio == T ~ contains_go_stepchild / n_stepchild,
                                        TRUE ~ contains_go_stepchild / n_children  )) %>%
      mutate(go_spec = .data$contains_go / .data$n_children) %>%
      filter(.data$contains_go > 0) %>%
      arrange(desc(.data$go_spec_stepchild))


    ## check for other terms that cover the same stepchild as the current top term

    stepchild_new$contains_go = F
    #child-stepchild of top parent term
    all_child_stepchild = stepchild_new$child_go_id[stepchild_new$go_term == parent_term_df_new.ordered$parents[1]]
    all_child_stepchild = all_child_stepchild[all_child_stepchild %in% all_go_terms]
    #all go terms (under investigation) covered by top parent term
    terms_invest = all_go_terms[all_go_terms %in% all_child_stepchild]
    # find all terms that are stepchild (to ANY parent) and set to T
    stepchild_new$contains_go[stepchild_new$child_go_id %in% terms_invest] = T
    # check all parents that are in the tested parent terms (ignoring any ordering so we can check if there is a better match)
    top_relatives = stepchild_new %>% filter(.data$contains_go == T) %>%
      filter(.data$go_term %in% parent_term_df_new.ordered$parents)
    top_relatives = table(top_relatives$go_term) %>% as.data.frame()

    colnames(top_relatives) = c("go_term", "top_stepchild")

    parent_term_df_new.ordered =
      parent_term_df_new %>%
      left_join(top_relatives, by= c('parents'='go_term')) %>%
      mutate(go_spec = .data$contains_go / .data$n_children) %>%
      mutate(go_spec_stepchild = case_when(stepchild_ratio == T ~ contains_go / n_stepchild,
                                        TRUE ~ contains_go / n_children  )) %>%
      arrange(desc(.data$top_stepchild), desc(.data$go_spec_stepchild))

    # parent_term_df.ordered =
    #   parent_term_df %>%
    #   left_join(top_stepchild[,-2], by= c('parents'='Var1')) %>% rename(top_stepchild=Freq) %>%
    #   mutate(go_spec = contains_go / n_children) %>%
    #   mutate(go_spec_stepchild = contains_go_stepchild / n_stepchild) %>%
    #   arrange(desc(top_stepchild),desc(go_spec))

    go_terms_categorised =
      all_go_terms[which(all_go_terms %in% stepchild_new$child_go_id[stepchild_new$go_term == parent_term_df_new.ordered$parents[1]] )]

    all_go_terms = all_go_terms[-which(all_go_terms %in% go_terms_categorised)]
    parent_terms_cat = rbind(parent_terms_cat, parent_term_df_new.ordered[1,])

  }

  ### check for overlapping categories and remove parents that are entriely covered by another 'top' parent term
  ## I.E Remove redundancy

  all_go_terms = go_terms
  parents = parent_terms_cat$parents

  parent_term_df_new = make_parent_df(all_go_terms, parents, children_new, stepchild_new)

  m = match(parent_term_df_new$parents[order(parent_term_df_new$contains_go_stepchild, decreasing = F)], parents)
  co_terms = lapply(m, function(x) all_go_terms[which(all_go_terms %in% stepchild_new$child_go_id[stepchild_new$go_term == parents[x]])])
  names(co_terms) = parent_term_df_new$parent_description[order(parent_term_df_new$contains_go_stepchild, decreasing = F)]
  redundant = NULL

  co_term_names = names(co_terms)
  for(c_n in seq_along(co_terms)){
    c_name = co_term_names[c_n]
    z = which(names(co_terms) == c_name)
    if(any(unlist(lapply(co_terms[-z], function(x) all(co_terms[[z]] %in% x))))){
      redundant = c(redundant, names(co_terms[z]))
      co_terms = co_terms[-z]
    }
  }

  parent_term_df_new = parent_term_df_new[parent_term_df_new$parent_description %in% names(co_terms),]

  ### make new name which covers all parent names
  new_parent = paste0(parent_term_df_new$parent_description, collapse = ";\n")
  return(new_parent)

}
