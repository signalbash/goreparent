#' plot the output from a goreparent add_go_group()
#'
#' @param go_results output from add_go_groups()
#' data.frame with original GO terms, new groups, and parent terms.
#' must have column named "parent_description"
#' and optionally "Description" for displaying individual GO ids (with \code{collapse=FALSE})
#' @param type \code{"dot"}. Currently nonfunctional, only dot available.
#' @param collapse collapse individual terms onto the same y-axis value for each group?. default \code{TRUE}
#' @param n_top number of groups to plot, ordered by \code{sort_top_by}.
#' default = \code{10}
#' @param sort_top_by how to sort groups for n_top. The first match to a character variable will be used.
#' default \code{c("FDR", "p.adjust", "pval", "p", "p.val")}
#' @param color which variable to color dots by. default \code{"FDR"}
#' @param size which variable to size dots by. default \code{NULL}
#' @param x.axis which variable to use as the x-axis. default \code{c("GeneRatio", "Odds.Ratio")}
#' @importFrom data.table rbindlist
#' @import dplyr
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom stats formula
#' @importFrom forcats fct_reorder
#' @importFrom forcats as_factor
#' @export
#' @examples
#' enrich_file = system.file("extdata","go_enrich_results.txt", package = "goreparent")
#' go_input = read.table(enrich_file, sep="\t", header=TRUE)
#' go_parents = add_go_groups(go_input, descriptive_parent=FALSE, max_from_top=3)
#' go_clusters = add_go_groups(go_input, descriptive_parent=TRUE, method = "cluster")
#'
#' plot_go_parents(go_parents)
#' plot_go_parents(go_clusters)
#'
#' # un-collapse parent terms and show all children for the top 2 categories
#' plot_go_parents(go_parents, collapse=FALSE, n_top=2)
#'
plot_go_parents = function(go_results, type="dot", collapse=TRUE, n_top=10,
                           sort_top_by = c("FDR", "p.adjust", "pval", "p", "p.val"),
                           color="FDR", size = NULL, x.axis = c("GeneRatio", "Odds.Ratio")){

  compare = FALSE
  # options not yet added for plottypes other than "dot"
  type="dot"

  if(class(go_results) == "list" & length(go_results) > 1){

    top = lapply(go_results, function(x) get_top_parent_descriptions(x, n_top = n_top))
    top = unique(unlist(top))

    go_results = data.table::rbindlist(go_results, idcol=TRUE) %>%
      as.data.frame() %>% rename(set = .data$.id)
    compare = T
  }else if(class(go_results) == "list" & length(go_results) == 1){
    message("you might be trying to plot a result which is formatted as a list of length 1...")
    message("converting to underlying data.frame for plotting...")
    go_results = go_results[[1]]
  }

  FDR_or_p = NA
  m = which(!is.na(match(sort_top_by, colnames(go_results))))
  if(length(m) > 0){
    FDR_or_p = sort_top_by[m][1]
  }else{
    stop("Can't find any of the sorting variables (sort_top_by = c(", paste(sort_top_by, collapse = ", "), ")) in go_results")
  }

  ## check x-axis, size, colour are variables
  if(!any(colnames(go_results) %in% c(color))){
    message(color, " is not found in the data.frame supplied, defaulting to ", FDR_or_p)
    color = FDR_or_p
  }
  if(!any(colnames(go_results) %in% c(size)) & !is.null(size)){
    message(size, " is not found in the data.frame supplied, set size to NULL if you don't want to plot size")
    size = FDR_or_p
  }
  if(!any(colnames(go_results) %in% x.axis)){
    message(paste(x.axis, collapse = " or "), " is not found in the data.frame supplied, defaulting to ", FDR_or_p)
    x.axis = FDR_or_p
  }else{
    x.axis = x.axis[x.axis %in% colnames(go_results)][1]
  }

  if(compare){


    if(collapse){
      go_results_summary = aggregate(formula(paste0(FDR_or_p,  " ~ parent_description + set")),
                                     data=go_results, FUN=function(x) length(x))

      colnames(go_results_summary)[ncol(go_results_summary)] = "GO_subcat_count"
      go_results_summary$min_FDR = aggregate(formula(paste0(FDR_or_p,  " ~ parent_description + set")),
                                             data=go_results, FUN=function(x) min(x))[,3]
      colnames(go_results_summary)[ncol(go_results_summary)] = FDR_or_p

      p = go_results_summary %>% filter(.data$parent_description %in% top) %>%
        mutate(parent_description = factor(.data$parent_description, levels = rev(top))) %>%
        ggplot(aes_string(x="set", y="parent_description", size="GO_subcat_count", col=sprintf("-log10(%s)", color))) +
        geom_point()+
        theme_minimal() + theme(panel.border = element_rect(fill=NA, size=0.1)) +
        scale_color_gradient(low="red", high="blue") +
        scale_size("GO child\ncount") +
        scale_y_discrete("GO parent term") + theme(axis.title.x = element_blank())
    }else{

      p = go_results %>%
        filter(.data$parent_description %in% top) %>%
        mutate(Description = as_factor(.data$Description)) %>%
        ggplot(aes_string(x="set", y = "Description",  col=sprintf("-log10(%s)", color), size=size)) +
        geom_point()+ scale_color_gradient(low="red", high="blue") +
        facet_grid(rows = vars(.data$parent_description), scales="free_y", space="free_y", switch="y")+
        theme_minimal()+ theme(panel.border = element_rect(fill=NA, size=0.1)) +
        theme(strip.placement = "top") + coord_cartesian(clip = 'off')+
        theme(strip.text.y.left = element_text(angle = 0, face="bold", hjust=1))+ theme(axis.title.x = element_blank())
    }


  }

  else{
    go_results = go_results[order(go_results[,FDR_or_p], decreasing = FALSE),]

    top = go_results %>%
      filter(!duplicated(.data$parent_description)) %>%
      head(n_top) %>%
      dplyr::select(.data$parent_description)

    if(collapse == TRUE & type == "dot"){
      p = go_results %>%
        filter(.data$parent_description %in% top$parent_description) %>%
        dplyr::rename(TEST = matches(x.axis)) %>%
        mutate(parent_description = fct_reorder(.data$parent_description, .data$TEST, .fun=max)) %>%
        dplyr::rename_with(function(x) x.axis,  matches("TEST")) %>%
        mutate(Description = as_factor(.data$Description)) %>%
        ggplot(aes_string(x=x.axis, y = "parent_description", col=sprintf("-log10(%s)", color), size=size)) +
        geom_point() + scale_color_gradient(low="red", high="blue") +
        theme_minimal() + theme(panel.border = element_rect(fill=NA, size=0.1)) + scale_y_discrete("GO parent term")
    }else if(collapse == FALSE & type == "dot"){
      p = go_results %>%
        filter(.data$parent_description %in% top$parent_description) %>%
        dplyr::rename(TEST = matches(x.axis)) %>%
        mutate(parent_description = fct_reorder(.data$parent_description, .data$TEST, .fun=max)) %>%
        dplyr::rename_with(function(x) x.axis,  matches("TEST")) %>%
        mutate(Description = as_factor(.data$Description)) %>%

        ggplot(aes_string(x=x.axis, y = "Description",  col=sprintf("-log10(%s)", color), size=size)) +
        geom_point()+ scale_color_gradient(low="red", high="blue") +
        facet_grid(rows = vars(.data$parent_description), scales="free_y", space="free_y", switch="y")+
        theme_minimal()+ theme(panel.border = element_rect(fill=NA, size=0.1)) +
        #ggforce::facet_col(facet = vars(.data$parent_description), scales= "free_y", space="free")
        theme(strip.placement = "top") + coord_cartesian(clip = 'off')+
        theme(strip.text.y.left = element_text(angle = 0, face="bold", hjust=1))
    }
  }
  print(p)
}
