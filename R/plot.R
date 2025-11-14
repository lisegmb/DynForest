#' Plot function in dynforest
#'
#' This function displays a plot of CIF for a given node and tree (for class \code{dynforest}),
#' the most predictive variables with the minimal depth (for class \code{dynforestvardepth}),
#' the variable importance (for class \code{dynforestvimp}), grouped variable importance (for class \code{dynforestgvimp}),
#' or ICE and PDP plots (for class \code{dynforestpdp}).
#'
#' @param x Object inheriting from classes \code{dynforest}, \code{dynforestvardepth}, \code{dynforestvimp}, \code{dynforestgvimp}, or \code{dynforestpdp}, to respectively plot the CIF, the minimal depth, the variable importance, grouped variable importance, or ICE/PDP plots.
#' @param tree For \code{dynforest} class, integer indicating the tree identifier.
#' @param nodes For \code{dynforest} class, identifiers for the selected nodes.
#' @param id For \code{dynforest} and \code{dynforestpred} classes, identifier for a given subject.
#' @param max_tree For \code{dynforest} class, integer indicating the number of trees to display when using \code{id} argument.
#' @param plot_level For \code{dynforestvardepth} class, compute the statistic at predictor (\code{plot_level} = "predictor") or feature (\code{plot_level} = "feature") level.
#' @param PCT For \code{dynforestvimp} or \code{dynforestgvimp} classes, logical to display VIMP statistic in percentage. Default is FALSE.
#' @param ordering For \code{dynforestvimp} class, logical to order predictors according to VIMP value. Default is TRUE.
#' @param type For \code{dynforestpdp} class, character indicating the plot type to display: "both" (default), "ice", or "pdp".
#' @param conf_band For \code{dynforestpdp} class, logical indicating whether to display confidence bands around the PDP. Default is FALSE.
#' @param alpha For \code{dynforestpdp} class, numeric specifying the transparency level of the confidence band. Default is 0.2.
#' @param ... Optional parameters to be passed to the low level function.
#'
#' @import ggplot2 viridis dplyr tidyr magrittr
#' @importFrom stringr str_order
#'
#' @seealso [dynforest()] [compute_ooberror()] [compute_vimp()] [compute_gvimp()] [compute_vardepth()]
#'
#' @return \code{plot()} function displays:
#' \tabular{ll}{
#'    With \code{dynforestvardepth} \tab the minimal depth for each predictor/feature \cr
#'    \cr
#'    With \code{dynforestvimp} \tab the VIMP for each predictor \cr
#'    \cr
#'    With \code{dynforestgvimp} \tab the grouped-VIMP for each given group \cr
#'    \cr
#'    With \code{dynforestpdp} \tab ICE (Individual Conditional Expectation) and/or PDP (Partial Dependence Plot) with optional confidence bands \cr
#' }

#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Get Gaussian distribution for longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
#'
#' # Sample 100 subjects
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id%in%id_sample)
#'
#' pbc2_train <- pbc2[id_row,]
#'
#  Build longitudinal data
#' timeData_train <- pbc2_train[,c("id","time",
#'                                 "serBilir","SGOT",
#'                                 "albumin","alkaline")]
#'
#' # Create object with longitudinal association for each predictor
#' timeVarModel <- list(serBilir = list(fixed = serBilir ~ time,
#'                                      random = ~ time),
#'                      SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                  random = ~ time + I(time^2)),
#'                      albumin = list(fixed = albumin ~ time,
#'                                     random = ~ time),
#'                      alkaline = list(fixed = alkaline ~ time,
#'                                      random = ~ time))
#'
#' # Build fixed data
#' fixedData_train <- unique(pbc2_train[,c("id","age","drug","sex")])
#'
#' # Build outcome data
#' Y <- list(type = "surv",
#'           Y = unique(pbc2_train[,c("id","years","event")]))
#'
#' # Run dynforest function
#' res_dyn <- dynforest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Plot estimated CIF at nodes 17 and 32
#' plot(x = res_dyn, tree = 1, nodes = c(17,32))
#'
#' # Run var_depth function
#' res_varDepth <- compute_vardepth(res_dyn)
#'
#' # Plot minimal depth
#' plot(x = res_varDepth, plot_level = "feature")
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_vimp(dynforest_obj = res_dyn, ncores = 2)
#'
#' # Plot VIMP
#' plot(x = res_dyn_VIMP, PCT = TRUE)
#'
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gvimp(dynforest_obj = res_dyn,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2)
#'
#' # Plot gVIMP
#' plot(x = res_dyn_gVIMP, PCT = TRUE)
#'
#' # Sample 5 subjects to predict the event
#' set.seed(123)
#' id_pred <- sample(id, 5)
#'
#' # Create predictors objects
#' pbc2_pred <- pbc2[which(pbc2$id%in%id_pred),]
#' timeData_pred <- pbc2_pred[,c("id", "time", "serBilir", "SGOT", "albumin", "alkaline")]
#' fixedData_pred <- unique(pbc2_pred[,c("id","age","drug","sex")])
#'
#' # Predict the CIF function for the new subjects with landmark time at 4 years
#' pred_dyn <- predict(object = res_dyn,
#'                     timeData = timeData_pred, fixedData = fixedData_pred,
#'                     idVar = "id", timeVar = "time",
#'                     t0 = 4)
#'
#' # Plot predicted CIF for subjects 26 and 110
#' plot(x = pred_dyn, id = c(26, 110))
#'
#' }
#'
#' @rdname plot.dynforest
#' @export
plot.dynforest <- function(x, tree = NULL, nodes = NULL, id = NULL, max_tree = NULL, ...){

  if (!methods::is(x,"dynforest")){
    cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You've supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }

  if (!is.null(tree)){

    if (!inherits(tree, "numeric")){
      cli_abort(c(
        "{.var tree} must be a numeric object containing the tree identifier",
        "x" = "You've supplied a {.cls {class(tree)}} object"
      ))
    }

    if (!any(tree==seq(x$param$ntree))){
      cli_abort(c(
        "{.var tree} must be chosen between 1 and {x$param$ntree}",
        "x" = "You've chosen {tree}"
      ))
    }

    if (all(!is.null(nodes))){
      if (!all(inherits(nodes, "numeric"))){
        cli_abort(c(
          "{.var nodes} must be a numeric vector containing the node identifiers",
          "x" = "You've supplied a {.cls {class(nodes)}} object"
        ))
      }
      if (!all(nodes%in%names(x$rf[,tree]$Y_pred))){
        cli_abort(c(
          "At least one selected node in {.var nodes} doesn't exist in {.var tree} {tree}"
        ))
      }
      if (any(sapply(nodes, FUN = function(node) is.null(x$rf[,tree]$Y_pred[[as.character(node)]])))){
        cli_abort(c(
          "At least one selected node in {.var nodes} doesn't exist in {.var tree} {tree}"
        ))
      }
    }else{
      nodes <- get_treenodes(dynforest_obj = x, tree = tree)
    }

    # data transformation for ggplot2
    CIFs_nodes_list <- lapply(nodes, FUN = function(node){

      CIFs_node <- x$rf[,tree]$Y_pred[[as.character(node)]]

      CIFs_node_list <- lapply(names(CIFs_node), FUN = function(y){

        CIF_node_cause <- CIFs_node[[y]]

        out <- data.frame(Node = rep(node, nrow(CIF_node_cause)),
                          Cause = rep(y, nrow(CIF_node_cause)),
                          Time = CIF_node_cause$times,
                          CIF = CIF_node_cause$traj)

        return(out)

      })

      return(do.call(rbind, CIFs_node_list))

    })

    data.CIF.plot <- do.call(rbind, CIFs_nodes_list)

    g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
      geom_step(aes_string(group = "Cause", color = "Cause")) +
      facet_wrap(~ Node) +
      ylim(0,1) +
      theme_bw()

    return(print(g))

  }

  if (!is.null(id)){

    nodes <- apply(x$rf, 2, FUN = function(y){
      leaf_tree <- y$leaves[which(y$idY==id)]
    })

    data.CIF.plot <- NULL

    for (tree_id in seq(length(nodes))){

      if (length(nodes[[tree_id]])>0){
        tree_node <- nodes[[tree_id]]
      }else{
        next()
      }

      CIFs_node <- x$rf[,tree_id]$Y_pred[[as.character(tree_node)]]

      CIFs_node_list <- lapply(names(CIFs_node), FUN = function(y){

        CIF_node_cause <- CIFs_node[[y]]

        out <- data.frame(Tree = rep(tree_id, nrow(CIF_node_cause)),
                          Node = rep(tree_node, nrow(CIF_node_cause)),
                          Cause = rep(y, nrow(CIF_node_cause)),
                          Time = CIF_node_cause$times,
                          CIF = CIF_node_cause$traj)

        return(out)

      })

      data.CIF.plot <- rbind(data.CIF.plot, do.call(rbind, CIFs_node_list))

    }

    data.CIF.plot$Tree_Node <- paste0("Tree ", data.CIF.plot$Tree, " / Node ", data.CIF.plot$Node)
    data.CIF.plot$Tree_Node <- factor(data.CIF.plot$Tree_Node, levels = unique(data.CIF.plot$Tree_Node))

    if (!is.null(max_tree)){
      max_tree_id <- unique(data.CIF.plot$Tree)[seq(max_tree)]
      data.CIF.plot <- data.CIF.plot[which(data.CIF.plot$Tree%in%max_tree_id),]
    }

    g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
      geom_step(aes_string(group = "Cause", color = "Cause")) +
      facet_wrap(~ Tree_Node) +
      ylim(0,1) +
      theme_bw()

    return(print(g))

  }

}


#' @name plot.dynforest
#' @export
plot.dynforestvardepth <- function(x, plot_level = c("predictor","feature"), ...){

  # checking
  if (!all(plot_level%in%c("predictor","feature"))){
    stop("Only 'predictor' and 'feature' options are allowed for plot_level argument!")
  }

  if (length(plot_level)>1){
    plot_level <- plot_level[1]
  }

  min_depth_all <- x$var_node_depth
  depth.df <- data.frame(var = rep(min_depth_all$var, ncol(min_depth_all)-1),
                         id_node = unlist(c(x$var_node_depth[,2:ncol(min_depth_all)])),
                         tree = rep(seq(ncol(min_depth_all)-1), each = nrow(min_depth_all)))
  depth.df$group <- sub("\\..*", "", depth.df$var)
  depth.df <- depth.df[order(depth.df$var),]

  if (plot_level=="feature"){

    depth.nbtree <- aggregate(id_node ~ var, data = depth.df, FUN = function(x){
      return(sum(!is.na(x)))
    })

    g <- ggplot(depth.df, aes_string(x = "var", y = "id_node")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depth.nbtree, aes_string(x = "var", label = "id_node"),
                y = max(depth.df$id_node, na.rm = T) - 1) +
      scale_x_discrete(limits=rev(unique(depth.df$var)[str_order(unique(depth.df$var))])) +
      xlab("Features") +
      ylab("Minimal depth") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold")) +
      coord_flip()

    return(print(g))

  }

  if (plot_level=="predictor"){

    depthVar.df <- aggregate(id_node ~ group + tree, data = depth.df, min, na.rm = T)

    depthVar.nbtree <- aggregate(id_node ~ group, data = depthVar.df, FUN = length)

    g <- ggplot(depthVar.df, aes_string(x = "group", y = "id_node")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depthVar.nbtree, aes_string(x = "group", label = "id_node"),
                y = max(depthVar.df$id_node, na.rm = T) - 1) +
      scale_x_discrete(limits=rev(unique(depthVar.df$group)[str_order(unique(depthVar.df$group))])) +
      xlab("Predictors") +
      ylab("Minimal depth") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold")) +
      coord_flip()

    return(print(g))

  }

}

#' @rdname plot.dynforest
#' @export
plot.dynforestvimp <- function(x, PCT = FALSE, ordering = TRUE, ...){

  vimp.df <- data.frame(var = unlist(x$Inputs),
                        vimp = unlist(x$Importance))

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(x$tree_oob_err, na.rm = T) # vimp relative
  }

  if (ordering){
    g <- ggplot(vimp.df) +
      geom_bar(aes_string("var", "vimp"), stat = "identity") +
      scale_x_discrete(limits=vimp.df$var[order(vimp.df$vimp)]) +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }else{
    g <- ggplot(vimp.df) +
      geom_bar(aes_string("var", "vimp"), stat = "identity") +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }

  return(print(g))
}

#' @rdname plot.dynforest
#' @export
plot.dynforestgvimp <- function(x, PCT = FALSE, ...){

  vimp.df <- data.frame(var = names(x$gVIMP),
                        vimp = x$gVIMP)

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(x$tree_oob_err, na.rm = T) # vimp relative
  }

  g <- ggplot(vimp.df) +
    geom_bar(aes_string("var", "vimp"), stat = "identity") +
    scale_x_discrete(limits=vimp.df$var[order(vimp.df$vimp)]) +
    xlab("Group of predictors") +
    ylab(ifelse(PCT,"% grouped-VIMP","grouped-VIMP")) +
    coord_flip() +
    theme_bw()

  return(print(g))

}

#' @rdname plot.dynforest
#' @export
plot.dynforestpred <- function(x, id = NULL, ...){

  if (!methods::is(x,"dynforestpred")){
    stop("'x' should be a 'dynforestpred' class!")
  }

  if (is.null(id)){
    stop("'id' cannot be NULL!")
  }

  if (!all(id%in%rownames(x$pred_indiv))){
    stop("Predictions are not available for some subjects. Please verify the subjects identifiers!")
  }

  data.CIF <- x$pred_indiv
  data.CIF <- data.CIF[which(rownames(data.CIF)%in%id),, drop = FALSE]

  times <- x$times
  n.times <- length(times)

  data.CIF.plot <- data.frame(id = as.factor(rep(id, each = n.times)),
                              Time = rep(times, length(id)),
                              CIF = c(t(data.CIF)))

  g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
    geom_step(aes(group = id, color = id)) +
    ylim(0,1) +
    geom_vline(xintercept = x$t0, linetype = "dashed") +
    theme_bw()

  print(g)

}

#' @rdname plot.dynforest
#' @method plot dynforestpdp
#' @export
plot.dynforestpdp <- function(x, x_label = NULL, y_label = NULL,
                              title = "ICE and PDP Plot", conf_band = FALSE, alpha = 0.2,
                              type = c("both", "ice", "pdp"),
                              target_class = NULL, ...) {
  type <- match.arg(type)
  df <- x

  # ---------------------------------------------------------------
  #  Detect the model type from the data structure
  #  - Survival     → columns: id, replicate_id, time, value (+1 parameter)
  #  - Classification → presence of 'target_class'
  #  - Regression    → numeric response
  # ---------------------------------------------------------------
  if (all(c("id","replicate_id","time","value") %in% colnames(df)) && ncol(df) == 5) {
    model_type <- "surv"
  } else if ("target_class" %in% colnames(df)) {
    model_type <- "factor"
  } else {
    model_type <- "numeric"
  }

  # ---------------------------------------------------------------
  #  Identify the variable of interest (the one varied in PDP/ICE)
  # ---------------------------------------------------------------
  exclude_cols <- c("id","replicate_id","value")
  if (model_type=="surv") exclude_cols <- c(exclude_cols,"time")
  if (model_type=="factor") exclude_cols <- c(exclude_cols,"target_class")

  var_cols <- setdiff(colnames(df), exclude_cols)
  color_var <- if(length(var_cols)>0) var_cols[1] else NULL

  # ---------------------------------------------------------------
  #  Determine the type of the variable of interest
  #  numeric  → continuous predictor
  #  factor   → categorical predictor
  #  longitudinal → dynamic time-varying predictor
  # ---------------------------------------------------------------
  var_type <- if (!is.null(color_var)) class(df[[color_var]])[1] else NA
  var_type <- ifelse(var_type %in% c("numeric","integer"), "numeric",
                     ifelse(var_type %in% c("factor"), "factor", "longitudinal"))

  # ---------------------------------------------------------------
  #  Automatic axis labels depending on the model type
  # ---------------------------------------------------------------
  if (is.null(y_label)) {
    y_label <- switch(model_type,
                      surv = "Predicted risk (CIF)",
                      factor = "Predicted probability",
                      numeric = "Predicted value")
  }
  if (is.null(x_label)) {
    if (model_type == "surv") {
      x_label <- "Time"
    } else {
      x_label <- if (!is.null(color_var)) color_var else ""
    }
  }

  # ---------------------------------------------------------------
  #  Filter for a specific class (classification only)
  # ---------------------------------------------------------------
  if (model_type=="factor" && !is.null(target_class)) {
    df <- df[df$target_class %in% target_class, , drop=FALSE]
  }

  p_ICE <- NULL
  p_PDP <- NULL

  # ===============================================================
  #                        SURVIVAL MODELS
  # ===============================================================
  if(model_type=="surv") {

    # Identify the variable used for color/facets
    extra_cols <- setdiff(colnames(df), c("id","replicate_id","time","value"))
    color_var_plot <- if ("transformation" %in% colnames(df)) "transformation"
    else if(length(extra_cols)>0) extra_cols[1] else NULL

    if(!is.null(color_var_plot))
      df[[color_var_plot]] <- factor(df[[color_var_plot]], levels=sort(unique(df[[color_var_plot]])))

    # -------------------------- ICE PLOT --------------------------
    p_ICE <- ggplot(df, aes(x=time, y=value,
                            group=interaction(id, .data[[color_var_plot]]),
                            color=.data[[color_var_plot]])) +
      geom_line(alpha=0.5, linewidth=0.7) +
      labs(x=x_label, y=y_label, title=paste("ICE -", title), color=color_var_plot) +
      theme_minimal()

    # -------------------------- PDP PLOT --------------------------
    df_mean <- df %>%
      group_by(time, .data[[color_var_plot]]) %>%
      summarise(mean_value = mean(value),
                var_mc = stats::var(value)/n(),
                n = n(),
                .groups="drop") %>%
      mutate(se_mc = sqrt(var_mc),
             lower = mean_value - if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else 0,
             upper = mean_value + if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else 0)

    p_PDP <- ggplot(df_mean, aes(x=time, y=mean_value,
                                 color=.data[[color_var_plot]],
                                 fill=.data[[color_var_plot]])) +
      geom_line(linewidth=1) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=alpha, color=NA) +
      labs(x=x_label, y=y_label, title=paste("PDP -", title),
           color=color_var_plot, fill=color_var_plot) +
      theme_minimal()
  }

  # ===============================================================
  #                REGRESSION + CLASSIFICATION MODELS
  # ===============================================================
  else if(model_type %in% c("numeric","factor")) {

    # x-axis variable
    x_val <- if(var_type=="longitudinal") "transformation" else color_var
    y_val <- "value"

    # ---------------------------------------------------------------
    #     Classification + numeric predictor
    # ---------------------------------------------------------------
    if(model_type=="factor" && var_type=="numeric") {

      # -------------------------- ICE --------------------------
      p_ICE <- ggplot(df, aes_string(
        x = x_val, y = y_val,
        group = "interaction(id, target_class)",
        color = "target_class"
      )) +
        geom_line(alpha = 0.4, linewidth = 0.7) +
        labs(x = x_label, y = y_label, title = paste("ICE -", title), color = "Class") +
        theme_minimal()

      # -------------------------- PDP --------------------------
      df_mean <- df %>%
        group_by(.data[[x_val]], target_class) %>%
        summarise(mean_value = mean(value),
                  var_mc = stats::var(value)/n(),
                  n = n(),
                  .groups="drop") %>%
        mutate(se_mc = sqrt(var_mc),
               lower = mean_value - if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else NA,
               upper = mean_value + if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else NA)

      p_PDP <- ggplot(df_mean, aes_string(x=x_val, y="mean_value",
                                          color="target_class", fill="target_class")) +
        geom_line(linewidth=1) +
        {if(conf_band) geom_ribbon(aes(ymin=lower, ymax=upper), alpha=alpha, color=NA) else NULL} +
        labs(x=x_label, y=y_label, title=paste("PDP -", title),
             color="Class", fill="Class") +
        theme_minimal()
    }

    # ---------------------------------------------------------------
    # Classification + categorical or longitudinal predictor
    # ---------------------------------------------------------------
    else if(model_type=="factor" && var_type %in% c("factor","longitudinal")) {

      # -------------------------- ICE --------------------------
      p_ICE <- ggplot(df, aes_string(x=x_val, y=y_val, fill=x_val)) +
        geom_boxplot(alpha=0.5) +
        labs(x=x_label, y=y_label, title=paste("ICE -", title), fill=x_val) +
        theme_minimal() +
        facet_wrap(~target_class,
                   labeller = labeller(target_class = function(x) paste0("Class : ", x)))

      # -------------------------- PDP --------------------------
      df_mean <- df %>%
        group_by(.data[[x_val]], target_class) %>%
        summarise(mean_value = mean(value),
                  sd_value = stats::sd(value),
                  n = n(),
                  .groups="drop") %>%
        mutate(lower = mean_value - sd_value,
               upper = mean_value + sd_value,
               lower_band = if(conf_band) mean_value - stats::qt(0.975, df=pmax(n-1,1))*sd_value/sqrt(n) else NA,
               upper_band = if(conf_band) mean_value + stats::qt(0.975, df=pmax(n-1,1))*sd_value/sqrt(n) else NA)

      p_PDP <- ggplot(df_mean, aes_string(x=x_val, y="mean_value")) +
        geom_point(size=3) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
        {if(conf_band) geom_errorbar(aes(ymin=lower_band, ymax=upper_band),
                                     width=0.5, color="blue", alpha=0.3) else NULL} +
        labs(x=x_label, y=y_label, title=paste("PDP -", title)) +
        theme_minimal() +
        facet_wrap(~target_class,
                   labeller = labeller(target_class = function(x) paste0("Class : ", x)))
    }

    # ---------------------------------------------------------------
    # Regression models
    # ---------------------------------------------------------------
    else if(model_type=="numeric") {

      # -------------------------- ICE --------------------------
      if(var_type=="numeric") {
        p_ICE <- ggplot(df, aes_string(x=x_val, y=y_val, group="id")) +
          geom_line(alpha=0.5) +
          labs(x=x_label, y=y_label, title=paste("ICE -", title)) +
          theme_minimal()
      } else {
        p_ICE <- ggplot(df, aes_string(x=x_val, y=y_val, fill=x_val)) +
          geom_boxplot(alpha=0.5) +
          labs(x=x_label, y=y_label, title=paste("ICE -", title), fill=x_val) +
          theme_minimal()
      }

      # -------------------------- PDP --------------------------
      if(var_type=="numeric") {
        df_mean <- df %>%
          group_by(.data[[x_val]]) %>%
          summarise(mean_value = mean(value),
                    var_mc = stats::var(value)/n(),
                    n = n(),
                    .groups="drop") %>%
          mutate(se_mc = sqrt(var_mc),
                 lower = mean_value - if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else NA,
                 upper = mean_value + if(conf_band) stats::qt(0.975, df=pmax(n-1,1))*se_mc else NA)

        p_PDP <- ggplot(df_mean, aes_string(x=x_val, y="mean_value")) +
          geom_line(linewidth=1, color="blue") +
          {if(conf_band) geom_ribbon(aes(ymin=lower, ymax=upper),
                                     alpha=alpha, fill="blue", color=NA) else NULL} +
          labs(x=x_label, y=y_label, title=paste("PDP -", title)) +
          theme_minimal()
      } else {
        df_mean <- df %>%
          group_by(.data[[x_val]]) %>%
          summarise(mean_value = mean(value),
                    sd_value = stats::sd(value),
                    n = n(),
                    .groups="drop") %>%
          mutate(lower = mean_value - sd_value,
                 upper = mean_value + sd_value,
                 lower_band = if(conf_band) mean_value - stats::qt(0.975, df=pmax(n-1,1))*sd_value/sqrt(n) else NA,
                 upper_band = if(conf_band) mean_value + stats::qt(0.975, df=pmax(n-1,1))*sd_value/sqrt(n) else NA)

        p_PDP <- ggplot(df_mean, aes_string(x=x_val, y="mean_value")) +
          geom_point(size=3) +
          geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
          {if(conf_band) geom_errorbar(aes(ymin=lower_band, ymax=upper_band),
                                       width=0.5, color="blue", alpha=0.3) else NULL} +
          labs(x=x_label, y=y_label, title=paste("PDP -", title)) +
          theme_minimal()
      }
    }
  }

  # ---------------------------------------------------------------
  #  Final output depending on selected type
  # ---------------------------------------------------------------
  if(type=="ice"){ print(p_ICE); return(invisible(p_ICE)) }
  else if(type=="pdp"){ print(p_PDP); return(invisible(p_PDP)) }
  else { print(p_ICE); print(p_PDP); return(invisible(list(ICE=p_ICE, PDP=p_PDP))) }
}
