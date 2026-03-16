#' Compute ICE and PDP for dynforest models with automatic data retrieval
#'
#' Compute Individual Conditional Expectation (ICE) and Partial Dependence Plot (PDP)
#' values for a specified predictor from a fitted \code{dynforest} model.
#' Supports numeric, categorical, and longitudinal predictors.
#'
#' @param timeData Optional data.frame with subject ID, time variable, and longitudinal predictors.
#' @param fixedData Optional data.frame with subject ID and fixed predictors.
#' @param model A fitted \code{dynforest} object.
#' @param var_name Name of the predictor variable for ICE/PDP.
#' @param values Numeric vector, factor levels, or named list of functions for longitudinal.
#' @param target_class Classes to retain for classification.
#' @param grid_type Grid type for numeric predictors: "regular", "observed", "quantile".
#' @param grid_size Number of points in grid.
#' @param idVar Subject ID variable (auto-detected if NULL).
#' @param timeVar Time variable (defaults to model$timeVar or "time").
#' @param t0 Landmark time for dynamic prediction.
#' @param ncores Number of cores for parallel computation.
#' @param conf_level Confidence level for PDP intervals.
#'
#' @return Object of class \code{dynforestpdp} with elements \code{ice}, \code{pdp}, \code{model_type}, \code{var_type}, \code{grid}.
#' @export
compute_ice_pdp <- function(
    timeData = NULL,
    fixedData = NULL,
    model,
    var_name,
    values = NULL,
    timeVar = NULL,
    idVar = NULL,
    t0 = NULL,
    target_class = NULL,
    grid_type = c("regular", "observed", "quantile"),
    grid_size = 50,
    ncores = 1,
    conf_level = 0.95
) {

  # Auto-detect idVar / timeVar if missing
  if (is.null(idVar)) {
    if (!is.null(model$data$Numeric) && length(model$data$Numeric) >= 3) {
      idVar <- names(model$data$Numeric)[3]
    } else if (!is.null(model$data$Factor) && length(model$data$Factor) >= 3) {
      idVar <- names(model$data$Factor)[3]
    } else if (!is.null(model$data$Longitudinal) && length(model$data$Longitudinal) >= 3) {
      idVar <- names(model$data$Longitudinal)[3]
    } else stop("Cannot auto-detect idVar. Please provide it explicitly.")
  }

  if (is.null(timeVar)) timeVar <- ifelse(!is.null(model$timeVar), model$timeVar, "time")

  # Determine model type
  model_type <- model$type
  if (!(model_type %in% c("surv","numeric","factor"))) stop("Unsupported model type.")

  # Determine variable type
  var_type <- NULL
  if (!is.null(model$data$Longitudinal$X) && var_name %in% colnames(model$data$Longitudinal$X)) var_type <- "longitudinal"
  else if (!is.null(model$data$Numeric$X) && var_name %in% colnames(model$data$Numeric$X)) var_type <- "numeric"
  else if (!is.null(model$data$Factor$X) && var_name %in% colnames(model$data$Factor$X)) var_type <- "factor"
  else stop("Variable '", var_name, "' not found in model data.")

  # Reconstruct timeData if missing
  if (is.null(timeData) && !is.null(model$data$Longitudinal$X)) {
    long_raw <- model$data$Longitudinal$X
    long_raw[[idVar]] <- model$data$Longitudinal$id
    long_raw[[timeVar]] <- model$data$Longitudinal$time
    value_cols <- setdiff(colnames(long_raw), c(idVar, timeVar))
    timeData <- long_raw %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(value_cols),
                          names_to = "Variable", values_to = "Value") %>%
      dplyr::select(dplyr::all_of(c(idVar, timeVar, "Variable", "Value"))) %>%
      tidyr::pivot_wider(names_from = "Variable", values_from = "Value")
  }

  # Reconstruct fixedData if missing
  if (is.null(fixedData)) {
    factor_part <- NULL
    numeric_part <- NULL
    if (!is.null(model$data$Factor$X)) {
      factor_part <- model$data$Factor$X
      factor_part[[idVar]] <- model$data$Factor$id
    }
    if (!is.null(model$data$Numeric$X)) {
      numeric_part <- model$data$Numeric$X
      numeric_part[[idVar]] <- model$data$Numeric$id
    }
    if (!is.null(factor_part) && !is.null(numeric_part)) {
      fixedData <- dplyr::full_join(factor_part, numeric_part, by = idVar)
    } else if (!is.null(factor_part)) {
      fixedData <- factor_part
    } else if (!is.null(numeric_part)) {
      fixedData <- numeric_part
    }
  }

  # Build grid of values
  if (!is.null(values)) grid_type_used <- "manual"
  else {
    grid_type <- match.arg(grid_type)
    grid_type_used <- grid_type
    if (var_type == "numeric") {
      x <- fixedData[[var_name]]
      values <- switch(
        grid_type,
        regular = seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out = grid_size),
        observed = sort(unique(x)),
        quantile = unique(as.numeric(stats::quantile(x, probs = seq(0,1,length.out=grid_size))))
      )
    } else if (var_type == "factor") values <- unique(fixedData[[var_name]])
  }

  # Single prediction function
  compute_single <- function(i, value_i=NULL, func_label=NULL, func=NULL) {
    fixedData_rep <- fixedData
    timeData_rep <- timeData
    fixedData_rep$replicate_id <- i
    if (!is.null(timeData_rep)) timeData_rep$replicate_id <- i

    if (var_type == "longitudinal") {
      if (is.null(func)) stop("For longitudinal variables, provide a function.")
      timeData_mod <- timeData_rep
      timeData_mod[[var_name]] <- unlist(lapply(split(timeData_mod[[var_name]], timeData_mod[[idVar]]), func))
      pred_dyn <- predict(model, fixedData = fixedData_rep, timeData = timeData_mod,
                          idVar = idVar, timeVar = timeVar, t0 = t0)
    } else {
      fixedData_rep[[var_name]] <- value_i
      pred_dyn <- predict(model, fixedData = fixedData_rep, timeData = timeData_rep,
                          idVar = idVar, timeVar = timeVar, t0 = t0)
    }

    # Format predictions
    if (model_type == "surv") {
      df <- as.data.frame(pred_dyn$pred_indiv)
      df$id <- rownames(df)
      df$replicate_id <- i
      if (var_type != "longitudinal") df[[var_name]] <- value_i
      if (var_type == "longitudinal") df$transformation <- func_label
      df_long <- df %>% tidyr::pivot_longer(cols = tidyselect::starts_with("V"),
                                            names_to = "time", values_to = "ice_value")
      df_long$time <- as.numeric(rep(pred_dyn$times, length.out = nrow(df_long)))
    } else if (model_type == "numeric") {
      df_long <- data.frame(
        id = names(pred_dyn$pred_indiv),
        ice_value = as.numeric(pred_dyn$pred_indiv),
        replicate_id = i
      )
      if (var_type != "longitudinal") df_long[[var_name]] <- value_i
      if (var_type == "longitudinal") df_long$transformation <- func_label
    } else if (model_type == "factor") {
      tree_mat <- pred_dyn$pred_indiv_tree
      n_indiv <- ncol(tree_mat)
      if (is.null(target_class)) {
        classes <- sort(unique(as.vector(tree_mat)))
        df_list <- lapply(classes, function(cls) {
          prob <- numeric(n_indiv)
          for (j in seq_len(n_indiv)) prob[j] <- mean(tree_mat[, j] == cls)
          df_cls <- data.frame(
            id = names(pred_dyn$pred_indiv),
            ice_value = prob,
            replicate_id = i,
            target_class = cls
          )
          if (var_type != "longitudinal") df_cls[[var_name]] <- value_i
          if (var_type == "longitudinal") df_cls$transformation <- func_label
          return(df_cls)
        })
        df_long <- dplyr::bind_rows(df_list)
      } else {
        prob <- numeric(n_indiv)
        for (j in seq_len(n_indiv)) prob[j] <- mean(tree_mat[, j] == target_class)
        df_long <- data.frame(
          id = names(pred_dyn$pred_indiv),
          ice_value = prob,
          replicate_id = i,
          target_class = target_class
        )
        if (var_type != "longitudinal") df_long[[var_name]] <- value_i
        if (var_type == "longitudinal") df_long$transformation <- func_label
      }
    }
    return(df_long)
  }

  # Parallel or sequential computation
  all_preds <- list()
  if (var_type == "longitudinal") {
    if (is.null(values)) stop("Provide a named list of functions for longitudinal variables.")
    for (i in seq_along(values)) {
      func_label <- names(values)[i]
      func <- values[[i]]
      all_preds[[func_label]] <- compute_single(i, func_label = func_label, func = func)
    }
  } else {
    if (.Platform$OS.type == "windows" && ncores > 1) {
      cl <- parallel::makeCluster(ncores)
      parallel::clusterExport(cl, varlist = c("model","fixedData","timeData","var_name",
                                              "t0","idVar","timeVar","values","target_class","compute_single"),
                              envir = environment())
      all_preds <- parallel::parLapply(cl, seq_along(values),
                                       function(i) compute_single(i, value_i = values[i]))
      parallel::stopCluster(cl)
    } else if (ncores > 1) {
      all_preds <- parallel::mclapply(seq_along(values),
                                      function(i) compute_single(i, value_i = values[i]),
                                      mc.cores = ncores)
    } else {
      all_preds <- lapply(seq_along(values), function(i) compute_single(i, value_i = values[i]))
    }
  }

  # Combine predictions
  df_long_all <- dplyr::bind_rows(all_preds)

  # Compute PDP
  group_vars <- c()
  if ("time" %in% colnames(df_long_all)) group_vars <- c(group_vars, "time")
  if ("target_class" %in% colnames(df_long_all)) group_vars <- c(group_vars, "target_class")
  if (var_type == "longitudinal") group_vars <- c(group_vars, "transformation") else group_vars <- c(group_vars, var_name)

  pdp <- df_long_all %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      pdp_value = mean(ice_value, na.rm = TRUE),
      pdp_sd = stats::sd(ice_value, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      alpha = 1 - conf_level,
      t_val = stats::qt(1 - alpha/2, df = pmax(n-1,1)),
      lower = pdp_value - t_val * pdp_sd / sqrt(n),
      upper = pdp_value + t_val * pdp_sd / sqrt(n)
    ) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::select(-n, -alpha, -t_val)

  # Return S3 object
  result <- list(
    ice = df_long_all,
    pdp = pdp,
    model_type = model_type,
    var_type = var_type,
    grid = list(
      grid_type = grid_type_used,
      grid_size = length(values)
    )
  )
  class(result) <- "dynforestpdp"
  return(result)
}
