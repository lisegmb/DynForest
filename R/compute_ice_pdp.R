#' Compute ICE and PDP for dynforest models with automatic data retrieval
#'
#' Compute Individual Conditional Expectation (ICE) and Partial Dependence Plot (PDP)
#' values for a specified predictor from a fitted \code{dynforest} model.
#' The function supports numeric, categorical, and longitudinal predictors. It can
#' reconstruct \code{timeData} and \code{fixedData} from the model object if they are not provided,
#' and attempts to auto-detect \code{idVar} and \code{timeVar} when not specified.
#'
#' @param timeData Optional \code{data.frame} with the subject ID, time variable, and time-dependent predictors.
#'                 If \code{NULL}, the function attempts to reconstruct it from \code{model$data$Longitudinal}.
#' @param fixedData Optional \code{data.frame} with the subject ID and time-fixed predictors (numeric or factor).
#'                  If \code{NULL}, it is reconstructed from \code{model$data$Factor} and \code{model$data$Numeric}.
#' @param model A fitted \code{dynforest} object (result of \code{dynforest()}).
#' @param var_name Character. Name of the predictor variable for which ICE/PDP are computed.
#' @param values For numeric predictors: a numeric vector of values to evaluate.
#'               For categorical predictors: a vector of factor levels (defaults to all observed levels).
#'               For longitudinal predictors: a **named list of functions** to summarize each subject's trajectory (required).
#' @param target_class Optional vector of classes to retain (used for classification models only).
#' @param grid_type Character. Type of grid for numeric variables (\code{"regular"}, \code{"observed"}, or \code{"quantile"}).
#' @param grid_size Integer. Number of points in the grid if \code{values} is not provided.
#' @param idVar Optional character. Name of the subject identifier variable. Auto-detected if \code{NULL}.
#' @param timeVar Optional character. Name of the time variable. Defaults to \code{model$timeVar} or \code{"time"}.
#' @param t0 Optional numeric. Landmark time for dynamic prediction; if \code{NULL}, predictions use all available information from time 0.
#' @param ncores Integer. Number of cores to use for parallel computation. Defaults to 1.
#' @param conf_level Numeric between 0 and 1. Confidence level for PDP intervals. Defaults to 0.95.
#'
#' @return An object of class \code{dynforestpdp}, which is a list containing:
#' \describe{
#'   \item{\code{ice}}{A \code{data.frame} of individual conditional expectation (ICE) values with columns such as
#'        \code{id}, \code{replicate_id}, \code{time}, \code{ice_value}, \code{transformation} (for longitudinal predictors),
#'        \code{var_name} (for fixed predictors), and \code{target_class} (for classification models).}
#'   \item{\code{pdp}}{A \code{data.frame} of partial dependence plot (PDP) values with columns
#'        \code{pdp_value}, \code{pdp_sd}, \code{lower}, \code{upper}, and optionally \code{time},
#'        \code{target_class}, or \code{transformation}.}
#'   \item{\code{model_type}}{Type of the model: \code{"numeric"}, \code{"factor"}, or \code{"surv"}.}
#'   \item{\code{var_type}}{Type of the predictor variable: \code{"numeric"}, \code{"factor"}, or \code{"longitudinal"}.}
#'   \item{\code{grid}}{A list containing \code{grid_type} and \code{grid_size} used for ICE/PDP computation.}
#' }
#'
#' @details
#' The function automatically handles three model types:
#' \itemize{
#'   \item \strong{surv} - longitudinal survival predictions (requires \code{id}, \code{time}, \code{ice_value}).
#'   \item \strong{factor} - classification with \code{target_class}.
#'   \item \strong{numeric} - regression.
#' }
#' For longitudinal predictors, \code{values} must be a named list of summary/transformation functions
#' (e.g., \code{list(mean = mean, slope = function(x) coef(lm(x ~ time))[2])}).
#'
#' @seealso \code{\link{dynforest}}, \code{\link{plot.dynforestpdp}} for plotting ICE/PDP curves.
#'
#' @importFrom dplyr bind_rows mutate select full_join group_by summarise
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect starts_with all_of
#' @importFrom stats quantile qt sd
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
      tidyr::pivot_longer(cols = tidyselect::all_of(value_cols), names_to = "Variable", values_to = "Value") %>%
      dplyr::select(dplyr::all_of(c(idVar, timeVar, "Variable", "Value"))) %>%
      tidyr::pivot_wider(names_from = "Variable", values_from = "Value")
  }

  # Reconstruct fixedData if missing
  if (is.null(fixedData)) {
    factor_part <- if (!is.null(model$data$Factor$X)) dplyr::mutate(model$data$Factor$X, !!idVar := model$data$Factor$id) else NULL
    numeric_part <- if (!is.null(model$data$Numeric$X)) dplyr::mutate(model$data$Numeric$X, !!idVar := model$data$Numeric$id) else NULL
    if (!is.null(factor_part) && !is.null(numeric_part)) fixedData <- dplyr::full_join(factor_part, numeric_part, by = idVar)
    else if (!is.null(factor_part)) fixedData <- factor_part
    else if (!is.null(numeric_part)) fixedData <- numeric_part
  }

  # Build grid of values for variable
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
      pred_dyn <- predict(model, fixedData = fixedData_rep, timeData = timeData_mod, idVar = idVar, timeVar = timeVar, t0 = t0)
    } else {
      fixedData_rep[[var_name]] <- value_i
      pred_dyn <- predict(model, fixedData = fixedData_rep, timeData = timeData_rep, idVar = idVar, timeVar = timeVar, t0 = t0)
    }

    # Format predictions
    if (model_type == "surv") {
      df <- as.data.frame(pred_dyn$pred_indiv)
      df$id <- rownames(df)
      df$replicate_id <- i
      if (var_type != "longitudinal") df[[var_name]] <- value_i
      if (var_type == "longitudinal") df$transformation <- func_label
      df_long <- df %>% tidyr::pivot_longer(cols = tidyselect::starts_with("V"), names_to = "time", values_to = "ice_value")
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

  # Parallel computation or sequential
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
      parallel::clusterExport(cl, varlist = c("model","fixedData","timeData","var_name","t0","idVar","timeVar","values","target_class","compute_single"), envir = environment())
      parallel::clusterEvalQ(cl, { library(DynForest); library(dplyr); library(tidyr) })
      all_preds <- parallel::parLapply(cl, seq_along(values), function(i) compute_single(i, value_i = values[i]))
      parallel::stopCluster(cl)
    } else if (ncores > 1) {
      all_preds <- parallel::mclapply(seq_along(values), function(i) compute_single(i, value_i = values[i]), mc.cores = ncores)
    } else {
      all_preds <- lapply(seq_along(values), function(i) compute_single(i, value_i = values[i]))
    }
  }

  # Combine all predictions
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
      pdp_sd = sd(ice_value, na.rm = TRUE),
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
