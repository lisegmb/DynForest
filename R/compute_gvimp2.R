#' Compute grouped variable importance (gVIMP) for dynforest models
#'
#' This function computes the importance of predefined variable groups.
#' Longitudinal markers can be permuted either by shuffling values or swapping trajectories between subjects.
#'
#' @param dynforest_obj A dynforest object from \code{dynforest()}.
#' @param IBS.min Minimal time for Integrated Brier Score computation (survival only). Default is 0.
#' @param IBS.max Maximal time for Integrated Brier Score computation (survival only). Default is the maximum observed time.
#' @param group A named list of variable groups. Each element contains variable names for that group.
#' @param ncores Number of cores for parallel computation. Default is available cores minus 1.
#' @param seed Random seed for reproducibility.
#' @param permute_trajectory If TRUE, longitudinal markers are permuted by swapping trajectories instead of shuffling values. Default is FALSE.
#'
#' @return A list containing Inputs, groups, gVIMP values, tree OOB errors, and IBS range.
#' @export
compute_gvimp2 <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                           group = NULL, ncores = NULL, seed = 1234,
                           permute_trajectory = FALSE) {

  # --- Validation checks ---
  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You provided a {.cls {class(dynforest_obj)}} object"
    ))
  }
  if (dynforest_obj$type == "surv" && is.null(IBS.max)) {
    IBS.max <- madynforest_obj$data$Y$Y[,1])
  }
  if (is.null(group)) stop("'group' argument cannot be NULL!")

  # --- Extract model data ---
  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  timeVar <- rf$timeVar
  idVar <- rf$idVar
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs[!sapply(rf$Inputs, is.null)])
  all_ids <- unique(Longitudinal$id)

  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
  pbapply::pboptions(type = "none")

  # --- Compute baseline OOB error ---
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(pbapply); library(doParallel) })

  tree_oob_err <- pbapply::pbsapply(1:ntree, function(i) {
    DynForest:::OOB.tree(rf$rf[, i], Longitudinal, Numeric, Factor, Y,
                         timeVar, IBS.min, IBS.max, cause = rf$cause)
  }, cl = cl)
  parallel::stopCluster(cl)

  # --- Initialize gVIMP vector ---
  gVIMP <- numeric(length(group))
  names(gVIMP) <- names(group)
  set.seed(seed)

  # --- Loop over each group ---
  for (g in seq_along(group)) {
    g_vars <- group[[g]]

    # Copy data structures to permute
    Longitudinal_perm <- Longitudinal
    Numeric_perm <- Numeric
    Factor_perm <- Factor

    # --- Permute longitudinal markers if any ---
    if ("Longitudinal" %in% Inputs) {
      markers_in_group <- g_vars[g_vars %in% colnames(Longitudinal$X)]

      if (permute_trajectory) {
        df_long_perm <- data.frame()
        for (idA in all_ids) {
          idB_choices <- setdiff(all_ids, idA)
          if (length(idB_choices) == 0) next
          idB <- sample(idB_choices, 1)

          df_long_perm <- rbind(
            df_long_perm,
            permute_longitudinal_group(
              Longitudinal = Longitudinal,
              id_var = idVar,
              time_var = timeVar,
              marker_indices = which(colnames(Longitudinal$X) %in% markers_in_group),
              idA = idA,
              idB = idB,
              seed = seed
            )
          )
        }
        df_long_perm <- df_long_perm[order(df_long_perm$id, df_long_perm$time), ]
        Longitudinal_perm <- list(
          type = "Longitudinal",
          X = df_long_perm[, colnames(Longitudinal$X), drop = FALSE],
          id = df_long_perm$id,
          time = df_long_perm$time,
          model = Longitudinal$model
        )
      } else {
        for (marker in markers_in_group) {
          Longitudinal_perm$X[, marker] <- sample(Longitudinal_perm$X[, marker])
        }
      }
    }

    # --- Permute numeric and factor variables using consistent mapping ---
    id_mapping <- setNames(sample(all_ids, length(all_ids), replace = FALSE), all_ids)

    if ("Numeric" %in% Inputs) {
      vars_numeric <- g_vars[g_vars %in% colnames(Numeric$X)]
      for (var in vars_numeric) {
        Numeric_perm$X[, var] <- Numeric$X[match(id_mapping[as.character(Longitudinal$id)], Longitudinal$id), var]
      }
    }

    if ("Factor" %in% Inputs) {
      vars_factor <- g_vars[g_vars %in% colnames(Factor$X)]
      for (var in vars_factor) {
        Factor_perm$X[, var] <- Factor$X[match(id_mapping[as.character(Longitudinal$id)], Longitudinal$id), var]
      }
    }

    # --- Compute OOB errors on permuted data ---
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, library(doRNG))

    errs <- foreach::foreach(k = 1:ntree, .combine = "c", .options.RNG = seed) %dorng% {
      DynForest:::OOB.tree(rf$rf[, k], Longitudinal_perm, Numeric_perm, Factor_perm,
                           Y, timeVar, IBS.min, IBS.max, cause = rf$cause)
    }
    parallel::stopCluster(cl)

    # --- Compute gVIMP for the current group ---
    gVIMP[g] <- mean(errs - tree_oob_err)
  }

  # --- Return results ---
  out <- list(
    Inputs = dynforest_obj$Inputs,
    group = group,
    gVIMP = gVIMP,
    tree_oob_err = tree_oob_err,
    IBS.range = c(IBS.min, IBS.max)
  )
  class(out) <- "dynforestgvimp"
  return(out)
}

