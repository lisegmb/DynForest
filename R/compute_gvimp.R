#' Compute the grouped variable importance (gVIMP) statistic for a dynforest object
#'
#' This function computes the importance of predefined variable groups for a dynforest object.
#' Longitudinal markers can be permuted either by shuffling values or swapping trajectories between subjects.
#'
#' @inheritParams compute_vimp
#' @param group A named list of variable groups. Each element contains variable names for that group.
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#'    \code{Inputs} \tab Names of predictors separated by type (Longitudinal, Numeric, Factor) \cr
#'    \code{group} \tab List of each group defined in the \code{group} argument \cr
#'    \code{gVIMP} \tab Numeric vector containing gVIMP for each group \cr
#'    \code{tree_oob_err} \tab OOB error for each tree used to compute the VIMP statistic \cr
#'    \code{IBS.range} \tab Vector containing the IBS min and max values
#' }
#'
#' @export
compute_gvimp <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                           group = NULL, ncores = NULL, seed = 1234,
                           permute_trajectory = FALSE) {

  # Check that input object is a dynforest
  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort("{.var dynforest_obj} must be a dynforest object")
  }

  # For survival objects, set IBS.max if missing
  if (dynforest_obj$type == "surv" && is.null(IBS.max)) {
    IBS.max <- max(dynforest_obj$data$Y$Y[,1])
  }

  if (is.null(group)) stop("'group' argument cannot be NULL!")

  # Extract data from dynforest object
  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  timeVar <- rf$timeVar
  idVar <- rf$idVar
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs[!sapply(rf$Inputs, is.null)])

  unique_ids <- unique(Longitudinal$id)  # unique subject IDs

  # Force all IDs to factors to prevent NA issues
  if (!is.null(Longitudinal)) Longitudinal$id <- as.factor(Longitudinal$id)
  if (!is.null(Numeric)) Numeric$id <- as.factor(Numeric$id)
  if (!is.null(Factor)) Factor$id <- as.factor(Factor$id)

  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
  pbapply::pboptions(type = "none")

  # Compute baseline OOB error for all trees
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(pbapply); library(doParallel) })

  tree_oob_err <- pbapply::pbsapply(1:ntree, function(i) {
    DynForest:::OOB.tree(rf$rf[, i], Longitudinal, Numeric, Factor, Y,
                         timeVar, IBS.min, IBS.max, cause = rf$cause)
  }, cl = cl)
  parallel::stopCluster(cl)

  gVIMP <- numeric(length(group))
  names(gVIMP) <- names(group)
  set.seed(seed)

  # Create a mapping for permutation of IDs
  id_mapping <- setNames(sample(unique_ids, length(unique_ids), replace = FALSE), unique_ids)

  # Loop over each group to compute gVIMP
  for (g in seq_along(group)) {

    g_vars <- group[[g]]

    Longitudinal_perm <- Longitudinal
    Numeric_perm <- Numeric
    Factor_perm <- Factor

    # Permute longitudinal markers
    if ("Longitudinal" %in% Inputs) {
      markers_in_group <- intersect(g_vars, colnames(Longitudinal$X))

      if (length(markers_in_group) > 0) {
        if (permute_trajectory) {
          # Permute entire trajectories between subjects
          df_long_perm <- data.frame()
          for (id in unique_ids) {
            mapped_id <- id_mapping[as.character(id)]
            df_long_perm <- rbind(
              df_long_perm,
              permute_longitudinal_group(
                Longitudinal = Longitudinal,
                id_var = idVar,
                time_var = timeVar,
                marker_indices = which(colnames(Longitudinal$X) %in% markers_in_group),
                idA = id,
                idB = mapped_id
              )
            )
          }
          df_long_perm <- df_long_perm[order(df_long_perm$id, df_long_perm$time), ]
          Longitudinal_perm <- list(
            type  = "Longitudinal",
            X     = df_long_perm[, colnames(Longitudinal$X), drop = FALSE],
            id    = df_long_perm$id,
            time  = df_long_perm$time,
            model = Longitudinal$model
          )
        } else {
          # Permute values individually
          for (marker in markers_in_group) {
            Longitudinal_perm$X[, marker] <- sample(Longitudinal$X[, marker])
          }
        }
      }
    }

    # Permute numeric variables (1 row = 1 patient)
    if ("Numeric" %in% Inputs) {
      vars_numeric <- g_vars[g_vars %in% colnames(Numeric$X)]
      perm_idx <- match(id_mapping[as.character(Numeric$id)], Numeric$id)
      if (any(is.na(perm_idx))) stop("Permutation Numeric: NA in perm_idx (id_mapping incorrect)")
      for (var in vars_numeric) {
        Numeric_perm$X[, var] <- Numeric$X[perm_idx, var]
      }
    }

    # Permute factor variables (1 row = 1 patient)
    if ("Factor" %in% Inputs) {
      vars_factor <- g_vars[g_vars %in% colnames(Factor$X)]
      perm_idx <- match(id_mapping[as.character(Factor$id)], Factor$id)
      if (any(is.na(perm_idx))) stop("Permutation Factor: NA in perm_idx (id_mapping incorrect)")
      for (var in vars_factor) {
        Factor_perm$X[, var] <- Factor$X[perm_idx, var]
      }
    }

    # Compute OOB error for permuted data
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, library(doRNG))

    errs <- foreach::foreach(k = 1:ntree, .combine = "c", .options.RNG = seed) %dorng% {
      DynForest:::OOB.tree(rf$rf[, k], Longitudinal_perm, Numeric_perm, Factor_perm,
                           Y, timeVar, IBS.min, IBS.max, cause = rf$cause)
    }
    parallel::stopCluster(cl)

    # Compute gVIMP for this group
    gVIMP[g] <- mean(errs - tree_oob_err)
  }

  # Return results as a dynforestgvimp object
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
