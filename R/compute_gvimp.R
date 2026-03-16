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

  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort("{.var dynforest_obj} must be a dynforest object")
  }

  if (dynforest_obj$type == "surv" && is.null(IBS.max)) {
    IBS.max <- max(dynforest_obj$data$Y$Y[,1])
  }

  if (is.null(group)) stop("'group' argument cannot be NULL!")

  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  timeVar <- rf$timeVar
  idVar <- rf$idVar
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs[!sapply(rf$Inputs, is.null)])
  unique_ids <- unique(Longitudinal$id)

  if (!is.null(Longitudinal)) Longitudinal$id <- as.factor(Longitudinal$id)
  if (!is.null(Numeric)) Numeric$id <- as.factor(Numeric$id)
  if (!is.null(Factor)) Factor$id <- as.factor(Factor$id)

  if (is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)

  pbapply::pboptions(type = "none")

  # Baseline OOB errors
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  tree_oob_err <- pbapply::pbsapply(1:ntree, function(i) {
    OOB.tree(rf$rf[, i], Longitudinal, Numeric, Factor, Y,
             timeVar, IBS.min, IBS.max, cause = rf$cause)
  }, cl = cl)

  parallel::stopCluster(cl)

  gVIMP <- numeric(length(group))
  names(gVIMP) <- names(group)

  set.seed(seed)

  id_mapping <- setNames(
    sample(unique_ids, length(unique_ids), replace = FALSE),
    unique_ids
  )

  for (g in seq_along(group)) {

    g_vars <- group[[g]]

    Longitudinal_perm <- Longitudinal
    Numeric_perm <- Numeric
    Factor_perm <- Factor

    # Permutation Longitudinal
    if ("Longitudinal" %in% Inputs) {

      markers_in_group <- intersect(g_vars, colnames(Longitudinal$X))

      if (length(markers_in_group) > 0) {

        if (permute_trajectory) {

          df_long_perm <- do.call(rbind, lapply(unique_ids, function(id) {

            mapped_id <- id_mapping[as.character(id)]

            permute_longitudinal_group(
              Longitudinal,
              idVar,
              timeVar,
              marker_indices = which(colnames(Longitudinal$X) %in% markers_in_group),
              idA = id,
              idB = mapped_id
            )

          }))

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
            Longitudinal_perm$X[, marker] <- sample(Longitudinal$X[, marker])
          }

        }

      }
    }

    # Permutation Numeric
    if ("Numeric" %in% Inputs) {

      vars_numeric <- g_vars[g_vars %in% colnames(Numeric$X)]

      if (length(vars_numeric) > 0) {

        perm_idx <- match(id_mapping[as.character(Numeric$id)], Numeric$id)

        for (var in vars_numeric) {
          Numeric_perm$X[, var] <- Numeric$X[perm_idx, var]
        }

      }
    }

    # Permutation Factor
    if ("Factor" %in% Inputs) {

      vars_factor <- g_vars[g_vars %in% colnames(Factor$X)]

      if (length(vars_factor) > 0) {

        perm_idx <- match(id_mapping[as.character(Factor$id)], Factor$id)

        for (var in vars_factor) {
          Factor_perm$X[, var] <- Factor$X[perm_idx, var]
        }

      }
    }

    # Compute OOB errors
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    errs <- foreach::foreach(k = 1:ntree, .combine = "c") %dopar% {

      OOB.tree(
        rf$rf[, k],
        Longitudinal_perm,
        Numeric_perm,
        Factor_perm,
        Y,
        timeVar,
        IBS.min,
        IBS.max,
        cause = rf$cause
      )

    }

    parallel::stopCluster(cl)

    gVIMP[g] <- mean(errs - tree_oob_err)

  }

  out <- list(
    Inputs = dynforest_obj$Inputs,
    group = group,
    gVIMP = gVIMP,
    tree_oob_err = tree_oob_err,
    IBS.range = c(IBS.min, IBS.max),
    verbose = TRUE
  )

  class(out) <- "dynforestgvimp"

  return(out)
}
