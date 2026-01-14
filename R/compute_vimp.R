#' Compute the variable importance (VIMP) statistic for a dynforest object
#'
#' This function computes the importance of predictors for a dynforest object.
#' Longitudinal markers can be permuted either by shuffling values or swapping entire trajectories between subjects.
#'
#' @param dynforest_obj A \code{dynforest} object
#' @param IBS.min (Survival only) Minimal time to compute Integrated Brier Score. Default is 0
#' @param IBS.max (Survival only) Maximal time to compute Integrated Brier Score. Default is maximal observed time
#' @param ncores Number of cores to use for parallel computation. Default: all cores minus 1
#' @param seed Random seed for reproducibility. Default: 1234
#' @param permute_trajectory Logical. If TRUE, longitudinal markers are permuted by swapping trajectories between subjects; if FALSE, markers are permuted independently. Default is FALSE
#'
#' @return A list containing:
#' \tabular{ll}{
#'   \code{Inputs} \tab Names of predictors separated by type (Longitudinal, Numeric, Factor) \cr
#'   \code{Importance} \tab List of VIMP values for Longitudinal, Numeric, and Factor predictors \cr
#'   \code{tree_oob_err} \tab OOB error for each tree \cr
#'   \code{IBS.range} \tab Vector with IBS min and max
#' }
#' @export
compute_vimp <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                         ncores = NULL, seed = 1234, permute_trajectory = FALSE) {

  # Check input object
  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = paste0("You supplied a {.cls {class(dynforest_obj)}} object")
    ))
  }

  # Set IBS.max for survival outcome if missing
  if (dynforest_obj$type == "surv" && is.null(IBS.max)) {
    IBS.max <- max(dynforest_obj$data$Y$Y[, 1])
  }

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

  # Force all IDs as factors
  Longitudinal$id <- as.factor(Longitudinal$id)
  if (!is.null(Numeric)) Numeric$id <- as.factor(Numeric$id)
  if (!is.null(Factor)) Factor$id <- as.factor(Factor$id)

  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
  pbapply::pboptions(type = "none")

  # Compute baseline OOB errors for all trees
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(pbapply)
    library(doParallel)
    library(doRNG)
  })

  tree_oob_err <- pbapply::pbsapply(1:ntree, function(i) {
    DynForest:::OOB.tree(rf$rf[, i], Longitudinal, Numeric, Factor, Y,
                         timeVar, IBS.min, IBS.max, cause = rf$cause)
  }, cl = cl)
  suppressWarnings(parallel::stopCluster(cl))

  # Initialize importance vectors
  Importance <- list(Longitudinal = NULL, Numeric = NULL, Factor = NULL)

  # Compute VIMP for longitudinal predictors
  if ("Longitudinal" %in% Inputs) {
    marker_names <- colnames(Longitudinal$X)
    Importance$Longitudinal <- numeric(length(marker_names))
    names(Importance$Longitudinal) <- marker_names
    all_ids <- unique(Longitudinal$id)

    for (p in seq_along(marker_names)) {
      set.seed(seed + p)

      # Shuffle values independently or swap full trajectories
      if (!permute_trajectory) {
        Longitudinal.perm <- Longitudinal
        Longitudinal.perm$X[, p] <- sample(na.omit(Longitudinal$X[, p]),
                                           size = nrow(Longitudinal$X),
                                           replace = TRUE)
      } else {
        df_long_perm <- data.frame()
        for (idA in all_ids) {
          idB_choices <- setdiff(all_ids, idA)
          if (length(idB_choices) == 0) next
          idB <- sample(idB_choices, 1)
          df_long_perm <- rbind(df_long_perm,
                                permute_longitudinal_group(
                                  Longitudinal = Longitudinal,
                                  id_var = idVar,
                                  time_var = timeVar,
                                  marker_indices = p,
                                  idA = idA,
                                  idB = idB
                                ))
        }

        remaining_ids <- setdiff(all_ids, unique(df_long_perm$id))
        if (length(remaining_ids) > 0) {
          df_rest <- data.frame(
            id = Longitudinal$id[Longitudinal$id %in% remaining_ids],
            time = Longitudinal$time[Longitudinal$id %in% remaining_ids],
            Longitudinal$X[Longitudinal$id %in% remaining_ids, , drop = FALSE]
          )
          df_long_perm <- rbind(df_long_perm, df_rest)
        }

        df_long_perm <- df_long_perm[order(df_long_perm$id, df_long_perm$time), ]
        Longitudinal.perm <- list(
          type = "Longitudinal",
          X = df_long_perm[, marker_names, drop = FALSE],
          id = df_long_perm$id,
          time = df_long_perm$time,
          model = Longitudinal$model
        )
      }

      errs <- numeric(ntree)
      for (k in 1:ntree) {
        errs[k] <- DynForest:::OOB.tree(rf$rf[, k], Longitudinal.perm, Numeric, Factor, Y,
                                        timeVar, IBS.min, IBS.max, cause = rf$cause)
      }
      Importance$Longitudinal[p] <- mean(errs - tree_oob_err)
    }
  }

  # Compute VIMP for numeric predictors
  if ("Numeric" %in% Inputs) {
    library(foreach)
    library(doRNG)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    set.seed(seed)

    Importance$Numeric <- foreach::foreach(
      p = 1:ncol(Numeric$X),
      .combine = "c",
      .options.RNG = seed
    ) %dorng% {
      Numeric.perm <- Numeric
      Numeric.perm$X[, p] <- sample(Numeric$X[, p])

      errs <- numeric(ntree)
      for (k in 1:ntree) {
        errs[k] <- DynForest:::OOB.tree(rf$rf[, k], Longitudinal, Numeric.perm, Factor, Y,
                                        timeVar, IBS.min, IBS.max, cause = rf$cause)
      }
      mean(errs - tree_oob_err)
    }
    suppressWarnings(parallel::stopCluster(cl))
  }

  # Compute VIMP for factor predictors
  if ("Factor" %in% Inputs) {
    library(foreach)
    library(doRNG)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    set.seed(seed)

    Importance$Factor <- foreach::foreach(
      p = 1:ncol(Factor$X),
      .combine = "c",
      .options.RNG = seed
    ) %dorng% {
      Factor.perm <- Factor
      Factor.perm$X[, p] <- sample(Factor$X[, p])

      errs <- numeric(ntree)
      for (k in 1:ntree) {
        errs[k] <- DynForest:::OOB.tree(rf$rf[, k], Longitudinal, Numeric, Factor.perm, Y,
                                        timeVar, IBS.min, IBS.max, cause = rf$cause)
      }
      mean(errs - tree_oob_err)
    }
    suppressWarnings(parallel::stopCluster(cl))
  }

  # Return results
  out <- list(
    Inputs = dynforest_obj$Inputs,
    Importance = Importance,
    tree_oob_err = tree_oob_err,
    IBS.range = c(IBS.min, IBS.max)
  )
  class(out) <- "dynforestvimp"
  return(out)
}
