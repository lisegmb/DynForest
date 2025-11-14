#' Compute the importance of variables (VIMP) statistic
#'
#' @param dynforest_obj dynforest object resulting from \code{dynforest()}
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default is 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default is the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default is the number of available cores minus 1.
#' @param seed Seed for reproducibility
#' @param permute_trajectory Logical. If TRUE, permutes longitudinal variables by exchanging patient trajectories instead of shuffling raw values. Default is FALSE.
#'
#' @importFrom methods is
#' @import doRNG
#' @export
#'
#' @return A list with VIMP results, including \code{Inputs}, \code{Importance}, \code{tree_oob_err}, and \code{IBS.range}.
#'
#' @seealso [dynforest()]
compute_vimp2 <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                          ncores = NULL, seed = 1234, permute_trajectory = FALSE) {

  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You've supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }

  if (dynforest_obj$type == "surv" && is.null(IBS.max)) {
    IBS.max <- max(dynforest_obj$data$Y$Y[, 1])
  }

  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  timeVar <- rf$timeVar
  idVar <- rf$idVar
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs[!sapply(rf$Inputs, is.null)])

  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 1
  }

  pbapply::pboptions(type = "none")

  suppressWarnings({
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, {
      library(pbapply)
      library(doParallel)
      library(doRNG)
    })
  })

  # OOB error d'origine
  tree_oob_err <- pbapply::pbsapply(1:ntree, function(i) {
    DynForest:::OOB.tree(rf$rf[, i], Longitudinal, Numeric, Factor, Y, timeVar,
                         IBS.min, IBS.max, cause = rf$cause)
  }, cl = cl)

  suppressWarnings(parallel::stopCluster(cl))

  Importance <- list(Longitudinal = NULL, Numeric = NULL, Factor = NULL)

  ### === LONGITUDINAL VARIABLES === ###
  if ("Longitudinal" %in% Inputs) {
    marker_names <- colnames(Longitudinal$X)
    Importance$Longitudinal <- numeric(length(marker_names))
    names(Importance$Longitudinal) <- marker_names
    all_ids <- unique(Longitudinal$id)

    for (p in seq_along(marker_names)) {
      marker <- marker_names[p]
      set.seed(seed + p)

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
                                permute_longitudinal_marker_patient(
                                  Longitudinal = Longitudinal,
                                  id_var = idVar,
                                  time_var = timeVar,
                                  marker_index = p,
                                  idA = idA,
                                  idB = idB,
                                  seed = seed + p + as.integer(idA)
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

  ### === NUMERIC VARIABLES === ###
  if ("Numeric" %in% Inputs) {
    suppressWarnings({
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      parallel::clusterEvalQ(cl, library(doRNG))
    })

    Importance$Numeric <- foreach::foreach(p = 1:ncol(Numeric$X),
                                           .combine = "c", .options.RNG = seed) %dorng% {
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

  ### === FACTOR VARIABLES === ###
  if ("Factor" %in% Inputs) {
    suppressWarnings({
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      parallel::clusterEvalQ(cl, library(doRNG))
    })

    Importance$Factor <- foreach::foreach(p = 1:ncol(Factor$X),
                                          .combine = "c", .options.RNG = seed) %dorng% {
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

  out <- list(
    Inputs = dynforest_obj$Inputs,
    Importance = Importance,
    tree_oob_err = tree_oob_err,
    IBS.range = c(IBS.min, IBS.max)
  )
  class(out) <- "dynforestvimp"
  return(out)
}
