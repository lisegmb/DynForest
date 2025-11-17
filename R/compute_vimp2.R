#' Compute the importance of variables (VIMP) statistic
#'
#' @param dynforest_obj dynforest_obj \code{dynforest} object.
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer minus 1.
#' @param seed Seed to replicate results. Default is 1234.
#' @param permute_trajectory Logical. If TRUE, longitudinal markers are permuted by swapping trajectories between subjects; if FALSE, markers are shuffled independently. Default is FALSE.
#'
#' @importFrom methods is
#' @import doRNG
#'
#' @return \code{compute_vimp2()} function returns a list with the following elements:\tabular{ll}{
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Importance} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains a numeric vector of VIMP statistic predictor in \code{Inputs} value \cr
#'    \tab \cr
#'    \code{tree_oob_err} \tab A numeric vector containing the OOB error for each tree needed to compute the VIMP statistic \cr
#'    \tab \cr
#'    \code{IBS.range} \tab A vector containing the IBS min and max \cr
#'    \tab \cr
#'    \code{permute_trajectory} \tab Logical flag indicating whether longitudinal trajectories were permuted as whole trajectories or shuffled independently \cr
#' }
#'
#' @export
#'
#' @seealso [dynforest()]
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Transform longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
#'
#' # Sample 100 subjects
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id %in% id_sample)
#' pbc2_train <- pbc2[id_row, ]
#'
#' # Build longitudinal data
#' timeData_train <- pbc2_train[, c("id", "time",
#'                                  "serBilir", "SGOT",
#'                                  "albumin", "alkaline")]
#'
#' # Longitudinal models
#' timeVarModel <- list(
#'   serBilir = list(fixed = serBilir ~ time, random = ~ time),
#'   SGOT = list(fixed = SGOT ~ time + I(time^2), random = ~ time + I(time^2)),
#'   albumin = list(fixed = albumin ~ time, random = ~ time),
#'   alkaline = list(fixed = alkaline ~ time, random = ~ time)
#' )
#'
#' # Fixed data
#' fixedData_train <- unique(pbc2_train[, c("id", "age", "drug", "sex")])
#'
#' # Outcome
#' Y <- list(type = "surv", Y = unique(pbc2_train[, c("id", "years", "event")]))
#'
#' # Run dynforest
#' res_dyn <- dynforest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Compute VIMP statistic with trajectory permutation
#' res_dyn_VIMP2 <- compute_vimp2(dynforest_obj = res_dyn, ncores = 2, seed = 1234,
#'                                permute_trajectory = TRUE)
#' }

compute_vimp2 <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                          ncores = NULL, seed = 1234, permute_trajectory = FALSE) {

  if (!methods::is(dynforest_obj, "dynforest")) {
    cli::cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }

  # For survival models, set IBS.max if NULL
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

  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
  pbapply::pboptions(type = "none")

  # --- Compute baseline OOB errors ---
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

  Importance <- list(Longitudinal = NULL, Numeric = NULL, Factor = NULL)

  # --- Longitudinal variables ---
  if ("Longitudinal" %in% Inputs) {
    marker_names <- colnames(Longitudinal$X)
    Importance$Longitudinal <- numeric(length(marker_names))
    names(Importance$Longitudinal) <- marker_names
    all_ids <- unique(Longitudinal$id)

    for (p in seq_along(marker_names)) {
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
                                permute_longitudinal_group(
                                  Longitudinal = Longitudinal,
                                  id_var = idVar,
                                  time_var = timeVar,
                                  marker_indices = p,
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

  # --- Numeric variables ---
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

  # --- Factor variables ---
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

  out <- list(
    Inputs = dynforest_obj$Inputs,
    Importance = Importance,
    tree_oob_err = tree_oob_err,
    IBS.range = c(IBS.min, IBS.max)
  )
  class(out) <- "dynforestvimp"
  return(out)
}
