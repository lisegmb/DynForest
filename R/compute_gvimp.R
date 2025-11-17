#' Compute the grouped importance of variables (gVIMP) statistic
#'
#' This function computes the importance of predefined variable groups for a dynforest object.
#' Longitudinal markers can be permuted either by shuffling values or swapping trajectories between subjects.
#'
#' @inheritParams compute_vimp
#' @param group A named list of variable groups. Each element contains variable names for that group.
#'
#' @importFrom methods is
#' @import doRNG
#'
#' @return \code{compute_gvimp()} function returns a list with the following elements:\tabular{ll}{
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{group} \tab A list of each group defined in \code{group} argument \cr
#'    \tab \cr
#'    \code{gVIMP} \tab A numeric vector containing the gVIMP for each group defined in \code{group} argument \cr
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
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id %in% id_sample)
#' pbc2_train <- pbc2[id_row, ]
#' timeData_train <- pbc2_train[, c("id","time","serBilir","SGOT","albumin","alkaline")]
#' timeVarModel <- list(
#'   serBilir = list(fixed = serBilir ~ time, random = ~ time),
#'   SGOT = list(fixed = SGOT ~ time + I(time^2), random = ~ time + I(time^2)),
#'   albumin = list(fixed = albumin ~ time, random = ~ time),
#'   alkaline = list(fixed = alkaline ~ time, random = ~ time)
#' )
#' fixedData_train <- unique(pbc2_train[, c("id","age","drug","sex")])
#' Y <- list(type = "surv", Y = unique(pbc2_train[, c("id","years","event")]))
#' res_dyn <- dynforest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#' res_dyn_gVIMP <- compute_gvimp(
#'   dynforest_obj = res_dyn,
#'   group = list(group1 = c("serBilir","SGOT"),
#'                group2 = c("albumin","alkaline")),
#'   ncores = 2, seed = 1234, permute_trajectory = TRUE
#' )
#' }

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
  all_ids <- unique(Longitudinal$id)

  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
  pbapply::pboptions(type = "none")

  # Compute baseline OOB error
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

  # Loop over groups
  for (g in seq_along(group)) {
    g_vars <- group[[g]]

    Longitudinal_perm <- Longitudinal
    Numeric_perm <- Numeric
    Factor_perm <- Factor

    # Permute longitudinal markers
    if ("Longitudinal" %in% Inputs) {
      markers_in_group <- g_vars[g_vars %in% colnames(Longitudinal$X)]
      if (permute_trajectory) {
        df_long_perm <- data.frame()
        for (idA in all_ids) {
          idB_choices <- setdiff(all_ids, idA)
          if (length(idB_choices) == 0) next
          idB <- sample(idB_choices, 1)
          df_long_perm <- rbind(df_long_perm,
                                permute_longitudinal_group(Longitudinal, idVar, timeVar,
                                                           which(colnames(Longitudinal$X) %in% markers_in_group),
                                                           idA, idB, seed))
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

    # Permute numeric/factor variables consistently
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

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, library(doRNG))

    errs <- foreach::foreach(k = 1:ntree, .combine = "c", .options.RNG = seed) %dorng% {
      DynForest:::OOB.tree(rf$rf[, k], Longitudinal_perm, Numeric_perm, Factor_perm,
                           Y, timeVar, IBS.min, IBS.max, cause = rf$cause)
    }
    parallel::stopCluster(cl)

    gVIMP[g] <- mean(errs - tree_oob_err)
  }

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
