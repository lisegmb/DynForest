#' Compute Out-Of-Bag Prediction Error for a Tree
#'
#' This function computes the predictive performance of a single tree using
#' the Out-Of-Bag (OOB) samples. For survival outcomes, it returns the
#' Integrated Brier Score (IBS) across time, **weighted using Inverse
#' Probability of Censoring (IPCW)** to account for censored observations.
#' For numeric and factor outcomes, it computes the standard prediction error
#' (mean squared error or misclassification rate).
#'
#' @param tree A tree object resulting from the \code{Rtmax_surv} function.
#' @param Longitudinal A list of longitudinal predictors with elements:
#'   \code{X} (data frame of repeated measurements), \code{id} (vector of subject IDs),
#'   \code{time} (measurement times), and optionally \code{model}.
#' @param Numeric A list of numeric predictors with elements: \code{X} (data frame) and \code{id}.
#' @param Factor A list of factor predictors with elements: \code{X} (data frame) and \code{id}.
#' @param Y A list describing the outcome with elements: \code{type} ("surv", "numeric", or "factor"),
#'   \code{Y} (outcome values), and \code{id} (subject identifiers matching input IDs).
#' @param timeVar Character specifying the name of the time variable (if applicable).
#' @param IBS.min Minimum time for computing IBS (only for survival outcomes, default 0).
#' @param IBS.max Maximum time for computing IBS (only for survival outcomes, default = max observed time).
#' @param cause Integer specifying the event of interest for competing risks (default 1).
#'
#' @return A numeric value representing the average OOB prediction error:
#'   - Survival outcomes: Integrated Brier Score (IBS) using IPCW.
#'   - Numeric outcomes: Mean squared error.
#'   - Factor outcomes: Misclassification rate.
#'
#' @import pec
#' @import prodlim
#' @keywords internal
#' @export
OOB.tree <- function(tree, Longitudinal = NULL, Numeric = NULL, Factor = NULL, Y,
                     timeVar = NULL, IBS.min = 0, IBS.max = NULL, cause = 1){

  Inputs <- c(Longitudinal$type, Numeric$type, Factor$type)

  BOOT <- tree$boot
  OOB <- setdiff(unique(Y$id), BOOT)

  Numeric_current <- NULL
  Factor_current <- NULL
  Longitudinal_current <- NULL

  if (Y$type=="surv"){ # survival outcome
    allTimes <- sort(unique(c(0,Y$Y[,1])))

    if (is.null(IBS.max)){
      IBS.max <- max(allTimes)
    }

    Y.surv <- data.frame(id = Y$id,
                         time.event = Y$Y[,1],
                         event = ifelse(Y$Y[,2]==cause,1,0))

    Y.surv <- Y.surv[order(Y.surv$time.event, -Y.surv$event),]

    #### NEW #####

    allTimes_IBS <- allTimes[which(allTimes>=IBS.min&allTimes<=IBS.max)]

    Y.surv <- Y.surv[which(Y.surv$time.event>=IBS.min),]

    # KM estimate of the survival function for the censoring
    G <- pec::ipcw(formula = Surv(time.event, event) ~ 1,
                   data = Y.surv,
                   method = "marginal",
                   times=allTimes_IBS,
                   subjectTimes=allTimes_IBS)

    OOB_IBS <- sort(Y.surv$id[which(Y.surv$id%in%OOB & Y.surv$time.event>=IBS.min)])

    xerror <- rep(NA,length(OOB_IBS))

    for (i in OOB_IBS){

      id_wY <- which(Y$id==i)
      if (is.element("Longitudinal",Inputs)==TRUE) {

        if (IBS.min==0){
          id_wXLongitudinal <- which(Longitudinal$id==i) # all measurements
        }else{
          id_wXLongitudinal <- which(Longitudinal$id==i&Longitudinal$time<=IBS.min) # only measurements until IBS.min
        }

        Longitudinal_current <- list(type="Longitudinal",X=Longitudinal$X[id_wXLongitudinal,,drop=FALSE],
                                     id=Longitudinal$id[id_wXLongitudinal],time=Longitudinal$time[id_wXLongitudinal],
                              model=Longitudinal$model)
      }

      if (is.element("Factor",Inputs)==TRUE){
        id_wXFactor <- which(Factor$id==i)
        Factor_current <- list(type="Factor",X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
      }

      if (is.element("Numeric",Inputs)==TRUE){
        id_wXNumeric <- which(Numeric$id==i)
        Numeric_current <- list(type="Numeric",X=Numeric$X[id_wXNumeric,,drop=FALSE], id=Numeric$id[id_wXNumeric])
      }

      pred_current <- tryCatch(pred.MMT(tree, Longitudinal = Longitudinal_current,
                                        Numeric = Numeric_current, Factor = Factor_current,
                                        timeVar = timeVar),
                               error = function(e) return(NA)) # handle permutation issue

      pred_current_chr <- as.character(pred_current)

      if (is.na(pred_current_chr)){

        xerror[which(OOB_IBS==i)] <- NA

      }else{

        # IPCW
        Wi_event <- (ifelse(Y$Y[id_wY,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[id_wY,2]!=0,1,0))/(G$IPCW.subjectTimes[which(Y.surv$id==i)])
        Wi_censored <- ifelse(Y$Y[id_wY,1] > allTimes_IBS, 1, 0)/(G$IPCW.times)
        Wi <- Wi_event + Wi_censored

        # CIF
        if (IBS.min == 0){
          pi_t <- tree$Y_pred[[pred_current_chr]][[as.character(cause)]]$traj[which(allTimes%in%allTimes_IBS)]
        }else{
          pi_t <- tree$Y_pred[[pred_current_chr]][[as.character(cause)]]$traj[which(allTimes%in%allTimes_IBS)] # pi(t)
          pi_s <- tree$Y_pred[[pred_current_chr]][[as.character(cause)]]$traj[sum(allTimes<IBS.min)] # pi(s)
          s_s <- 1 - sum(unlist(lapply(tree$Y_pred[[pred_current_chr]], FUN = function(x){
            return(x$traj[sum(allTimes<IBS.min)])
          }))) # s(s)
          pi_t <- (pi_t - pi_s)/s_s # P(S<T<S+t|T>S)
        }

        # Individual Brier Score
        Di <- ifelse(Y$Y[id_wY,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[id_wY,2]==cause,1,0) # D(t) = 1(s<Ti<s+t, event = cause)
        pec.res <- list()
        pec.res$AppErr$matrix <- Wi*(Di-pi_t)^2 # BS(t)
        pec.res$models <- "matrix"
        pec.res$time <- allTimes_IBS
        class(pec.res) <- "pec"

        #xerror[which(OOB==i)] <- pec::ibs(pec.res, start = IBS.min, times = IBS.max)[1] # IBS
        xerror[which(OOB_IBS==i)] <- pec::ibs(pec.res, start = IBS.min, times = max(allTimes_IBS))[1] # IBS

      }

    }

    ##############

  }

  if (Y$type == "factor" || Y$type == "numeric"){
    w_XLongitudinal <- NULL
    w_XNumeric <- NULL
    w_XFactor <- NULL
    w_y <- NULL

    if (is.element("Longitudinal",Inputs)==TRUE) w_XLongitudinal <- which(Longitudinal$id%in%OOB)
    if (is.element("Numeric",Inputs)==TRUE) w_XNumeric <- which(Numeric$id%in%OOB)
    if (is.element("Factor",Inputs)==TRUE) w_XFactor <- which(Factor$id%in%OOB)

    w_y <- which(Y$id%in%OOB)

    if (is.element("Longitudinal",Inputs)==TRUE) Longitudinal_current <- list(type="Longitudinal",X=Longitudinal$X[w_XLongitudinal,,drop=FALSE],
                                                                              id=Longitudinal$id[w_XLongitudinal],time=Longitudinal$time[w_XLongitudinal],
                                                                model=Longitudinal$model)
    if (is.element("Numeric",Inputs)==TRUE) Numeric_current  <- list(type="Numeric", X=Numeric$X[w_XNumeric,, drop=FALSE], id=Numeric$id[w_XNumeric])
    if (is.element("Factor",Inputs)==TRUE) Factor_current  <- list(type="Factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])

    pred_current <- pred.MMT(tree, Longitudinal = Longitudinal_current, Numeric = Numeric_current,
                     Factor = Factor_current, timeVar = timeVar)

    pred_current_chr <- as.character(pred_current)

    pred <- unlist(sapply(pred_current_chr, FUN = function(x) {
      ifelse(!is.null(tree$Y_pred[[x]]), tree$Y_pred[[x]], NA)
    }))

    if (Y$type=="numeric"){xerror <- (Y$Y[w_y]-pred)^2}
    if (Y$type=="factor"){xerror <- 1*(pred!=Y$Y[w_y])}

  }

  return(mean(xerror, na.rm = T))
}
