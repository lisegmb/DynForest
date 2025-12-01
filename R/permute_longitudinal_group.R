#' Permute a group of longitudinal markers between two subjects
#'
#' This function swaps the trajectories of a set of longitudinal markers
#' from subject B into subject A while keeping the other markers for subject A unchanged.
#'
#' @param Longitudinal A list containing longitudinal data (as used in dynforest).
#' @param id_var Name of the subject identifier variable (e.g., "id").
#' @param time_var Name of the time variable (e.g., "time").
#' @param marker_indices Indices of the markers to permute.
#' @param idA The ID of the subject whose trajectories will be replaced.
#' @param idB The ID of the subject whose trajectories will be used for replacement.
#' @param seed Optional seed for reproducibility.
#'
#' @return A data.frame containing the new trajectory for subject A.
#' @export
permute_longitudinal_group <- function(Longitudinal, id_var, time_var,
                                       marker_indices, idA, idB) {

  marker_names <- colnames(Longitudinal$X)
  selected_markers <- marker_names[marker_indices]
  other_markers <- setdiff(marker_names, selected_markers)

  # Extract data for subjects A and B
  df_A <- Longitudinal$X[Longitudinal$id == idA, , drop = FALSE]
  t_A  <- Longitudinal$time[Longitudinal$id == idA]

  df_B <- Longitudinal$X[Longitudinal$id == idB, , drop = FALSE]
  t_B  <- Longitudinal$time[Longitudinal$id == idB]

  # Combine all unique time points
  new_times <- sort(unique(c(t_A, t_B)))
  n_new <- length(new_times)

  # Create final data frame for subject A
  df_final <- data.frame(id = rep(idA, n_new), time = new_times)

  # Fill the selected markers with values from subject B
  for (mk in selected_markers) {
    x_B <- rep(NA_real_, n_new)
    idx_B <- which(new_times %in% t_B)
    x_B[idx_B] <- df_B[[mk]][match(new_times[idx_B], t_B)]
    df_final[[mk]] <- x_B
  }

  # Fill the remaining markers with original values from subject A
  for (mk in other_markers) {
    x_A <- rep(NA_real_, n_new)
    idx_A <- which(new_times %in% t_A)
    x_A[idx_A] <- df_A[[mk]][match(new_times[idx_A], t_A)]
    df_final[[mk]] <- x_A
  }

  # Return the final data frame sorted by time
  df_final <- df_final[order(df_final$time), ]
  return(df_final)
}
