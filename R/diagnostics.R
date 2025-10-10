
# Diagnostics builders & plotting (ggplot2) ------------------------------

.build_ess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  if (requireNamespace("coda", quietly = TRUE)) {
    return(as.numeric(coda::effectiveSize(x)))
  }
  ac <- stats::acf(x, plot = FALSE, lag.max = min(100, length(x)-1))$acf[-1]
  pos <- which(ac > 0)
  tau <- 1 + 2*sum(ac[pos])
  n <- length(x); n / tau
}

.collect_cluster_series <- function(clist, cluster_label) {
  if (length(clist) == 0) return(NULL)
  rows <- list()
  keep <- 0L
  for (s in seq_along(clist)) {
    cs <- clist[[s]]
    if (is.null(cs)) next
    kp <- cs$kern_pars
    rows[[length(rows)+1]] <- data.frame(save=s, param="tau_B", value=cs$tau_B, cluster=cluster_label, stringsAsFactors = FALSE)
    rows[[length(rows)+1]] <- data.frame(save=s, param="L11", value=cs$Lsel["L11"], cluster=cluster_label, stringsAsFactors = FALSE)
    if (!is.na(cs$Lsel["L21"])) rows[[length(rows)+1]] <- data.frame(save=s, param="L21", value=cs$Lsel["L21"], cluster=cluster_label, stringsAsFactors = FALSE)
    if (!is.na(cs$Lsel["L31"])) rows[[length(rows)+1]] <- data.frame(save=s, param="L31", value=cs$Lsel["L31"], cluster=cluster_label, stringsAsFactors = FALSE)
    if (length(kp)) for (nm in names(kp)) {
      rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("kern:",nm), value=kp[[nm]], cluster=cluster_label, stringsAsFactors = FALSE)
    }
    for (m in seq_along(cs$eta)) {
      rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("eta[ch",m,"]"), value=cs$eta[m], cluster=cluster_label, stringsAsFactors = FALSE)
    }
    for (m in seq_along(cs$sigma2)) {
      rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("sigma2[ch",m,"]"), value=cs$sigma2[m], cluster=cluster_label, stringsAsFactors = FALSE)
    }
    rows[[length(rows)+1]] <- data.frame(save=s, param="tau_sigma", value=cs$tau_sigma, cluster=cluster_label, stringsAsFactors = FALSE)
    if (!is.null(cs$pi_g) && nrow(cs$pi_g) > 0) {
      for (lev in rownames(cs$pi_g)) {
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("pi_",lev), value=cs$pi_g[lev,"pi"], cluster=cluster_label, stringsAsFactors = FALSE)
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("g_",lev),  value=cs$pi_g[lev,"g"],  cluster=cluster_label, stringsAsFactors = FALSE)
      }
    }
    keep <- keep + 1L
  }
  if (!keep || !length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Build diagnostics list from a wicmad() result
#' @export
build_diagnostics <- function(res, Y=NULL, t=NULL) {
  if (is.null(res$diagnostics)) stop("This run didn't store diagnostics. Re-run with diagnostics=TRUE.")
  diag <- res$diagnostics
  diag$global$ess <- list(
    K_occ = .build_ess(diag$global$K_occ),
    alpha = .build_ess(diag$global$alpha),
    loglik = .build_ess(diag$global$loglik)
  )
  dfA <- .collect_cluster_series(diag$clusters$A, "A")
  dfB <- .collect_cluster_series(diag$clusters$B, "B")
  dfAB <- rbind(dfA, dfB)
  if (!is.null(dfAB) && nrow(dfAB)) {
    splits <- split(dfAB, list(dfAB$cluster, dfAB$param), drop = TRUE)
    ess_rows <- vector("list", length(splits))
    keep <- 0L
    for (i in seq_along(splits)) {
      sub <- splits[[i]]
      if (is.null(sub) || !nrow(sub)) next
      ord <- order(sub$save)
      x <- sub$value[ord]
      x <- x[is.finite(x)]
      if (!length(x)) {
        ess_val <- NA_real_
      } else {
        ess_val <- .build_ess(x)
      }
      keep <- keep + 1L
      ess_rows[[i]] <- data.frame(cluster=sub$cluster[1], param=sub$param[1], ess=ess_val)
    }
    if (keep) {
      diag$clusters$ess <- do.call(rbind, ess_rows[!vapply(ess_rows, is.null, logical(1))])
    }
  }
  diag
}

#' Plot diagnostics (returns a named list of ggplot objects)
#' @export
plot_diagnostics <- function(diag) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install.packages('ggplot2').")
  }
  library(ggplot2)
  out <- list()
  df_g <- data.frame(iter=seq_along(diag$global$alpha),
                     K=diag$global$K_occ, alpha=diag$global$alpha,
                     loglik=diag$global$loglik)
  out$trace_K <- ggplot(df_g, aes(iter, K)) + geom_line() + labs(title="Trace: K_occ", x="save index", y="K")
  if (!is.null(diag$global$K_occ_all) && length(diag$global$K_occ_all)) {
    df_g_all <- data.frame(iter = seq_along(diag$global$K_occ_all),
                           K = diag$global$K_occ_all)
    out$trace_K_all <- ggplot(df_g_all, aes(iter, K)) +
      geom_line() + labs(title = "Trace: K_occ (full chain)",
                         x = "iteration", y = "K")
  }
  out$trace_alpha <- ggplot(df_g, aes(iter, alpha)) + geom_line() + labs(title="Trace: alpha", x="save index", y="alpha")
  out$trace_ll <- ggplot(df_g, aes(iter, loglik)) + geom_line() + labs(title="Trace: log-likelihood", x="save index", y="loglik")

  if (length(diag$ari) >= 1) {
    df_ari <- data.frame(iter=seq_along(diag$ari)+1, ARI=diag$ari)
    out$ari <- ggplot(df_ari, aes(iter, ARI)) + geom_line() + labs(title="ARI between successive saved draws", x="save index", y="ARI")
  }

  if (!is.null(diag$ari_all) && any(is.finite(diag$ari_all))) {
    df_ari_all <- data.frame(iter=seq_along(diag$ari_all), ARI=diag$ari_all)
    out$ari_all <- ggplot(df_ari_all, aes(iter, ARI)) +
      geom_line() + labs(title="ARI between successive iterations (full chain)", x="iteration", y="ARI")
  }

  dfA <- .collect_cluster_series(diag$clusters$A, "A")
  dfB <- .collect_cluster_series(diag$clusters$B, "B")
  dfAB <- rbind(dfA, dfB)
  if (!is.null(dfAB) && nrow(dfAB)) {
    out$cluster_traces <- ggplot(dfAB, aes(save, value, color=cluster)) + geom_line() +
      facet_wrap(~param, scales="free_y") + labs(title="Cluster parameter traces (A: largest, B: second)",
                                                 x="save index", y="value", color="cluster")

    build_cluster_acf <- function(df) {
      if (is.null(df) || !nrow(df)) return(NULL)
      splits <- split(df, list(df$param, df$cluster), drop = TRUE)
      rows <- vector("list", length(splits))
      keep <- 0L
      for (i in seq_along(splits)) {
        sub <- splits[[i]]
        ord <- order(sub$save)
        x <- sub$value[ord]
        ok <- is.finite(x)
        if (sum(ok) <= 1) next
        x <- x[ok]
        max_lag <- min(30L, length(x) - 1L)
        if (max_lag < 1) next
        ac <- stats::acf(x, plot = FALSE, lag.max = max_lag)$acf[-1]
        if (all(!is.finite(ac))) next
        rows[[i]] <- data.frame(
          lag = seq_len(length(ac)),
          acf = ac,
          param = sub$param[1],
          cluster = sub$cluster[1]
        )
        keep <- keep + 1L
      }
      if (!keep) return(NULL)
      do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    }

    df_acf <- build_cluster_acf(dfAB)
    if (!is.null(df_acf) && nrow(df_acf)) {
      out$cluster_acfs <- ggplot(df_acf) +
        geom_hline(yintercept = 0, color = "grey60") +
        geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf, color = cluster)) +
        facet_wrap(~param, scales = "free_y") +
        labs(title = "Cluster parameter autocorrelations (A: largest, B: second)",
             x = "lag", y = "ACF", color = "cluster")
    }
  }

  if (length(diag$clusters$A)) {
    lastA <- diag$clusters$A[[length(diag$clusters$A)]]
    if (!is.null(lastA$sparsity)) {
      suppressWarnings({
        dfs <- as.data.frame(as.table(lastA$sparsity))
        names(dfs) <- c("level","channel","p_gamma")
        out$sparsity_A <- ggplot(dfs, aes(level, p_gamma, group=channel, color=channel)) +
          geom_line() + geom_point() + labs(title="Sparsity profile by level (cluster A, last save)",
                                            x="wavelet level", y="Pr(γ=1)")
      })
    }
  }
  out
}

#' Save plots to PNG files
#' @export
save_plots_png <- function(plots, out_dir = "diagnostics",
                           dpi = 150, width = 9, height = 6, units = "in") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list()
  nms <- names(plots)
  if (is.null(nms)) nms <- sprintf("plot_%02d", seq_along(plots))
  for (i in seq_along(plots)) {
    nm <- nms[i]
    p  <- plots[[i]]
    fp <- file.path(out_dir, paste0(nm, ".png"))
    try(ggplot2::ggsave(filename = fp, plot = p, dpi = dpi,
                        width = width, height = height, units = units), silent = TRUE)
    paths[[nm]] <- fp
  }
  invisible(paths)
}

#' One-liner: run model then save diagnostics to disk
#' @export
save_diagnostics_from_res <- function(res, Y=NULL, t=NULL, out_dir = "diagnostics",
                                      dpi = 150, width = 9, height = 6, units = "in") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  diag  <- build_diagnostics(res, Y = Y, t = t)
  plots <- plot_diagnostics(diag)
  paths <- save_plots_png(plots, out_dir = out_dir, dpi = dpi,
                          width = width, height = height, units = units)
  ess_rows <- list()
  if (!is.null(diag$global$ess) && length(diag$global$ess)) {
    eg <- diag$global$ess
    ess_rows[[length(ess_rows)+1]] <- data.frame(
      category = "global",
      cluster = NA_character_,
      parameter = names(eg),
      ess = as.numeric(eg),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(diag$clusters$ess) && nrow(diag$clusters$ess) > 0) {
    ec <- diag$clusters$ess
    ess_rows[[length(ess_rows)+1]] <- data.frame(
      category = "cluster",
      cluster = ec$cluster,
      parameter = ec$param,
      ess = ec$ess,
      stringsAsFactors = FALSE
    )
  }
  if (length(ess_rows)) {
    ess_df <- do.call(rbind, ess_rows)
    ess_path <- file.path(out_dir, "effective_sample_sizes.txt")
    utils::write.table(ess_df, file = ess_path, sep = "\t",
                       row.names = FALSE, quote = FALSE, na = "NA")
    paths$ess <- ess_path
  }
  invisible(list(diag = diag, paths = paths))
}
