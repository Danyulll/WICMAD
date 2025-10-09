
# Diagnostics builders & plotting (ggplot2) ------------------------------

.build_ess <- function(x) {
  if (length(x) < 10 || all(!is.finite(x))) return(NA_real_)
  if (requireNamespace("coda", quietly = TRUE)) {
    return(as.numeric(coda::effectiveSize(x)))
  }
  ac <- stats::acf(x, plot = FALSE, lag.max = min(100, length(x)-1))$acf[-1]
  pos <- which(ac > 0)
  tau <- 1 + 2*sum(ac[pos])
  n <- length(x); n / tau
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
  out$trace_alpha <- ggplot(df_g, aes(iter, alpha)) + geom_line() + labs(title="Trace: alpha", x="save index", y="alpha")
  out$trace_ll <- ggplot(df_g, aes(iter, loglik)) + geom_line() + labs(title="Trace: log-likelihood", x="save index", y="loglik")

  if (length(diag$ari) >= 1) {
    df_ari <- data.frame(iter=seq_along(diag$ari)+1, ARI=diag$ari)
    out$ari <- ggplot(df_ari, aes(iter, ARI)) + geom_line() + labs(title="ARI between successive saved draws", x="save index", y="ARI")
  }

  build_cluster_df <- function(clist, lab) {
    if (length(clist) == 0) return(NULL)
    S <- length(clist)
    rows <- list()
    for (s in 1:S) {
      cs <- clist[[s]]
      if (is.null(cs)) next
      kp <- cs$kern_pars
      rows[[length(rows)+1]] <- data.frame(save=s, param="tau_B", value=cs$tau_B, cluster=lab)
      rows[[length(rows)+1]] <- data.frame(save=s, param="L11", value=cs$Lsel["L11"], cluster=lab)
      if (!is.na(cs$Lsel["L21"])) rows[[length(rows)+1]] <- data.frame(save=s, param="L21", value=cs$Lsel["L21"], cluster=lab)
      if (!is.na(cs$Lsel["L31"])) rows[[length(rows)+1]] <- data.frame(save=s, param="L31", value=cs$Lsel["L31"], cluster=lab)
      if (length(kp)) for (nm in names(kp)) {
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("kern:",nm), value=kp[[nm]], cluster=lab)
      }
      for (m in seq_along(cs$eta)) {
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("eta[ch",m,"]"), value=cs$eta[m], cluster=lab)
      }
      for (m in seq_along(cs$sigma2)) {
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("sigma2[ch",m,"]"), value=cs$sigma2[m], cluster=lab)
      }
      rows[[length(rows)+1]] <- data.frame(save=s, param="tau_sigma", value=cs$tau_sigma, cluster=lab)
      for (lev in rownames(cs$pi_g)) {
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("pi_",lev), value=cs$pi_g[lev,"pi"], cluster=lab)
        rows[[length(rows)+1]] <- data.frame(save=s, param=paste0("g_",lev),  value=cs$pi_g[lev,"g"],  cluster=lab)
      }
    }
    do.call(rbind, rows)
  }
  dfA <- build_cluster_df(diag$clusters$A, "A")
  dfB <- build_cluster_df(diag$clusters$B, "B")
  dfAB <- rbind(dfA, dfB)
  if (!is.null(dfAB)) {
    out$cluster_traces <- ggplot(dfAB, aes(save, value, color=cluster)) + geom_line() +
      facet_wrap(~param, scales="free_y") + labs(title="Cluster parameter traces (A: largest, B: second)",
                                                 x="save index", y="value", color="cluster")
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
  diag  <- build_diagnostics(res, Y = Y, t = t)
  plots <- plot_diagnostics(diag)
  paths <- save_plots_png(plots, out_dir = out_dir, dpi = dpi,
                          width = width, height = height, units = units)
  invisible(list(diag = diag, paths = paths))
}
