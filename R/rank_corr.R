rank_corr <- function(data, variable = "Drug", value = "AUC", type = "pearson",
                      min.per.corr = 3, plots = TRUE, FDR = 0.05,
                      xlab = colnames(data)[2], ylab = value,
                      position.x = "mid", position.y = "max", se = TRUE) {
  message("Running correlations and regressions...")

  # check class of input is as expected
  testthat::expect_is(data, "data.frame")

  # use parallel computing if possible
  if (requireNamespace("parallel") &
    requireNamespace("snow") &
    requireNamespace("doSNOW")) {
    cores <- parallel::detectCores() # number of cores available
    if (cores[1] > 1) {
      cl <- snow::makeCluster(cores[1] - 1) # cluster using all but 1 core
      doSNOW::registerDoSNOW(cl) # register cluster
    }
  }

  # run correlations and regressions
  not_all_na <- function(x) any(!is.na(x))
  all.data.corr <- data
  all.data.corr <- all.data.corr %>% select_if(not_all_na)
  variable.set <- colnames(all.data.corr[, 3:ncol(all.data.corr)])
  corr <- data.frame(variable.set)
  colnames(corr)[1] <- variable
  corr[, c(
    "Pearson.est", "Pearson.p", "Pearson.q",
    "Spearman.est", "Spearman.p", "Spearman.q",
    "Slope", "Intercept", "R.squared",
    "Rank.slope", "Rank.intercept", "Rank.R.squared", "N"
  )] <- NA
  data.corr <- list()
  for (j in 3:ncol(all.data.corr)) {
    data.corr[[j - 2]] <- na.omit(all.data.corr[, c(2, j)])
    a <- nrow(data.corr[[j - 2]]) # number of samples in the correlation (N)
    if (a >= min.per.corr) {
      corr$N[j - 2] <- a

      # list of rank metrics
      x <- as.numeric(data.corr[[j - 2]][, 1])

      # list of gene expression values for each gene j
      y <- as.numeric(data.corr[[j - 2]][, 2])

      # run linear regression (value-based)
      Regression <- stats::lm(y ~ x)
      corr$Slope[j - 2] <- Regression$coeff[[2]]
      corr$Intercept[j - 2] <- Regression$coeff[[1]]
      corr$R.squared[j - 2] <- summary(Regression)$r.squared

      # run linear regression (rank-based)
      Rank.regression <- stats::lm(rank(y) ~ rank(x))
      corr$Rank.slope[j - 2] <- Rank.regression$coeff[[2]]
      corr$Rank.intercept[j - 2] <- Rank.regression$coeff[[1]]
      corr$Rank.R.squared[j - 2] <- summary(Rank.regression)$r.squared

      # run Pearson correlation
      Pearson <- stats::cor.test(x, y, method = "pearson")
      corr$Pearson.est[j - 2] <- Pearson$estimate
      corr$Pearson.p[j - 2] <- Pearson$p.value

      # run Spearman correlation
      Spearman <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
      corr$Spearman.est[j - 2] <- Spearman$estimate
      corr$Spearman.p[j - 2] <- Spearman$p.value
    }
  }
  corr.no.na <- corr[!is.na(corr$N), ]

  # perform multiple hypothesis testing (Benjamini-Hochberg)
  corr.no.na$Pearson.q <-
    qvalue::qvalue(p = corr.no.na$Pearson.p, pi0 = 1)$qvalues
  corr.no.na$Spearman.q <-
    qvalue::qvalue(p = corr.no.na$Spearman.p, pi0 = 1)$qvalues

  # check class of output is as expected
  testthat::expect_is(corr.no.na, "data.frame")

  # scatter plots for significant results
  if (plots) {
    # load themes for plots
    ng.theme <- ggplot2::theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 8, colour = "black")
    )

    bg.theme <- ggplot2::theme(
      legend.background = element_rect(), legend.position = "top",
      legend.text = element_text(size = 14), legend.key = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(lineheight = .8, face = "bold", size = 36)
    )

    # prepare data frame for plots
    id.var <- colnames(data)[1]
    rank.var <- colnames(data)[2]
    plot.data <-
      reshape2::melt(all.data.corr, id = c(id.var, rank.var),
                     variable.name = variable, value.name = value, na.rm = TRUE)
    if (type == "pearson") {
      results <- unique(corr.no.na[which(corr.no.na$Pearson.q <= FDR), 1])
      if (length(results) > 0) {
        a <- list()
        for (i in seq_len(length(results))) {
          # identify significant correlations
          sig.data <- plot.data[plot.data[, 3] == results[i], ]

          # set plot parameters
          min.x <- min(sig.data[, c(rank.var)])
          max.x <- max(sig.data[, c(rank.var)])
          mid.x <- 0.5 * (min.x + max.x)
          min.y <- min(sig.data[, c(value)])
          max.y <- max(sig.data[, c(value)])
          mid.y <- 0.5 * (min.y + max.y)
          if (position.x == "min") {
            pos.x <- min.x
          } else if (position.x == "mid") {
            pos.x <- mid.x
          } else if (position.x == "max") {
            pos.x <- max.x
          } else if (is.numeric(position.x)) {
            pos.x <- position.x
          }
          if (position.y == "min") {
            pos.y <- min.y
          } else if (position.y == "mid") {
            pos.y <- mid.y
          } else if (position.y == "max") {
            pos.y <- max.y
          } else if (is.numeric(position.y)) {
            pos.y <- position.y
          }

          # get correlation parameters for plot
          Pearson.est <- corr.no.na[corr.no.na[, 1] == results[i], ]$Pearson.est
          Pearson.p <- corr.no.na[corr.no.na[, 1] == results[i], ]$Pearson.p
          stats_pearson <- substitute(
            r == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = format(Pearson.est, digits = 3),
              p = format(Pearson.p, digits = 3)
            )
          )

          # generate plot for each significant correlation
          a[[i]] <- ggplot2::ggplot(data = sig.data,
                                    aes_string(x = rank.var, y = value)) +
            ggplot2::geom_point() +
            ggplot2::labs(x = xlab, y = ylab) +
            ggplot2::ggtitle(results[i]) +
            ggplot2::geom_smooth(method = "lm", size = 1.5,
                                 linetype = "solid", color = "blue",
                                 se = se, na.rm = TRUE) +
            ggplot2::geom_text(
              x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
              colour = "blue", parse = TRUE,
              label = as.character(as.expression(stats_pearson)), size = 8
            ) +
            ng.theme +
            bg.theme
        }
        scatter.plots <- gridExtra::marrangeGrob(a, nrow = 1, ncol = 1)
      } else {
        scatter.plots <- NA
        warning("No correlations met the FDR cut-off to produce scatter plots")
      }
    } else if (type == "spearman") {
      results <- unique(corr.no.na[which(corr.no.na$Spearman.q <= FDR), 1])
      if (length(results) > 0) {
        a <- list()
        for (i in seq_len(length(results))) {
          # identify significant correlations
          sig.data <- plot.data[plot.data[, 3] == results[i], ]

          # set plot parameters
          min.x <- 1
          mid.x <- 0.5 * length(sig.data[, c(rank.var)])
          max.x <- length(sig.data[, c(rank.var)])
          min.y <- 1
          mid.y <- 0.5 * length(sig.data[, c(value)])
          max.y <- length(sig.data[, c(value)])
          if (position.x == "min") {
            pos.x <- min.x
          } else if (position.x == "mid") {
            pos.x <- mid.x
          } else if (position.x == "max") {
            pos.x <- max.x
          } else if (is.numeric(position.x)) {
            pos.x <- position.x
          }
          if (position.y == "min") {
            pos.y <- min.y
          } else if (position.y == "mid") {
            pos.y <- mid.y
          } else if (position.y == "max") {
            pos.y <- max.y
          } else if (is.numeric(position.y)) {
            pos.y <- position.y
          }

          # get correlation parameters for plot
          Spearman.est <-
            corr.no.na[corr.no.na[, 1] == results[i], ]$Spearman.est
          Spearman.p <- corr.no.na[corr.no.na[, 1] == results[i], ]$Spearman.p
          stats_spearman <- substitute(
            rho == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = format(Spearman.est, digits = 3),
              p = format(Spearman.p, digits = 3)
            )
          )

          # generate plot for each significant correlation
          a[[i]] <- ggplot2::ggplot(data = sig.data, aes_string(
            x = rank(sig.data[, c(rank.var)]),
            y = rank(sig.data[, c(value)])
          )) +
            ggplot2::geom_point() +
            ggplot2::labs(x = paste0(xlab, " Rank"), y = paste0(ylab, " Rank")) +
            ggplot2::ggtitle(results[i]) +
            ggplot2::geom_smooth(
              method = "lm", size = 1.5, linetype = "solid", color = "blue",
              se = se, na.rm = TRUE
            ) +
            ggplot2::geom_text(
              x = pos.x, y = pos.y,
              label = as.character(as.expression(stats_spearman)),
              colour = "blue", size = 8, parse = TRUE
            ) +
            ng.theme +
            bg.theme
        }
        scatter.plots <- gridExtra::marrangeGrob(a, nrow = 1, ncol = 1)
      } else {
        warning("No correlations met the FDR cut-off to produce scatter plots")
        scatter.plots <- NULL
      }
    } else {
      warning(paste("Type must be specified as either spearman or",
                    "pearson to produce scatter plots"))
      scatter.plots <- NULL
    }
  } else {
    scatter.plots <- NULL
  }

  if (requireNamespace("parallel") &
    requireNamespace("snow") &
    requireNamespace("doSNOW")) {
    if (cores[1] > 1) {
      snow::stopCluster(cl) # stop cluster
      rm(cl)
    }
  }

  return(list(result = corr.no.na,
              scatter.plots = scatter.plots))
}
