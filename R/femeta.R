#' Generic Fixed-effects Meta-analysis with Several 2x2 Tables
#'
#' Generic function for obtaining the confidence interval of logarithm of Cochran-Mantel-Haenszel estimator,
#' Woolf's estimator of MLE of an assumed-common log odds ratio without assuming
#' homogeneity
#'
#' @param x an object used to select a method. Currently x can be either a 2x2xK array, which links to \link{femeta.array} method,
#' or a integer vector, which leads to the use of \link{femeta.default}.
#' @param ... further arguments passed to or from other methods.
#' @export
femeta <- function(x, ...) {
  UseMethod("femeta", x)
}


#' Array Method of Improved Inference for Fixed-effects Meta-analysis
#'
#' @param x vector to specify the 2x2 table frequencies (upper left cell).
#' @param correction_factor A scalar used for zero-cell correction. Intuitively it is the total number of
#' subjects added to each group of a study. A small correction factor induces less bias. The default is 0.01.
#' See 'Details'.
#' @param drop00 logical indicating whether studies with no cases/events (or only cases) in both groups
#' should be dropped when calculating the observed outcomes of the individual studies. See 'Details'.
#' @param estimator a character string indicating which estimator is used to summarize to overall association (either "log-CMH",
#' "Woolf" or "CEMLE"). See 'Details'.
#' @param vtype a character string indicating whether to assume homogeneity when calculating the variance (either "Het" or "Hom).
#' @param alternative a character string indicating the alternative hypothesis and must be one of "two.sided", "greater"
#' or "less". You can specify just the initial letter.
#' @param conf.level a scalar indicating the confidence level for the returned confidence interval
#' @param null.value scalar indicating the value of the parameter of interest (some summary of the true effect sizes)
#' @param ... other arguments passed to \link{femeta.default}.
#' @importFrom stats glm pchisq pnorm qnorm
#'
#' @examples
#' X <- array(sample(1:20, size = 12, replace = TRUE), dim = c(2, 2, 3))
#' femeta(X)

#' @export
femeta.array <- function(x, correction_factor = 0.01,
                         drop00 = FALSE, estimator = "log-CMH", vtype = "Het",
                         alternative = "two.sided", conf.level = 0.95,
                         null.value = 0, ...) {
  ### check if x is a 3-dimensional array
  if (is.array(x)) {
    if (length(dim(x)) == 3L) {
      if (anyNA(x))
        stop("NAs are not allowed")
      if (any(dim(x) < 2L))
        stop("each dimension in table must be >= 2")
    }
    else stop("'x' must be a 3-dimensional array")
  }

  ### check if each stratum is non-empty
  if (any(apply(x, 3L, sum) < 2))
    stop("sample size in each stratum must be > 1")

  ### check if each table is 2x2
  if (dim(x)[1L] != 2 | dim(x)[2L] != 2)
    stop("each subtable must be 2x2.")

  ### use the default function
  femeta.default(ai = x[1, 1, ], bi = x[1, 2, ], ci = x[2, 1, ], di = x[2, 2, ],
                 correction_factor = 0.01, drop00 = FALSE, estimator = "log-CMH", vtype = "Het",
                 alternative = "two.sided", conf.level = 0.95,
                 null.value = 0, ...)
}


#' Default Method of Improved Inference for Fixed-effects meta-analysis
#'
#' @param ai vector to specify the 2x2 table frequencies (upper left cell).
#' @param bi vector to specify the 2x2 table frequencies (upper right cell).
#' @param ci vector to specify the 2x2 table frequencies (lower left cell).
#' @param di vector to specify the 2x2 table frequencies (lower right cell).
#' @param n1i vector to specify the group sizes or column totals (first group/column).
#' @param n2i vector to specify the groups sizes or column totals (second group/column).
#' @param x1i vector to specify the number of events (first group).
#' @param x2i vector to specify the number of evetns (second group).
#' @param data optional data frame containing the variables given to the arguments above.
#' @param slab optional vector with labels for the studies.
#' @param subset optional vector indicating the subset of studies that should be used.
#' This can be a logical vector or a numeric vector indicating the indices of the studies to include.
#' @param correction_factor A scalar used for zero-cell correction. Intuitively it is the total number of
#' subjects added to each group of a study. A small correction factor induces less bias. The default is 0.01.
#' See 'Details'.
#' @param drop00 logical indicating whether studies with no cases/events (or only cases) in both groups
#' should be dropped when calculating the observed outcomes of the individual studies. See 'Details'.
#' @param estimator a character string indicating which estimator is used to summarize to overall association (either "log-CMH",
#' "Woolf" or "CEMLE"). See 'Details'.
#' @param vtype a character string indicating whether to assume homogeneity when calculating the variance (either "Het" or "Hom).
#' @param alternative a character string indicating the alternative hypothesis and must be one of "two.sided", "greater"
#' or "less". You can specify just the initial letter.
#' @param conf.level a scalar indicating the confidence level for the returned confidence interval.of the
#' @param null.value scalar indicating the value of the parameter of interest (some summary of the true effect sizes)

#' @importFrom stats glm pchisq pnorm qnorm
#'
#' @export
femeta.default <- function(ai, bi, ci, di, n1i, n2i, x1i, x2i,
                   data, slab, subset, correction_factor = 0.01,
                   drop00 = FALSE, estimator = "log-CMH", vtype = "Het",
                   alternative = "two.sided", conf.level = 0.95,
                   null.value = 0) {

  ### check argument specifications



  if (!is.element(estimator, c("log-CMH",
                            "Woolf",
                            "CEMLE")))

  stop("Unknown 'estimator' argument specified. ")

  if (!is.element(vtype, c("Het", "Hom"))) ### vtype can be an entire vector, so use any() and na.rm=TRUE
    stop("Unknown 'vtype' argument specified.")

  alternative <- char.expand(alternative, c("two.sided",
                                            "less", "greater"))

  if (!is.element(alternative, c("two.sided", "less", "greater"))) {
    stop("Invalid 'alternative' argument.")
  }
  ### check if data argument has been specified

  if (missing(data))
    data <- NULL

  ### need this at the end to check if append=TRUE can actually be done

  no.data <- is.null(data)

  ### check if data argument has been specified

  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data))
      data <- data.frame(data)
  }

  mf <- match.call()

  ### get slab and subset arguments (will be NULL when unspecified)

  mf.slab   <- mf[[match("slab",   names(mf))]]
  mf.subset <- mf[[match("subset", names(mf))]]
  slab      <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
  subset    <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))





  #########################################################################
  #########################################################################
  #########################################################################



      mf.ai  <- mf[[match("ai",  names(mf))]]
      mf.bi  <- mf[[match("bi",  names(mf))]]
      mf.ci  <- mf[[match("ci",  names(mf))]]
      mf.di  <- mf[[match("di",  names(mf))]]
      mf.n1i <- mf[[match("n1i", names(mf))]]
      mf.n2i <- mf[[match("n2i", names(mf))]]
      ai     <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
      bi     <- eval(mf.bi,  data, enclos=sys.frame(sys.parent()))
      ci     <- eval(mf.ci,  data, enclos=sys.frame(sys.parent()))
      di     <- eval(mf.di,  data, enclos=sys.frame(sys.parent()))
      n1i    <- eval(mf.n1i, data, enclos=sys.frame(sys.parent()))
      n2i    <- eval(mf.n2i, data, enclos=sys.frame(sys.parent()))
      if (is.null(bi)) bi <- n1i - ai
      if (is.null(di)) di <- n2i - ci

      if (!is.null(subset)) {
        ai <- ai[subset]
        bi <- bi[subset]
        ci <- ci[subset]
        di <- di[subset]
        n1i <- n1i[subset]
        n2i <- n2i[subset]
      }

      if (length(ai)==0L || length(bi)==0L || length(ci)==0L || length(di)==0L)
        stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

      if (!all(length(ai) == c(length(ai),length(bi),length(ci),length(di))))
        stop("Supplied data vectors are not all of the same length.")

      if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
        stop("One or more event counts are larger than the corresponding group sizes.")

      if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
        stop("One or more counts are negative.")

      ni.u <- ai + bi + ci + di ### unadjusted total sample sizes

      k <- length(ai) ### number of studies before adjusting

      ### if drop00=TRUE, drop the studies with zero counts in both arms
      if (drop00) {
        id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
        ai <- ai[!id00]
        bi <- bi[!id00]
        ci <- ci[!id00]
        di <- di[!id00]
      }

      ### save unadjusted counts
      ai.u <- ai
      bi.u <- bi
      ci.u <- ci
      di.u <- di
      n1i.u <- ai + ci
      n2i.u <- bi + di
      ni.u <- n1i.u + n2i.u


      ### number of studies
      K <- length(ai)

      ### zero-cell correction: add the reciprocal of the sample size of the opposite arm to each cell in tables with zeroes
      any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

      for (jj in 1:K) {
        if (any0[jj]) {
          ai[jj] <- ai[jj] + correction_factor * n1i.u[jj] / ni.u[jj]
          bi[jj] <- bi[jj] + correction_factor * n2i.u[jj] / ni.u[jj]
          ci[jj] <- ci[jj] + correction_factor * n1i.u[jj] / ni.u[jj]
          di[jj] <- di[jj] + correction_factor * n2i.u[jj] / ni.u[jj]
        }
      }



      ### recompute group and total sample sizes (after add/to adjustment)

      n1i <- ai + ci
      n2i <- bi + di
      ni <- n1i + n2i ### ni.u computed earlier is always the 'unadjusted' total sample size
      nn <- sum(ni) ### total sample size

      ### compute proportions for the two groups (unadjusted and adjusted)
      ESTIMATOR <- switch(estimator,
                          "log-CMH" = "Logarithm of Cochrane-Mantel-Haenszel estimator",
                          "Woolf" = "Woolf's inverse-variance estimator",
                          "CEMLE" = "Maximum Likelihood Estimator of the assumed common log odds ratio")
      TYPE <- switch(vtype,
                     "Het" = "without assuming homogeneity",
                     "Hom" = "assuming homogeneity")
      if (estimator == "log-CMH") {
        ## logarithm of mantel-haenszel statistics
        theta_MH <- sum(ai * di / nn) / sum(bi * ci / nn)
        EST <- log(theta_MH)

        if (vtype == "Het" ) {
          VAR <- sum((ai * ci / n1i * (di + bi * theta_MH) ^ 2 + bi * di / n2i * (ai + ci * theta_MH) ^ 2) / nn ^ 2) /
            (sum(bi * ci / nn) * theta_MH) ^ 2
        } else if (vtype == "Hom") {
          VAR <- sum((bi * ci / nn) ^ 2 * (1 / ai + 1 / bi + 1 / ci + 1 / di)) /
            sum(bi * ci / nn) ^ 2
        }



      } else if (estimator == "Woolf") {
        ## subtable-specific log odds ratio
        ff <- log((ai * di) / (bi * ci))

        ## weights
        ww <- 1 / (1 / ai + 1 / bi + 1 / ci + 1 / di)

        ## Woolf's estimator
        EST <- sum(ww * ff) / sum(ww)
        if (vtype == "Hom") {

          ### Under homogeneity, compute the variance through inverse-variance method
          VAR <- 1 / sum(ww)
        } else if (vtype == "Het") {

          ### Under heterogeneity, compute the variance through methods stated in Li 2020
          ### weights
          ww <- 1 / (1 / ai + 1 / bi + 1 / ci + 1 / di)

          ### inflation factors
          Delta1 <- ((2 * ai / n1i - 1) * (bi * di / n2i) ^ 2 - (2 * bi / n2i - 1) * (ai * ci / n1i) ^ 2 ) /
            (ai * ci / n1i + bi * di / n2i) ^ 2

          Delta2 <- ((2 * ai / n1i - 1) ^ 2 * (bi * di / n2i) ^ 3 + (2 * bi / n2i - 1) ^ 2 * (ai * ci / n1i) ^ 3 ) /
            (ai * ci / n1i + bi * di / n2i) ^ 3


          VAR <- 1 / sum(ww) + (-2 * sum(ww * ff * Delta1) + sum(ww * ff ^ 2 * Delta2)) / sum(ww) ^ 2 +
            2 * (sum(ww * ff) * sum(ww * Delta1) - sum(ww * ff) * sum(ww * ff * Delta2)) / sum(ww) ^ 3 +
            sum(ww * ff) ^ 2 * sum(ww * Delta2) / sum(ww) ^ 4

        }
      } else if (estimator == "CEMLE") {

        ## formulate the data for running logistic regression
        dat_reshape <- rbind(data.frame(rr = ai, nr = ci, xx = 1, t = 1:K),
                             data.frame(rr = bi, nr = di, xx = 0, t = 1:K))
        suppressWarnings(
          logit_mod <- glm(cbind(rr, nr) ~ 0 + as.factor(t) + xx, data = dat_reshape, family = "binomial")
        )
        ## MLE <- logit_mod$coefficients["xx"]


        MLE <- logit_mod$coefficients

        ## compute the asymptotic variance-covariance matrix of beta_hat
        aa <- as.numeric(MLE[1:K])
        bb <- as.numeric(MLE[K + 1])

        ## estimator: MLE of the assumed common log odds ratio
        EST <- bb

        ## make sandwich estimator of variance
        ## make matrix J.hat as in the document (negtive Fisher Information)
        uu <- - n1i * exp(aa + bb) / ((1 + exp(aa + bb)) ^ 2 * nn)
        vv <- - n2i * exp(aa) / ((1 + exp(aa)) ^ 2 * nn)

        J <- matrix(0, nrow = K + 1, ncol = K + 1)
        diag(J) <- c(uu + vv, sum(uu))
        J[K + 1, 1:K] <- J[1:K, K + 1] <- uu

        if (vtype == "Hom") {
          OMEGA <- solve(-J) ### covariance matrix of the MLE
        } else if (vtype == "Het") {
          ## make matrix U.hat as in the paper -- the meat in the sandwich
          ss <- ai * ci / (n1i * nn)
          tt <- bi * di / (n2i * nn)

          U <- matrix(0, nrow = K + 1, ncol = K + 1)
          diag(U) <- c(ss + tt, sum(ss))
          U[K + 1, 1:K] <- U[1:K, K + 1] <- ss

          OMEGA <- solve(-J) %*% U %*% solve(-J)
        }

        ### approximate variance of CEMLE
        VAR <- OMEGA[K + 1, K + 1] / nn

      }

      ### approximate standard error of the estimator
      SE <- sqrt(VAR)

      ### confidence interval for the estimators
      alpha <- 1 - conf.level
      if (alternative == "less") {
        CI <- c(EST + qnorm(alpha) * SE, Inf)
        PVAL <- pnorm(q = EST, mean = null.value, sd = SE, lower.tail = T)
      } else if (alternative == "greater") {
        CI <- c(-Inf, EST + qnorm(conf.level) * SE)
        PVAL <- pnorm(q = EST, mean = null.value, sd = SE, lower.tail = F)
      } else if (alternative == "two.sided") {
        CI <- EST + qnorm(c(alpha / 2, 1 - alpha / 2)) * SE
        PVAL <- pchisq(q = (EST - null.value) ^ 2 / VAR, df = 1, lower.tail = F)
      }

      RESULTS <- structure(list(estimator = ESTIMATOR,
                                type = TYPE,
                                estimate = EST,
                                variance = VAR,
                                standard_error = SE,
                                alternative = alternative,
                                conf.int = CI,
                                null.value = null.value,
                                p.value = PVAL),
                           class = "femeta")
      return(RESULTS)
}



#' Printing results of fixed-effects meta-analysis
#'
#' @param x object of class "femeta"
#' @param digits number of significant digits to be used
#' @param prefix string, passed to strwrap for displaying the method component of the \code{femeta} object.
#'
#' @export
print.femeta <- function(x, digits = getOption("digits"), prefix = "\t") {
  cat("\n")
  cat(x$estimator, x$type, "\n", sep = " ")
  cat("\n")
  out <- character()
  out <- c(out, paste("Estimate =", format(x$estimate, digits = max(1L, digits - 2L))))

    if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (substr(fp, 1L,
                                              1L) == "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to",
                           less = "less than", greater = "greater than")
        cat("true overall log OR is ",
            alt.char, " ", x$null.value, "\n",
            sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n",
            sep = "")
        print(x$null.value, digits = digits)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")),
        " percent confidence interval:\n", " ",
        paste(format(x$conf.int[1:2], digits = digits), collapse = " "),
        "\n", sep = "")
  }

  invisible(x)

}


