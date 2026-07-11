#This code extends MorphoGAM to add a by.t into the GAM 

morphoGAM_byt <- function(Y,
                           curve.fit,
                           design,
                           by.t,
                           shrinkage = FALSE,
                           min.count.per.gene = 10,
                           return.fx = TRUE,
                           offset = NULL,
                           knots.t = NULL,
                           knots.r = NULL) {

  if (is.null(offset)) {
    l.o <- log(colSums(Y))
  } else {
    l.o <- offset
  }

  median.depth <- median(colSums(Y))
  n <- ncol(Y)
  p <- nrow(Y)

  if (length(by.t) != n) {
    stop("by.t must have length equal to ncol(Y).")
  }

  data <- data.frame(
    y = Y[1, ],
    t = curve.fit$xyt$t,
    r = curve.fit$xyt$r,
    by.t = by.t,
    l.o = l.o
  )

  fit.setup <- mgcv::gam(
    formula = design,
    family = mgcv::nb(),
    data = data,
    fit = FALSE
  )

  H <- diag(ncol(fit.setup$X))
  H[1, 1] <- 0

  results <- matrix(0, nrow = p, ncol = 7)
  colnames(results) <- c(
    "peak.t", "range.t", "pv.t",
    "peak.r", "range.r", "pv.r",
    "intercept"
  )

  if (return.fx) {
    fxs.t <- matrix(0, nrow = p, ncol = n)
    fxs.r <- matrix(0, nrow = p, ncol = n)
  }

  bar <- txtProgressBar(min = 0, max = nrow(Y) - 1, initial = 0)

  for (i in 1:nrow(Y)) {
    setTxtProgressBar(bar, i)

    if (rowSums(Y)[i] < min.count.per.gene) {
      results[i, ] <- c(0, 0, 1, 0, 0, 1, 0)
      next
    }

    data$y <- Y[i, ]

    knot.list <- list()
    if (!is.null(knots.t)) knot.list$t <- knots.t
    if (!is.null(knots.r)) knot.list$r <- knots.r

    fit <- mgcv::gam(
      formula = design,
      family = mgcv::nb(),
      offset = l.o,
      data = data,
      H = H,
      knots = if (length(knot.list) > 0) knot.list else NULL
    )

    beta_g0 <- fit$coefficients[1]
    results[i, "intercept"] <- beta_g0

    se_beta <- diag(vcov(fit))[-1]

    if (shrinkage) {
      my.ash <- ashr::ash(fit$coefficients[-1], se_beta)
      beta.shrink <- apply(
        ashr::get_post_sample(my.ash, 1000),
        2,
        median
      )
    } else {
      beta.shrink <- fit$coefficients[-1]
    }

    basis.functions <- mgcv::predict.gam(fit, type = "lpmatrix")[, -1, drop = FALSE]

    pred <- predict(
      fit,
      newdata = data.frame(
        t = data$t,
        r = data$r,
        by.t = data$by.t,
        l.o = data$l.o
      ),
      type = "terms",
      se.fit = TRUE
    )

    term.names <- colnames(pred$fit)

    t.term <- grep("s\\(t\\)", term.names, value = TRUE)
    r.term <- grep("s\\(r\\)", term.names, value = TRUE)

    if (length(t.term) > 0) {
      t.idx <- match(t.term[1], term.names)

      fx.t <- pred$fit[, t.idx] / (1 + pred$se.fit[, t.idx]^2)

      if (return.fx) {
        fxs.t[i, ] <- fx.t
      }

      fx.t.sig <- ifelse(
        pred$fit[, t.idx] > 1.96 * pred$se.fit[, t.idx],
        pred$fit[, t.idx] - 1.96 * pred$se.fit[, t.idx],
        ifelse(
          pred$fit[, t.idx] < -1.96 * pred$se.fit[, t.idx],
          pred$fit[, t.idx] + 1.96 * pred$se.fit[, t.idx],
          0
        )
      )

      peak.t <- max(abs(fx.t.sig))
      range.t <- max(exp(beta_g0 + fx.t.sig)) - min(exp(beta_g0 + fx.t.sig))
      range.t <- median.depth * range.t

      pv.t <- summary(fit)$s.table[t.term[1], 4]

      results[i, c("peak.t", "range.t", "pv.t")] <- c(peak.t, range.t, pv.t)
    }

    if (length(r.term) > 0) {
      r.idx <- match(r.term[1], term.names)

      fx.r <- pred$fit[, r.idx] / (1 + pred$se.fit[, r.idx]^2)

      if (return.fx) {
        fxs.r[i, ] <- fx.r
      }

      fx.r.sig <- ifelse(
        pred$fit[, r.idx] > 1.96 * pred$se.fit[, r.idx],
        pred$fit[, r.idx] - 1.96 * pred$se.fit[, r.idx],
        ifelse(
          pred$fit[, r.idx] < -1.96 * pred$se.fit[, r.idx],
          pred$fit[, r.idx] + 1.96 * pred$se.fit[, r.idx],
          0
        )
      )

      peak.r <- max(abs(fx.r.sig))
      range.r <- max(exp(beta_g0 + fx.r.sig)) - min(exp(beta_g0 + fx.r.sig))
      range.r <- median.depth * range.r

      pv.r <- summary(fit)$s.table[r.term[1], 4]

      results[i, c("peak.r", "range.r", "pv.r")] <- c(peak.r, range.r, pv.r)
    }
  }

  close(bar)

  out <- list()
  out$results <- as.data.frame(results)
  rownames(out$results) <- rownames(Y)

  if (return.fx) {
    rownames(fxs.t) <- rownames(Y)
    rownames(fxs.r) <- rownames(Y)

    try({
      out$fpca.t <- irlba::irlba(fxs.t)
      out$fpca.r <- irlba::irlba(fxs.r)
    })

    out$fxs.t <- fxs.t
    out$fxs.r <- fxs.r
  }

  return(out)
}