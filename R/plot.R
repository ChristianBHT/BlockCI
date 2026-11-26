#' Plot the null distribution of a block CCI test
#'
#' Histogram of the null statistics `T*` (from resampled copies of X) with
#' the observed statistic `T_obs` marked by a vertical line (one line per
#' observed split when `npval > 1`). For a parametric p-value the fitted normal
#' curve is added.
#'
#' @param x a `BlockCCI` object from [Block.CCI.test()]
#' @param breaks number of histogram bins (passed to [graphics::hist()])
#' @param main plot title; the default shows the formula, method, resampling and
#'   p-value
#' @param xlab x-axis label
#' @param col_null,col_obs colours of the null bars and the observed-statistic line
#' @param ... further arguments passed to [graphics::plot()]
#' @return `x`, invisibly
#' @examples
#' \dontrun{
#' r <- Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)
#' plot(r)
#' }
#' @export
plot.BlockCCI <- function(x, breaks = 30, main = NULL, xlab = expression(T == M[1] - M[0]),
                          col_null = "grey80", col_obs = "#1b9e77", ...) {
  T0 <- x$null.distribution
  Tobs <- if (is.null(x$T.splits)) unname(x$statistic) else x$T.splits   # one per split
  if (is.null(main)) {
    main <- sprintf("%s  (%s, %s, R = %d)\np = %.3f (%s)", paste(deparse(x$formula), collapse = ""),
                    x$method, x$resampling, x$nperm, x$p.value,
                    if (is.null(x$pvalue)) "empirical" else x$pvalue)
  }
  h <- graphics::hist(c(T0, Tobs), breaks = breaks, plot = FALSE)
  graphics::plot(h, col = col_null, border = "white", main = main, xlab = xlab,
                 ylab = "Null samples", ...)
  usr <- graphics::par("usr")
  param <- identical(x$pvalue, "parametric") && is.finite(stats::sd(T0)) && stats::sd(T0) > 0
  if (param) {   # fitted normal null, scaled to counts
    xs <- seq(usr[1], usr[2], length.out = 200)
    graphics::lines(xs, stats::dnorm(xs, mean(T0), stats::sd(T0)) * length(T0) * diff(h$breaks)[1],
                    col = "grey30", lwd = 2, lty = 2)
  }
  graphics::abline(v = Tobs, col = col_obs, lwd = 3)
  graphics::legend("topright", bty = "n",
                   lwd = c(NA, 3, if (param) 2), lty = c(NA, 1, if (param) 2),
                   pch = c(15, NA, if (param) NA),
                   col = c(col_null, col_obs, if (param) "grey30"),
                   legend = c(expression(paste("Null ", T^"*")),
                              if (length(Tobs) == 1) bquote(T[obs] == .(signif(Tobs, 3)))
                              else bquote(T[obs] ~ "(" * .(length(Tobs)) ~ "splits)"),
                              if (param) "Normal fit"))
  invisible(x)
}
