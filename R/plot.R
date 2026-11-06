#' Plot the null distribution of an SO-PLS CCI test
#'
#' Histogram of the null statistics `T*` (from rotation or permutation of X) with
#' the observed statistic `T_obs` marked. The p-value is the share of `T*` at or
#' below `T_obs` (left tail), so that tail is shaded.
#'
#' @param x a `SOPLSCCI` object from [SOPLS.CCI.test()]
#' @param breaks number of histogram bins (passed to [graphics::hist()])
#' @param main plot title; the default shows the formula, resampling method and p-value
#' @param xlab x-axis label
#' @param col_null,col_tail,col_obs colours of the null bars, the shaded
#'   rejection tail and the observed-statistic line
#' @param ... further arguments passed to [graphics::plot()]
#' @return `x`, invisibly
#' @examples
#' \dontrun{
#' r <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)
#' plot(r)
#' }
#' @export
plot.SOPLSCCI <- function(x, breaks = 30, main = NULL, xlab = expression(T == M[1] - M[0]),
                          col_null = "grey80", col_tail = "#d95f02", col_obs = "#1b9e77", ...) {
  T0 <- x$null.distribution
  Tobs <- unname(x$statistic)
  if (is.null(main)) {
    main <- sprintf("%s  (%s, R = %d)\np = %.3f", paste(deparse(x$formula), collapse = ""),
                    x$method, x$nperm, x$p.value)
  }
  h <- graphics::hist(c(T0, Tobs), breaks = breaks, plot = FALSE)
  graphics::plot(h, col = col_null, border = "white", main = main, xlab = xlab,
                 ylab = "Null samples", ...)
  usr <- graphics::par("usr")
  # shade the region at or below T_obs: the tail whose share of T* is the p-value
  graphics::rect(usr[1], usr[3], Tobs, usr[4], border = NA,
                 col = grDevices::adjustcolor(col_tail, alpha.f = 0.2))
  graphics::abline(v = Tobs, col = col_obs, lwd = 3)
  graphics::legend("topright", bty = "n", lwd = c(NA, NA, 3), pch = c(15, 15, NA),
                   col = c(col_null, grDevices::adjustcolor(col_tail, alpha.f = 0.4), col_obs),
                   legend = c(expression(paste("Null ", T^"*")), "Rejection tail (p-value)",
                              bquote(T[obs] == .(signif(Tobs, 3)))))
  invisible(x)
}
