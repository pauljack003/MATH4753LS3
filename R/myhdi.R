#' @title Highest Density Interval
#'
#' @param funame q R function
#' @param credwidth internval width desired
#' @param tol this is the tolerance for the q function
#' @param ... parameters for the q function
#'
#' @return list of two components, L and U
#' @examples
#' # Example with Beta distribution
#' myhdi(funame = qbeta, credwidth = 0.95, tol = 1e-8, shape1 = 5, shape2 = 7)
#'
#' # Example with Gamma distribution
#' myhdi(funame = qgamma, credwidth = 0.95, tol = 1e-8, rate = 2, shape = 4)
#'
#' @importFrom stats dbeta dchisq dgamma dnorm optimize
#' @importFrom graphics abline curve
#' @export
myhdi <- function(funame, credwidth = 0.95, tol = 1e-8, ...) {

  incredwidth <- 1.0 - credwidth
  intervalw <- function(lowertail, funame, credwidth, ...) {
    funame(credwidth + lowertail, ...) - funame(lowertail, ...)
  }

  LI <- stats::optimize(f = intervalw, interval = c(0, incredwidth),
                 funame = funame, credwidth = credwidth,
                 tol = tol,
                 ...)

  hdilTail <- LI$minimum
  BCI <- list(L = funame(hdilTail, ...),
              U = funame(credwidth + hdilTail, ...))

  # Map the quantile function (funame) to the corresponding density function
  density_fun <- switch(
    as.character(substitute(funame)),
    qbeta = dbeta,
    qgamma = dgamma,
    qnorm = dnorm,
    qchisq = dchisq,
    stop("Unsupported distribution")
  )

  # Generate values for the x-axis
 x_vals <- seq(min(BCI$L, 0), max(BCI$U, 3), length.out = 200)

  # Plot the dynamically chosen density function using an anonymous function to pass additional parameters
  x <- NULL
  curve(density_fun(x, ...),
        from = min(x_vals), to = max(x_vals),
        main = "Density Curve with HDI",
        xlab = "x", ylab = "Density")

  # Add vertical lines for HDI bounds
  abline(v = c(BCI$L, BCI$U), col = "green", lwd = 2)

  # Return the HDI bounds
  return(BCI)
}
