#' @export
print.vasicekfit <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nVasicek model parameters:\n")
  cat("  p (PD)  =", format(x$p, digits = 6), "\n")
  cat("  rho (AVC) =", format(x$rho, digits = 6), "\n")
  if (length(x$kappa) > 0) {
    cat("  kappa:\n")
    print(x$kappa, digits = 6)
  }
  invisible(x)
}

#' @export
summary.vasicekfit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)

  cat("\nVasicek model parameters:\n")
  cat("  p (PD)    =", format(object$p, digits = 6), "\n")
  cat("  rho (AVC) =", format(object$rho, digits = 6), "\n")
  cat("  sigma2    =", format(object$sigma2, digits = 6), "\n")
  if (length(object$kappa) > 0) {
    cat("\n  Macro-factor sensitivities (kappa):\n")
    print(object$kappa, digits = 6)
  }
  cat("\nBias correction:", object$bias_correct, "\n")
  if (!is.null(object$portfolio_size)) {
    cat("Portfolio size:", object$portfolio_size, "\n")
  }

  cat("\n--- Underlying OLS regression (probit space) ---\n")
  print(summary(object$lm_fit))

  invisible(object)
}

#' @export
coef.vasicekfit <- function(object, ...) {
  c(p = object$p, rho = object$rho, object$kappa)
}

#' @export
fitted.vasicekfit <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.vasicekfit <- function(object, ...) {
  object$residuals
}
