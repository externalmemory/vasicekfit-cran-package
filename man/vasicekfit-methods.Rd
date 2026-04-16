\name{vasicekfit-methods}
\alias{print.vasicekfit}
\alias{summary.vasicekfit}
\alias{coef.vasicekfit}
\alias{fitted.vasicekfit}
\alias{residuals.vasicekfit}
\alias{vcov.vasicekfit}
\alias{confint.vasicekfit}
\title{Methods for vasicekfit Objects}
\description{
Print, summary, coefficient extraction, covariance, confidence interval,
and accessor methods for \code{"vasicekfit"} objects.
}
\usage{
\method{print}{vasicekfit}(x, \dots)
\method{summary}{vasicekfit}(object, \dots)
\method{coef}{vasicekfit}(object, \dots)
\method{fitted}{vasicekfit}(object, \dots)
\method{residuals}{vasicekfit}(object, \dots)
\method{vcov}{vasicekfit}(object, \dots)
\method{confint}{vasicekfit}(object, parm, level = 0.95, \dots)
}
\arguments{
  \item{x, object}{a \code{\link{vasicekfit}} object.}
  \item{parm}{a specification of which parameters to compute confidence
    intervals for, either a vector of names or indices. If missing, all
    parameters are considered.}
  \item{level}{the confidence level required.}
  \item{\dots}{additional arguments (currently unused).}
}
\details{
\code{coef} returns the recovered Vasicek parameters (p, rho, and kappas)
in the original parameter space.

\code{fitted} and \code{residuals} return values in probit (Y) space,
consistent with the underlying OLS regression.

\code{vcov} returns the delta-method covariance matrix for the
original-space parameters \eqn{(p, \rho, \kappa_1, \ldots, \kappa_m)},
obtained by propagating the joint covariance of \eqn{(\hat\beta, \hat\sigma^2)}
through the closed-form recovery transform. The \eqn{(\hat\beta, \hat\sigma^2)}
block is treated as block-diagonal by the independence of sample mean and
variance under normality. The variance of \eqn{\hat\sigma^2} uses
\eqn{2 \sigma^4 (N - m - 1)/N^2} for the MLE and \eqn{2 \sigma^4 / (N - m - 1)}
when \code{bias_correct = TRUE}.

\code{confint} returns Wald intervals based on \code{vcov}.
}
\value{
\code{print} and \code{summary} return \code{object} invisibly.

\code{coef} returns a named numeric vector.

\code{fitted} and \code{residuals} return numeric vectors.

\code{vcov} returns a \eqn{(m+2) \times (m+2)} covariance matrix.

\code{confint} returns a two-column matrix of lower and upper bounds.
}
\seealso{
\code{\link{vasicekfit}}, \code{\link{predict.vasicekfit}}
}
