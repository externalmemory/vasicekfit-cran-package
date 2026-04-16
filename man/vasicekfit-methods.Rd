\name{vasicekfit-methods}
\alias{print.vasicekfit}
\alias{summary.vasicekfit}
\alias{coef.vasicekfit}
\alias{fitted.vasicekfit}
\alias{residuals.vasicekfit}
\title{Methods for vasicekfit Objects}
\description{
Print, summary, coefficient extraction, and accessor methods for
\code{"vasicekfit"} objects.
}
\usage{
\method{print}{vasicekfit}(x, \dots)
\method{summary}{vasicekfit}(object, \dots)
\method{coef}{vasicekfit}(object, \dots)
\method{fitted}{vasicekfit}(object, \dots)
\method{residuals}{vasicekfit}(object, \dots)
}
\arguments{
  \item{x, object}{a \code{\link{vasicekfit}} object.}
  \item{\dots}{additional arguments (currently unused).}
}
\details{
\code{coef} returns the recovered Vasicek parameters (p, rho, and kappas)
in the original parameter space.

\code{fitted} and \code{residuals} return values in probit (Y) space,
consistent with the underlying OLS regression.
}
\value{
\code{print} and \code{summary} return \code{object} invisibly.

\code{coef} returns a named numeric vector.

\code{fitted} and \code{residuals} return numeric vectors.
}
\seealso{
\code{\link{vasicekfit}}, \code{\link{predict.vasicekfit}}
}
