# vasicekfit

An R package implementing the extended Vasicek single-factor credit loss
model, where the probability of default depends on macroeconomic covariates.
All parameters, including the asset value correlation, are estimated in
closed form via probit-transformed OLS.

## Installation

From CRAN (once published):

```r
install.packages("vasicekfit")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("externalmemory/vasicekfit-cran-package")
```

## Example

```r
library(vasicekfit)

set.seed(1)
n <- 500

unemp <- rnorm(n)
hpi   <- rnorm(n)
z     <- rnorm(n)

y <- pnorm((qnorm(0.03) + 0.13 * unemp - 0.07 * hpi + sqrt(0.05) * z) / sqrt(0.95))
d <- data.frame(default_rate = y, unemp = unemp, hpi = hpi)

fit <- vasicekfit(default_rate ~ unemp + hpi, data = d)

coef(fit)
confint(fit)
```

<code>summary(fit)<code> prints a coefficient table with delta-method standard errors, z-values, and p-values for the recovered parameters (p, rho, kappa).

## References

Vasicek, O. A. (2002), The distribution of loan portfolio value. Risk, 15(12), 160-162.

Yang, Bill Huajian (2014), Estimating Long-Run PD, Asset Correlation, and Portfolio Level PD by Vasicek Models, MPRA Paper No. 57244 https://mpra.ub.uni-muenchen.de/57244/1/MPRA_paper_57244.pdf

Mayorov, Dmitriy, Macroeconomic Sensitivity in the Vasicek Credit Loss Model: Closed-Form Maximum Likelihood Estimation via OLS (April 01, 2026). Available at SSRN: https://ssrn.com/abstract=6506378 or http://dx.doi.org/10.2139/ssrn.6506378 

## License

MIT
