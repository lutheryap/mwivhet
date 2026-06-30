# Inference with Many Weak Instruments and Heterogeneity

## Introduction

The package `mwivhet` implements a procedure for valid inference for the
many instrumental variable (IV) model with many weak instruments and
heterogeneous treatment effects, based on Yap (2026). A leading
application of this setting is the examiner design, where judge or
caseworker indicators are instruments. In this setting, the number of
instruments can be large relative to the sample size.

Existing inference procedures, such as those of Mikusheva and Sun (2022)
and Matsushita and Otsu (2022), are robust to many weak instruments but
assume homogeneity in treatment effects. Yap (2026) shows that these
procedures need not control type I error, and therefore have incorrect
size, under heterogeneity. This package implements the procedure from
Yap (2026), which uses a novel leave-three-out (L3O) variance estimator
and has both correct size and is asymptotically uniformly most powerful
unbiased under both many weak instruments and heterogeneous treatment
effects.

The point estimator is the unbiased jackknife IV estimator (UJIVE) of
Kolesár (2013), which removes own-observation dependence by applying the
leave-one-out construction so that both the instrument coefficients and
the covariate adjustment are estimated excluding observation $`i`$ and
the fitted value is constructed entirely from other observations. Under
the monotonicity condition of Imbens and Angrist (1994), UJIVE can be
interpreted as a weighted average of causal effects for individuals
whose treatment is affected by the instrument. When monotonicity fails,
the estimand remains well-defined. The confidence interval is obtained
by inverting the Lagrange Multiplier (LM) test, where the test statistic
uses a novel leave-three-out (L3O) variance estimator that is consistent
even when reduced-form coefficients are not consistently estimable.

## The Suffolk County Application

We demonstrate the package using the `suffolk` dataset, which contains
administrative records from the Suffolk County District Attorney’s
Office. The dataset is obtained from Agan, Doleac, and Harvey (2023),
who study the effect of misdemeanor prosecution on recidivism, which is
measured by criminal complaint in two years. Following Agan, Doleac, and
Harvey (2023), the paper frames prosecution as the treatment; here the
treatment variable $`X_i`$ is coded as an indicator for immediate
non-prosecution (`ng_immed_all`), the complement, and the outcome
$`Y_i`$ is an indicator for criminal complaint in two years
(`anyr_twoyears_arrest2`). Instruments $`Z_i`$ are prosecutor indicators
(`first_pros`), and covariates $`W_i`$ are constructed from
court-by-year and court-by-day-of-week combinations.

### Data Preparation

``` r
library(mwivhet)
library(dplyr)
library(fixest)
```

``` r
# 0. Setup ----------------------------------------------------------------
df <- suffolk

# 1. Preparation ----------------------------------------------------------
# 1.1. Variable Assignment
df <- df %>%
  mutate(
    X = ng_immed_all,
    Y = anyr_twoyears_arrest2,
    groupZ = first_pros
  )

# 1.2. Construct Covariate Groups (groupW)
df <- df %>%
  mutate(
    monthID = as.integer(factor(court_month2)),
    dowID   = as.integer(factor(court_dow2)),
    groupW  = as.integer(factor(paste(dowID, monthID, sep = "_")))
  )

# 1.3. Construct Interaction Groups (groupQ / group)
df <- df %>%
  mutate(
    groupQ = as.integer(factor(paste(groupW, groupZ, sep = "_"))),
    group  = groupQ
  )

# 1.4. Filter for Group Size
df <- df %>%
  group_by(groupQ) %>%
  filter(n() >= 4) %>%
  ungroup()

# 1.5. Get residualized objects
lmX <- fixest::feols(X ~ 1 | group, data = df)
df$MX <- lmX$residuals

lmY <- fixest::feols(Y ~ 1 | group, data = df)
df$MY <- lmY$residuals
```

### Inference

The UJIVE point estimate is computed via `GetLM`, which computes the
weighted cross-product $`\sum_i\sum_{j\neq i} G_{ij} A_i B_j`$ for two
input variables $`A`$ and $`B`$. The ratio of
$`\sum_{i}\sum_{j\neq i}G_{ij}X_iY_j`$ to
$`\sum_i\sum_{j\neq i}G_{ij}X_iX_j`$ gives
$`\hat\beta_{\mathrm{UJIVE}}`$.

The confidence interval is obtained by inverting the LM test (with L3O
variance estimator) via `GetCIcoef` and `GetCItypebd`. The test rejects
if $`KT_{LM}^2/\hat{V}_{LM} \geq \Phi^{-1}(1-\alpha/2)^2`$, where
$`\Phi(\cdot)`$ is the standard normal CDF. Since $`T_{LM}^2`$ and
$`\hat{V}_{LM}`$ are both quadratic in $`\beta_0`$, the confidence set
is obtained by solving a quadratic inequality. `GetCIcoef` returns the
coefficients $`(C_0, C_1, C_2)`$ of this quadratic, and `GetCItypebd`
solves the inequality and returns the confidence interval type and
bounds.

``` r
# 2. Inference ------------------------------------------------------------
# UJIVE point estimate
S_hat <- GetLM(df, X, X, groupW, group, noisy = FALSE)
UJIVE <- GetLM(df, X, Y, groupW, group, noisy = FALSE) / S_hat
UJIVE
#> [1] -0.1440279

# L3O confidence interval
L3OCIcoef <- GetCIcoef(df, groupW, group, X, Y, MX, MY, noisy = FALSE)
L3OCI <- GetCItypebd(L3OCIcoef)[2:3]
L3OCI
#> [1] -0.22183974 -0.06643859
```

### Results

``` r
# 3. Table ----------------------------------------------------------------
# 3.1. Create the matrix using calculated L3O variables
L3O_data <- matrix(c(
  L3OCI[1],                # Lower Bound
  L3OCI[2],                # Upper Bound
  UJIVE,                    # Point Estimate
  L3OCI[2] - L3OCI[1]      # CI Length
), ncol = 1)

# 3.2. Add labels for row/col
rownames(L3O_data) <- c("LB", "UB", "Estimate", "CIlength")
colnames(L3O_data) <- "L3O"

# 3.3 Round and display
CItab_L3O_round <- round(L3O_data, 3)
CItab_L3O_round
#>             L3O
#> LB       -0.222
#> UB       -0.066
#> Estimate -0.144
#> CIlength  0.155
```

The UJIVE estimate of $`-0.144`$ shows how non-prosecution reduces the
probability of criminal complaint within two years by about 14
percentage points. The L3O 95% confidence interval $`[-0.222, -0.066]`$
excludes zero. These results show that the findings in Agan, Doleac, and
Harvey (2023) are robust to many weak instruments and heterogeneous
treatment effects.

## References

Agan, Amanda, Jennifer L. Doleac, and Anna Harvey. 2023. “Misdemeanor
Prosecution.” *The Quarterly Journal of Economics* 138.

Imbens, Guido W., and Joshua D. Angrist. 1994. “Identification and
Estimation of Local Average Treatment Effects.” *Econometrica* 62.

Kolesár, Michal. 2013. “Estimation in an Instrumental Variables Model
with Treatment Effect Heterogeneity.” Working Papers 2013-2. Princeton
University, Economics Department.

Yap, Luther. 2026. “Inference with Many Weak Instruments and
Heterogeneity.”
