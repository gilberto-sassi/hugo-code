---
title: "GeoFourierFDA: R package for Kriging Methods to Function-valued Data"
date: 2020-07-31T14:34:08-03:00
draft: false
tags: ["Pacotes R"]
---

In this package, we have observed *n* curves
*χ*<sub>**s**<sub>1</sub></sub>(*t*), ⋯, *χ*<sub>**s**<sub>*n*</sub></sub>(*t*)
at locations **s**<sub>1</sub>, ⋯, **s**<sub>*n*</sub> in a specified
region, where
**s**<sub>*i*</sub> = (*θ*<sub>*i*</sub>, *η*<sub>*i*</sub>), *i* = 1, …, *n*,
*θ*<sub>*i*</sub> is the latitute and *η*<sub>*i*</sub> is the
longitude. The aim is to estimate the unobserved curve
*χ*<sub>**s**<sub>0</sub></sub>(*t*) where
**s**<sub>0</sub> ∉ {**s**<sub>1</sub>, ⋯**s**<sub>*n*</sub>}. The ideia
proposed by Giraldo, Delicado, and Mateu (2011) was uncomplicated: the
curve *χ*<sub>**s**<sub>0</sub></sub>(*t*) at
**s**<sub>0</sub> ∉ {**s**<sub>1</sub>, ⋯**s**<sub>*n*</sub>} is a
linear combination of all the curves
*χ*<sub>**s**<sub>1</sub></sub>(*t*), ⋯, *χ*<sub>**s**<sub>*n*</sub></sub>(*t*),
i.e., *χ*<sub>**s**<sub>0</sub></sub>(*t*)=
*λ*<sub>1</sub>*χ*<sub>**s**<sub>1</sub></sub>(*t*) + ⋯ + *λ*<sub>1</sub>*χ*<sub>**s**<sub>*n*</sub></sub>(*t*),
where *λ*<sub>1</sub>, ⋯, *λ*<sub>*n*</sub> is the solution of the
linear system

![Sistema Linear](/files/geoFourierFDA/linear_system.png)

where *μ* is a constant and the function
*γ*(*h*) = ∫*γ*(*h*; *t*)*d**t*, (1)
is called trace-variogram and *γ*(*h*, *t*) is the semivariogram. More
precisely, for each *t* a weakly stationary and isotropic spatial
process is assumed and the semivariogram can be computed. Furthermore,
the trace-variogram is an integration of *γ*(*h*; *t*) over *t*.
Usually, the integral in the equation (1) is approximated using an
modified version of empirical semivariogram (see Giraldo, Delicado, and
Mateu 2011 for more details). In this paper, we propose an approach
using Legendre-Gauss quadrature, which is intuitive since it explicitly
uses the definition of trace-variogram.

Installing
==========

    devtools::install_github('gilberto-sassi/geoFourierFDA')

How to interpolate a curve at an unmonitored location
=====================================================

```r
# interpolating curve at Halifax using all remaining curves in the functional dataset
data(canada)

# removing the data from Halifax and getting its coordinate
l_item <- 2
log_item <- !(1 : ncol(m_data) %in% l_item)
y_true <- m_data[, l_item]
t_true <- seq(from = -pi, to = pi, length.out = nrow(m_data))
new_coord <- m_coord[l_item, ] %>% matrix(ncol = 2)
m_data <- m_data[,log_item]
m_coord <- m_coord[log_item, ]

# Fourier series polynomial
m <- ceiling(1 + log2(nrow(m_data))) # Sturge's rule

# Estimating the curve at halifax using all remained curves
y_s0_mine <- geo_fda(m_data, m_coord, new_coord,
                     t = t_true, m = m)
```

Coefficients of smoothing using Fourier series polynomial
=========================================================

```r
# Coefficients of smoothing using Fourier series polynomial
# Coefficients of smoothing at Halifax
data(canada)

# the data from Halifax and getting
l_item <- 2
log_item <- !(1 : ncol(m_data) %in% l_item)
f <- m_data[, l_item]

# order of Fourier series polynomial
m <- (1 + length(f) %>% log2()) %>% ceiling() 

coefs <- coef_fourier(f, m)
```

Smoothed curve using Fourier series
===================================

```r
# Coefficients of smoothing using Fourier series polynomial
# Coefficients of smoothing at Halifax
data(canada)

# the data from Halifax and getting
l_item <- 2
log_item <- !(1 : ncol(m_data) %in% l_item)
y <- m_data[, l_item]

# order of Fourier series polynomial
m <- (1 + length(f) %>% log2()) %>% ceiling() 

# coefficients of Fourier series polynomial
coefs <- coef_fourier(f, m)

# points to evaluate curve at interval [-pi, pi]
x <- seq(from = -pi, to = pi, by = 0.01)

# smoothed curve at some points
y_est <- fourier_b(coefs, x)
```

Installation
============

To install this package in R, it can be used the R package `devtools`: 

```r
devtools::install_github('gilberto-sassi/geoFourierFDA')
```

Referências
===========

Giraldo, R, P Delicado, and J Mateu. 2011. “Ordinary Kriging for
Function-Valued Spatial Data.” *Environmental and Ecological Statistics*
18 (3): 411–26.
