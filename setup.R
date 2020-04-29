library(knitr)
opts_chunk$set(fig.align = "center",
               fig.width = 6, fig.height = 5.5,
               dev.args = list(pointsize=10),
               out.width = "90%", dpi = 300,
               cache = FALSE,
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & ouput code in chunks
               warning = FALSE, message = FALSE)

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
tab_theme <- ttheme_default(
  core = list(# bg_params = list(fill = NA, col=NA),
              fg_params=list(fontface=1, fontsize = 10)),
  colhead = list(fg_params=list(col="black", fontface=2, fontsize = 10)),
  rowhead = list(fg_params=list(col="black", fontface=2, hjust=0, x=0,
                                fontsize = 10),
                 bg_params=list(fill = c("grey95", "grey90"))))

source("misc/predictMBBnls.R")  # experimental

logpos <- function(x) ifelse(x > 0, log(x), 1)

AICc <- function(object, ..., k = 2) 
{
  n      <- nobs(object)
  loglik <- logLik(object)
  df     <- attr(loglik, "df")
  AIC <- -2 * c(loglik) + k * df
  AICc <- AIC + k * df * (df +1)/(n - df - 1)
  return(AICc)
}

cols <- c("Exponential" = "red3", "Logistic" = "dodgerblue3",
          "Gompertz" = "forestgreen", "Richards" = "azure4")

peaks <- function (x, span = 3) 
{
  z <- embed(x, span)
  s <- span%/%2
  v <- max.col(z, ties.method = "first") == (1 + s)
  result <- c(rep(FALSE, s), v, rep(FALSE, s))
  result
}


library(numDeriv)

exponential <- function(x, th1, th2) th1 * exp(th2 * x)

exponentialGrad <- function(object, x, ...)
{
  th <- coef(object)
  f <- function(x) exponential(x, th[1], th[2])
  numDeriv::grad(f, x)
}

logisticGrad <- function(object, x, ...)
{
  th <- coef(object)
  f <- function(x) SSlogis(x, th[1], th[2], th[3])
  numDeriv::grad(f, x)
}

gompertzGrad <- function(object, x, ...)
{
  th <- coef(object)
  f <- function(x) SSgompertz(x, th[1], th[2], th[3])
  numDeriv::grad(f, x)
}

richardsGrad <- function(object, x, ...)
{
  th <- coef(object)
  f <- function(x) richards(x, th[1], th[2], th[3])
  numDeriv::grad(f, x)
}

reldiff <- function(x) c(NA, diff(x)/x[-length(x)])
