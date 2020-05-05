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

day2date <- function(x, origin) as.Date(as.numeric(x), origin = origin)

cols = c("black", tableau_color_pal(palette = "Classic 10")(10))

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


SSexponential <- structure(
  function(x, th1, th2) 
  {
    .value <- th1 * exp(th2 * x)
    .actualArgs <- as.list(match.call()[c("th1", "th2")])
    if(all(vapply(.actualArgs, is.name, NA))) 
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("th1", "th2")))
      .grad[, "th1"] <- exp(th2*x)
      .grad[, "th2"] <- x * th1 * exp(th2 * x)
      dimnames(.grad) <- list(NULL, .actualArgs)
      attr(.value, "gradient") <- .grad
    }
    .value
  }, 
  initial = function (mCall, data, LHS) 
  {
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 3) 
      stop("too few distinct input values to fit the Exponential model")
    xyL <- xy
    pars <- coef(lm(log(y) ~ x, data = xyL))
    pars <- c(exp(pars[1]), pars[2])
    setNames(pars, mCall[c("th1", "th2")])
  }, 
  pnames = c("th1", "th2"), 
  class = "selfStart"
)

Exponential_gradient <- function(x, pars)
{
  pars[1] * pars[2] * exp(pars[2]*x)
}

Exponential_growthrate <- function(x, pars)
{
# growth rate = d^2 GrothCurve(x)/x^2 / d GrothCurve(x)/x
  pars[2]
}

Exponential_R0 <- function(x, pars, ET, type = 2)
{
# Estimate reproduction number R_0 from Exponential growth
# pars = parameters of nonlinear Exponential growth model
# ET = Expected incubation time, 7.5 days for COVID-19
# type = algorithm used for computing R0
  t = floor(ET)
  if(type == 2)
  {
    r = Exponential_growthrate(x, pars)
    R = exp(r*ET)
  } else
  {
    dyhat = Exponential_gradient(x, pars)
    beta = dyhat/frollsum(shift(dyhat), t)
    R = ET * beta
  }
  return(R)
}


GompertzCurve <- function(x, pars)
{
  asym <- pars[1]
  rate <- pars[2]
  infl <- pars[3]
  asym * exp(-exp(-rate*(x - infl)))
}

Gompertz_initial <- function(data)
{
  xy <- sortedXyData("x", "y", data)
  xy$y <- cumsum(xy$y)
  if(nrow(xy) < 4) 
    stop("too few distinct input values to fit the Gompertz model")
  xyL <- xy
  xyL$y <- log(abs(xyL$y))
  pars <- NLSstAsymptotic(xyL)
  asym <- max(exp(pars["b0"] + pars["b1"]*(1-exp(-exp(pars["lrc"]) * max(xyL["x"])))),
              1.1*max(xy$y, na.rm=TRUE))
  xy$y <- log(log(asym) - xyL$y)
  pars <- coef(lm(y ~ x, data = xy))
  pars <- c("asym" = asym, "rate" = -pars[[2]], "infl" = -pars[[1]]/pars[[2]])
  return(pars)
}

Gompertz <- function(x, pars)
{
  asym <- pars[1]
  rate <- pars[2]
  infl <- pars[3]
  asym * rate * exp(-rate*(x - infl) - exp(-rate*(x - infl)))
}

Gompertz_growthrate <- function(x, pars)
{
# growth rate = d^2 GrothCurve(x)/x^2 / d GrothCurve(x)/x
  asym <- pars[1]
  rate <- pars[2]
  infl <- pars[3]
  # d1 = asym * rate * exp(-rate*(x - infl) - exp(-rate*(x - infl)))
  # d2 = asym * rate * exp(-rate*(x - infl) - exp(-rate*(x - infl))) * (rate*exp(-rate*(x - infl)) - rate)
  # d2/d1
  rate*exp(-rate*(x - infl)) - rate
}

Gompertz_R0 <- function(x, pars, ET, type = 2)
{
# Estimate reproduction number R_0 from Gompertz growth
# x = time variable
# pars = parameters of nonlinear Gompertz growth model
# ET = Expected incubation time, 7.5 days for COVID-19
# type = algorithm used for computing R0
  t = floor(ET)
  if(type == 2)
  {
    r = Gompertz_growthrate(x, pars)
    R = exp(r*ET)
  } else
  {
    dyhat = Gompertz_gradient(x, pars)
    beta = dyhat/frollsum(shift(dyhat), t)
    R = ET * beta
  }
  return(unname(R))
}

# Moving Block Bootstrap (MBB) forecast for PoissonGompertz model

# this version uses the current fitted object to estimate the trend. 
# Residuals are then bootstrapped on moving blocks to create bootstrap 
# samples

MBB <- function(x, block.size) 
{
  n <- length(x)
  nw <- floor(n/block.size)
  bx <- rep(0.0, (nw + 2) * block.size)
  for(i in 1:(nw + 2))
  {
    s <- sample(1:(n - block.size + 1), 1)
    bx[((i - 1) * block.size + 1):(i * block.size)] <- 
      x[s:(s + block.size - 1)]
  }
  start_from <- sample(0:(block.size - 1), 1) + 1
  return(bx[start_from:(start_from + n - 1)])
}

predictMBB.nls <- function(object, newdata, 
                           nboot = 999, 
                           conf.level = 0.95, 
                           levels = c((1 - conf.level)/2, 
                                      (1 + conf.level)/2),
                           block.size = NULL,
                           verbose = interactive())
{

  stopifnot(inherits(object, "nls"))
  stopifnot(!missing(newdata))
  #
  # browser()
  cstart <- eval(object$call$start)
  cdata <- eval(object$call$data)
  cfmla <- object$call$formula
  y <- cdata[[as.character(cfmla)[[2]]]]
  n <- length(y)
  df <- n - length(coef(object))
  
  #####################################################################
  # generate bootstrap samples
  if(is.null(block.size)) 
    block.size <- round(3*n^(1/3))
  trend <- fitted(object)
  residuals <- residuals(object)
  meanres <- mean(residuals)
  residuals <- residuals - meanres
  yboot <- vector(mode = "list", length = n)
  # browser()
  for(i in 1:nboot) 
  {
    yboot[[i]] <- trend + meanres + MBB(residuals, block.size)
    # plot(y, type = "b")
    # points(yboot[[i]], type = "b", col = 2)
    # pause()
  }
  
  #####################################################################

  if(verbose) 
  { 
    flush.console()
    cat("Moving block bootstrap ...\n")
    pbar <- txtProgressBar(min = 0, max = nboot, style = 3)
    on.exit(close(pbar))
  }
  
  # browser()
  newdata <- as.data.frame(newdata)
  boot <- matrix(as.double(NA), nrow = nrow(newdata), ncol = nboot)
  for(b in 1:nboot)
  {
    if(verbose) 
      setTxtProgressBar(pbar, b) 
    
    call <- object$call
    cdata[,as.character(cfmla)[[2]]] <- yboot[[b]]
    call$data <- cdata
    mod <- try(eval(call), silent = TRUE)
    if(inherits(mod, "try-error"))
    { 
      call$start <- coef(object)
      mod <- try(eval(call), silent = TRUE)
    }
    if(inherits(mod, "try-error")) next()
    #
    ypred <- predict(mod, newdata)
    rstdev <- sqrt(deviance(mod)/df)
    ypred <- ypred + rnorm(length(ypred), mean = 0, sd = rstdev)
    boot[,b] <- ypred
    # plot(data$y, xlim = c(0,max(newdata$x)),
    #      ylim = range(0,predict(object, newdata)*1.2))
    # lines(1:max(newdata$x), 
    #       predict(object, data.frame(x = 1:max(newdata$x))))
    # lines(1:max(newdata$x), 
    #       predict(mod, data.frame(x = 1:max(newdata$x))), col = 2)
    # points(yboot[[b]], pch = 15, col = "grey")
    # points(newdata$x, ypred, pch = 20, col = 2)
    # pause()
  }
  #
  ci <- t(apply(boot, 1, quantile, probs = levels, na.rm = TRUE))
  out <- data.frame(fit = predict(object, newdata),
                    lwr = ci[,1], upr = ci[,2])
  return(out)
}


