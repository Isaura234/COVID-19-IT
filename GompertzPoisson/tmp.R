library(bbmle)
source("setup.R")

set.seed(1)
n = 50
x = seq(n)
mu = Gompertz(x, c(asym = 100, rate = 0.2, infl = 20))
y = rpois(n, mu)
Y = cumsum(y)
x = x[Y>0]
y = y[Y>0]
mu = mu[Y > 0]
Y = Y[Y>0]
plot(x, Y, type = "b")
points(x, y, type = "b", pch = 8)
lines(x, mu)

loglik <- function(pars, data)
{
  mu <- Gompertz(data$x, pars)
  sum(dpois(data$y, mu, log = TRUE))
}

loglik(c(asym = 100, rate = 0.2, infl = 20), data = list(x = x, y = y))

start = Gompertz_initial(data = list(x = x, y = y))
loglik(start, data = list(x = x, y = y))

MLE = optim(loglik, par = start, data = list(x = x, y = y),
            method = "BFGS", hessian = TRUE,
            control = list(trace = 2, fnscale = -1))
MLE$par
sqrt(diag(solve(-MLE$hessian)))

yfit = Gompertz(data$x, MLE$par)
plot(yfit, y - yfit) 
                  
plot(x, y, type = "b", pch = 8)
lines(x, mu)
lines(x, Gompertz(data$x, MLE$par), col = 3)

plot(x, Y, type = "b")
lines(x, GompertzCurve(data$x, MLE$par), col = 3)






PoissonGompertz <- function(y, x, wts)
{
  
  loglik <- function(pars, data)
  {
    mu <- Gompertz(data$x, pars)
    sum(data$wts*dpois(data$y, mu, log = TRUE))
  }
  
  data <- data.frame(y = as.vector(y),
                     x = as.vector(x),
                     wts = as.vector(wts/mean(wts)))
  data <- subset(data, y > 0)
  n <- nrow(data)
  
  start <- Gompertz_initial(data)
  MLE <- optim(loglik, par = start, data = data, 
               method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1))
  
  df <- length(MLE$par)
  out <- list(x = data$x, y = data$y, wts = data$wts, 
              estimate = MLE$par, 
              se = sqrt(diag(solve(-MLE$hessian))),
              loglik = MLE$value, df = length(MLE$par),
              AIC = 2*MLE$value - 2*df,
              BIC = 2*MLE$value - log(n)*df,
              fitted = Gompertz(data$x, MLE$par),
              curvefit = GompertzCurve(data$x, MLE$par))
  return(out)
}
              

mod = PoissonGompertz(y = COVID19$new_casi, x = COVID19$day, wts = COVID19$new_tamponi)
unlist(mod[c("loglik", "df", "AIC", "BIC")])
as.data.frame(mod[c("estimate","se")])

with(mod, plot(x, y, cex = bubble2Size(wts)$cex, type = "b"))
with(mod, lines(x, fitted, col = 3))
abline(v = mod$estimate[3], lty = 2)

with(mod, plot(x, cumsum(y), cex = bubble2Size(wts)$cex, type = "b", ylim = c(0,mod$estimate[1])))
with(mod, lines(x, curvefit, col = 3))
abline(h = mod$estimate[1], lty = 2)

with(mod, plot(fitted, y - fitted)); abline(h = 0, lty = 2) 
  
plot(mod$x, Gompertz_R0(mod$x, mod$estimate, ET = 7.5, type = 2), type = "b", ylim = c(0,3))
abline(h = 1, lty = 2)  
abline(v = mod$estimate[3], lty = 2)
Gompertz_R0(max(mod$x), mod$estimate, ET = 7.5, type = 2)
