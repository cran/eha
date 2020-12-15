## ----setup, include = FALSE---------------------------------------------------
library(eha)
options(digits = 4, fig.width=8)

## ----6hazs,fig.cap = "Selected hazard functions.", echo=FALSE,fig.height=5,message=FALSE----
library(eha)
oldpar <- par(mfrow = c(2, 2))
x <- seq(0.001, 10.001, length = 1000)
y <- hweibull(x, shape = 3, scale = 6)
plot(x, y, main = "Weibull, shape = 3", type = "l", 
     ylab = "", xlab = "Time", lwd = 2, ylim = c(0, 1.5), las = 1)
abline(h = 0, v = 0)
y <- hweibull(x, shape = 1/2, scale = 2)
plot(x, y, main = "Weibull, shape = 1/3", type = "l",
     ylab = "", xlab = "Time", ylim = c(0, 1.5), lwd = 2, las = 1)
abline(h = 0, v = 0)
y <- hgompertz(x, shape = 1/6, rate = 0.2, param = "rate")
plot(x, y, main = "Gompertz, rate = 2", type = "l",
     ylab = "", xlab = "Time", ylim = c(0, 1.5), lwd = 2, las = 1)
abline(h = 0, v = 0)
y <- hgompertz(x, shape = 1.5, rate = -0.2, param = "rate")
plot(x, y, main = "Gompertz, rate = -2", type = "l",
     ylab = "", xlab = "Time", ylim = c(0, 1.5), lwd = 2, las = 1)
abline(h = 0, v = 0)
par(oldpar)

## ----gowecof8, echo = TRUE----------------------------------------------------
fit.g <- phreg(Surv(enter - 60, exit - 60, event) ~ sex + civ + region, 
             dist = "gompertz", data = oldmort)
summary(fit.g)	     

## ----coxold8, echo = TRUE-----------------------------------------------------
fit.c <- coxreg(Surv(enter - 60, exit - 60, event) ~ sex + civ + region, data = oldmort)
summary(fit.c)

## ----check, echo = TRUE-------------------------------------------------------
check.dist(fit.c, fit.g)

## ----aftph6,fig.cap="Proportional hazards (left) and accelerated failure time model (right). The baseline distribution is Loglogistic with shape 5 (dashed).",echo=FALSE----
op <- par(las = 1)
x <- seq(0, 3, length = 1000)
par(mfrow = c(1, 2))
plot(x, 2 * hllogis(x, shape = 5), type = "l", ylab = "", main = "PH", xlab = "Time")
lines(x, hllogis(x, shape = 5), lty = 2)
plot(x, 2 * hllogis(2 * x, shape = 5), type = "l", ylab = "", main = "AFT", xlab = "Time")
lines(x, hllogis(x, shape = 5), lty = 2)
par(op)

## ----oldmort6.aft-------------------------------------------------------------
fit.g1 <- aftreg(Surv(enter - 60, exit - 60, event) ~ sex + civ + region, data = oldmort, dist = "gompertz")
summary(fit.g1)

