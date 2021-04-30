## ----setup, include = FALSE---------------------------------------------------
library(eha)
options(digits = 4)

## ----gowecof8, fig.cap = "Baseline cumulative hazards for Cox and Gompertz regressions."----
fit.c <- coxreg(Surv(enter - 60, exit - 60, event) ~ sex + region,
                data = oldmort)
fit.g <- phreg(Surv(enter - 60, exit - 60, event) ~ sex + region, 
             dist = "gompertz", param = "rate", data = oldmort)
plot(fit.c, xlab = "Years above age 60.")
haz.g <- hazards(fit.g, cum = TRUE)
lines(haz.g$x, haz.g$y, lty = 2)
legend("topleft", legend = c("Cox regression", "Gompertz regression"), lty = 1:2)

## ----gocofph8-----------------------------------------------------------------
summary(fit.g)

## ----cocofph8-----------------------------------------------------------------
summary(fit.c)

## ----gocofaft8----------------------------------------------------------------
fit.gaft <- aftreg(Surv(enter - 60, exit - 60, event) ~ sex + region, 
                   id = id, param = "lifeExp", dist = "gompertz", 
                   data = oldmort)
summary(fit.gaft)

