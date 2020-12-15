## ----setup, include = FALSE---------------------------------------------------
library(eha)
options(digits = 4)

## ----coxph--------------------------------------------------------------------
fit <- survival::coxph(Surv(enter, exit, event) ~ sex + civ + imr.birth, 
                       data = oldmort)
summary(fit)

## ----coxreg-------------------------------------------------------------------
fit <- coxreg(Surv(enter, exit, event) ~ sex + civ + imr.birth, data = oldmort)
summary(fit)

## ----boot, eval = FALSE-------------------------------------------------------
#  library(eha)
#  fit <- coxreg(Surv(enter, exit, event) ~ sex, data = oldmort, boot = 100,
#                control = list(trace = TRUE))

## ----mppl, cache = FALSE, eval = FALSE----------------------------------------
#  library(eha)
#  fit <- coxreg(Surv(enter, exit, event) ~ sex, data = oldmort, method = "ml")
#  plot(c(60, fit$hazards[[1]][, 1]), c(0, cumsum(fit$hazards[[1]][, 2])), type = "l")

