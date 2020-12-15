## ----setup, include = FALSE---------------------------------------------------
library(eha)
options(digits = 4, fig.width=8)

## ----swe----------------------------------------------------------------------
head(swepop)
head(swedeaths)

## ----cleantab-----------------------------------------------------------------
dat <- swepop[, c("age", "sex", "year", "pop")]
dat$deaths <- swedeaths$deaths
rownames(dat) <- 1:NROW(dat) # Simplify rownames.
head(dat)
tail(dat)

## ----glm----------------------------------------------------------------------
fit.glm <- glm(deaths ~ offset(log(pop)) + I(year - 2000) + sex + 
                 factor(age), data = dat, family = poisson)
summary(fit.glm)$coefficients[2:3, ]

## ----glmbase------------------------------------------------------------------

lhaz <- coefficients(fit.glm)[-(2:3)]
n <- length(lhaz)
lhaz[-1] <- lhaz[-1] + lhaz[1]
haz <- exp(lhaz)

## ----glmbasefig,  fig.cap = "Age-specific mortality, Sweden 1968-2019. Poisson regression."----
oldpar <- par(las = 1, lwd = 1.5, mfrow = c(1, 2))
plot(0:(n-1), haz, type = "s", main = "log(hazards)", 
     xlab = "Age", ylab = "", log = "y")
plot(0:(n-1), haz, type = "s", main = "hazards", 
     xlab = "Age", ylab = "Deaths / Year")

## ----restore1, echo = FALSE---------------------------------------------------
par(oldpar) # Restore the default

## ----tpchreg1-----------------------------------------------------------------
fit <- tpchreg(oe(deaths, pop) ~ I(year - 2000) + sex, 
               time = age, last = 101, data = dat)

## ----showf--------------------------------------------------------------------
summary(fit)

## ----tpplot, fig.cap = "Age-specific mortality, Sweden 1968-2019. 'tpch' regression. Baseline refers to women and the year 2000."----
oldpar <- par(mfrow = c(1, 2), las = 1, lwd = 1.5)
plot(fit, fn = "haz", log = "y", main = "log(hazards)", 
     xlab = "Age")
plot(fit, fn = "haz", log = "", main = "hazards", 
     xlab = "Age", ylab = "Deaths / Year")

## ----tpTpchex-----------------------------------------------------------------
head(oldmort[, c("enter", "exit", "event", "sex", "civ", "birthdate")])
oldmort$birthyear <- floor(oldmort$birthdate) - 1800
om <- toTpch(Surv(enter, exit, event) ~ sex + civ + birthyear, 
             cuts = seq(60, 100, by = 2), data = oldmort)
head(om)

## ----rbef---------------------------------------------------------------------
fit3 <- tpchreg(oe(event, exposure) ~ sex + civ + 
                  birthyear, time = age, data = om)
summary(fit3)


## ----plotom, fig.cap = "Old age mortality, Skellefteå 1860-1880."-------------
oldpar <- par(mfrow = c(1, 2), las = 1, lwd = 1.5)
plot(fit3, fn = "haz", log = "y", main = "log(hazards)", 
     xlab = "Age", ylab = "log(Deaths / Year)", col = "blue")
plot(fit3, fn = "haz", log = "", main = "hazards", 
     xlab = "Age", ylab = "Deaths / Year", col = "blue")

## ----plotsurcum, fig.cap = "Old age mortality, Skellefteå 1860-1880. Cumulative hazards and survivor functions."----
oldpar <- par(mfrow = c(1, 2), las = 1, lwd = 1.5)
plot(fit3, fn = "cum", log = "y", main = "Cum. hazards", 
     xlab = "Age", col = "blue")
plot(fit3, fn = "sur", log = "", main = "Survivor function", 
     xlab = "Age", col = "blue")
par(oldpar)

## ----coxorig------------------------------------------------------------------
fit4 <- coxreg(Surv(enter, exit, event) ~ sex + civ + I(birthdate - 1800), 
               data = oldmort)
summary(fit4)

## ----coxgraphs, fig.cap = "Old age mortality, Skellefteå 1860-1880. Cox regression with original data."----
oldpar <- par(mfrow = c(1, 2), lwd = 1.5, las = 1)
plot(fit4, main = "Cumulative hazards", xlab = "Age", 
     col = "blue")
plot(fit4, main = "Survivor function", xlab = "Age", 
     fn = "surv", col = "blue")

## ----echo = FALSE-------------------------------------------------------------
par(oldpar)

