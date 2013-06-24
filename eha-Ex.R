pkgname <- "eha"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('eha')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BCoxreg")
### * BCoxreg

flush(stderr()); flush(stdout())

### Name: Coxreg
### Title: Cox regression
### Aliases: Coxreg
### Keywords: survival regression

### ** Examples


 dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
                status=c(1,1,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
 Coxreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
 


cleanEx()
nameEx("SurvSplit")
### * SurvSplit

flush(stderr()); flush(stdout())

### Name: SurvSplit
### Title: Split a survival object at specified durations.
### Aliases: SurvSplit
### Keywords: manip

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Y, cuts){
    if (NCOL(Y) == 2) Y <- cbind(rep(0, NROW(Y)), Y)
    indat <- cbind(Y, 1:NROW(Y), rep(-1, NROW(Y)))
    colnames(indat) <- c("enter", "exit", "event", "idx", "ivl")
    n <- length(cuts)
    cuts <- sort(cuts)
    if ((cuts[1] <= 0) || (cuts[n] == Inf))
        stop("'cuts' must be positive and finite.")
    cuts <- c(0, cuts, Inf)
    n <- n + 1
    out <- list()
    indat <- as.data.frame(indat)
    for (i in 1:n){
        out[[i]] <- age.window(indat, cuts[i:(i+1)])
        out[[i]]$ivl <- i
        out[[i]] <- t(out[[i]])
    }
    Y <- matrix(unlist(out), ncol = 5, byrow = TRUE)
    colnames(Y) <- colnames(indat)
    list(Y = Y[, 1:3],
         ivl = Y[, 5],
         idx = Y[, 4]
         )
  }



cleanEx()
nameEx("aftreg")
### * aftreg

flush(stderr()); flush(stdout())

### Name: aftreg
### Title: Accelerated Failure Time Regression
### Aliases: aftreg
### Keywords: survival regression

### ** Examples

data(mort)
aftreg(Surv(enter, exit, event) ~ ses, data = mort)



cleanEx()
nameEx("age.window")
### * age.window

flush(stderr()); flush(stdout())

### Name: age.window
### Title: Age cut of survival data
### Aliases: age.window
### Keywords: survival

### ** Examples

dat <- data.frame(enter = 0, exit = 5.731, event = 1, x = 2)
window <- c(2, 5.3)
dat.trim <- age.window(dat, window)  



cleanEx()
nameEx("cal.window")
### * cal.window

flush(stderr()); flush(stdout())

### Name: cal.window
### Title: Calendar time cut of survival data
### Aliases: cal.window
### Keywords: survival

### ** Examples

dat <- data.frame(enter = 0, exit = 5.731, event = 1,
birthdate = 1962.505, x = 2)
window <- c(1963, 1965)
dat.trim <- cal.window(dat, window)  



cleanEx()
nameEx("check.dist")
### * check.dist

flush(stderr()); flush(stdout())

### Name: check.dist
### Title: Graphical goodness-of-fit test
### Aliases: check.dist
### Keywords: distribution

### ** Examples

data(mort)
oldpar <- par(mfrow = c(2, 2))
fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
fit.w <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
fit.g <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "gompertz")
fit.ln <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "lognormal")
fit.ev <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "ev")
check.dist(fit.cr, fit.w)
check.dist(fit.cr, fit.g)
check.dist(fit.cr, fit.ln)
check.dist(fit.cr, fit.ev)
par(oldpar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("check.surv")
### * check.surv

flush(stderr()); flush(stdout())

### Name: check.surv
### Title: Check the integrity of survival data.
### Aliases: check.surv
### Keywords: manip survival

### ** Examples

xx <- data.frame(enter = c(0, 1), exit = c(1.5, 3), event = c(0, 1), id =
c(1,1))
check.surv(xx$enter, xx$exit, xx$event, xx$id)



cleanEx()
nameEx("coxreg")
### * coxreg

flush(stderr()); flush(stdout())

### Name: coxreg
### Title: Cox regression
### Aliases: coxreg
### Keywords: survival regression

### ** Examples


 dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
                status=c(1,1,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
 coxreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
 # Same as:
 rs <- risksets(Surv(dat$time, dat$status), strata = dat$sex)
 coxreg( Surv(time, status) ~ x, data = dat, rs = rs) #stratified model
 


cleanEx()
nameEx("coxreg.fit")
### * coxreg.fit

flush(stderr()); flush(stdout())

### Name: coxreg.fit
### Title: Cox regression
### Aliases: coxreg.fit
### Keywords: survival regression

### ** Examples

 X <- as.matrix(data.frame(
                x=     c(0, 2,1,4,1,0,3),
                sex=   c(1, 0,0,0,1,1,1)))
 time <- c(1,2,3,4,5,6,7)
 status <- c(1,1,1,0,1,1,0)
 stratum <- rep(1, length(time))

 coxreg.fit(X, Surv(time, status), strats = stratum, max.survs = 6,
     control = list(eps=1.e-4, maxiter = 10, trace = FALSE))



cleanEx()
nameEx("cro")
### * cro

flush(stderr()); flush(stdout())

### Name: cro
### Title: Creates a minimal representation of a data frame.
### Aliases: cro
### Keywords: manip

### ** Examples

dat <- data.frame(y = c(1.1, 2.3, 0.7), x1 = c(1, 0, 1), x2 = c(0, 1, 0))
cro(dat)



cleanEx()
nameEx("fert")
### * fert

flush(stderr()); flush(stdout())

### Name: fert
### Title: Marital fertility nineteenth century
### Aliases: fert
### Keywords: datasets

### ** Examples

data(fert)
fit <- coxreg(Surv(next.ivl, event) ~ ses + prev.ivl, data = fert, subset =
(parity == 1))
drop1(fit, test = "Chisq")



cleanEx()
nameEx("ghq")
### * ghq

flush(stderr()); flush(stdout())

### Name: ghq
### Title: Gauss-Hermite
### Aliases: ghq
### Keywords: math

### ** Examples

ghq(15, FALSE)



cleanEx()
nameEx("glmmML")
### * glmmML

flush(stderr()); flush(stdout())

### Name: glmmML
### Title: Generalized Linear Models with random intercept
### Aliases: glmmML
### Keywords: regression

### ** Examples

id <- factor(rep(1:20, rep(5, 20)))
y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
x <- rnorm(100)
dat <- data.frame(y = y, x = x, id = id)
glmmML(y ~ x, data = dat, cluster = id)



cleanEx()
nameEx("glmmML.fit")
### * glmmML.fit

flush(stderr()); flush(stdout())

### Name: glmmML.fit
### Title: Generalized Linear Model with random intercept
### Aliases: glmmML.fit
### Keywords: regression

### ** Examples

x <- cbind(rep(1, 14), rnorm(14))
y <- rbinom(14, prob = 0.5, size = 1)
id <- rep(1:7, 2)

glmmML.fit(x, y, cluster = id)





cleanEx()
nameEx("glmmboot")
### * glmmboot

flush(stderr()); flush(stdout())

### Name: glmmboot
### Title: Generalized Linear Models with fixed effects grouping
### Aliases: glmmboot
### Keywords: regression nonlinear

### ** Examples

## Not run:
id <- factor(rep(1:20, rep(5, 20)))
y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
x <- rnorm(100)
dat <- data.frame(y = y, x = x, id = id)
res <- glmmboot(y ~ x, cluster = id, data = dat, boot = 500)
## End(Not run)
##system.time(res.glm <- glm(y ~ x + id, family = binomial))



cleanEx()
nameEx("glmmbootFit")
### * glmmbootFit

flush(stderr()); flush(stdout())

### Name: glmmbootFit
### Title: Generalized Linear Models with fixed effects grouping
### Aliases: glmmbootFit
### Keywords: regression nonlinear

### ** Examples

## Not run
x <- matrix(rnorm(1000), ncol = 1)
id <- rep(1:100, rep(10, 100))
y <- rbinom(1000, size = 1, prob = 0.4)
fit <- glmmbootFit(x, y, cluster = id, boot = 200)
summary(fit)
## End(Not run)
## Should show no effects.



cleanEx()
nameEx("hweibull")
### * hweibull

flush(stderr()); flush(stdout())

### Name: hweibull
### Title: The (Cumulative) Hazard Function of a Weibull Distribution
### Aliases: hweibull Hweibull
### Keywords: survival

### ** Examples

hweibull(3, 2, 1)
dweibull(3, 2, 1) / pweibull(3, 2, 1, lower.tail = FALSE)
Hweibull(3, 2, 1)
-pweibull(3, 2, 1, log.p = TRUE, lower.tail = FALSE)



cleanEx()
nameEx("infants")
### * infants

flush(stderr()); flush(stdout())

### Name: infants
### Title: Infant mortality and maternal death, Sweeden 1821-1894.
### Aliases: infants
### Keywords: datasets

### ** Examples

data(infants)
fit <- coxreg(Surv(enter, exit, event) ~ strata(stratum) + mother, data
= infants)
fit
fit.w <- phreg(Surv(enter, exit, event) ~ mother + parish + ses, data =
infants)
fit.w ## Weibull proportional hazards model.



cleanEx()
nameEx("logrye")
### * logrye

flush(stderr()); flush(stdout())

### Name: logrye
### Title: Rye prices, Scania, southern Sweden, 1801-1894.
### Aliases: logrye
### Keywords: datasets

### ** Examples

data(logrye)
summary(logrye)



cleanEx()
nameEx("ltx")
### * ltx

flush(stderr()); flush(stdout())

### Name: ltx
### Title: LaTeX printing of coxreg results.
### Aliases: ltx
### Keywords: printing

### ** Examples

data(oldmort)
fit <- coxreg(Surv(enter, exit, event) ~ civ + sex, data = oldmort)
dr <- drop1(fit, test = "Chisq")
ltx(fit, dr = dr, caption = "A test example.", label = "tab:test1") 



cleanEx()
nameEx("make.communal")
### * make.communal

flush(stderr()); flush(stdout())

### Name: make.communal
### Title: Put communals in "fixed" data frame
### Aliases: make.communal
### Keywords: survival

### ** Examples

dat <- data.frame(enter = 0, exit = 5.731, event = 1,
birthdate = 1962.505, x = 2)
## Birth date: July 2, 1962 (approximately).
com.dat <- data.frame(price = c(12, 3, -5, 6, -8, -9, 1, 7))
dat.com <- make.communal(dat, com.dat, start = 1962.000) 



cleanEx()
nameEx("male.mortality")
### * male.mortality

flush(stderr()); flush(stdout())

### Name: male.mortality
### Title: Male mortality in ages 40-60, nineteenth century
### Aliases: male.mortality
### Keywords: datasets

### ** Examples

data(male.mortality)
coxreg(Surv(enter, exit, event) ~ ses, data = male.mortality)



cleanEx()
nameEx("mlreg")
### * mlreg

flush(stderr()); flush(stdout())

### Name: mlreg
### Title: ML proportional hazards regression
### Aliases: mlreg
### Keywords: survival regression

### ** Examples


 dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
                status=c(1,1,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
 mlreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
 # Same as:
 rs <- risksets(Surv(dat$time, dat$status), strata = dat$sex)
 mlreg( Surv(time, status) ~ x, data = dat, rs = rs) #stratified model
 


cleanEx()
nameEx("mort")
### * mort

flush(stderr()); flush(stdout())

### Name: mort
### Title: Male mortality in ages 40-60, nineteenth century
### Aliases: mort
### Keywords: datasets

### ** Examples

data(mort)
coxreg(Surv(enter, exit, event) ~ ses, data = mort)



cleanEx()
nameEx("oldmort")
### * oldmort

flush(stderr()); flush(stdout())

### Name: oldmort
### Title: Old age mortality, Sundsvall, Sweden, 1860-1880.
### Aliases: oldmort
### Keywords: datasets

### ** Examples

data(oldmort)
summary(oldmort)
## maybe str(oldmort) ; plot(oldmort) ...



cleanEx()
nameEx("phreg")
### * phreg

flush(stderr()); flush(stdout())

### Name: phreg
### Title: Parametric Proportional Hazards Regression
### Aliases: phreg
### Keywords: survival regression

### ** Examples

data(mort)
fit <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
fit
plot(fit)
fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
check.dist(fit.cr, fit)



cleanEx()
nameEx("plot.Surv")
### * plot.Surv

flush(stderr()); flush(stdout())

### Name: plot.Surv
### Title: Plots of survivor functions.
### Aliases: plot.Surv
### Keywords: survival

### ** Examples

time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50)
plot.Surv(Surv(time0, time1, event), strata = group)



cleanEx()
nameEx("plot.aftreg")
### * plot.aftreg

flush(stderr()); flush(stdout())

### Name: plot.aftreg
### Title: Plots output from a Weibull regression
### Aliases: plot.aftreg
### Keywords: dplot survival

### ** Examples

y <- rllogis(40, shape = 1, scale = 1)
x <- rep(c(1,1,2,2), 10)
fit <- aftreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
plot(fit)



cleanEx()
nameEx("plot.coxreg")
### * plot.coxreg

flush(stderr()); flush(stdout())

### Name: plot.coxreg
### Title: Plots of survivor functions.
### Aliases: plot.coxreg
### Keywords: survival

### ** Examples

time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50) 
fit <- coxreg(Surv(time0, time1, event) ~ strata(group))
plot.coxreg(fit)



cleanEx()
nameEx("plot.hazdata")
### * plot.hazdata

flush(stderr()); flush(stdout())

### Name: plot.hazdata
### Title: Plots of survivor functions.
### Aliases: plot.hazdata
### Keywords: survival

### ** Examples

time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50)
fit <- coxreg(Surv(time0, time1, event) ~ strata(group))
plot(fit$hazards)



cleanEx()
nameEx("plot.phreg")
### * plot.phreg

flush(stderr()); flush(stdout())

### Name: plot.phreg
### Title: Plots output from a Weibull regression
### Aliases: plot.phreg
### Keywords: dplot survival

### ** Examples

y <- rllogis(40, shape = 1, scale = 1)
x <- rep(c(1,1,2,2), 10)
fit <- phreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
plot(fit)



cleanEx()
nameEx("plot.weibreg")
### * plot.weibreg

flush(stderr()); flush(stdout())

### Name: plot.weibreg
### Title: Plots output from a Weibull regression
### Aliases: plot.weibreg
### Keywords: dplot survival

### ** Examples

y <- rweibull(4, shape = 1, scale = 1)
x <- c(1,1,2,2)
fit <- weibreg(Surv(y, c(1,1,1,1)) ~ x)
plot(fit)



cleanEx()
nameEx("risksets")
### * risksets

flush(stderr()); flush(stdout())

### Name: risksets
### Title: Finds the compositions and sizes of risk sets
### Aliases: risksets
### Keywords: survival

### ** Examples

 enter <- c(0, 1, 0, 0)
 exit <- c(1, 2, 3, 4)
 event <- c(1, 1, 1, 0)
 risksets(Surv(enter, exit, event))



cleanEx()
nameEx("scania")
### * scania

flush(stderr()); flush(stdout())

### Name: scania
### Title: Old age mortality, Scania, southern Sweden, 1813-1894.
### Aliases: scania
### Keywords: datasets

### ** Examples

data(scania)
summary(scania)



cleanEx()
nameEx("summary.aftreg")
### * summary.aftreg

flush(stderr()); flush(stdout())

### Name: summary.aftreg
### Title: Prints aftreg objects
### Aliases: summary.aftreg
### Keywords: survival print

### ** Examples

## The function is currently defined as
function (object, ...) 
print(object)



cleanEx()
nameEx("summary.coxreg")
### * summary.coxreg

flush(stderr()); flush(stdout())

### Name: summary.coxreg
### Title: Prints coxreg objects
### Aliases: summary.coxreg
### Keywords: survival print

### ** Examples

## The function is currently defined as
function (object, ...) 
print(object)



cleanEx()
nameEx("summary.phreg")
### * summary.phreg

flush(stderr()); flush(stdout())

### Name: summary.phreg
### Title: Prints phreg objects
### Aliases: summary.phreg
### Keywords: survival print

### ** Examples

## The function is currently defined as
function (object, ...) 
print(object)



cleanEx()
nameEx("summary.weibreg")
### * summary.weibreg

flush(stderr()); flush(stdout())

### Name: summary.weibreg
### Title: Prints a weibreg object
### Aliases: summary.weibreg
### Keywords: survival print

### ** Examples

## The function is currently defined as
function (object, ...) 
print(object)



cleanEx()
nameEx("swe07")
### * swe07

flush(stderr()); flush(stdout())

### Name: swe07
### Title: Swedish population and deaths in ages 61-80, 2007.
### Aliases: swe07
### Keywords: datasets

### ** Examples

data(swe07)
fit <- glm(deaths ~ offset(log.pop) + sex * as.factor(age), family = poisson, data = swe07)
drop1(fit, test = "Chisq")  ## Proportional hazards?



cleanEx()
nameEx("table.events")
### * table.events

flush(stderr()); flush(stdout())

### Name: table.events
### Title: Calculating failure times, risk set sizes and No. of events in
###   each risk set
### Aliases: table.events
### Keywords: survival

### ** Examples

exit = c(1,2,3,4,5)
event = c(1,1,0,1,1)
table.events(exit = exit, event = event)



cleanEx()
nameEx("toBinary")
### * toBinary

flush(stderr()); flush(stdout())

### Name: toBinary
### Title: Transforms a "survival" data frame into a data frame suitable
###   for binary (logistic) regression
### Aliases: toBinary
### Keywords: survival cluster

### ** Examples

enter <- rep(0, 4)
exit <- 1:4
event <- rep(1, 4)
z <- rep(c(-1, 1), 2)
dat <- data.frame(enter, exit, event, z)
binDat <- toBinary(dat)
dat
binDat
coxreg(Surv(enter, exit, event) ~ z, method = "ml", data = dat)
## Same as:
summary(glm(event ~ z + riskset, data = binDat, family = binomial(link = cloglog)))



cleanEx()
nameEx("toDate")
### * toDate

flush(stderr()); flush(stdout())

### Name: toDate
### Title: Convert time in years since "0000-01-01" to a date.
### Aliases: toDate
### Keywords: survival

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
toDate(1897.357)



cleanEx()
nameEx("toTime")
### * toTime

flush(stderr()); flush(stdout())

### Name: toTime
### Title: Calculate duration in years from "0000-01-01" to a given date
### Aliases: toTime
### Keywords: survival

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
toTime(c("1897-05-16", "1901-11-21"))



cleanEx()
nameEx("weibreg")
### * weibreg

flush(stderr()); flush(stdout())

### Name: weibreg
### Title: Weibull Regression
### Aliases: weibreg
### Keywords: survival regression

### ** Examples

 dat <- data.frame(time = c(4, 3, 1, 1, 2, 2, 3),
                status = c(1, 1, 1, 0, 1, 1, 0),
                x = c(0, 2, 1, 1, 1, 0, 0),
                sex = c(0, 0, 0, 0, 1, 1, 1))
 weibreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
