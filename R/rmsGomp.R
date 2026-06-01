#' Restricted mean survival calculation for the Gompertz distribution by
#' numerical integration of the survival function.
#' 
#' @param from Start time for the integration
#' @param to Stop time of the integration.
#' @param scale The scale parameter.
#' @param shape The shape parameter.
#' @param param Which parametrization? Defaults to 'canonical', the other choice
#' is 'rate' (not developed yet, avoid until next version of eha!).
#' @param fit The result of a call to 'aftGompreg'. If not missing, scale and 
#' shape will be taken from the fit, and corresponding argumens in the call will 
#' be ignored.
#' @examples
#' res <- aftGompreg(Surv(enter-60, exit-60, event) ~ ses.50 + strata(sex), data = oldmort)
#' rmsGomp(from = 0, to = 40, fit = res)
#' @export 

rmsGomp <- function(from = 0, to = Inf, scale = 1, shape = 1,
                    param = c("canonical", "rate"), fit){
### For "canonical" parametrization (mandatory in AFT Gompertz regression)
### or "rate" parametrization (best in PH Gompertz regression).
    ##
    ## NOTE: if 'fit' is mot missing, all parameter values are taken
    ## from 'fit'. param == "lifeAcc" or "lifeExp" implies
    ## param = "canonical".  
    if (!is.numeric(scale)) stop("scale non-numeric")
    if (!missing(fit)){
        scale <- fit$scale
        shape <- fit$shape
        param <- fit$param
        if (param %in% c("lifeExp", "lifeAcc")) param <- "canonical"
        if (param != "canonical") stop("Only canonical for the moment...")
        if (!is.na(fit$n.strata) & (fit$n.strata >= 2)){
            names(scale) <- fit$strata
            names(shape) <- fit$strata
        }
    }else{
        param <- param[1]
    }
### NOTE:
    if (param == "rate") stop ("param = 'rate' not implemented yet.")
###    
    n <- length(scale)
    if (!(length(shape) == n)) stop("scale and shape have different lengths.")
    if (n == 1){    
        fun <- function(x) pgompertz(x, scale = scale, shape = shape,
                                     param = param,
                                     lower.tail = FALSE)
        out <- integrate(fun, from, to)$value
    }else{
        out <- numeric(n)
        for (i in 1:n){
            fun <- function(x) pgompertz(x, scale = scale[i], shape = shape[i],
                                         param = param,
                                         lower.tail = FALSE)
            out[i] <- integrate(fun, from, to)$value
        }
        names(out) <- fit$strata
    }
    out
}
