#' Summary for tpchreg objects
#' 
#' @param object A \code{tpchreg} object.
#' @param ci Logical. If TRUE, confidence limits are given instead of se's.
#' @param level Confidence level, used if ci.
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{tpchreg}}
#' @keywords survival summary
#' @examples
#' 
#' ## The function is currently defined as
#' ## function (object, ...) 
#' 
#' @export
summary.tpchreg <- function(object, ci = FALSE, level = 0.95, ...){
    haz <- object$hazards
    n <- length(object$cuts)
    
    shift <- object$cuts[1]
    new.cuts <- object$cuts[-1] - shift
    
    new.last.value <- new.cuts[n - 1]
    new.cuts <- new.cuts[-(n - 1)]
    
    if (object$n.strata == 1){
        rmean <- integrate(ppch, 0, new.last.value, cuts = new.cuts, 
                           levels = object$hazards, lower.tail = FALSE)$value
        psurv <- ppch(new.last.value, cuts = new.cuts, 
                      levels = object$hazards, lower.tail = FALSE)
    }else{
        rmean <- rep(0, object$n.strata)
        psurv <- rep(0, object$n.strata)
        for (i in 1:object$n.strata){
            rmean[i] <- integrate(ppch, 0, new.last.value, cuts = new.cuts, 
                                 levels = object$hazards[i, ], 
                                 lower.tail = FALSE)$value
            psurv[i] <- ppch(new.last.value, cuts = new.cuts, 
                          levels = object$hazards[i, ], lower.tail = FALSE)
        }
        names(rmean) <- object$sstrata
        names(psurv) <- names(rmean)
    }
    if (!object$nullModel){
        if (is.null(object$var)){
            sd <- rep(NA, length(object$coefficients))
        }else{
            sd <- sqrt(diag(object$var))
        }
        object$ci <- ci
        if (ci){
            coef <- object$coefficients
            clim <- qnorm(1 - (1 - level / 2))
            low <- exp(coef - clim * sd)
            high <- exp(coef + clim * sd)
            coefficients <- cbind(exp(coef), low, high)
        }else{
            coefficients <- cbind(object$coefficients, 
                                  exp(object$coefficients),
                                  sd)
            zval <- coefficients[, 1] / coefficients[, 3]
            pval <- pchisq(zval^2, df = 1, lower.tail = FALSE )
            coefficients <- cbind(coefficients, zval, pval)
            colnames(coefficients) <- c("coef", "exp(coef)", "se(coef)", "z", "Wald p")
            rownames(coefficients) <- names(object$coefficients)
        }
        object$coefficients <- coefficients
    
        object$dr <- drop1(object, test = "Chisq")
    }
    
    object$rmean <- rmean
    object$psurv <- psurv
    
    class(object) <- "summary.tpchreg"
    object
}
