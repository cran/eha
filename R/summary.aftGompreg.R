#' Creates a summary of aftGompreg objects
#' 
#' 
#' @param object A \code{aftGompreg} object
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{summary.aftreg}}
#' @keywords survival summary
#' @examples
#' 
#' ## The function is currently defined as
#' function (object, ...) 
#' print(object)
#' 
#' @export
summary.aftGompreg <- function(object, ...){
    dr <- drop1(object, test = "Chisq")
    object$dr <- dr
    ncoef <- object$df
    ## Split coefficients into regression parameters and hazard ditto.
    
    hazards <- object$coefficients[-(1:ncoef)]
    
    coefficients <- object$coefficients[1:ncoef]
    
    ## Regression parameters:
    rawnames <- names(coefficients)
    varcoef <- diag(object$var[1:ncoef, 1:ncoef, drop = FALSE])
    varhaz <- diag(object$var[-(1:ncoef), -(1:ncoef), drop = FALSE])
    class(object) <- "summary.aftGompreg"
    coefficients <- cbind(coefficients, 
                          exp(coefficients),
                          sqrt(varcoef))
    zval <- coefficients[, 1] / coefficients[, 3]
    pval <- pchisq(zval^2, df = 1, lower.tail = FALSE )
    coefficients <- cbind(coefficients, zval, pval)
    colnames(coefficients) <- c("coef", "exp(coef)", "se(coef)", "z", "Wald p"
    )
    rownames(coefficients) <- rawnames
    
    ## Hazard parameters:
    haznames <- names(hazards)
    hazards <- cbind(hazards, sqrt(varhaz))
    colnames(hazards) <- c("par", "se(par)")
    rownames(hazards) <- haznames
    
    ##list(fit = object, coefficients = coefficients)
    object$coefficients <- coefficients
    object$hazards <- hazards
    class(object) <- c("summary.aftGompreg")
    object
}
